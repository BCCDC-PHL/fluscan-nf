#!/usr/bin/env python3
"""Merge flumut amino-acid marker mutations into a Nextclade freebayes-style VCF.

Reads an existing per-sample/per-segment VCF (produced upstream by
nextclade_to_freebayes_vcf.py -> fix_mnv_codons.py), the flumut->Nextclade coordinate
map JSON (build_reference_coordinate_map.py), and flumut's per-sample melted mutations
TSV. It labels every existing AA change with which tool(s) found it and appends new rows
for flumut markers that Nextclade/fluscan did not call.

Each variant row gets an INFO field:
    DETECTED_BY=fluscan   variant/AA change called by the Nextclade (fluscan) branch only
    DETECTED_BY=flumut     AA marker called by flumut only (a newly inserted row)
    DETECTED_BY=both       the same exact AA substitution called by both tools

flumut only reports amino-acid substitutions, so newly inserted flumut-only rows carry
filler nucleotide alleles (REF=N, ALT=N); their POS is derived from the coordinate map
(flumut AA -> Nextclade AA) and the target GFF (Nextclade AA -> nucleotide codon start).
The amino-acid change itself is carried in the ANN= annotation, matching the snpEff-style
format used by the existing converter.

This script does not modify any existing file; it writes a new VCF.

Pipeline position: after post_process_vcf_batch (nextclade workflow), consuming
coordinate_conversion + flumut outputs.
"""

import argparse
import json
import logging
import os
import re
import sys

import pandas as pd

# Reuse helpers from the existing (untouched) converter and the shared VCF I/O in tools.py.
# Both live in bin/, so adding this script's directory to sys.path makes them importable.
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from nextclade_to_freebayes_vcf import (  # noqa: E402
    _build_ann_entry,
    _cds_to_genomic,
    _format_hgvs_p_sub,
    variant_row,
)
from tools import read_vcf, write_vcf  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

QUAL = 60

DETECTED_BY_HEADER = (
    '##INFO=<ID=DETECTED_BY,Number=1,Type=String,'
    'Description="Tool(s) detecting this AA change: flumut, fluscan, or both">'
)

# flumut protein names carry a trailing "-<digits>" that is either a FluMutDB reference
# version (HA1-5 -> HA1) or, for NS, a protein number (NS-1 -> NS1, NS-2 -> NEP). Rather
# than guess, gene resolution generates candidate names and keeps whichever exists as a
# gene key in the coordinate map for the segment being processed. This alias table covers
# names that differ between flumut and the Nextclade annotation.
GENE_ALIASES = {
    "NS2": ["NEP"],
    "NEP": ["NS2"],
    "M": ["MP"],
    "MP": ["M"],
}

HGVS_P_RE = re.compile(r"p\.([A-Za-z*])(\d+)([A-Za-z*])")
MUTATION_RE = re.compile(r"^([A-Za-z*])(\d+)([A-Za-z*])$")


# ---------------------------------------------------------------------------
# Gene-name resolution
# ---------------------------------------------------------------------------

def resolve_gene(raw: str, available_genes: set) -> str | None:
    """Resolve a flumut gene name to a coordinate-map gene key for this segment.

    Returns the matching key, or None if the gene does not belong to this segment
    (or cannot be reconciled with the annotation).
    """
    candidates: list[str] = [raw]
    m = re.match(r"^(.*?)-(\d+)$", raw)
    if m:
        base, suffix = m.group(1), m.group(2)
        candidates.append(base)            # HA1-5 -> HA1 ; NS-1 -> NS
        candidates.append(base + suffix)   # NS-1 -> NS1 ; NS-2 -> NS2
    for c in list(candidates):
        candidates.extend(GENE_ALIASES.get(c, []))
    for c in candidates:
        if c in available_genes:
            return c
    return None


# ---------------------------------------------------------------------------
# VCF annotation helpers
# ---------------------------------------------------------------------------

def ann_entries(info: str) -> list[list[str]]:
    """Extract ANN entries (each a list of fields) from a VCF INFO string."""
    for token in info.split(";"):
        if token.startswith("ANN="):
            return [entry.split("|") for entry in token[4:].split(",") if entry]
    return []


def row_aa_changes(info: str) -> set:
    """Set of (gene, aa_pos, alt_aa) AA substitutions annotated on a VCF row."""
    changes = set()
    for fields in ann_entries(info):
        if len(fields) < 11:
            continue
        gene = fields[3]
        hp = HGVS_P_RE.search(fields[10])
        if not hp:
            continue
        ref_aa, pos, alt_aa = hp.group(1), int(hp.group(2)), hp.group(3)
        if ref_aa == alt_aa:  # synonymous
            continue
        changes.add((gene, pos, alt_aa))
    return changes


def contig_chrom(header_lines: list) -> str | None:
    """CHROM id from the ##contig header line (fallback when the VCF has no data rows)."""
    for line in header_lines:
        if line.startswith("##contig="):
            m = re.search(r"ID=([^,>]+)", line)
            if m:
                return m.group(1)
    return None


# ---------------------------------------------------------------------------
# Inputs (coordinate map, target GFF, flumut mutations)
# ---------------------------------------------------------------------------

def load_segment_genes(coord_map_path: str, segment: str) -> dict:
    """Per-gene aa/nt maps for one segment from the coordinate map JSON ({} if absent)."""
    with open(coord_map_path) as fh:
        coord = json.load(fh)
    seg_map = coord.get("segments", {}).get(segment)
    if not seg_map:
        log.warning("Segment %s not present in coordinate map; copying VCF with fluscan labels only", segment)
        return {}
    return seg_map.get("genes", {})


def read_flumut_mutations(path: str, sample: str) -> pd.DataFrame:
    """Read the melted flumut mutations TSV, restricted to one sample.

    Returns a DataFrame with columns sample/gene/mutation/allele (as strings); empty when
    the sample has no mutations.
    """
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = {"sample", "gene", "mutation"} - set(df.columns)
    if missing:
        log.error("flumut mutations TSV %s missing column(s): %s", path, ", ".join(sorted(missing)))
        sys.exit(1)
    return df[df["sample"] == sample].reset_index(drop=True)


# ---------------------------------------------------------------------------
# flumut marker -> VCF rows
# ---------------------------------------------------------------------------

def flumut_only_row(chrom: str, gene_key: str, ref_aa: str, target_pos: int,
                    alt_aa: str, target_cds: dict | None) -> list[str] | None:
    """Build a flumut-only VCF row (column list) with filler nucleotides and an AA ANN.

    POS is the target CDS codon start, derived from the target_cds placement carried in the
    coordinate map (no GFF needed). REF/ALT are filler (N). Returns None if the marker
    cannot be placed on the target CDS.
    """
    if not target_cds:
        log.warning("No target CDS placement for gene %s; cannot place %s%d%s; skipping",
                    gene_key, ref_aa, target_pos, alt_aa)
        return None
    segments = [tuple(seg) for seg in target_cds["segments"]]
    genomic0 = _cds_to_genomic((target_pos - 1) * 3, segments, target_cds["strand"])
    if genomic0 is None:
        log.warning("Target AA %d out of CDS range for gene %s; skipping", target_pos, gene_key)
        return None
    feature_id = target_cds.get("feature_id", gene_key)
    ann_entry = _build_ann_entry(
        allele="N",
        annotation="missense_variant",
        impact="MODERATE",
        gene_name=gene_key,
        gene_id=feature_id,
        feature_id=feature_id,
        transcript_biotype="protein_coding",
        hgvs_c="",
        hgvs_p=_format_hgvs_p_sub(ref_aa, target_pos, alt_aa),
        cds_pos=(target_pos - 1) * 3 + 1,
        cds_len=target_cds.get("cds_length", 0),
        aa_pos=target_pos,
        aa_len=target_cds.get("aa_length", 0),
        warnings="INFO_FLUMUT_INSERTED_NO_NUCLEOTIDE_EVIDENCE",
    )
    # Reuse variant_row for the freebayes INFO/FORMAT shape, then split to columns and
    # append the DETECTED_BY tag to the INFO field (column index 7).
    cols = variant_row(chrom, genomic0 + 1, "N", "N", QUAL, ann="ANN=" + ann_entry).split("\t")
    cols[7] = cols[7] + ";DETECTED_BY=flumut"
    return cols


def build_flumut_inserts(flumut_df: pd.DataFrame, genes_map: dict, available_genes: set,
                         fluscan_changes: set, chrom: str) -> tuple[set, list, dict]:
    """Map flumut markers into target coordinates.

    Returns (flumut_keys, new_rows, stats):
      flumut_keys : {(gene, target_pos, alt_aa)} markers that resolved to this segment
      new_rows    : VCF column lists for flumut-only markers (not already called by fluscan)
      stats       : {'unmapped': int, 'other_segment': int}
    """
    flumut_keys: set = set()
    new_rows: list[list[str]] = []
    stats = {"unmapped": 0, "other_segment": 0}

    for gene_raw, mutation in zip(flumut_df["gene"], flumut_df["mutation"]):
        gene_key = resolve_gene(gene_raw, available_genes)
        if gene_key is None:
            stats["other_segment"] += 1
            continue
        mm = MUTATION_RE.match(mutation)
        if not mm:
            log.warning("Unparseable flumut mutation '%s' (gene %s)", mutation, gene_raw)
            continue
        ref_aa, src_pos, alt_aa = mm.group(1), int(mm.group(2)), mm.group(3)

        target_pos = genes_map[gene_key].get("aa_source_to_target", {}).get(str(src_pos))
        if target_pos is None:
            stats["unmapped"] += 1
            log.warning("flumut %s:%s%d%s has no target mapping (gapped/unaligned); skipping",
                        gene_key, ref_aa, src_pos, alt_aa)
            continue

        key = (gene_key, int(target_pos), alt_aa)
        flumut_keys.add(key)
        if key in fluscan_changes:
            continue  # detected by both -> existing row retagged by caller, no insert

        row = flumut_only_row(chrom, gene_key, ref_aa, int(target_pos), alt_aa,
                              genes_map[gene_key].get("target_cds"))
        if row is not None:
            new_rows.append(row)

    return flumut_keys, new_rows, stats


# ---------------------------------------------------------------------------
# Output assembly
# ---------------------------------------------------------------------------

def assemble_output(df: pd.DataFrame, new_rows: list, flumut_keys: set,
                    row_changes: list) -> pd.DataFrame:
    """Label existing rows (both/fluscan), append flumut-only rows, and sort by POS.

    row_changes holds the pre-parsed AA-change set for each existing row (index-aligned
    with df), so a row's ANN is parsed once rather than twice.
    """
    labels = ["both" if (rc & flumut_keys) else "fluscan" for rc in row_changes]
    df = df.copy()
    df["INFO"] = df["INFO"].astype(str) + ";DETECTED_BY=" + pd.Series(labels, index=df.index)

    if new_rows:
        new_df = pd.DataFrame(new_rows, columns=list(df.columns))
        out_df = pd.concat([df, new_df], ignore_index=True)
    else:
        out_df = df
    out_df["POS"] = out_df["POS"].astype(int)
    return out_df.sort_values("POS", kind="stable").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def get_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="Existing per-sample/per-segment VCF to augment")
    p.add_argument("--coord-map", required=True, help="coordinate_map.json (flumut -> Nextclade)")
    p.add_argument("--flumut-mutations", required=True, help="Melted flumut mutations TSV (sample/gene/mutation/allele)")
    p.add_argument("--consensus", required=False, help="Cleaned consensus FASTA (context; not required for placement)")
    p.add_argument("--segment", required=True, help="Segment label, e.g. HA, M, NS")
    p.add_argument("--sample", required=True, help="Sample id (matches flumut TSV 'sample' column)")
    p.add_argument("-o", "--output", required=True, help="Output VCF path")
    return p.parse_args()


def main() -> None:
    args = get_args()
    segment = args.segment.upper()

    df, header_lines = read_vcf(args.vcf)
    chrom = str(df["CHROM"].iloc[0]) if len(df) else contig_chrom(header_lines)
    if chrom is None:
        log.error("Could not determine CHROM from %s (no data rows and no ##contig header)", args.vcf)
        sys.exit(1)

    genes_map = load_segment_genes(args.coord_map, segment)
    available_genes = set(genes_map.keys())

    # AA changes already called by the fluscan (Nextclade) branch. Parse each row's ANN
    # once and reuse the per-row sets for both the union here and the labelling later.
    row_changes = [row_aa_changes(info) for info in df["INFO"].astype(str)]
    fluscan_changes: set = set().union(*row_changes) if row_changes else set()

    # flumut markers mapped into Nextclade coordinates (placement comes from the map's
    # per-gene target_cds block, so no GFF is needed here).
    flumut_df = read_flumut_mutations(args.flumut_mutations, args.sample)
    flumut_keys, new_rows, stats = build_flumut_inserts(
        flumut_df, genes_map, available_genes, fluscan_changes, chrom,
    )

    # Label, append, sort, and write (original *_final.vcf is left untouched).
    out_df = assemble_output(df, new_rows, flumut_keys, row_changes)
    out_header = header_lines[:-1] + [DETECTED_BY_HEADER, header_lines[-1]]
    write_vcf(out_df, out_header, args.output)

    log.info(
        "Wrote %s: %d existing rows retagged, %d flumut-only rows inserted "
        "(%d flumut marker(s) matched fluscan/both, %d unmapped, %d other-segment) for %s/%s",
        args.output, len(df), len(new_rows),
        len(flumut_keys & fluscan_changes), stats["unmapped"], stats["other_segment"],
        args.sample, segment,
    )


if __name__ == "__main__":
    main()
