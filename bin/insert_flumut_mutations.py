#!/usr/bin/env python3
"""Insert flumut amino-acid marker mutations into a Nextclade freebayes-style VCF.

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
from pathlib import Path

# Reuse helpers from the existing (untouched) converter. Both scripts are staged from
# bin/, so adding this script's directory to sys.path makes the import resolvable.
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from nextclade_to_freebayes_vcf import (  # noqa: E402
    _build_ann_entry,
    _cds_to_genomic,
    _format_hgvs_p_sub,
    load_gff3,
    variant_row,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

QUAL = 60

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
# VCF parsing / editing
# ---------------------------------------------------------------------------

def parse_vcf(path: str) -> tuple[list, list, str | None]:
    """Return (header_lines, data_rows, chrom).

    header_lines includes the '#CHROM' column line. data_rows are raw (unsplit) strings.
    chrom is taken from the first data row, else from a ##contig header line.
    """
    header_lines: list[str] = []
    data_rows: list[str] = []
    chrom: str | None = None
    contig_chrom: str | None = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header_lines.append(line)
                if line.startswith("##contig=") and contig_chrom is None:
                    mm = re.search(r"ID=([^,>]+)", line)
                    if mm:
                        contig_chrom = mm.group(1)
                continue
            if not line.strip():
                continue
            data_rows.append(line)
            if chrom is None:
                chrom = line.split("\t", 1)[0]
    if chrom is None:
        chrom = contig_chrom
    return header_lines, data_rows, chrom


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


def add_info(row: str, key_value: str) -> str:
    """Append ';key=value' to the INFO column (index 7) of a raw VCF row."""
    cols = row.split("\t")
    cols[7] = f"{cols[7]};{key_value}"
    return "\t".join(cols)


# ---------------------------------------------------------------------------
# Target GFF transcript lookup
# ---------------------------------------------------------------------------

def index_transcripts(transcripts: list) -> dict:
    """Map gene identifiers -> transcript, preferring main (non-mature) CDS transcripts."""
    by_name: dict = {}
    for tx in transcripts:
        # Prefer a real CDS transcript over a mature-peptide entry for the same key.
        for key in (tx.get("gene_name"), tx.get("feature_id"), tx.get("protein_id")):
            if not key:
                continue
            existing = by_name.get(key)
            if existing is None or (existing.get("is_mature") and not tx.get("is_mature")):
                by_name[key] = tx
    return by_name


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def get_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="Existing per-sample/per-segment VCF to augment")
    p.add_argument("--coord-map", required=True, help="coordinate_map.json (flumut -> Nextclade)")
    p.add_argument("--flumut-mutations", required=True, help="Melted flumut mutations TSV (sample/gene/mutation/allele)")
    p.add_argument("--consensus", required=False, help="Cleaned consensus FASTA (context; not required for placement)")
    p.add_argument("--dataset-dir", required=True, help="Target Nextclade dataset dir (for genome_annotation.gff3)")
    p.add_argument("--segment", required=True, help="Segment label, e.g. HA, M, NS")
    p.add_argument("--sample", required=True, help="Sample id (matches flumut TSV 'sample' column)")
    p.add_argument("-o", "--output", required=True, help="Output VCF path")
    return p.parse_args()


def main() -> None:
    args = get_args()
    segment = args.segment.upper()

    header_lines, data_rows, chrom = parse_vcf(args.vcf)
    if chrom is None:
        log.error("Could not determine CHROM from %s (no data rows and no ##contig header)", args.vcf)
        sys.exit(1)

    with open(args.coord_map) as fh:
        coord = json.load(fh)
    seg_map = coord.get("segments", {}).get(segment)
    if not seg_map:
        log.warning("Segment %s not present in coordinate map; copying VCF with fluscan labels only", segment)
        genes_map = {}
    else:
        genes_map = seg_map.get("genes", {})
    available_genes = set(genes_map.keys())

    gff_path = Path(args.dataset_dir) / "genome_annotation.gff3"
    tx_by_key: dict = {}
    if gff_path.exists():
        _, transcripts = load_gff3(str(gff_path))
        tx_by_key = index_transcripts(transcripts)
    else:
        log.warning("No genome_annotation.gff3 in %s; cannot place flumut-only rows", args.dataset_dir)

    # --- fluscan AA changes already present in the VCF ---
    fluscan_changes: set = set()
    for row in data_rows:
        info = row.split("\t")[7]
        fluscan_changes |= row_aa_changes(info)

    # --- flumut mutations mapped into target (Nextclade) coordinates ---
    flumut_keys: set = set()          # (gene, target_pos, alt_aa) both-eligible keys
    new_rows: list[tuple[int, str]] = []  # (pos, raw_vcf_row) for flumut-only inserts
    n_unmapped = 0
    n_skipped_segment = 0

    with open(args.flumut_mutations) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_sample = header.index("sample")
            i_gene = header.index("gene")
            i_mut = header.index("mutation")
        except ValueError:
            log.error("Unexpected flumut mutations header: %s", header)
            sys.exit(1)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_sample, i_gene, i_mut):
                continue
            if parts[i_sample] != args.sample:
                continue
            gene_raw = parts[i_gene]
            gene_key = resolve_gene(gene_raw, available_genes)
            if gene_key is None:
                n_skipped_segment += 1
                continue
            mm = MUTATION_RE.match(parts[i_mut])
            if not mm:
                log.warning("Unparseable flumut mutation '%s' (gene %s)", parts[i_mut], gene_raw)
                continue
            ref_aa, src_pos, alt_aa = mm.group(1), int(mm.group(2)), mm.group(3)

            aa_map = genes_map[gene_key].get("aa_source_to_target", {})
            target_pos = aa_map.get(str(src_pos))
            if target_pos is None:
                n_unmapped += 1
                log.warning("flumut %s:%s%d%s has no target mapping (gapped/unaligned); skipping",
                            gene_key, ref_aa, src_pos, alt_aa)
                continue

            key = (gene_key, int(target_pos), alt_aa)
            flumut_keys.add(key)

            if key in fluscan_changes:
                continue  # detected by both -> existing row is retagged below, no insert

            # flumut-only: build a new row, deriving a real POS from the target CDS.
            tx = tx_by_key.get(gene_key)
            if tx is None:
                log.warning("No target transcript for gene %s; cannot place %s%d%s; skipping",
                            gene_key, ref_aa, target_pos, alt_aa)
                continue
            genomic0 = _cds_to_genomic((int(target_pos) - 1) * 3, tx["segments"], tx["strand"])
            if genomic0 is None:
                log.warning("Target AA %d out of CDS range for gene %s; skipping", target_pos, gene_key)
                continue
            pos = genomic0 + 1
            cds_pos = (int(target_pos) - 1) * 3 + 1
            ann_entry = _build_ann_entry(
                allele="N",
                annotation="missense_variant",
                impact="MODERATE",
                gene_name=gene_key,
                gene_id=tx.get("feature_id", gene_key),
                feature_id=tx.get("feature_id", gene_key),
                transcript_biotype="protein_coding",
                hgvs_c="",
                hgvs_p=_format_hgvs_p_sub(ref_aa, int(target_pos), alt_aa),
                cds_pos=cds_pos,
                cds_len=tx.get("cds_length", 0),
                aa_pos=int(target_pos),
                aa_len=tx.get("aa_length", 0),
                warnings="INFO_FLUMUT_INSERTED_NO_NUCLEOTIDE_EVIDENCE",
            )
            row = variant_row(chrom, pos, "N", "N", QUAL, ann="ANN=" + ann_entry)
            row = add_info(row, "DETECTED_BY=flumut")
            new_rows.append((pos, row))

    # --- retag existing fluscan rows (both vs fluscan) ---
    tagged_rows: list[tuple[int, str]] = []
    for row in data_rows:
        cols = row.split("\t")
        pos = int(cols[1])
        info = cols[7]
        detected = "both" if (row_aa_changes(info) & flumut_keys) else "fluscan"
        tagged_rows.append((pos, add_info(row, f"DETECTED_BY={detected}")))

    all_rows = sorted(tagged_rows + new_rows, key=lambda t: (t[0], t[1]))

    # --- write output: original header + DETECTED_BY INFO line, then sorted rows ---
    detected_by_line = (
        '##INFO=<ID=DETECTED_BY,Number=1,Type=String,'
        'Description="Tool(s) detecting this AA change: flumut, fluscan, or both">'
    )
    out_header: list[str] = []
    inserted = False
    for line in header_lines:
        if line.startswith("#CHROM") and not inserted:
            out_header.append(detected_by_line)
            inserted = True
        out_header.append(line)
    if not inserted:  # no #CHROM line found (degenerate); append at end of headers
        out_header.append(detected_by_line)

    with open(args.output, "w") as out:
        out.write("\n".join(out_header) + "\n")
        for _, row in all_rows:
            out.write(row + "\n")

    log.info(
        "Wrote %s: %d existing rows retagged, %d flumut-only rows inserted "
        "(%d flumut marker(s) matched fluscan/both, %d unmapped, %d other-segment) for %s/%s",
        args.output, len(tagged_rows), len(new_rows),
        len(flumut_keys & fluscan_changes), n_unmapped, n_skipped_segment,
        args.sample, segment,
    )


if __name__ == "__main__":
    main()
