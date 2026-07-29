#!/usr/bin/env python3
"""Scan a final annotated VCF for the substitutions on a custom watchlist.

This reads the VCF the pipeline actually publishes, not a derived TSV. The chain it
replaces went VCF -> summarize_vcf.py -> *_mutations.tsv -> here, which meant a second
ANN= parser, a second HGVS.p parser and a second idea of what a gene is called -- and an
intermediate file whose publishDir was commented out, so nothing about the result was
reviewable. Reading the VCF directly also means this tool runs post hoc against any
published *_final.vcf without a pipeline run.

    mutation_watcher.py --vcf sample_HA_final.vcf --watchlist nai.csv \\
                        --sample SAMPLE-ID --segment HA -o hits.tsv

Watchlist hits are deliberately NOT written back into the VCF. FLUMUT= is a curated,
always-on annotation and belongs there permanently; a watchlist is a list we set on our
own terms, and writing it into the VCF would make the same consensus produce different
VCFs depending on the contents of a CSV. What the two share is the mechanics of reading
ANN=, and those live once, in lib/aa_markers.py.

Matching is on (gene, position, alt_aa). The watchlist's REF residue is reported and
checked, never matched on -- see bin/validate_watchlist.py.

Gene resolution goes through the canonical alias table, then falls back to the ANN entry's
Feature_ID. That fallback is what makes NCBI-style HA annotations work: a mature peptide's
ANN entry names the *parent* gene (HA) and only its Feature_ID names the subunit (HA1). It
never runs in reverse. A watchlist saying HA must not match an HA1 annotation, because HA1
numbering starts after HA's 16-residue signal peptide -- keying an HA1 marker on a bare HA
CDS returns the wrong residue at 93% of positions, silently (ISSUES.md I-4).

Rows flagged MNV_UNRESOLVED by fix_mnv_codons.py are skipped: the pipeline has already
determined their HGVS.p to be wrong, so a hit there would be a false positive.

Pipeline position: after post_process_vcf / post_process_vcf_batch, on the final VCF.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
from pathlib import Path

# bin/ holds entry points and is on PATH for every task; the modules they import live in
# bin/lib/. Python already puts a script's own directory on sys.path, but this insert uses
# realpath so the import still resolves when bin/ is reached through a symlink.
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from lib.aa_markers import (  # noqa: E402
    aa_changes_from_info,
    ann_entries,
    has_info_flag,
    parse_hgvs_p,
)
from lib.flu_gene_names import canonical_gene  # noqa: E402
from lib.tools import read_vcf  # noqa: E402
from validate_watchlist import (  # noqa: E402
    WatchlistEntry,
    WatchlistError,
    check_reference_residues,
    load_reference_residues,
    load_watchlist,
    report_conflicts,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

SCRIPT_VERSION = "2.0.0"

OUTPUT_COLUMNS = (
    "sample_id",
    "segment",
    "gene",
    "position",
    "ref_aa",
    "alt_aa",
    "mutation",
    "chrom",
    "pos",
    "hgvs_p",
    "note",
)

COVERAGE_COLUMNS = ("sample_id", "segment", "gene", "in_vcf", "watchlist_entries", "hits")


def resolve_gene(gene: str, feature_id: str, watched_genes: set[str]) -> str | None:
    """The watchlist gene key an ANN entry should be looked up under, or None.

    ``gene`` (ANN Gene_Name) is tried first, then ``feature_id`` (ANN Feature_ID). The
    fallback exists for mature peptides, whose ANN entry carries the parent gene's name --
    an NCBI-style HA annotation emits Gene_Name=HA, Feature_ID=HA1 for an HA1 change, and
    only the second names the protein the numbering belongs to.

    The fallback is strictly subunit-ward. There is no path from a parent gene to one of
    its subunits, by design: HA and HA1 disagree by a 16-residue signal peptide, so letting
    a watchlist's HA entry match an HA1 annotation would silently report the wrong residue
    at 93% of positions rather than report nothing.
    """
    for candidate in (gene, feature_id):
        if not candidate:
            continue
        key = canonical_gene(candidate)
        if key in watched_genes:
            return key
    return None


def scan_vcf(
    vcf_path: Path,
    entries: list[WatchlistEntry],
    sample_id: str,
    segment: str,
) -> tuple[list[dict], dict[str, int], int]:
    """Match watchlist entries against a VCF's annotations.

    Returns (hit rows, {gene: annotated changes seen}, MNV_UNRESOLVED rows skipped). The
    per-gene count is what backs the coverage report: a watchlist gene that never appears
    in any annotation is the failure this tool exists to make visible, and it is
    indistinguishable from "present but unmutated" unless the absence is counted.
    """
    by_key: dict[tuple[str, int, str], list[WatchlistEntry]] = {}
    for entry in entries:
        by_key.setdefault(entry.key, []).append(entry)
    watched_genes = {entry.gene for entry in entries}

    df, _ = read_vcf(str(vcf_path))
    genes_seen: dict[str, int] = {gene: 0 for gene in watched_genes}
    hits: list[dict] = []
    mnv_unresolved = 0

    for row in df.itertuples(index=False):
        info = str(row.INFO)
        if has_info_flag(info, "MNV_UNRESOLVED"):
            mnv_unresolved += 1
            continue
        for gene, feature_id, position, alt_aa in aa_changes_from_info(info):
            key_gene = resolve_gene(gene, feature_id, watched_genes)
            if key_gene is None:
                continue
            genes_seen[key_gene] += 1
            matched = by_key.get((key_gene, position, alt_aa))
            if not matched:
                continue
            # The observed reference residue comes from the row, not the watchlist: when the
            # two disagree the sample's own reference is the fact and the watchlist's REF is
            # the claim, and the output has to show which is which.
            hgvs_p = source_hgvs_p(info, gene, feature_id, position, alt_aa)
            parsed = parse_hgvs_p(hgvs_p)
            observed_ref = parsed[0] if parsed else ""
            for entry in matched:
                if observed_ref and observed_ref != entry.ref_aa:
                    log.warning(
                        "%s %s: reference residue at %s:%d is %s, not %s as the watchlist "
                        "says. Reported as a hit -- matching is on position and ALT.",
                        sample_id, segment, key_gene, position, observed_ref, entry.ref_aa,
                    )
                hits.append({
                    "sample_id": sample_id,
                    "segment": segment,
                    "gene": key_gene,
                    "position": position,
                    "ref_aa": observed_ref or entry.ref_aa,
                    "alt_aa": alt_aa,
                    "mutation": f"{observed_ref or entry.ref_aa}{position}{alt_aa}",
                    "chrom": str(row.CHROM),
                    "pos": str(row.POS),
                    "hgvs_p": hgvs_p,
                    "note": entry.note,
                })

    return hits, genes_seen, mnv_unresolved


def source_hgvs_p(info: str, gene: str, feature_id: str, position: int, alt_aa: str) -> str:
    """The HGVS.p string of the ANN entry a change came from, or ''.

    aa_changes_from_info de-duplicates and does not carry the source string. Re-deriving it
    here keeps that function's return shape shared with the merge rather than widened for
    one caller's reporting needs, and it only runs on rows that already matched.
    """
    for fields in ann_entries(info):
        if len(fields) < 11 or fields[3] != gene or fields[6] != feature_id:
            continue
        parsed = parse_hgvs_p(fields[10])
        if parsed and parsed[1] == position and parsed[2] == alt_aa:
            return fields[10]
    return ""


def report_coverage(
    entries: list[WatchlistEntry],
    genes_seen: dict[str, int],
    hits: list[dict],
) -> list[dict]:
    """One row per watchlist gene, and a warning for any that the VCF never annotated.

    This is the safeguard. An earlier watchlist keyed 88 HA rows on the gene name 'HA'
    while the production datasets annotate HA1 and HA2, so a quarter of the list matched
    nothing -- and nothing said so, because zero hits is what a clean sample looks like
    too. A gene the VCF has never heard of is a different statement from a gene with no
    mutations, and the two must not print the same.
    """
    per_gene_entries: dict[str, int] = {}
    for entry in entries:
        per_gene_entries[entry.gene] = per_gene_entries.get(entry.gene, 0) + 1
    per_gene_hits: dict[str, int] = {}
    for hit in hits:
        per_gene_hits[hit["gene"]] = per_gene_hits.get(hit["gene"], 0) + 1

    rows = []
    for gene in sorted(per_gene_entries):
        rows.append({
            "gene": gene,
            "in_vcf": "yes" if genes_seen.get(gene, 0) else "no",
            "watchlist_entries": per_gene_entries[gene],
            "hits": per_gene_hits.get(gene, 0),
        })
    return rows


def warn_unmatched_genes(coverage: list[dict], annotated_genes: set[str]) -> None:
    for row in coverage:
        if row["in_vcf"] == "yes":
            continue
        present = ", ".join(sorted(annotated_genes)) or "none"
        suggestion = ""
        # The HA1/HA2 case is worth naming outright: it is the one that has actually bitten,
        # and the corrected key is sitting right there in the list of genes present.
        near = sorted(g for g in annotated_genes if g.startswith(row["gene"]) and g != row["gene"])
        if near:
            suggestion = f" Did you mean {' or '.join(near)}?"
        log.warning(
            "watchlist gene %r matched no annotated gene in this VCF (%d entr%s unused). "
            "Genes present: %s.%s",
            row["gene"], row["watchlist_entries"],
            "y" if row["watchlist_entries"] == 1 else "ies", present, suggestion,
        )


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict]) -> None:
    """Always writes the header, so 'scanned, found nothing' is not 'did not run'."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def get_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # The process's provenance block reads this rather than repeating the literal, so the
    # recorded tool_version cannot drift from SCRIPT_VERSION.
    p.add_argument("--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}")
    p.add_argument("--vcf", required=True, type=Path, help="Final annotated VCF to scan")
    p.add_argument("--watchlist", required=True, type=Path, help="Watchlist CSV (gene,mutation[,note])")
    # Taken as an argument rather than parsed out of the filename: the previous version used
    # args.input.split('_')[0], which truncated every sample id at its first underscore.
    p.add_argument("--sample", required=True, help="Sample id, reported verbatim")
    p.add_argument("--segment", default="", help="Segment label, reported verbatim")
    p.add_argument("--reference-aa", type=Path,
                   help="coordinate_map.aa.tsv to validate watchlist REF residues against")
    p.add_argument("-o", "--output", required=True, type=Path, help="Output hits TSV")
    p.add_argument("--coverage-output", type=Path,
                   help="Coverage sidecar TSV (default: <output> with _hits -> _coverage)")
    return p.parse_args()


def main() -> int:
    args = get_args()

    try:
        entries = load_watchlist(args.watchlist)
    except WatchlistError as exc:
        # Fatal, unlike the previous implementation's `except Exception: return dict()`,
        # under which an unreadable watchlist reported no mutations of concern.
        log.error("%s", exc)
        return 1

    for message in report_conflicts(entries):
        log.warning("%s", message)

    if args.reference_aa:
        try:
            residues = load_reference_residues(args.reference_aa)
        except WatchlistError as exc:
            log.error("%s", exc)
            return 1
        for message in check_reference_residues(entries, residues):
            log.warning("%s", message)

    hits, genes_seen, mnv_unresolved = scan_vcf(args.vcf, entries, args.sample, args.segment)

    if mnv_unresolved:
        log.warning(
            "%d row(s) carry the MNV_UNRESOLVED flag; their HGVS.p is known to be wrong, "
            "so they were skipped rather than scanned", mnv_unresolved,
        )

    coverage = report_coverage(entries, genes_seen, hits)
    warn_unmatched_genes(coverage, {g for g, n in genes_seen.items() if n})

    write_tsv(args.output, OUTPUT_COLUMNS, hits)

    coverage_path = args.coverage_output
    if coverage_path is None:
        name = args.output.name
        coverage_path = args.output.with_name(
            name.replace("_hits.tsv", "_coverage.tsv") if "_hits.tsv" in name
            else f"{args.output.stem}_coverage.tsv"
        )
    write_tsv(
        coverage_path, COVERAGE_COLUMNS,
        [dict(sample_id=args.sample, segment=args.segment, **row) for row in coverage],
    )

    log.info(
        "%s %s: %d hit(s) from %d watchlist entr%s -> %s",
        args.sample, args.segment or "-", len(hits), len(entries),
        "y" if len(entries) == 1 else "ies", args.output,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
