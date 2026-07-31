#!/usr/bin/env python3
"""Scan an annotated VCF for the amino-acid substitutions on a watchlist.

    mutation_watcher.py --vcf sample_HA_final.vcf --watchlist nai.csv \\
                        --sample SAMPLE-ID --segment HA -o hits.tsv

A watchlist is a CSV of `gene,mutation[,note]`. See docs/watchlists.md.

Matching keys on (gene, position, ALT). The watchlist's REF residue is not part of the key:
what makes a sample resistant is ending up with Y at 275, whatever residue it started from.
REF is instead compared against what the VCF reports at that position, and a disagreement
is warned about.

Hits are never written back into the VCF. FLUMUT= is a curated annotation and belongs
there; a watchlist is a file we edit, and merging it in would make the same consensus
produce different VCFs depending on the contents of a CSV.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
from dataclasses import dataclass
from pathlib import Path

# bin/ holds entry points and is on PATH for every task; the modules they import live in
# bin/lib/. Python already puts a script's own directory on sys.path, but this insert uses
# realpath so the import still resolves when bin/ is reached through a symlink.
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from lib.aa_markers import annotated_changes, has_info_flag, parse_hgvs_p  # noqa: E402
from lib.flu_gene_names import KNOWN_GENE_KEYS, canonical_gene  # noqa: E402
from lib.tools import read_vcf  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

SCRIPT_VERSION = "3.0.0"

COLUMNS = (
    "sample_id", "segment", "gene", "position", "ref_aa", "alt_aa",
    "mutation", "chrom", "pos", "hgvs_p", "note",
)

REQUIRED_COLUMNS = ("gene", "mutation")


class WatchlistError(Exception):
    """A watchlist that cannot be scanned with. Always fatal, never swallowed."""


@dataclass(frozen=True)
class Entry:
    """One watchlist row. ``line`` is the physical line in the file, for error messages."""

    gene: str
    mutation: str
    ref_aa: str
    position: int
    alt_aa: str
    note: str
    line: int

    @property
    def key(self) -> tuple[str, int, str]:
        """The matching key. REF is deliberately absent."""
        return (self.gene, self.position, self.alt_aa)


def load_watchlist(path: Path) -> dict[tuple[str, int, str], Entry]:
    """Watchlist entries keyed on (gene, position, ALT). Raises WatchlistError.

    Anything malformed is fatal: a watchlist that silently failed to load would report no
    mutations of concern, which is exactly what a clean sample looks like.

    Lines beginning with '#' are comments; physical line numbers are tracked so an error
    points at the line to edit. Duplicate keys collapse to the first entry.
    """
    numbered = [
        (number, line)
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1)
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not numbered:
        raise WatchlistError(f"{path}: no data rows (file is empty or all comments)")

    reader = csv.DictReader([line for _, line in numbered])
    columns = [c.strip() for c in (reader.fieldnames or [])]
    missing = [c for c in REQUIRED_COLUMNS if c not in columns]
    if missing:
        raise WatchlistError(
            f"{path}: missing required column(s) {', '.join(missing)}. Expected "
            f"{','.join(REQUIRED_COLUMNS)}[,note]; found {', '.join(columns) or 'nothing'}."
        )

    # numbered[0] is the header, so reader row i came from numbered[i + 1].
    entries = [
        _entry(path, line, row)
        for (line, _), row in zip(numbered[1:], reader)
    ]
    if not entries:
        raise WatchlistError(f"{path}: no data rows (file is empty or all comments)")

    # setdefault rather than a comprehension: this must keep the first of a repeated key
    # *and* the file's order, which decides the order hits are reported in.
    watchlist: dict[tuple[str, int, str], Entry] = {}
    for entry in entries:
        watchlist.setdefault(entry.key, entry)
    return watchlist


def _entry(path: Path, line: int, row: dict) -> Entry:
    """One validated Entry from a CSV row, or WatchlistError naming the line."""
    gene = canonical_gene((row.get("gene") or "").strip())
    mutation = (row.get("mutation") or "").strip()

    if not gene or not mutation:
        raise WatchlistError(
            f"{path}:{line}: both 'gene' and 'mutation' are required "
            f"(gene={gene!r}, mutation={mutation!r})"
        )

    parsed = parse_hgvs_p(f"p.{mutation}")
    if parsed is None:
        raise WatchlistError(
            f"{path}:{line}: cannot parse mutation {mutation!r}. Expected a single "
            "amino-acid substitution as [REF][POS][ALT], 1-letter (E627K) or 3-letter "
            "(Glu627Lys). Insertions, deletions and frameshifts are not supported."
        )
    ref_aa, position, alt_aa = parsed

    if ref_aa == alt_aa:
        raise WatchlistError(
            f"{path}:{line}: {mutation!r} is synonymous (REF and ALT are both {ref_aa!r}), "
            "so it can never match a non-synonymous annotation."
        )
    if gene not in KNOWN_GENE_KEYS:
        raise WatchlistError(
            f"{path}:{line}: unknown gene {gene!r}. "
            f"Known genes: {', '.join(sorted(KNOWN_GENE_KEYS))}."
        )

    return Entry(gene, mutation, ref_aa, position, alt_aa, (row.get("note") or "").strip(), line)


def resolve_gene(gene: str, feature_id: str, watched: set[str]) -> str | None:
    """The watchlist gene key an ANN entry should be looked up under, or None.

    Feature_ID is tried before Gene_Name because the feature is what HGVS.p numbering
    belongs to: an NCBI-style HA annotation emits Gene_Name=HA with Feature_ID=HA1 for an
    HA1 change, and only the second names the protein the position is numbered in.

    There is no path from a parent gene to one of its subunits. HA1 numbering starts after
    HA's 16-residue signal peptide, so matching a bare-HA watchlist entry against an HA1
    annotation would report the wrong residue at 93% of positions.
    """
    return next(
        (key for key in map(canonical_gene, (feature_id, gene)) if key in watched),
        None,
    )


def index_vcf(
    path: Path, watched: set[str]
) -> tuple[dict[tuple[str, int, str], tuple], dict[tuple[str, int], str], int]:
    """Index a VCF's annotations for lookup, and count the rows skipped.

    Returns:
        changes       (gene, position, ALT) -> (change, chrom, pos), which drives hits.
        observed_ref  (gene, position) -> the reference residue the VCF reports there,
                      which drives the REF-disagreement warning.
        skipped       rows flagged MNV_UNRESOLVED. fix_mnv_codons.py has already determined
                      their HGVS.p to be wrong, so a hit on one would be a false positive.

    Where a position is annotated more than once the first observation wins; duplicate
    annotations of one position report the same reference residue.
    """
    df, _ = read_vcf(str(path))
    rows = [(str(r.CHROM), str(r.POS), str(r.INFO)) for r in df.itertuples(index=False)]
    scanned = [row for row in rows if not has_info_flag(row[2], "MNV_UNRESOLVED")]

    observations = [
        (gene, chrom, pos, change)
        for chrom, pos, info in scanned
        for change in annotated_changes(info)
        if (gene := resolve_gene(change.gene, change.feature_id, watched)) is not None
    ]
    changes = {
        (gene, c.position, c.alt_aa): (c, chrom, pos)
        for gene, chrom, pos, c in reversed(observations)
    }
    observed_ref = {
        (gene, c.position): c.ref_aa for gene, _, _, c in reversed(observations)
    }
    return changes, observed_ref, len(rows) - len(scanned)


def find_hits(
    watchlist: dict[tuple[str, int, str], Entry],
    changes: dict[tuple[str, int, str], tuple],
    sample_id: str,
    segment: str,
) -> list[dict]:
    """One output row per watchlist entry the VCF annotates, in watchlist order.

    The reference residue reported is the one observed, not the one the watchlist claims:
    where they differ, the sample's own reference is the fact and the watchlist is a claim.
    """
    return [
        _hit_row(entry, changes[key], sample_id, segment)
        for key, entry in watchlist.items()
        if key in changes
    ]


def _hit_row(entry: Entry, observation: tuple, sample_id: str, segment: str) -> dict:
    change, chrom, pos = observation
    return {
        "sample_id": sample_id,
        "segment": segment,
        "gene": entry.gene,
        "position": entry.position,
        "ref_aa": change.ref_aa,
        "alt_aa": entry.alt_aa,
        "mutation": f"{change.ref_aa}{entry.position}{entry.alt_aa}",
        "chrom": chrom,
        "pos": pos,
        "hgvs_p": change.hgvs_p,
        "note": entry.note,
    }


def warn_ref_mismatches(
    watchlist: dict[tuple[str, int, str], Entry],
    observed_ref: dict[tuple[str, int], str],
    path: Path,
) -> None:
    """Warn where an entry's REF disagrees with the residue the VCF reports there.

    Only entries at a position the VCF annotates can be checked, because a VCF carries
    variants rather than every residue.

    A mismatch has two causes the data cannot separate: the entry uses a numbering
    convention this dataset does not, or this reference genuinely differs there. A
    systematic numbering error lights up many entries at once; a real difference shows up
    as one or two. Either way the entry is still scanned, since matching is on position
    and ALT.
    """
    mismatched = [
        (entry, seen)
        for entry in watchlist.values()
        if (seen := observed_ref.get((entry.gene, entry.position))) and seen != entry.ref_aa
    ]
    for entry, seen in mismatched:
        log.warning(
            "%s %s (%s:%d): VCF reports reference %s at %s:%d, not %s",
            entry.gene, entry.mutation, path.name, entry.line,
            seen, entry.gene, entry.position, entry.ref_aa,
        )


def write_tsv(path: Path, rows: list[dict]) -> None:
    """Always writes the header, so 'scanned, found nothing' is not 'did not run'."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def get_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # The process's provenance block reads this rather than repeating the literal, so the
    # recorded tool_version cannot drift from SCRIPT_VERSION.
    p.add_argument("--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}")
    p.add_argument("--vcf", required=True, type=Path, help="Annotated VCF to scan")
    p.add_argument("--watchlist", required=True, type=Path,
                   help="Watchlist CSV (gene,mutation[,note])")
    p.add_argument("--sample", required=True, help="Sample id, reported verbatim")
    p.add_argument("--segment", default="", help="Segment label, reported verbatim")
    p.add_argument("-o", "--output", required=True, type=Path, help="Output hits TSV")
    return p.parse_args()


def main() -> int:
    args = get_args()

    try:
        watchlist = load_watchlist(args.watchlist)
    except WatchlistError as exc:
        log.error("%s", exc)
        return 1

    watched = {entry.gene for entry in watchlist.values()}
    changes, observed_ref, skipped = index_vcf(args.vcf, watched)

    if skipped:
        log.warning(
            "%d row(s) carry the MNV_UNRESOLVED flag; their HGVS.p is known to be wrong, "
            "so they were skipped rather than scanned", skipped,
        )

    warn_ref_mismatches(watchlist, observed_ref, args.watchlist)
    hits = find_hits(watchlist, changes, args.sample, args.segment)
    write_tsv(args.output, hits)

    log.info(
        "%s %s: %d hit(s) from %d watchlist entr%s -> %s",
        args.sample, args.segment or "-", len(hits), len(watchlist),
        "y" if len(watchlist) == 1 else "ies", args.output,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
