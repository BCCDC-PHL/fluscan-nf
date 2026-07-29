#!/usr/bin/env python3
"""Load and validate a custom mutation watchlist, and report what cannot be trusted.

A watchlist is a CSV of amino-acid substitutions to scan every sample for. Unlike the
FluMutDB markers that `merge_flumut_mutations.py` writes into the VCF as FLUMUT=, this
list is set on our own terms, so nothing curates it but us -- which makes silent failure
the thing to design against. The previous implementation caught every load error with
`except Exception: return dict()`, so a malformed file scanned cleanly and reported no
mutations of concern.

    gene,mutation,note
    PB2,E627K,mammalian adaptation
    HA1,T195I,
    NA,H275Y,oseltamivir resistance

Lines beginning with '#' are comments and are ignored, so a list can record where it came
from without a second file.

`gene` is canonicalised through lib/flu_gene_names.py -- the same function the coordinate
map builder and the FluMut merge use -- so NS-1/NS1 and NS-2/NEP are interchangeable.

The column is `gene`, not `segment`: it is matched against a VCF ANN entry's gene name,
and HA1 is not a segment. A file still carrying a `segment` column is rejected with the
rename spelled out rather than silently read, because the two are not the same question
and reading one as the other is how 88 HA rows of an earlier list matched nothing at all.

REF/POS/ALT and matching
------------------------
`mutation` is [REF][POS][ALT] in either the 1-letter (E627K) or 3-letter (Glu627Lys)
spelling, parsed by the same `lib/aa_markers.parse_hgvs_p` the merge uses.

Matching is on (gene, position, alt) only -- the REF letter is *not* part of the key. What
makes a sample resistant is ending up with Y at 275, not which residue it started from; if
a dataset's reference carries N there and a sample reads p.N275Y, that sample is a hit.
This mirrors the merge, whose MARKER_QUERY omits ref_seq for the same stated reason.

REF's job here is to sanity-check the list. With --reference-aa it is compared against the
reference's actual residue at that gene and position, which is what catches a wrong
numbering convention: R292K checked against N1 coordinates finds the reference residue at
292 is not R, and the entry is named at load time instead of silently mis-scanning every
sample. A mismatch warns rather than fails, because it has two possible causes that the
data cannot distinguish -- a wrong entry, or a reference that genuinely differs there. In
practice a systematic numbering error lights up many entries at once and is unmistakable,
while a real reference difference shows up as one or two.

Exit codes: 0 valid, 1 invalid. With --strict, warnings also exit 1, for curation and CI.
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
from lib.aa_markers import parse_hgvs_p  # noqa: E402
from lib.flu_gene_names import KNOWN_GENE_KEYS, canonical_gene  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

SCRIPT_VERSION = "1.0.0"

REQUIRED_COLUMNS = ("gene", "mutation")
OPTIONAL_COLUMNS = ("note",)


class WatchlistError(Exception):
    """A watchlist that cannot be scanned with. Always fatal -- never swallowed."""


@dataclass(frozen=True)
class WatchlistEntry:
    """One validated watchlist row.

    ``gene`` is canonical; ``raw_gene`` is what the file said, so messages can quote it.
    ``line`` is the physical line number in the source file, for the same reason.
    """

    gene: str
    raw_gene: str
    mutation: str
    ref_aa: str
    position: int
    alt_aa: str
    note: str
    line: int

    @property
    def key(self) -> tuple[str, int, str]:
        """The matching key: gene, position, alt. REF is deliberately absent."""
        return (self.gene, self.position, self.alt_aa)


def _decommented(handle) -> tuple[list[str], list[int]]:
    """File lines with '#' comments and blank lines removed, line numbers preserved.

    csv.DictReader has no comment support, so comments are stripped first. Physical line
    numbers are tracked separately and attached to each entry, so an error message points
    at the line the user has to edit rather than at a post-filter row index.
    """
    kept: list[str] = []
    numbers: list[int] = []
    for number, line in enumerate(handle, start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        kept.append(line)
        numbers.append(number)
    return kept, numbers


def _non_ascii(text: str) -> str | None:
    """The first non-ASCII character in ``text``, or None.

    Worth its own check because the failure it catches is invisible on screen: an en dash
    (U+2013) in '245-248del' is indistinguishable from a hyphen in most editors, and the
    resulting parse failure would otherwise be reported as a malformed mutation with no
    hint as to why.
    """
    for char in text:
        if ord(char) > 127:
            return char
    return None


def load_watchlist(path: Path) -> list[WatchlistEntry]:
    """Parse and validate a watchlist file. Raises WatchlistError on anything fatal.

    Fatal: a missing required column, a `segment` column, a non-ASCII or unparseable
    mutation, a synonymous mutation, an unknown gene. Exact duplicate rows are dropped
    with a warning; everything else that is merely suspicious is left to
    ``report_conflicts`` and ``check_reference_residues``.
    """
    with path.open(newline="", encoding="utf-8") as handle:
        lines, numbers = _decommented(handle)

    if not lines:
        raise WatchlistError(f"{path}: no data rows (file is empty or all comments)")

    reader = csv.DictReader(lines)
    columns = [c.strip() for c in (reader.fieldnames or [])]

    if "segment" in columns and "gene" not in columns:
        raise WatchlistError(
            f"{path}: the column is named 'gene', not 'segment'. Rename the header column "
            "and check its values: they are matched against a VCF annotation's gene name, "
            "so HA1 and HA2 are the keys for HA -- not 'HA', whose numbering differs by a "
            "16-residue signal peptide."
        )

    missing = [c for c in REQUIRED_COLUMNS if c not in columns]
    if missing:
        raise WatchlistError(
            f"{path}: missing required column(s) {', '.join(missing)}. "
            f"Expected {', '.join(REQUIRED_COLUMNS)}[,{','.join(OPTIONAL_COLUMNS)}]; "
            f"found {', '.join(columns) or 'nothing'}."
        )

    entries: list[WatchlistEntry] = []
    seen: dict[tuple[str, str], int] = {}
    # numbers[0] is the header's physical line; row i of the reader is numbers[i + 1].
    for index, row in enumerate(reader):
        line = numbers[index + 1] if index + 1 < len(numbers) else -1
        raw_gene = (row.get("gene") or "").strip()
        mutation = (row.get("mutation") or "").strip()
        note = (row.get("note") or "").strip()

        if not raw_gene or not mutation:
            raise WatchlistError(
                f"{path}:{line}: both 'gene' and 'mutation' are required "
                f"(gene={raw_gene!r}, mutation={mutation!r})"
            )

        for field, value in (("gene", raw_gene), ("mutation", mutation)):
            bad = _non_ascii(value)
            if bad is not None:
                raise WatchlistError(
                    f"{path}:{line}: non-ASCII character {bad!r} (U+{ord(bad):04X}) in "
                    f"{field} {value!r}. Editors render an en dash and a hyphen almost "
                    "identically; retype the value in plain ASCII."
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
                f"{path}:{line}: {mutation!r} is synonymous (REF and ALT are both "
                f"{ref_aa!r}), so it can never match a non-synonymous annotation."
            )

        gene = canonical_gene(raw_gene)
        if gene not in KNOWN_GENE_KEYS:
            raise WatchlistError(
                f"{path}:{line}: unknown gene {raw_gene!r}"
                + (f" (canonicalised to {gene!r})" if gene != raw_gene else "")
                + f". Known genes: {', '.join(sorted(KNOWN_GENE_KEYS))}."
            )

        dedup_key = (gene, f"{ref_aa}{position}{alt_aa}")
        if dedup_key in seen:
            log.warning(
                "%s:%d: duplicate of line %d (%s %s) -- keeping the first",
                path, line, seen[dedup_key], raw_gene, mutation,
            )
            continue
        seen[dedup_key] = line

        entries.append(
            WatchlistEntry(
                gene=gene,
                raw_gene=raw_gene,
                mutation=mutation,
                ref_aa=ref_aa,
                position=position,
                alt_aa=alt_aa,
                note=note,
                line=line,
            )
        )

    if not entries:
        raise WatchlistError(f"{path}: no usable entries after validation")

    return entries


def report_conflicts(entries: list[WatchlistEntry]) -> list[str]:
    """Internal contradictions that usually mean two numbering conventions in one file.

    Two shapes, both seen in the real NA list:

    1. One (gene, position) given two different REF residues. A position has one reference
       residue; two means two coordinate systems.
    2. The same substitution one position apart -- I222L and I223L. That is one marker
       written under N2 and N1 numbering, not two markers.

    Reported, never auto-resolved: which convention is wanted depends on what this lab
    reports in and where each entry came from, and a tool that guessed would silently
    delete a marker somebody deliberately added.
    """
    messages: list[str] = []

    by_position: dict[tuple[str, int], list[WatchlistEntry]] = {}
    for entry in entries:
        by_position.setdefault((entry.gene, entry.position), []).append(entry)
    for (gene, position), group in sorted(by_position.items()):
        refs = sorted({e.ref_aa for e in group})
        if len(refs) > 1:
            messages.append(
                f"{gene}:{position} is given more than one reference residue "
                f"({', '.join(refs)}) by lines "
                f"{', '.join(str(e.line) for e in sorted(group, key=lambda e: e.line))}. "
                "A position has one reference residue; two means two numbering conventions."
            )

    by_substitution: dict[tuple[str, str, str], list[WatchlistEntry]] = {}
    for entry in entries:
        by_substitution.setdefault((entry.gene, entry.ref_aa, entry.alt_aa), []).append(entry)
    for (gene, ref_aa, alt_aa), group in sorted(by_substitution.items()):
        ordered = sorted(group, key=lambda e: e.position)
        for first, second in zip(ordered, ordered[1:]):
            if second.position - first.position == 1:
                messages.append(
                    f"{gene}: {first.mutation} (line {first.line}) and {second.mutation} "
                    f"(line {second.line}) are the same substitution one position apart, "
                    "which is what one marker looks like under two numbering conventions. "
                    "Keep whichever matches the reference in use."
                )

    return messages


def load_reference_residues(path: Path) -> dict[tuple[str, int], str]:
    """{(gene, target_aa_pos): residue} from a coordinate_map.aa.tsv.

    The *target* columns are the ones to read: watchlist positions are in the coordinates
    the VCF reports, which is the coordinate map's target reference.
    """
    residues: dict[tuple[str, int], str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene", "target_aa_pos", "target_residue"}
        if not required.issubset(reader.fieldnames or []):
            raise WatchlistError(
                f"{path}: not a coordinate_map.aa.tsv "
                f"(missing {', '.join(sorted(required - set(reader.fieldnames or [])))})"
            )
        for row in reader:
            target_pos = (row.get("target_aa_pos") or "").strip()
            residue = (row.get("target_residue") or "").strip()
            if not target_pos or not residue:
                continue
            residues[(canonical_gene(row["gene"].strip()), int(target_pos))] = residue
    return residues


def check_reference_residues(
    entries: list[WatchlistEntry], residues: dict[tuple[str, int], str]
) -> list[str]:
    """Compare each entry's REF against the reference's actual residue at that position.

    Warnings, not errors: a mismatch means either a wrong entry or a reference that
    genuinely differs there, and nothing in the data separates the two. Positions the map
    does not cover are reported separately -- a position past the end of a protein is a
    stronger signal than a residue disagreement.
    """
    messages: list[str] = []
    for entry in entries:
        actual = residues.get((entry.gene, entry.position))
        if actual is None:
            messages.append(
                f"{entry.gene}:{entry.mutation} (line {entry.line}): position "
                f"{entry.position} is not present in the reference map for {entry.gene}."
            )
        elif actual != entry.ref_aa:
            messages.append(
                f"{entry.gene}:{entry.mutation} (line {entry.line}): reference residue at "
                f"{entry.position} is {actual!r}, not {entry.ref_aa!r}. Either the entry "
                "uses a different numbering convention, or this reference differs there. "
                "Matching is on position and ALT, so the entry is still used."
            )
    return messages


def main() -> int:
    p = argparse.ArgumentParser(
        description="Validate a mutation watchlist CSV.",
        epilog=(
            "Exits 0 when the list is usable, 1 when it is not. --strict also exits 1 on "
            "warnings, for curation and CI."
        ),
    )
    p.add_argument("watchlist", type=Path, help="watchlist CSV (gene,mutation[,note])")
    p.add_argument(
        "--reference-aa",
        type=Path,
        help=(
            "coordinate_map.aa.tsv to check each entry's REF residue against. Without it, "
            "numbering errors are only visible at hit time."
        ),
    )
    p.add_argument(
        "--strict",
        action="store_true",
        help="treat warnings (duplicates, conflicts, REF mismatches) as failure",
    )
    # Shells out from the provenance block, so the recorded tool_version cannot drift.
    p.add_argument("--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}")
    args = p.parse_args()

    try:
        entries = load_watchlist(args.watchlist)
    except WatchlistError as exc:
        log.error("%s", exc)
        return 1

    warnings = report_conflicts(entries)

    if args.reference_aa:
        try:
            residues = load_reference_residues(args.reference_aa)
        except WatchlistError as exc:
            log.error("%s", exc)
            return 1
        warnings += check_reference_residues(entries, residues)

    for message in warnings:
        log.warning("%s", message)

    log.info(
        "%s: %d entries over %d gene(s): %s",
        args.watchlist,
        len(entries),
        len({e.gene for e in entries}),
        ", ".join(sorted({e.gene for e in entries})),
    )

    if warnings:
        log.info("%d warning(s). See above.", len(warnings))
        if args.strict:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
