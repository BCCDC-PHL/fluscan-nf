#!/usr/bin/env python3
"""Turn the WHO GISRS neuraminidase marker table into source watchlists.

    parse_nai_table.py --pdf human-nai-marker-table_for-publication_final_20240918.pdf \\
                       --out-dir assets/watchlists/source

Needs pdfplumber, which is not a pipeline dependency -- run it with an interpreter that
has it. This is a one-off producer for one document, kept so the derivation from the PDF
is reproducible; the general path into the watchlist system is a hand-written source CSV,
which `resolve_watchlist.py` reads directly. Nothing here is imported by anything.

The table's rows carry each position twice, in the subtype's own numbering and in N2
numbering. Transcribing it by hand mixed the two, which is what this replaces: the
native column is kept as the position, and the N2 column is carried alongside as an
independent cross-check.

One file is written per subtype section that has a local reference to resolve against.
Every row that is not written is counted and reported, so the totals reconcile against
the document.
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import pdfplumber

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

SCRIPT_VERSION = "1.0.0"

SOURCE_CITATION = (
    "human-nai-marker-table_for-publication_final_20240918.pdf "
    "(WHO GISRS, 18 September 2024)"
)

# The document is a neuraminidase table throughout, so every row is an NA entry.
GENE = "NA"

# Tables run to page 8; the rest is the reference list.
TABLE_PAGES = 8

# A row's first cell is one physical column holding "Type/subtype" and "Comments", so a
# subtype label and the first data row of its block share a cell.
SUBTYPE_RE = re.compile(r"^(?:A\(H\d+N\d+\)(?:pdm09|v)?|B)$")

# [REF][POS][ALT], ALT occasionally given as a choice (D151E/N), with trailing footnote
# letters. ALT is always upper case, so a trailing lower-case run is unambiguously a
# footnote and never part of the substitution.
SUBSTITUTION_RE = re.compile(r"^([A-Z])(\d+)([A-Z](?:/[A-Z])*)([a-z]*)$")

FIELDS = ("mutation", "n2", "oseltamivir", "zanamivir", "peramivir", "laninamivir",
          "source_of_viruses", "references")

COLUMNS = ("gene", "mutation", "subtype", "n2_position", "footnotes", "combination_of",
           "oseltamivir", "zanamivir", "peramivir", "laninamivir", "source_of_viruses",
           "references", "comment")

# Sections with a local reference bundle to resolve against, and the numbering each is
# written in. A(H3N2) is the definition of N2 numbering (the table's footnote b), so its
# native column and its N2 column are the same thing.
#
# The A(H5N1) section mixes two native frames: most rows follow A/turkey/Turkey/1/2005,
# which carries a 20-residue stalk deletion, while footnote-f rows follow a clade 2.3.4.4b
# strain that does not. Declaring both candidate offsets lets the resolver pick per entry
# by residue agreement instead of this script guessing.
SECTIONS = {
    "A(H1N1)pdm09": ("nai_h1n1pdm", "n1-pdm09", None),
    "A(H5N1)": ("nai_h5n1", "n1-h5-who", "0, 20"),
    "A(H3N2)": ("nai_h3n2", "n2", None),
}


@dataclass
class Tally:
    written: int = 0
    other_section: int = 0
    not_a_substitution: list[str] = field(default_factory=list)
    mismatched_n2: list[str] = field(default_factory=list)
    non_data_rows: int = 0


def cell(value: str | None) -> str:
    return (value or "").replace("\n", " ").strip()


def read_rows(pdf_path: Path) -> tuple[list[tuple[str, str, dict]], int]:
    """(subtype, comment, fields) per data row, plus the count of non-data rows.

    Cells are compacted rather than indexed. The number of raw columns varies by page
    (9 where cells are not split, 27 where merged cells spread each column over three)
    and the spread is not uniform -- one row on page 8 puts its zanamivir value a cell
    late -- so an index-based read silently mis-assigns. Every data row in the document
    compacts to exactly the eight fields below; anything else is not a data row.
    """
    rows: list[tuple[str, str, dict]] = []
    non_data = 0
    subtype = ""
    with pdfplumber.open(pdf_path) as pdf:
        for page in pdf.pages[:TABLE_PAGES]:
            table = page.extract_table()
            if table is None:
                continue
            for raw in table:
                head = cell(raw[0])
                values = [value for value in map(cell, raw[1:]) if value]
                if SUBTYPE_RE.match(head):
                    subtype, head = head, ""
                if len(values) != len(FIELDS):
                    non_data += 1
                    continue
                rows.append((subtype, head, dict(zip(FIELDS, values))))
    return rows, non_data


def expand(row: dict) -> list[dict] | str:
    """One table row as its constituent substitutions, or why it yielded none.

    A combination (`H275Y+I436N`) is split into its constituents, each paired with the
    matching part of the N2 column and tagged with the combination it came from. That
    means a constituent can flag on its own where the document claimed an effect only
    for the full set; `combination_of` is what records it.
    """
    parts = [part.strip() for part in row["mutation"].split("+")]
    n2_parts = [part.strip() for part in row["n2"].split("+")]
    if len(parts) != len(n2_parts):
        return f"{row['mutation']}: {len(parts)} substitution(s) but {len(n2_parts)} N2 position(s)"

    expanded = []
    for part, n2 in zip(parts, n2_parts):
        matched = SUBSTITUTION_RE.match(part)
        if not matched:
            # Deletions and ranges (`Del 247-250`) are not single substitutions and the
            # watcher cannot match them at all.
            return f"{row['mutation']}: {part!r} is not a single substitution"
        ref, position, alts, footnotes = matched.groups()
        expanded.extend(
            {
                "mutation": f"{ref}{position}{alt}",
                "n2_position": n2,
                "footnotes": footnotes,
                "combination_of": row["mutation"] if len(parts) > 1 else "",
            }
            for alt in alts.split("/")
        )
    return expanded


def write_section(path: Path, numbering: str, offsets: str | None, rows: list[dict]) -> None:
    header = [f"numbering: {numbering}"]
    if offsets:
        header.append(f"offsets: {offsets}")
    header.append(f"source: {SOURCE_CITATION}")
    header.append(f"generated: by scripts/parse_nai_table.py {SCRIPT_VERSION}")

    with path.open("w", newline="", encoding="utf-8") as handle:
        handle.writelines(f"# {line}\n" for line in header)
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def get_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}")
    p.add_argument("--pdf", required=True, type=Path, help="The WHO marker table PDF")
    p.add_argument("--out-dir", required=True, type=Path, help="Where to write source CSVs")
    return p.parse_args()


def main() -> int:
    args = get_args()
    rows, non_data = read_rows(args.pdf)
    tally = Tally(non_data_rows=non_data)
    collected: dict[str, list[dict]] = {name: [] for name in SECTIONS}

    for subtype, comment, row in rows:
        if subtype not in SECTIONS:
            tally.other_section += 1
            continue
        expanded = expand(row)
        if isinstance(expanded, str):
            tally.not_a_substitution.append(f"{subtype} {expanded}")
            continue
        collected[subtype].extend(
            {
                "gene": GENE,
                "subtype": subtype,
                "oseltamivir": row["oseltamivir"],
                "zanamivir": row["zanamivir"],
                "peramivir": row["peramivir"],
                "laninamivir": row["laninamivir"],
                "source_of_viruses": row["source_of_viruses"],
                "references": row["references"],
                "comment": comment,
                **entry,
            }
            for entry in expanded
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for subtype, (stem, numbering, offsets) in SECTIONS.items():
        entries = collected[subtype]
        tally.written += len(entries)
        path = args.out_dir / f"{stem}_source.csv"
        write_section(path, numbering, offsets, entries)
        log.info("%s: %d entr(ies) in %s numbering -> %s", subtype, len(entries), numbering, path)

    log.info("%d data row(s) read, %d non-data row(s) skipped", len(rows), tally.non_data_rows)
    log.info("%d row(s) in sections with no local reference", tally.other_section)
    for reason in tally.not_a_substitution:
        log.warning("not written: %s", reason)
    for reason in tally.mismatched_n2:
        log.warning("N2 mismatch: %s", reason)
    log.info("%d entr(ies) written from %d row(s)", tally.written, len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
