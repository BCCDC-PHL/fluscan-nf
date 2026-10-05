#!/usr/bin/env python3
"""Split a multi-segment influenza FASTA into per-segment files.

split_fasta.py split in.fasta  ->  in_PB2.s.fasta, in_HA.s.fasta, …
"""

import argparse
import re
from pathlib import Path
from Bio import SeqIO

# The 8 canonical influenza segment names (order matches standard numbering).
SEGMENTS = ("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")


def parse_segment(header: str) -> str | None:
    """Return the segment named by the header's last `_`/`|` field (e.g. "HA"), or None."""
    for seg in SEGMENTS:
        if re.search(rf"[_\|]{seg}$", header, re.IGNORECASE):
            return seg.upper()
    return None


def split(input_fasta: Path) -> None:
    """Read *input_fasta* and write one ``{stem}_{segment}.s.fasta`` per record."""
    stem = input_fasta.stem
    out_dir = input_fasta.parent

    records = list(SeqIO.parse(input_fasta, "fasta"))
    written: dict[str, int] = {}

    for rec in records:
        seg = parse_segment(rec.description)
        if seg is None:
            print(f"WARNING: could not determine segment for '{rec.id}' — skipping")
            continue
        out_path = out_dir / f"{stem}_{seg}.s.fasta"
        with open(out_path, "a") as fh:  # append so multiple records per segment work
            SeqIO.write(rec, fh, "fasta")
        written[seg] = written.get(seg, 0) + 1

    for seg, n in written.items():
        print(f"  {seg}: {n} record(s) → {stem}_{seg}.s.fasta")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="mode", required=True)

    sp = sub.add_parser("split", help="Split one FASTA into per-segment files")
    sp.add_argument("input", type=Path, help="Input multi-record FASTA")

    split(p.parse_args().input)


if __name__ == "__main__":
    main()
