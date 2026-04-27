#!/usr/bin/env python3
"""Split a multi-segment influenza FASTA into per-segment files, or combine them back.

Split:  fasta_splitter.py split  in.fasta          → in_PB2.s.fasta, in_HA.s.fasta, …
Combine: fasta_splitter.py combine in_*.s.fasta    → in.fasta
"""

import argparse
import re
from pathlib import Path
from Bio import SeqIO

# The 8 canonical influenza segment names (order matches standard numbering).
SEGMENTS = ("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")


def parse_segment(header: str) -> str | None:
    """Extract the segment name from a FASTA header.

    TODO: Customise this regex/logic to match your header format.
    Currently looks for one of the 8 segment names as a word in the header.
    Returns the segment string (e.g. "HA") or None if not found.
    """
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


def combine(input_fastas: list[Path], out_dir: Path) -> None:
    """Merge ``*_{segment}.s.fasta`` files back into a single FASTA.

    The output filename is derived by stripping the ``_{segment}.s`` suffix
    from the first file's name.
    """
    # Derive output name from the first file: drop _{SEGMENT}.s
    first_stem = input_fastas[0].stem                          # e.g. "sample_HA.s"
    base = re.sub(r"_(" + "|".join(SEGMENTS) + r")\.s$", "", first_stem, flags=re.IGNORECASE)
    out_path = out_dir / f"{base}.fasta"

    records = [rec for f in input_fastas for rec in SeqIO.parse(f, "fasta")]
    SeqIO.write(records, out_path, "fasta")
    print(f"  {len(records)} record(s) → {out_path}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="mode", required=True)

    sp = sub.add_parser("split", help="Split one FASTA into per-segment files")
    sp.add_argument("input", type=Path, help="Input multi-record FASTA")

    cp = sub.add_parser("combine", help="Combine *.s.fasta files into one FASTA")
    cp.add_argument("inputs", type=Path, nargs="+", help="Input *.s.fasta files")
    cp.add_argument("-o", "--outdir", type=Path, default=None,
                    help="Output directory (default: same as first input file)")

    args = p.parse_args()

    if args.mode == "split":
        split(args.input)
    else:
        out_dir = args.outdir or args.inputs[0].parent
        combine(args.inputs, out_dir)


if __name__ == "__main__":
    main()

