#!/usr/bin/env python3
"""
prepare_consensus.py

Validate and prepare an influenza consensus FASTA for the consensus workflows.

The sample name is supplied by the caller as --sample-id and is authoritative:
input headers contribute the segment label and nothing else. Every output header
becomes {sample_id}_{SEGMENT}, and the caller names both output files.

Steps performed:
  a) Read the segment label from each header (see below) and map it to a record.
  b) Uppercase sequences, optionally mask partially degenerate IUPAC nucleotides
     (R,Y,S,W,K,M,B,D,H,V) with N, drop alignment gaps (-, .), and coerce any
     other non-ACGTN character to N, so the emitted FASTA is always clean.
  c) Emit a QC row for every expected segment (present or MISSING) with a
     PASS/WARN/FAIL/MISSING status driven by segment length ranges and an
     ambiguity threshold.

Segment labels are read by splitting the header on | and _ and taking the
rightmost token that names a segment, which covers every header form in use:
  >WIN-AH-2024-FAV-0415-1-OS_segment1_PB2      (underscore, wild west)
  >WIN-AH-2024-FAV-0415-1-OS_PB2               (underscore, already canonical)
  >EPI_ISL_20399898|A_/_H5N1|NA|6|EPI5261179|2.3.4.4b_NA   (GISAID-style pipes)

One file must hold one sample, and two things are fatal because they mean it does
not: two records carrying the same segment label, and headers whose leading field
names more than one sample.
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from collections import Counter
from dataclasses import asdict, dataclass, fields
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --------------------------------------------------------------------------- #
# Configuration / constants
# --------------------------------------------------------------------------- #

EXPECTED_SEGMENTS = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS']

# Every label accepted in a header, mapped to its canonical name (MP is legacy M).
SEGMENT_LOOKUP = {segment: segment for segment in EXPECTED_SEGMENTS} | {'MP': 'M'}

# Header fields are separated by pipes, tokens within them by underscores.
TOKEN_SPLIT = re.compile(r'[|_]')

# Partially degenerate IUPAC nucleotide codes (excluding N, which is already "any").
DEGENERATE_BASES = set('RYSWKMBDHV')

# Unambiguous bases plus the fully-ambiguous N.
VALID_BASES = set('ACGTN')

# Alignment gap / pad characters that must never reach downstream tools.
GAP_CHARS = set('-.')

# Everything a sequence may contain before cleaning; the rest is coerced to N.
KNOWN_CHARS = VALID_BASES | DEGENERATE_BASES | GAP_CHARS

# Acceptable sequence length ranges (min_len, max_len) per segment.
SEGMENT_LENGTH_RANGES = {
    'PB2': (1700, 2400),
    'PB1': (1650, 2400),
    'PA':  (2100, 2280),
    'HA':  (1600, 1800),
    'NP':  (1450, 1600),
    'NA':  (1300, 1550),
    'M':   (900,  1050),
    'NS':  (800,   950),
}

# Segment QC status values.
STATUS_PASS = 'PASS'
STATUS_WARN = 'WARN'
STATUS_FAIL = 'FAIL'
STATUS_MISSING = 'MISSING'

# Exit codes.
EXIT_OK = 0
EXIT_INPUT_ERROR = 2

log = logging.getLogger('prepare_consensus')


@dataclass
class SegmentQC:
    """One row of the QC report; the defaults describe an absent segment."""
    sample_id: str
    segment: str
    seq_id: str = ''
    status: str = STATUS_MISSING
    length: int = 0
    length_ok: bool = False
    n_prop_before: float = 0.0
    n_replaced: int = 0
    n_gaps_removed: int = 0
    n_invalid: int = 0
    n_prop_after: float = 0.0
    notes: str = ''


# --------------------------------------------------------------------------- #
# Parsing and cleaning
# --------------------------------------------------------------------------- #

def segment_of(header: str) -> str | None:
    """Return the canonical segment label for a header, or None if it has none."""
    for token in reversed(TOKEN_SPLIT.split(header)):
        segment = SEGMENT_LOOKUP.get(token.upper())
        if segment:
            return segment
    return None


def sample_stem(header: str) -> str:
    """Return the leading header field, which names the sample in every accepted form.

    Only ever compared across the records of one file, never used to name anything.
    """
    return header.split('|')[0] if '|' in header else header.split('_')[0]


def n_proportion(sequence: str) -> float:
    """Return the proportion of N bases in the sequence (0.0 for empty)."""
    return sequence.count('N') / len(sequence) if sequence else 0.0


def process_segment(sample_id: str, segment: str, raw: str,
                    *, mask_degenerate: bool, max_n_frac: float) -> tuple[SeqRecord, SegmentQC]:
    """Clean one segment, returning its canonical record and its QC row."""
    upper = raw.upper()
    counts = Counter(upper)
    keep = VALID_BASES if mask_degenerate else VALID_BASES | DEGENERATE_BASES
    sequence = ''.join(b if b in keep else 'N' for b in upper if b not in GAP_CHARS)

    n_gaps_removed = sum(counts[char] for char in GAP_CHARS)
    n_replaced = sum(counts[base] for base in DEGENERATE_BASES) if mask_degenerate else 0
    n_invalid = len(upper) - sum(counts[char] for char in KNOWN_CHARS)

    length = len(sequence)
    n_prop_after = n_proportion(sequence)
    low, high = SEGMENT_LENGTH_RANGES[segment]
    length_ok = low <= length <= high

    notes: list[str] = []
    if length == 0 or n_prop_after >= 1.0:
        status = STATUS_FAIL
        notes.append('empty or all-N sequence')
    elif not length_ok:
        status = STATUS_FAIL
        notes.append(f'length {length} outside expected range [{low}, {high}]')
    else:
        status = STATUS_PASS
        if n_prop_after > max_n_frac:
            status = STATUS_WARN
            notes.append(f'N fraction {n_prop_after:.3f} exceeds {max_n_frac:g}')
        if n_invalid:
            status = STATUS_WARN
            notes.append(f'{n_invalid} invalid base(s) coerced to N')
        if n_gaps_removed:
            status = STATUS_WARN
            notes.append(f'{n_gaps_removed} gap char(s) stripped')

    if status != STATUS_PASS:
        log.warning("[%s] Segment %s: %s", sample_id, segment, '; '.join(notes))

    seq_id = f'{sample_id}_{segment}'
    return (
        SeqRecord(Seq(sequence), id=seq_id, description=''),
        SegmentQC(
            sample_id=sample_id,
            segment=segment,
            seq_id=seq_id,
            status=status,
            length=length,
            length_ok=length_ok,
            n_prop_before=round(n_proportion(upper), 6),
            n_replaced=n_replaced,
            n_gaps_removed=n_gaps_removed,
            n_invalid=n_invalid,
            n_prop_after=round(n_prop_after, 6),
            notes='; '.join(notes),
        ),
    )


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Validate and clean an influenza consensus FASTA for the consensus workflows.'
    )
    parser.add_argument('--input', required=True, type=Path,
                        help='Input consensus FASTA file')
    parser.add_argument('--sample-id', required=True,
                        help='Authoritative sample identifier; every output header becomes '
                             '{sample_id}_{SEGMENT}')
    parser.add_argument('--output-fasta', required=True, type=Path,
                        help='Path to write the cleaned FASTA')
    parser.add_argument('--output-qc', required=True, type=Path,
                        help='Path to write the per-segment QC report TSV')
    parser.add_argument('--mask-degenerate', action=argparse.BooleanOptionalAction, default=True,
                        help='Mask partially degenerate IUPAC bases (R,Y,S,W,K,M,B,D,H,V) with N '
                             '(default: enabled; use --no-mask-degenerate to keep them)')
    parser.add_argument('--max-n-frac', type=float, default=1.0,
                        help='WARN when the post-clean N fraction of a segment exceeds this value '
                             '(default: 1.0, i.e. never warns)')
    parser.add_argument('--log-level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help='Logging verbosity (default: INFO)')
    return parser.parse_args(argv)


def run(args: argparse.Namespace) -> int:
    sample_id = args.sample_id

    by_segment: dict[str, SeqRecord] = {}
    stems: set[str] = set()
    for record in SeqIO.parse(str(args.input), 'fasta'):
        segment = segment_of(record.id)
        if segment is None:
            log.warning("[%s] Could not determine segment for '%s'; skipping.",
                        sample_id, record.id)
            continue
        if segment in by_segment:
            log.error("[%s] Duplicate record for segment %s ('%s' and '%s'); one file must "
                      "hold one sample.", sample_id, segment, by_segment[segment].id, record.id)
            return EXIT_INPUT_ERROR
        by_segment[segment] = record
        stems.add(sample_stem(record.id))

    if not by_segment:
        log.error("[%s] No recognised influenza segments in %s (expected a label from: %s).",
                  sample_id, args.input, ', '.join(EXPECTED_SEGMENTS))
        return EXIT_INPUT_ERROR

    if len(stems) > 1:
        log.error("[%s] Headers name more than one sample (%s); one file must hold one sample.",
                  sample_id, ', '.join(sorted(stems)))
        return EXIT_INPUT_ERROR

    records: list[SeqRecord] = []
    qc_rows: list[SegmentQC] = []

    for segment in EXPECTED_SEGMENTS:
        record = by_segment.get(segment)
        if record is None:
            log.warning("[%s] Missing expected segment: %s", sample_id, segment)
            qc_rows.append(SegmentQC(sample_id, segment, notes='segment not found in input'))
            continue
        cleaned, qc_row = process_segment(sample_id, segment, str(record.seq),
                                          mask_degenerate=args.mask_degenerate,
                                          max_n_frac=args.max_n_frac)
        records.append(cleaned)
        qc_rows.append(qc_row)

    SeqIO.write(records, str(args.output_fasta), 'fasta')

    with open(args.output_qc, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=[f.name for f in fields(SegmentQC)], delimiter='\t')
        writer.writeheader()
        writer.writerows(asdict(row) for row in qc_rows)

    log.info("[%s] Processed %d segment(s). Cleaned FASTA: %s. QC report: %s.",
             sample_id, len(records), args.output_fasta, args.output_qc)
    return EXIT_OK


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(levelname)s [prepare_consensus]: %(message)s',
        stream=sys.stderr,
    )
    return run(args)


if __name__ == '__main__':
    sys.exit(main())
