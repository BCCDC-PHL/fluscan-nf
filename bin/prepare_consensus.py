#!/usr/bin/env python3
"""
prepare_consensus.py

Validate and prepare influenza consensus FASTA files for the consensus workflow.

This is the upstream gatekeeper for the consensus workflows: it is the single
authoritative source for the isolate name (derived from the header, NOT the
filename, which is unreliable), it normalises headers to a canonical form, it
optionally masks ambiguous/degenerate bases, and it writes a per-isolate QC
report describing every expected segment.

Steps performed per file:
  a) Derive the authoritative isolate name from the headers (all records must
     agree; the script exits non-zero if they disagree or none can be parsed).
  b) Enforce the canonical M label for the matrix segment (legacy MP -> M).
  c) Uppercase sequences and, optionally, mask partially degenerate IUPAC
     nucleotides (R,Y,S,W,K,M,B,D,H,V) with N.
  d) Strip alignment gap characters (-, .) and coerce any other non-ACGTN
     character to N so the emitted FASTA is always clean.
  e) Rewrite headers to the canonical {isolate}_{SEGMENT} form.
  f) Emit a QC report row for every expected segment (present or MISSING) with a
     PASS/WARN/FAIL/MISSING status driven by segment length ranges and an
     ambiguity threshold.

Accepted input header formats:
  Underscore-delimited with _segmentN token (the wild west):
    >WIN-AH-2024-FAV-0415-1-OS_segment1_PB2
    >R3980461859-3370-A-D05_segment1_PB2
    >WIN-AH-2022-FAV-1687-1_OS_2_segment1_PB2   (isolate carries a copy number)
  Underscore-delimited, already canonical:
    >WIN-AH-2024-FAV-0415-1-OS_PB2
  Pipe-delimited (GISAID-style):
    >EPI_ISL_20399898|A_/_H5N1|NA|6|EPI5261179|2.3.4.4b_NA

Isolate derivation (underscores in the isolate are always replaced with dashes,
since underscores are not allowed in isolate names):
  - Pipe-delimited:               isolate = first pipe field.
  - Underscore with _segmentN:    isolate = everything BEFORE the first _segmentN
                                  (e.g. WIN-AH-2022-FAV-1687-1_OS_2_segment1_PB2
                                   -> WIN-AH-2022-FAV-1687-1-OS-2).
  - Underscore, already canonical: isolate = header with the trailing _SEGMENT
                                  stripped.

Outputs (written to --outdir):
  {isolate}.fasta            Cleaned, canonically-headed consensus FASTA.
  {isolate}_qc_report.tsv    One row per expected segment, with status.
  isolate_id.txt             The derived isolate name, for the caller to capture.
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from dataclasses import dataclass, fields
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --------------------------------------------------------------------------- #
# Configuration / constants
# --------------------------------------------------------------------------- #

EXPECTED_SEGMENTS = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS']

# Legacy/alternate segment labels mapped to their canonical name.
SEGMENT_ALIASES = {'MP': 'M'}

# Partially degenerate IUPAC nucleotide codes (excluding N, which is already "any").
DEGENERATE_BASES = set('RYSWKMBDHV')

# Unambiguous bases plus the fully-ambiguous N.
VALID_BASES = set('ACGTN')

# Alignment gap / pad characters that must never reach downstream tools.
GAP_CHARS = set('-.')

# Anchor for underscore-delimited "wild west" headers: everything before the
# first _segmentN token is the isolate name.
SEGMENT_TOKEN_RE = re.compile(r'_segment\d+', re.IGNORECASE)

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


# --------------------------------------------------------------------------- #
# Data structures
# --------------------------------------------------------------------------- #

@dataclass
class CleanResult:
    """Outcome of cleaning a single raw sequence string."""
    sequence: str
    n_replaced: int      # degenerate IUPAC bases masked to N
    n_gaps_removed: int  # gap/pad characters stripped
    n_invalid: int       # other non-ACGTN characters coerced to N
    n_prop_before: float
    n_prop_after: float


@dataclass
class SegmentQC:
    """One row of the per-isolate QC report."""
    sample_id: str
    segment: str
    seq_id: str
    status: str
    length: int
    length_ok: bool
    n_prop_before: float
    n_replaced: int
    n_gaps_removed: int
    n_invalid: int
    n_prop_after: float
    notes: str = ''

    @classmethod
    def fieldnames(cls) -> list[str]:
        return [f.name for f in fields(cls)]

    def as_row(self) -> dict:
        return {
            'sample_id':      self.sample_id,
            'segment':        self.segment,
            'seq_id':         self.seq_id,
            'status':         self.status,
            'length':         self.length,
            'length_ok':      self.length_ok,
            'n_prop_before':  round(self.n_prop_before, 6),
            'n_replaced':     self.n_replaced,
            'n_gaps_removed': self.n_gaps_removed,
            'n_invalid':      self.n_invalid,
            'n_prop_after':   round(self.n_prop_after, 6),
            'notes':          self.notes,
        }


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Validate and clean influenza consensus FASTA for the consensus workflow.'
    )
    parser.add_argument('--input', required=True, type=Path,
                        help='Input consensus FASTA file')
    parser.add_argument('--sample_id', '--sample-id', dest='sample_id', default=None,
                        help='Provisional sample identifier, used only for log messages '
                             '(the authoritative isolate is derived from the header)')
    parser.add_argument('--outdir', default=Path('.'), type=Path,
                        help='Output directory for {isolate}.fasta, {isolate}_qc_report.tsv '
                             'and isolate_id.txt (default: current directory)')
    parser.add_argument('--mask-degenerate', action=argparse.BooleanOptionalAction, default=True,
                        help='Mask partially degenerate IUPAC bases (R,Y,S,W,K,M,B,D,H,V) with N '
                             '(default: enabled; use --no-mask-degenerate to keep them)')
    parser.add_argument('--strip-gaps', action=argparse.BooleanOptionalAction, default=True,
                        help='Strip alignment gap characters (-, .) from sequences '
                             '(default: enabled)')
    parser.add_argument('--max-n-frac', type=float, default=1.0,
                        help='WARN when the post-clean N fraction of a segment exceeds this value '
                             '(default: 1.0, i.e. never warns)')
    parser.add_argument('--expected-segments', default=','.join(EXPECTED_SEGMENTS),
                        help='Comma-separated ordered list of expected segment labels '
                             f'(default: {",".join(EXPECTED_SEGMENTS)})')
    parser.add_argument('--log-level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help='Logging verbosity (default: INFO)')
    return parser.parse_args(argv)


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format='%(levelname)s [prepare_consensus]: %(message)s',
        stream=sys.stderr,
    )


# --------------------------------------------------------------------------- #
# Parsing / isolate resolution
# --------------------------------------------------------------------------- #

def load_records(path: Path) -> list[SeqRecord]:
    """Parse a FASTA file into a list of records (may be empty)."""
    return list(SeqIO.parse(str(path), 'fasta'))


def detect_delimiter(records: list[SeqRecord]) -> str:
    """Detect the header field delimiter: '|' if any header contains a pipe, else '_'."""
    for rec in records:
        if '|' in rec.id:
            return '|'
    return '_'


def get_segment_label(record_id: str, delimiter: str = '_') -> str | None:
    """Extract the canonical segment label (e.g. PB2, HA, M) from a sequence ID.

    For pipe-delimited headers, each field is scanned for an underscore-separated
    suffix that matches a known segment name (e.g. '2.3.4.4b_NA' -> 'NA').
    For underscore-delimited headers, the last field is used directly.
    Legacy aliases (e.g. MP) are normalised via SEGMENT_ALIASES.
    Returns None if no recognised segment label can be found.
    """
    known = set(EXPECTED_SEGMENTS) | set(SEGMENT_ALIASES)
    if delimiter == '|':
        for field_str in record_id.split('|'):
            if '_' in field_str:
                candidate = field_str.split('_')[-1].upper()
                if candidate in known:
                    return SEGMENT_ALIASES.get(candidate, candidate)
        return None
    label = record_id.split(delimiter)[-1].upper()
    return SEGMENT_ALIASES.get(label, label)


def derive_isolate(record_id: str, delimiter: str, segment: str | None) -> str:
    """Derive the isolate name from a single FASTA header.

    Underscores are always replaced with dashes (underscores are not allowed in
    isolate names). See the module docstring for the per-format rules.
    """
    if delimiter == '|':
        isolate = record_id.split('|')[0]
    else:
        m = SEGMENT_TOKEN_RE.search(record_id)
        if m:
            isolate = record_id[:m.start()]
        elif segment and record_id.upper().endswith(f"_{segment}"):
            isolate = record_id[:-(len(segment) + 1)]
        else:
            isolate = record_id
    return isolate.replace('_', '-')


def resolve_isolate(records: list[SeqRecord], delimiter: str,
                    expected_segments: list[str]) -> tuple[str, dict[str, SeqRecord]]:
    """Map records to segments and derive the single authoritative isolate.

    Raises ValueError on fatal input problems: no parseable segments, or records
    that resolve to more than one isolate. Unlabelled records and duplicate
    segments are warned about (the first occurrence of each segment wins).
    """
    seg_to_record: dict[str, SeqRecord] = {}
    isolates: set[str] = set()

    for rec in records:
        segment = get_segment_label(rec.id, delimiter)
        if segment is None:
            log.warning("Could not determine segment for '%s'; skipping.", rec.id)
            continue

        isolate = derive_isolate(rec.id, delimiter, segment)
        if not isolate:
            raise ValueError(f"Could not derive an isolate name from header '{rec.id}'.")
        isolates.add(isolate)

        if segment in seg_to_record:
            log.warning("Duplicate record for segment %s ('%s'); keeping the first "
                        "occurrence ('%s') and ignoring this one.",
                        segment, rec.id, seg_to_record[segment].id)
            continue
        seg_to_record[segment] = rec

    if not seg_to_record:
        raise ValueError(
            "No valid influenza segments could be extracted. Check that headers "
            f"contain a recognised segment label (one of: {', '.join(expected_segments)})."
        )

    if len(isolates) > 1:
        raise ValueError(
            f"Records resolved to multiple isolate names ({', '.join(sorted(isolates))}); "
            "refusing to guess. Check the input headers."
        )

    return next(iter(isolates)), seg_to_record


# --------------------------------------------------------------------------- #
# Sequence cleaning / QC
# --------------------------------------------------------------------------- #

def n_proportion(sequence: str) -> float:
    """Return the proportion of N bases in the sequence (0.0 for empty)."""
    if not sequence:
        return 0.0
    return sequence.count('N') / len(sequence)


def clean_sequence(sequence: str, *, mask_degenerate: bool, strip_gaps: bool) -> CleanResult:
    """Uppercase and clean a raw sequence string.

    Optionally strips gap characters, optionally masks degenerate IUPAC bases to
    N, and always coerces any remaining non-ACGTN character to N so the emitted
    FASTA is guaranteed clean. N proportions are reported before and after.
    """
    upper = sequence.upper()
    n_prop_before = n_proportion(upper)

    cleaned: list[str] = []
    n_replaced = n_gaps_removed = n_invalid = 0

    for base in upper:
        if base in GAP_CHARS:
            if strip_gaps:
                n_gaps_removed += 1
                continue
            # Not stripping: an in-sequence gap is not a valid base -> mask.
            cleaned.append('N')
            n_invalid += 1
        elif base in DEGENERATE_BASES:
            if mask_degenerate:
                cleaned.append('N')
                n_replaced += 1
            else:
                cleaned.append(base)
        elif base in VALID_BASES:
            cleaned.append(base)
        else:
            cleaned.append('N')
            n_invalid += 1

    result_seq = ''.join(cleaned)
    return CleanResult(
        sequence=result_seq,
        n_replaced=n_replaced,
        n_gaps_removed=n_gaps_removed,
        n_invalid=n_invalid,
        n_prop_before=n_prop_before,
        n_prop_after=n_proportion(result_seq),
    )


def evaluate_segment(isolate: str, segment: str, seq_id: str, clean: CleanResult,
                     *, max_n_frac: float) -> SegmentQC:
    """Build the QC row for one present segment, assigning a PASS/WARN/FAIL status."""
    length = len(clean.sequence)

    length_ok = True
    lo = hi = None
    if segment in SEGMENT_LENGTH_RANGES:
        lo, hi = SEGMENT_LENGTH_RANGES[segment]
        length_ok = lo <= length <= hi

    notes: list[str] = []
    status = STATUS_PASS

    # FAIL conditions.
    if length == 0 or clean.n_prop_after >= 1.0:
        status = STATUS_FAIL
        notes.append('empty or all-N sequence')
    elif not length_ok:
        status = STATUS_FAIL
        notes.append(f'length {length} outside expected range [{lo}, {hi}]')

    # WARN conditions (only downgrade from PASS, never override FAIL).
    if status == STATUS_PASS:
        if clean.n_prop_after > max_n_frac:
            status = STATUS_WARN
            notes.append(f'N fraction {clean.n_prop_after:.3f} exceeds {max_n_frac:g}')
        if clean.n_invalid:
            status = STATUS_WARN
            notes.append(f'{clean.n_invalid} invalid base(s) coerced to N')
        if clean.n_gaps_removed:
            status = STATUS_WARN
            notes.append(f'{clean.n_gaps_removed} gap char(s) stripped')

    if status != STATUS_PASS:
        log.warning("Segment %s [%s]: %s", segment, isolate, '; '.join(notes) or status)

    return SegmentQC(
        sample_id=isolate,
        segment=segment,
        seq_id=seq_id,
        status=status,
        length=length,
        length_ok=length_ok,
        n_prop_before=clean.n_prop_before,
        n_replaced=clean.n_replaced,
        n_gaps_removed=clean.n_gaps_removed,
        n_invalid=clean.n_invalid,
        n_prop_after=clean.n_prop_after,
        notes='; '.join(notes),
    )


def missing_segment_qc(isolate: str, segment: str) -> SegmentQC:
    """Build a placeholder QC row for an expected segment that is absent."""
    return SegmentQC(
        sample_id=isolate,
        segment=segment,
        seq_id='',
        status=STATUS_MISSING,
        length=0,
        length_ok=False,
        n_prop_before=0.0,
        n_replaced=0,
        n_gaps_removed=0,
        n_invalid=0,
        n_prop_after=0.0,
        notes='segment not found in input',
    )


def build_cleaned_record(isolate: str, segment: str, sequence: str) -> SeqRecord:
    """Build a SeqRecord with the canonical {isolate}_{SEGMENT} header."""
    canonical_id = f"{isolate}_{segment}"
    return SeqRecord(Seq(sequence), id=canonical_id, description='')


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #

def write_fasta(records: list[SeqRecord], path: Path) -> None:
    SeqIO.write(records, str(path), 'fasta')


def write_qc_tsv(rows: list[SegmentQC], path: Path) -> None:
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=SegmentQC.fieldnames(), delimiter='\t')
        writer.writeheader()
        writer.writerows(row.as_row() for row in rows)


def write_isolate_id(isolate: str, path: Path) -> None:
    with open(path, 'w') as fh:
        fh.write(isolate + '\n')


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #

def run(args: argparse.Namespace) -> int:
    tag = args.sample_id or args.input.name
    expected_segments = [s.strip().upper() for s in args.expected_segments.split(',') if s.strip()]

    args.outdir.mkdir(parents=True, exist_ok=True)

    records = load_records(args.input)
    if not records:
        log.error("[%s] No sequences found in %s", tag, args.input)
        return EXIT_INPUT_ERROR

    delimiter = detect_delimiter(records)

    try:
        isolate, seg_to_record = resolve_isolate(records, delimiter, expected_segments)
    except ValueError as exc:
        log.error("[%s] %s", tag, exc)
        return EXIT_INPUT_ERROR

    # Authoritative isolate resolved: write the sidecar immediately so the caller
    # can always capture it, even if later steps only produce warnings.
    write_isolate_id(isolate, args.outdir / 'isolate_id.txt')

    found = set(seg_to_record)
    unexpected = sorted(found - set(expected_segments))
    if unexpected:
        log.warning("[%s] Unexpected segment labels found: %s", isolate, ', '.join(unexpected))

    qc_rows: list[SegmentQC] = []
    cleaned_records: list[SeqRecord] = []

    for segment in expected_segments:
        rec = seg_to_record.get(segment)
        if rec is None:
            log.warning("[%s] Missing expected segment: %s", isolate, segment)
            qc_rows.append(missing_segment_qc(isolate, segment))
            continue

        clean = clean_sequence(str(rec.seq),
                               mask_degenerate=args.mask_degenerate,
                               strip_gaps=args.strip_gaps)
        cleaned_rec = build_cleaned_record(isolate, segment, clean.sequence)
        cleaned_records.append(cleaned_rec)
        qc_rows.append(evaluate_segment(isolate, segment, cleaned_rec.id, clean,
                                        max_n_frac=args.max_n_frac))

    # Any unexpected-but-present segments are still cleaned and emitted so no data
    # is silently dropped from the output FASTA.
    for segment in unexpected:
        rec = seg_to_record[segment]
        clean = clean_sequence(str(rec.seq),
                               mask_degenerate=args.mask_degenerate,
                               strip_gaps=args.strip_gaps)
        cleaned_rec = build_cleaned_record(isolate, segment, clean.sequence)
        cleaned_records.append(cleaned_rec)
        qc_rows.append(evaluate_segment(isolate, segment, cleaned_rec.id, clean,
                                        max_n_frac=args.max_n_frac))

    output_fasta = args.outdir / f"{isolate}.fasta"
    output_qc = args.outdir / f"{isolate}_qc_report.tsv"
    write_fasta(cleaned_records, output_fasta)
    write_qc_tsv(qc_rows, output_qc)

    log.info("[%s] Processed %d segment(s). Cleaned FASTA: %s. QC report: %s.",
             isolate, len(cleaned_records), output_fasta, output_qc)
    return EXIT_OK


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    configure_logging(args.log_level)
    return run(args)


if __name__ == '__main__':
    sys.exit(main())
