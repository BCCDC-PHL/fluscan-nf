#!/usr/bin/env python3
"""
prepare_consensus.py

Validates and prepares influenza consensus FASTA files for the consensus workflow.

Steps performed:
  a) Check for presence of 8 expected influenza segments
    b) Enforce M label for the matrix segment by normalising legacy labels
  c) Replace partially degenerate IUPAC nucleotides (R,Y,S,W,K,M,B,D,H,V) with N
  d) Record the proportion of ambiguous (N) bases per sequence
  e) Enforce uppercase
  f) Check that each segment length falls within predefined ranges

Accepted input header formats:
  Underscore-delimited (legacy):
    >WIN-AH-2024-FAV-0415-1-OS_segment1_PB2
    >R3980461859-3370-A-D05_segment1_PB2
  Pipe-delimited (GISAID-style):
    >EPI_ISL_20399898|A_/_H5N1|NA|6|EPI5261179|2.3.4.4b_NA

Cleaned output headers are normalised to the canonical form:
  {sample_id}_{SEGMENT}
e.g. EPI-ISL-20399898_NA or WIN-AH-2024-FAV-0415-1-OS_PB2

For pipe-delimited headers, sample ID is taken from the first pipe field with
underscores replaced by dashes. Segment is extracted from the underscore suffix
of the field that contains a recognised influenza segment name.
"""

import argparse
import sys
import csv
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

EXPECTED_SEGMENTS = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS']

# All recognised segment labels, including the legacy MP alias for M
SEGMENT_NAMES = set(EXPECTED_SEGMENTS) | {'MP'}

# Partially degenerate IUPAC nucleotide codes (excluding N, which is already "any")
DEGENERATE_BASES = set('RYSWKMBDHV')

# Acceptable sequence length ranges (min_len, max_len) per segment
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


def parse_args():
    parser = argparse.ArgumentParser(
        description='Validate and clean influenza consensus FASTA for the consensus workflow.'
    )
    parser.add_argument('--input',      required=True,  help='Input consensus FASTA file')
    parser.add_argument('--sample_id',  required=True,  help='Sample identifier')
    parser.add_argument('--output',     default=None,   help='Output cleaned FASTA (default: {sample_id}_cleaned.fasta)')
    parser.add_argument('--output_qc',  default=None,   help='Output QC report TSV (default: {sample_id}_qc_report.tsv)')
    return parser.parse_args()


def detect_delimiter(records):
    """Detect the field delimiter used in FASTA headers.
    Defaults to '|' if any header contains a pipe, otherwise '_'.
    """
    for rec in records:
        if '|' in rec.id:
            return '|'
    return '_'


def get_segment_label(record_id, delimiter='_'):
    """Extract the segment label (e.g., PB2, HA, M) from the sequence ID.

    For pipe-delimited headers, each field is scanned for an underscore-separated
    suffix that matches a known influenza segment name, e.g.:
      '2.3.4.4b_NA' -> 'NA'
    For underscore-delimited headers, the last field is used directly.
    Legacy matrix label 'MP' is normalised to 'M'.
    Returns None if no recognised segment label can be found.
    """
    if delimiter == '|':
        for field in record_id.split('|'):
            if '_' in field:
                candidate = field.split('_')[-1].upper()
                if candidate in SEGMENT_NAMES:
                    return 'M' if candidate == 'MP' else candidate
        return None
    else:
        label = record_id.split(delimiter)[-1].upper()
        return 'M' if label == 'MP' else label


def make_canonical_header(record_id, delimiter, segment):
    """Build the canonical cleaned FASTA header: {normalized_sample_id}_{segment}.

    For pipe-delimited headers, the sample ID is the first pipe-delimited field
    with underscores replaced by dashes, e.g.:
      EPI_ISL_20399898|A_/_H5N1|NA|... -> EPI-ISL-20399898_NA

    For underscore-delimited headers, the sample ID is derived by stripping the
    trailing _segmentN_SEGMENT and/or _segmentN suffixes, e.g.:
      WIN-AH-2024-FAV-0415-1-OS_segment1_PB2 -> WIN-AH-2024-FAV-0415-1-OS_PB2
    """
    if delimiter == '|':
        sample_id = record_id.split('|')[0].replace('_', '-')
    else:
        sample_id = re.sub(r'_segment\d+_[A-Za-z0-9]+$', '', record_id)
        sample_id = re.sub(r'_segment\d+$', '', sample_id)
    return f"{sample_id}_{segment}"


def replace_degenerate(sequence):
    """Replace partially degenerate IUPAC bases with N. Returns (cleaned_seq, count_replaced)."""
    cleaned = []
    replaced = 0
    for base in sequence:
        if base in DEGENERATE_BASES:
            cleaned.append('N')
            replaced += 1
        else:
            cleaned.append(base)
    return ''.join(cleaned), replaced


def n_proportion(sequence):
    """Return the proportion of N bases in the sequence."""
    if len(sequence) == 0:
        return 0.0
    return sequence.count('N') / len(sequence)


def main():
    args = parse_args()

    output_fasta = args.output   or f"{args.sample_id}_cleaned.fasta"
    output_qc    = args.output_qc or f"{args.sample_id}_qc_report.tsv"

    records = list(SeqIO.parse(args.input, 'fasta'))

    if not records:
        print(f"ERROR: No sequences found in {args.input}", file=sys.stderr)
        sys.exit(1)

    delimiter = detect_delimiter(records)

    # Index records by segment label
    seg_to_record = {}
    for rec in records:
        segment = get_segment_label(rec.id, delimiter)
        if segment is None:
            print(f"WARNING [{args.sample_id}]: Could not determine segment for '{rec.id}'; skipping.",
                  file=sys.stderr)
            continue
        canonical_id = make_canonical_header(rec.id, delimiter, segment)
        rec.id = canonical_id
        rec.description = canonical_id
        rec.name = canonical_id
        seg_to_record[segment] = rec

    if not seg_to_record:
        print(f"ERROR [{args.sample_id}]: No valid influenza segments could be extracted from {args.input}. "
              f"Check that headers contain a recognised segment label (one of: {', '.join(EXPECTED_SEGMENTS)}).",
              file=sys.stderr)
        sys.exit(1)

    # --- a) Check for 8 expected segments ---
    found_segments = set(seg_to_record.keys())
    missing = [s for s in EXPECTED_SEGMENTS if s not in found_segments]
    unexpected = [s for s in found_segments if s not in EXPECTED_SEGMENTS]

    if missing:
        print(f"WARNING [{args.sample_id}]: Missing expected segments: {', '.join(missing)}", file=sys.stderr)
    if unexpected:
        print(f"WARNING [{args.sample_id}]: Unexpected segment labels found: {', '.join(unexpected)}", file=sys.stderr)

    qc_rows = []
    cleaned_records = []

    for seg in EXPECTED_SEGMENTS:
        if seg not in seg_to_record:
            continue

        rec = seg_to_record[seg]

        # --- e) Enforce uppercase ---
        seq_upper = str(rec.seq).upper()

        # --- d) N proportion before degenerate replacement ---
        n_prop_before = n_proportion(seq_upper)

        # --- c) Replace partially degenerate nucleotides with N ---
        seq_clean, n_replaced = replace_degenerate(seq_upper)

        # --- d) N proportion after degenerate replacement ---
        n_prop_after = n_proportion(seq_clean)

        # --- b) Canonical header already set at indexing time (MP normalised to M) ---
        new_id = rec.id

        # --- f) Check sequence length ---
        seq_len = len(seq_clean)
        length_ok = True
        if seg in SEGMENT_LENGTH_RANGES:
            lo, hi = SEGMENT_LENGTH_RANGES[seg]
            if not (lo <= seq_len <= hi):
                print(
                    f"WARNING [{args.sample_id}]: Segment {seg} length {seq_len} "
                    f"is outside expected range [{lo}, {hi}]",
                    file=sys.stderr
                )
                length_ok = False

        qc_rows.append({
            'sample_id':       args.sample_id,
            'segment':         seg,
            'seq_id':          new_id,
            'length':          seq_len,
            'length_ok':       length_ok,
            'n_prop_before':   round(n_prop_before, 6),
            'n_replaced':      n_replaced,
            'n_prop_after':    round(n_prop_after, 6),
        })

        cleaned_rec = SeqRecord(Seq(seq_clean), id=new_id, description='')
        cleaned_records.append(cleaned_rec)

    # Write cleaned FASTA
    SeqIO.write(cleaned_records, output_fasta, 'fasta')

    # Write QC report
    fieldnames = ['sample_id', 'segment', 'seq_id', 'length', 'length_ok',
                  'n_prop_before', 'n_replaced', 'n_prop_after']
    with open(output_qc, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(qc_rows)

    print(f"[{args.sample_id}] Processed {len(cleaned_records)} segments. "
          f"QC report: {output_qc}. Cleaned FASTA: {output_fasta}.")


if __name__ == '__main__':
    main()
