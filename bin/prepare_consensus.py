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

Expected input header format examples:
  >WIN-AH-2024-FAV-0415-1-OS_segment1_PB2
  >R3980461859-3370-A-D05_segment1_PB2
The segment name is the last underscore-delimited field of the sequence ID.
"""

import argparse
import sys
import csv
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

EXPECTED_SEGMENTS = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS']

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


def clean_header(record_id):
    """Clean the sequence ID by removing the trailing segment ordinal."""
    new_id = re.sub(r"_segment\d", "", record_id)
    return new_id

def get_segment_label(record_id, delimiter='_'):
    """Extract the segment label (e.g., PB2, HA, M) from the sequence ID."""
    label = record_id.split(delimiter)[-1].upper()
    # Normalise legacy matrix labels to M.
    if label == 'MP':
        label = 'M'
    return label


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
        new_id = clean_header(rec.id)  # Clean the header for output

        rec.id = new_id  # Update the record 
        rec.description = new_id
        rec.name = new_id

        seg_to_record[segment] = rec

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

        # --- b) Enforce M label in the header ---
        new_id = rec.id
        label_in_id = rec.id.split(delimiter)[-1].upper()
        if label_in_id == 'MP':
            # Replace the trailing legacy matrix label with M in the sequence ID.
            new_id = delimiter.join(rec.id.split(delimiter)[:-1]) + delimiter + 'M'

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
