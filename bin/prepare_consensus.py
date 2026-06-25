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
  Underscore-delimited with _segmentN token (the wild west):
    >WIN-AH-2024-FAV-0415-1-OS_segment1_PB2
    >R3980461859-3370-A-D05_segment1_PB2
    >WIN-AH-2022-FAV-1687-1_OS_2_segment1_PB2   (isolate carries a copy number)
  Underscore-delimited, already canonical:
    >WIN-AH-2024-FAV-0415-1-OS_PB2
  Pipe-delimited (GISAID-style):
    >EPI_ISL_20399898|A_/_H5N1|NA|6|EPI5261179|2.3.4.4b_NA

This script is the single authoritative source for the isolate name. The isolate
is derived from the header (NOT the filename, which is unreliable) and used for
both the rewritten headers and the output filenames:

  - Pipe-delimited:               isolate = first pipe field.
  - Underscore with _segmentN:    isolate = everything BEFORE the first _segmentN
                                  (e.g. WIN-AH-2022-FAV-1687-1_OS_2_segment1_PB2
                                   -> WIN-AH-2022-FAV-1687-1_OS_2).
  - Underscore, already canonical: isolate = header with the trailing _SEGMENT
                                  stripped.
  In every case underscores in the isolate are replaced with dashes (underscores
  are not allowed in isolate names).

Cleaned output headers are normalised to the canonical form {isolate}_{SEGMENT},
e.g. WIN-AH-2022-FAV-1687-1-OS-2_NA or EPI-ISL-20399898_NA.

All records in one file must resolve to the same isolate; the script exits
non-zero (strict) if they disagree or if no isolate can be parsed. The derived
isolate is written to {outdir}/isolate_id.txt so the caller can capture it.
"""

import argparse
import os
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

# Anchor for underscore-delimited "wild west" headers: everything before the
# first _segmentN token is the isolate name.
SEGMENT_TOKEN_RE = re.compile(r'_segment\d+', re.IGNORECASE)

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
    parser.add_argument('--sample_id',  default=None,   help='Provisional sample identifier, used only for log messages '
                                                              '(the authoritative isolate is derived from the header)')
    parser.add_argument('--outdir',     default='.',    help='Output directory for {isolate}.fasta, {isolate}_qc_report.tsv '
                                                              'and isolate_id.txt (default: current directory)')
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


def derive_isolate(record_id, delimiter, segment):
    """Derive the isolate name from a single FASTA header.

    The isolate is the authoritative sample identifier; underscores are always
    replaced with dashes since underscores are not allowed in isolate names.

    Pipe-delimited:
      EPI_ISL_20399898|A_/_H5N1|NA|... -> EPI-ISL-20399898

    Underscore-delimited with a _segmentN token (everything before it is the
    isolate, including any copy number):
      WIN-AH-2022-FAV-1687-1_OS_2_segment1_PB2 -> WIN-AH-2022-FAV-1687-1-OS-2
      WIN-AH-2024-FAV-0415-1-OS_segment1_PB2   -> WIN-AH-2024-FAV-0415-1-OS

    Underscore-delimited, already canonical (no _segmentN token) — strip the
    trailing recognised _SEGMENT:
      WIN-AH-2024-FAV-0415-1-OS_PB2 -> WIN-AH-2024-FAV-0415-1-OS
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


def make_canonical_header(isolate, segment):
    """Build the canonical cleaned FASTA header: {isolate}_{segment}."""
    return f"{isolate}_{segment}"


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

    # Provisional tag used only for log messages until the authoritative isolate
    # is derived from the headers.
    tag = args.sample_id or os.path.basename(args.input)

    os.makedirs(args.outdir, exist_ok=True)

    records = list(SeqIO.parse(args.input, 'fasta'))

    if not records:
        print(f"ERROR [{tag}]: No sequences found in {args.input}", file=sys.stderr)
        sys.exit(1)

    delimiter = detect_delimiter(records)

    # Index records by segment label, deriving the isolate from each header.
    seg_to_record = {}
    isolates = set()
    for rec in records:
        segment = get_segment_label(rec.id, delimiter)
        if segment is None:
            print(f"WARNING [{tag}]: Could not determine segment for '{rec.id}'; skipping.",
                  file=sys.stderr)
            continue
        isolate = derive_isolate(rec.id, delimiter, segment)
        if not isolate:
            print(f"ERROR [{tag}]: Could not derive an isolate name from header '{rec.id}'.",
                  file=sys.stderr)
            sys.exit(1)
        isolates.add(isolate)
        seg_to_record[segment] = (isolate, rec)

    if not seg_to_record:
        print(f"ERROR [{tag}]: No valid influenza segments could be extracted from {args.input}. "
              f"Check that headers contain a recognised segment label (one of: {', '.join(EXPECTED_SEGMENTS)}).",
              file=sys.stderr)
        sys.exit(1)

    # --- Strict: all records in a file must resolve to the same isolate ---
    if len(isolates) > 1:
        print(f"ERROR [{tag}]: Records resolved to multiple isolate names ({', '.join(sorted(isolates))}); "
              f"refusing to guess. Check the input headers in {args.input}.",
              file=sys.stderr)
        sys.exit(1)

    isolate = next(iter(isolates))

    # Rewrite each record's header to the canonical {isolate}_{SEGMENT} form.
    for segment, (_iso, rec) in seg_to_record.items():
        canonical_id = make_canonical_header(isolate, segment)
        rec.id = canonical_id
        rec.description = canonical_id
        rec.name = canonical_id
    seg_to_record = {seg: rec for seg, (_iso, rec) in seg_to_record.items()}

    output_fasta = os.path.join(args.outdir, f"{isolate}.fasta")
    output_qc    = os.path.join(args.outdir, f"{isolate}_qc_report.tsv")

    # --- a) Check for 8 expected segments ---
    found_segments = set(seg_to_record.keys())
    missing = [s for s in EXPECTED_SEGMENTS if s not in found_segments]
    unexpected = [s for s in found_segments if s not in EXPECTED_SEGMENTS]

    if missing:
        print(f"WARNING [{isolate}]: Missing expected segments: {', '.join(missing)}", file=sys.stderr)
    if unexpected:
        print(f"WARNING [{isolate}]: Unexpected segment labels found: {', '.join(unexpected)}", file=sys.stderr)

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
                    f"WARNING [{isolate}]: Segment {seg} length {seq_len} "
                    f"is outside expected range [{lo}, {hi}]",
                    file=sys.stderr
                )
                length_ok = False

        qc_rows.append({
            'sample_id':       isolate,
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

    # Sidecar so the caller (e.g. the Nextflow process) can capture the
    # authoritative isolate name derived from the header.
    with open(os.path.join(args.outdir, 'isolate_id.txt'), 'w') as f:
        f.write(isolate + '\n')

    print(f"[{isolate}] Processed {len(cleaned_records)} segments. "
          f"QC report: {output_qc}. Cleaned FASTA: {output_fasta}.")


if __name__ == '__main__':
    main()
