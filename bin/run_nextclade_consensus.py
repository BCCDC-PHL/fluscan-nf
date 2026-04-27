#!/usr/bin/env python3
"""
run_nextclade_consensus.py

For each segment in a consensus FASTA, determine the appropriate Nextclade
dataset from BLAST results and run Nextclade.

The BLAST filter TSV (produced by filter_alignments.py) contains a 'sseqid'
column with the format:
  accession|subtype|segment
  e.g. NC_026438.1|H1N1|PB2

The Nextclade datasets directory is expected to be organised as:
  <datasets_dir>/<subtype_lower>/<segment_lower>/
  e.g. /path/to/nextclade_datasets/h1n1/pb2/
"""

import argparse
import csv
import os
import subprocess
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description='Run Nextclade for each segment of a consensus FASTA.'
    )
    parser.add_argument('--fasta',        required=True, help='Cleaned consensus FASTA file')
    parser.add_argument('--blast_filter', required=True, help='Filtered BLAST TSV from filter_alignments.py')
    parser.add_argument('--datasets_dir', required=True, help='Root directory of Nextclade datasets')
    parser.add_argument('--sample_id',    required=True, help='Sample identifier')
    return parser.parse_args()


def load_segment_subtype_map(blast_filter_path):
    """
    Parse the filtered BLAST TSV to build a mapping of
    segment_lower -> subtype_lower from the best-hit reference headers.
    """
    mapping = {}
    with open(blast_filter_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sseqid = row.get('sseqid', '')
            parts = sseqid.split('|')
            if len(parts) < 2:
                print(f"WARNING: Unexpected sseqid format '{sseqid}'; skipping.", file=sys.stderr)
                continue
            subtype = parts[-2].lower()
            segment = parts[-1].lower()
            mapping[segment] = subtype
    return mapping


def load_consensus_by_segment(fasta_path):
    """
    Parse the consensus FASTA and index sequences by their segment label
    (last underscore-delimited field of the sequence ID, upper-cased).
    M is normalised to MP.
    """
    seqs = {}
    for rec in SeqIO.parse(fasta_path, 'fasta'):
        seg_label = rec.id.split('_')[-1].upper()
        if seg_label == 'M':
            seg_label = 'MP'
        seqs[seg_label] = rec
    return seqs


def main():
    args = parse_args()

    seg_subtype = load_segment_subtype_map(args.blast_filter)
    if not seg_subtype:
        print(f"ERROR: No segment-subtype mappings found in {args.blast_filter}", file=sys.stderr)
        sys.exit(1)

    seqs = load_consensus_by_segment(args.fasta)
    if not seqs:
        print(f"ERROR: No sequences found in {args.fasta}", file=sys.stderr)
        sys.exit(1)

    for segment, subtype in seg_subtype.items():
        seg_upper = segment.upper()
        if seg_upper not in seqs:
            print(f"WARNING [{args.sample_id}]: Segment {seg_upper} not found in consensus FASTA; skipping.", file=sys.stderr)
            continue

        # Write individual segment FASTA
        seg_fasta = f"{args.sample_id}_{segment}_tmp.fasta"
        SeqIO.write([seqs[seg_upper]], seg_fasta, 'fasta')

        dataset = os.path.join(args.datasets_dir, subtype, segment)
        if not os.path.isdir(dataset):
            print(f"WARNING [{args.sample_id}]: Nextclade dataset not found for {subtype}/{segment}: {dataset}", file=sys.stderr)
            continue

        output_tsv = f"{args.sample_id}_nextclade_{segment}.tsv"
        cmd = [
            'nextclade', 'run',
            '--input-dataset', dataset,
            '--output-tsv', output_tsv,
            seg_fasta,
        ]
        print(f"[{args.sample_id}] Running nextclade for {subtype}/{segment}: {' '.join(cmd)}")
        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"WARNING [{args.sample_id}]: Nextclade returned non-zero exit code for {subtype}/{segment}.", file=sys.stderr)


if __name__ == '__main__':
    main()
