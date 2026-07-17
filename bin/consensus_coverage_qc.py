#!/usr/bin/env python3
"""
consensus_coverage_qc.py

Compute reference-coverage statistics for a minimap2 PAF alignment.
Emits a single-row TSV per (sample, segment) and writes a WARNING to stderr
when coverage falls below threshold.

Usage:
    consensus_coverage_qc.py --paf <sorted.paf> --reference <ref.fasta>
                             --sample-id <id> --segment <name>
                             --pass-pct 90 --warn-pct 50 --out-tsv <output.tsv>
"""

import argparse
import sys
from pathlib import Path


def parse_tags(fields):
    """Return a dict of {tag_name: (type, value)} from PAF optional fields."""
    tags = {}
    for f in fields[12:]:
        parts = f.split(':', 2)
        if len(parts) == 3:
            tags[parts[0]] = (parts[1], parts[2])
    return tags


def union_coverage(intervals):
    """Total reference bases covered by a list of (start, end) half-open intervals."""
    if not intervals:
        return 0
    merged = sorted(intervals)
    total = 0
    cur_start, cur_end = merged[0]
    for start, end in merged[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
    total += cur_end - cur_start
    return total


def reference_info_from_fasta(fasta_path):
    """Return ``(record_id, length)`` for the first sequence in a FASTA file."""
    reference_id = None
    seq_len = 0
    started = False
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if started:
                    break
                started = True
                reference_id = line[1:].split()[0]
            else:
                seq_len += len(line.strip())
    return reference_id, seq_len


def main():
    parser = argparse.ArgumentParser(
        description='Compute reference coverage statistics from a PAF alignment.'
    )
    parser.add_argument('--paf', required=True, type=Path, help='Sorted PAF alignment file')
    parser.add_argument('--reference', required=True, type=Path, help='Single-segment reference FASTA')
    parser.add_argument('--sample-id', required=True, help='Sample/isolate identifier')
    parser.add_argument('--segment', required=True, help='Influenza segment label (e.g. HA)')
    parser.add_argument('--pass-pct', required=True, type=float, help='PASS coverage threshold (percent)')
    parser.add_argument('--warn-pct', required=True, type=float, help='WARN coverage threshold (percent)')
    parser.add_argument('--out-tsv', required=True, type=Path, help='Output single-row TSV')
    args = parser.parse_args()

    intervals = []
    ref_len = 0
    ref_name = ''
    n_alignments = 0
    weighted_dv_sum = 0.0
    total_aln_len = 0

    with open(args.paf) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 12:
                continue

            tags = parse_tags(fields)

            # Keep only primary alignments; if no tp tag present treat as primary
            tp = tags.get('tp', ('A', 'P'))[1]
            if tp not in ('P', 'p'):
                continue

            tlen   = int(fields[6])
            tstart = int(fields[7])
            tend   = int(fields[8])
            blen   = int(fields[10])

            if ref_len == 0:
                ref_len  = tlen
                ref_name = fields[5]

            intervals.append((tstart, tend))
            n_alignments += 1

            # Prefer dv:f: (per-base divergence), fall back to NM/blen
            if 'dv' in tags:
                dv = float(tags['dv'][1])
            elif 'NM' in tags and blen > 0:
                dv = int(tags['NM'][1]) / blen
            else:
                dv = 0.0

            weighted_dv_sum += dv * blen
            total_aln_len   += blen

    # Fall back to FASTA when PAF has no alignments.
    if ref_len == 0:
        ref_name, ref_len = reference_info_from_fasta(args.reference)

    aligned_bases = union_coverage(intervals)
    first_pos = min(s for s, _ in intervals) if intervals else 0
    last_pos  = max(e for _, e in intervals) if intervals else 0

    coverage_pct = aligned_bases / ref_len * 100 if ref_len > 0 else 0.0
    mean_identity_pct = (1.0 - weighted_dv_sum / total_aln_len) * 100 if total_aln_len > 0 else 0.0

    if n_alignments == 0:
        status = 'FAIL'
    elif coverage_pct >= args.pass_pct:
        status = 'PASS'
    elif coverage_pct >= args.warn_pct:
        status = 'WARN'
    else:
        status = 'FAIL'

    header = [
        'sample_id', 'segment', 'status', 'coverage_pct',
        'reference_id', 'reference_length', 'aligned_ref_bases', 'first_pos',
        'last_pos', 'n_alignments', 'mean_identity_pct',
    ]
    row = [
        args.sample_id,
        args.segment,
        status,
        f'{coverage_pct:.3f}',
        ref_name or 'NA',
        ref_len,
        aligned_bases,
        first_pos,
        last_pos,
        n_alignments,
        f'{mean_identity_pct:.3f}',
    ]

    with open(args.out_tsv, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        fh.write('\t'.join(str(x) for x in row) + '\n')

    if status != 'PASS':
        print(
            f'WARNING [{args.sample_id}]: Consensus coverage {status} for segment {args.segment}: '
            f'coverage={coverage_pct:.1f}% ({aligned_bases}/{ref_len} ref bp), '
            f'first_pos={first_pos}, last_pos={last_pos}, n_alignments={n_alignments}',
            file=sys.stderr,
        )


if __name__ == '__main__':
    main()
