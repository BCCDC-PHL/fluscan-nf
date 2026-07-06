#!/usr/bin/env python3
"""
paftools_qc_check.py

Compute reference-coverage statistics for a minimap2 PAF alignment.
Emits a single-row TSV per (sample, segment) and writes a WARNING to stderr
when coverage falls below threshold.

Usage:
    paftools_qc_check.py --paf <sorted.paf> --sample <id> --segment <name>
                         --reference <ref.fasta>
                         --warn-thresh 0.90 --fail-thresh 0.50
                         --out-tsv <output.tsv>
"""

import argparse
import sys


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


def ref_len_from_fasta(fasta_path):
    """Return length of the first sequence in a FASTA file."""
    seq_len = 0
    started = False
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if started:
                    break
                started = True
            else:
                seq_len += len(line.strip())
    return seq_len


def main():
    parser = argparse.ArgumentParser(
        description='Compute reference coverage statistics from a PAF alignment.'
    )
    parser.add_argument('--paf',         required=True,  help='Sorted PAF alignment file')
    parser.add_argument('--sample',      required=True,  help='Sample identifier')
    parser.add_argument('--segment',     required=True,  help='Segment name')
    parser.add_argument('--reference',   default='',     help='Reference FASTA (used as fallback for ref length when PAF is empty)')
    parser.add_argument('--warn-thresh', type=float, default=0.90, help='Coverage fraction at/above which status is PASS; below it is WARN')
    parser.add_argument('--fail-thresh', type=float, default=0.50, help='Coverage fraction below which status is FAIL')
    parser.add_argument('--out-tsv',     required=True,  help='Output TSV path')
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

    # Fall back to FASTA for ref_len when PAF has no alignments
    if ref_len == 0 and args.reference:
        ref_len = ref_len_from_fasta(args.reference)

    aligned_bases = union_coverage(intervals)
    first_pos = min(s for s, _ in intervals) if intervals else 0
    last_pos  = max(e for _, e in intervals) if intervals else 0

    cov_frac      = aligned_bases / ref_len if ref_len > 0 else 0.0
    mean_identity = (1.0 - weighted_dv_sum / total_aln_len) * 100 if total_aln_len > 0 else 0.0

    if n_alignments == 0:
        status = 'FAIL'
    elif cov_frac >= args.warn_thresh:
        status = 'PASS'
    elif cov_frac >= args.fail_thresh:
        status = 'WARN'
    else:
        status = 'FAIL'

    header = [
        'sample_id', 'segment', 'reference', 'ref_length', 'aligned_ref_bases',
        'coverage_fraction', 'first_pos', 'last_pos', 'n_alignments',
        'mean_identity', 'status',
    ]
    row = [
        args.sample,
        args.segment,
        ref_name or args.reference,
        ref_len,
        aligned_bases,
        f'{cov_frac:.4f}',
        first_pos,
        last_pos,
        n_alignments,
        f'{mean_identity:.2f}',
        status,
    ]

    with open(args.out_tsv, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        fh.write('\t'.join(str(x) for x in row) + '\n')

    if status != 'PASS':
        print(
            f'WARNING [{args.sample}]: Alignment coverage {status} for segment {args.segment}: '
            f'coverage={cov_frac:.1%} ({aligned_bases}/{ref_len} ref bp), '
            f'first_pos={first_pos}, last_pos={last_pos}, n_alignments={n_alignments}',
            file=sys.stderr,
        )


if __name__ == '__main__':
    main()
