#!/usr/bin/env python3
"""
aggregate_coverage_report.py

Concatenate per-segment coverage TSVs (from paf_coverage_check.py) into a
run-level summary TSV and a short Markdown report highlighting WARN/FAIL rows.

Usage:
    aggregate_coverage_report.py [tsv ...] --out-tsv run_coverage_summary.tsv
                                           --out-md  run_coverage_summary.md
"""

import argparse
import csv
import os
import sys


STATUS_ORDER = {'FAIL_NO_ALIGN': 0, 'WARN_PARTIAL': 1, 'PASS': 2}


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate per-segment alignment-coverage TSVs into a run summary.'
    )
    parser.add_argument('inputs', nargs='*', help='Per-segment coverage TSV files (positional)')
    parser.add_argument('--out-tsv', required=True, help='Output run-level TSV path')
    parser.add_argument('--out-md',  required=True, help='Output Markdown summary path')
    args = parser.parse_args()

    # Accept input files from positional args or glob current dir as fallback
    input_files = args.inputs
    if not input_files:
        import glob
        input_files = glob.glob('*_coverage.tsv')

    if not input_files:
        print('WARNING: No coverage TSV files found; writing empty summary.', file=sys.stderr)

    rows = []
    header = None

    for path in sorted(input_files):
        try:
            with open(path, newline='') as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                if header is None:
                    header = reader.fieldnames
                for row in reader:
                    rows.append(row)
        except OSError as e:
            print(f'WARNING: Could not read {path}: {e}', file=sys.stderr)

    # Sort: FAIL first, then WARN, then PASS; within status by sample_id / segment
    rows.sort(key=lambda r: (
        STATUS_ORDER.get(r.get('status', 'PASS'), 3),
        r.get('sample_id', ''),
        r.get('segment', ''),
    ))

    # ------------------------------------------------------------------
    # Write TSV
    # ------------------------------------------------------------------
    if header is None:
        header = [
            'sample_id', 'segment', 'reference', 'ref_length', 'aligned_ref_bases',
            'coverage_fraction', 'first_pos', 'last_pos', 'n_alignments',
            'mean_identity', 'status',
        ]

    with open(args.out_tsv, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

    # ------------------------------------------------------------------
    # Write Markdown
    # ------------------------------------------------------------------
    n_total = len(rows)
    n_pass  = sum(1 for r in rows if r.get('status') == 'PASS')
    n_warn  = sum(1 for r in rows if r.get('status') == 'WARN_PARTIAL')
    n_fail  = sum(1 for r in rows if r.get('status') == 'FAIL_NO_ALIGN')

    flagged = [r for r in rows if r.get('status') != 'PASS']

    with open(args.out_md, 'w') as fh:
        fh.write('# Alignment Coverage Summary\n\n')
        fh.write(f'| Status | Count |\n|---|---|\n')
        fh.write(f'| PASS | {n_pass} |\n')
        fh.write(f'| WARN_PARTIAL | {n_warn} |\n')
        fh.write(f'| FAIL_NO_ALIGN | {n_fail} |\n')
        fh.write(f'| **Total** | **{n_total}** |\n\n')

        if not flagged:
            fh.write('All segments passed alignment coverage thresholds.\n')
        else:
            fh.write('## Flagged segments (WARN / FAIL)\n\n')
            fh.write('| Sample | Segment | Coverage | Aligned bp | Ref bp | first_pos | last_pos | Status |\n')
            fh.write('|---|---|---|---|---|---|---|---|\n')
            for r in flagged:
                cov_pct = f"{float(r.get('coverage_fraction', 0)) * 100:.1f}%"
                fh.write(
                    f"| {r.get('sample_id','')} "
                    f"| {r.get('segment','')} "
                    f"| {cov_pct} "
                    f"| {r.get('aligned_ref_bases','')} "
                    f"| {r.get('ref_length','')} "
                    f"| {r.get('first_pos','')} "
                    f"| {r.get('last_pos','')} "
                    f"| **{r.get('status','')}** |\n"
                )


if __name__ == '__main__':
    main()
