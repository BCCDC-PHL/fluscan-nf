#!/usr/bin/env python3
"""Choose a sample's reference bundle from a BLAST of its HA and NA against a bundle pool.

Pool subjects are headed bundle|accession|subtype|segment. For HA and NA separately, each
bundle's hit is its highest-bitscore alignment to the same segment. A bundle is assigned
only when it is the top bundle for HA and for NA, and both hits pass the identity and
query-coverage cutoffs. Otherwise the sample is unassigned and the bundle file is empty.
"""

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from lib.blast import BLAST_COLUMNS

SELECT_SEGMENTS = ['HA', 'NA']

OUTPUT_COLUMNS = [
    'sample_id', 'status', 'db', 'subtype',
    *[f'{segment}_{field}' for segment in SELECT_SEGMENTS
      for field in ('db', 'ref', 'pident', 'qcov', 'bitscore')],
    'runner_up_db', 'runner_up_bitscore',
]


def get_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('blastn', type=Path, help='blastn -outfmt 6 TSV of the sample against the pool')
    parser.add_argument('--query', type=Path, required=True, help='Cleaned consensus FASTA ({sample_id}_{SEGMENT} headers)')
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--min-identity', type=float, required=True, help='Minimum percent identity of each HA/NA hit')
    parser.add_argument('--min-coverage', type=float, required=True, help='Minimum fraction of each HA/NA query aligned')
    parser.add_argument('--out-tsv', type=Path, required=True, help='One-row selection report')
    parser.add_argument('--out-db', type=Path, required=True, help='File holding the chosen bundle name, empty when unassigned')
    return parser.parse_args(argv)


def read_hits(path: Path) -> pd.DataFrame:
    """BLAST rows of same-segment HA/NA hits, one per (segment, bundle): the best one."""
    columns = ['db', 'accession', 'subtype', 'segment']
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS + columns + ['qcov'])
    hits = pd.read_csv(path, sep='\t', names=BLAST_COLUMNS)
    hits['query_segment'] = hits['qseqid'].str.rsplit('_', n=1).str[-1]
    hits[columns] = hits['sseqid'].str.split('|', expand=True)
    hits = hits[(hits['query_segment'] == hits['segment']) & hits['segment'].isin(SELECT_SEGMENTS)].copy()
    hits['qcov'] = ((hits['qend'] - hits['qstart']).abs() + 1) / hits['qlen']
    # Highest bitscore first, ties broken by bundle name so the choice is deterministic.
    hits = hits.sort_values(['segment', 'bitscore', 'db'], ascending=[True, False, True])
    return hits.drop_duplicates(['segment', 'db'])


def select(sample_id: str, segments: set[str], hits: pd.DataFrame,
           min_identity: float, min_coverage: float) -> dict:
    """The selection report row; `db` is empty unless `status` is assigned."""
    row = dict.fromkeys(OUTPUT_COLUMNS, '')
    row['sample_id'] = sample_id

    best = {}
    for segment in SELECT_SEGMENTS:
        segment_hits = hits[hits['segment'] == segment]
        if segment_hits.empty:
            continue
        top = best[segment] = segment_hits.iloc[0]
        row[f'{segment}_db'] = top['db']
        row[f'{segment}_ref'] = top['accession']
        row[f'{segment}_pident'] = f"{top['pident']:.2f}"
        row[f'{segment}_qcov'] = f"{top['qcov']:.3f}"
        row[f'{segment}_bitscore'] = f"{top['bitscore']:g}"

    # Combined HA+NA bitscore per bundle, a missing segment counting zero.
    combined = (hits.groupby('db')['bitscore'].sum().reset_index()
                .sort_values(['bitscore', 'db'], ascending=[False, True]))
    winner = best['HA']['db'] if 'HA' in best else None
    runners_up = combined[combined['db'] != winner]
    if not runners_up.empty:
        row['runner_up_db'] = runners_up.iloc[0]['db']
        row['runner_up_bitscore'] = f"{runners_up.iloc[0]['bitscore']:g}"

    if not set(SELECT_SEGMENTS) <= segments:
        row['status'] = 'missing_segment'
    elif len(best) < len(SELECT_SEGMENTS):
        row['status'] = 'no_hit'
    elif best['HA']['db'] != best['NA']['db']:
        row['status'] = 'disagree'
    elif any(hit['pident'] < min_identity or hit['qcov'] < min_coverage for hit in best.values()):
        row['status'] = 'below_threshold'
    else:
        row['status'] = 'assigned'
        row['db'] = winner
        row['subtype'] = best['HA']['subtype']
    return row


def main(argv=None) -> int:
    args = get_args(argv)
    segments = {record.id.rsplit('_', 1)[-1] for record in SeqIO.parse(args.query, 'fasta')}
    hits = read_hits(args.blastn)
    row = select(args.sample_id, segments, hits, args.min_identity, args.min_coverage)

    with args.out_tsv.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerow(row)
    args.out_db.write_text(row['db'])

    print(f"{args.sample_id}: {row['status']} {row['db']}".rstrip(), file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
