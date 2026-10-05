#!/usr/bin/env python3
import os
import sys
import pandas as pd
import argparse
from Bio import SeqIO

# Must match the -outfmt the blastn process passes
BLAST_COLUMNS = 'qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore score'.split()

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('blastn', help='BlastN results file in TSV format')
	parser.add_argument('-s', '--seqs', required=True, help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	parser.add_argument('-o','--out_tsv', required=True, help='Filtered results output')
	parser.add_argument('-O','--out_fasta', required=True, help='Filtered FASTA output')
	return parser.parse_args()

def parse_blast(filepath):
	blast_df = pd.read_csv(filepath, sep='\t', names=BLAST_COLUMNS)

	blast_df.insert(0, 'sample_name', os.path.basename(filepath).split("_")[0])

	blast_df.insert(3, 'segment', blast_df['sseqid'].str.split("|").str[-1])
	return blast_df


def filter_alignments(blast_results):
	'''Return the highest-bitscore hit for each reference segment, best first.'''
	blast_results = blast_results.drop_duplicates()
	filtered = blast_results.loc[blast_results.groupby('segment')['bitscore'].idxmax()]

	if filtered.shape[0] != 8:
		print(f"WARNING: expected 8 segments, found {filtered.shape[0]}", file=sys.stderr)

	return filtered.sort_values('bitscore', ascending=False).reset_index(drop=True)

def write_best_sequence(blast_results, database_path, fasta_out):
	'''
	Looks up the single top hit in the DB FASTA file and writes them to their own FASTA file.
	'''
	db_seqs = SeqIO.to_dict(SeqIO.parse(database_path, 'fasta'))

	final_references = [db_seqs[ref_name] for ref_name in blast_results['sseqid']]

	SeqIO.write(final_references, fasta_out, 'fasta')


def main():
	args = parse_args()

	blast_df = parse_blast(args.blastn)

	if blast_df.shape[0] == 0:
		print("WARNING: No data found in blast input file. Exiting.", file=sys.stderr)
		sys.exit(1)

	blast_df = filter_alignments(blast_df)

	blast_df.to_csv(args.out_tsv, index=False, sep='\t')

	write_best_sequence(blast_df, args.seqs, args.out_fasta)


if __name__ == '__main__':
	main()
