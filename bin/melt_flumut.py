#!/usr/bin/env python3
import os
import pandas as pd
from pathlib import Path
import argparse 

def melt_flumut(df):
	melt_df = df.melt(id_vars='Sample', var_name='temp', value_name='allele')
	cols = melt_df['temp'].str.split(":", expand=True)
	cols = cols.rename({0:'gene', 1:'mutation'}, axis=1)

	melt_df = melt_df.drop('temp', axis=1).rename({'Sample':'sample'},axis=1)

	melt_df.insert(1, 'mutation', cols['mutation'])
	melt_df.insert(1, 'gene', cols['gene'])

	return melt_df

def parse_args():
	parser = argparse.ArgumentParser(description="Melt flumut output")
	parser.add_argument("-i", "--input", required=True, type=Path, help="Path to input flumut file")
	parser.add_argument("-o", "--output", required=True, type=Path, help="Path to output melted file")
	return parser.parse_args()

if __name__ == "__main__":
	args = parse_args()

	df = pd.read_csv(args.input, sep="\t")
	melted_df = melt_flumut(df)
	melted_df.to_csv(args.output, sep="\t", index=False)