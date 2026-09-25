#!/usr/bin/env python3
import os
import pandas as pd
from pathlib import Path
import argparse 

def melt_flumut(df):
	# A sample with no FluMut hits has only the 'Sample' column; this yields a header-only table.
	melt_df = df.melt(id_vars='Sample', var_name='temp', value_name='allele').rename({'Sample': 'sample'}, axis=1)
	parts = melt_df['temp'].str.split(':')
	melt_df['gene'] = parts.str[0]
	melt_df['mutation'] = parts.str[1]
	return melt_df[['sample', 'gene', 'mutation', 'allele']]

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