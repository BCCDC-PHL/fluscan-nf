#!/usr/bin/env python3
import argparse
import logging
import os
import tempfile
import re
import pandas as pd
from subprocess import run

from lib.tools import read_vcf, write_vcf

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)


def rename_chrom(df):
    """Strip pipe-delimited prefix from CHROM, keeping only the accession."""
    df = df.copy()
    df['CHROM'] = df['CHROM'].str.split("|").str[0]
    return df


def filter_contig_headers(header_lines, accession):
    """Return header_lines with only the ##contig line matching *accession*."""
    filtered = []
    for line in header_lines:
        if line.startswith("##contig"):
            if accession in line:
                line = re.sub(r'ID=[^,]+', f"ID={accession}", line)
                filtered.append(line)
        else:
            filtered.append(line)
    return filtered


def run_snpeff(input_path, genome, output_path, config_path):
    """Run snpEff annotation and write the result to *output_path*."""
    command = ['snpEff', '-hgvs1LetterAa', '-config', config_path, genome, input_path]
    process = run(command, capture_output=True, check=True, timeout=600)

    with open(output_path, 'w') as outfile:
        outfile.write(process.stdout.decode('utf-8'))


def get_args():
    parser = argparse.ArgumentParser(description="Run snpEff on a multi-segment influenza VCF.")
    parser.add_argument('-i', '--input', required=True, help="Input VCF")
    parser.add_argument('-c', '--config', required=True, help="Path to the snpEff.config for the custom database")
    parser.add_argument('-o', '--outname', required=True, help="Output VCF name prefix")
    return parser.parse_args()


def main():
    args = get_args()

    vcf_df, header_lines = read_vcf(args.input)

    annotated_dfs = []
    ann_header_lines = None

    # Split the input VCF by flu segment, annotate each independently
    for ref_name, seg_df in vcf_df.groupby('CHROM'):
        log.info("Processing segment: %s", ref_name)

        # Parse out segment and subtype information from the CHROM field
        ref_name_fields = ref_name.split("|")
        segment = ref_name_fields[-1]
        subtype = ref_name_fields[-2]

        # Rename CHROM entries to match the snpEff database entry names
        renamed_df = rename_chrom(seg_df)
        accession = renamed_df['CHROM'].iloc[0]
        seg_header_lines = filter_contig_headers(header_lines, accession)

        # Use a temporary file for the single-segment VCF
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_vcf = os.path.join(tmpdir, 'segment.vcf')
            output_vcf = f'{args.outname}__{accession}__{subtype}__{segment}.vcf'

            write_vcf(renamed_df, seg_header_lines, tmp_vcf)
            run_snpeff(tmp_vcf, f'influenza_{subtype.lower()}', output_vcf, args.config)

        # Read back the annotated VCF and accumulate
        ann_df, ann_header_lines = read_vcf(output_vcf)
        annotated_dfs.append(ann_df)

    if not annotated_dfs:
        log.warning("No segments found in input VCF; nothing to concatenate.")
        return

    # Concatenate all segment DataFrames in canonical segment order
    final_vcf_df = pd.concat(annotated_dfs, ignore_index=True)

    # Write the final concatenated VCF
    final_vcf_path = f"{args.outname}_all_anno.vcf"
    log.info("Writing combined VCF to %s", final_vcf_path)
    write_vcf(final_vcf_df, ann_header_lines, final_vcf_path)


if __name__ == '__main__':
    main()    