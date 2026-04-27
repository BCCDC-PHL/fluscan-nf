#!/usr/bin/env python3

import argparse
import logging

import numpy as np
import pandas as pd

from tools import read_vcf

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

SNPEFF_ANN_FIELDS = [
    "allele",
    "annotation",
    "putative_impact",
    "gene_name",
    "gene_id",
    "feature_type",
    "feature_id",
    "transcript_biotype",
    "rank_total",
    "hgvs_c",
    "hgvs_p",
    "cdna_position_cdna_len",
    "cds_position_cds_len",
    "protein_position_protein_len",
    "distance_to_feature",
    "errors_warnings_info"
]

AA_CODES = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Asx': 'B',
    'Cys': 'C',
    'Glu': 'E',
    'Gln': 'Q',
    'Glx': 'Z',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V',
    'del': 'del',
    'Ter': '*'
}

INFO_FIELDS = "DP AF SB DP4 ANN".split()
DROP_COLS = "ID allele FILTER ANN INFO unknown transcript_biotype rank_total distance_to_feature errors_warnings_info".split()


def convert_aa_three_one(aa_three_str):
    """Convert a three-letter amino acid string to one-letter codes."""
    aa_three_str = str(aa_three_str)
    if aa_three_str in ('', 'nan'):
        return ''
    if aa_three_str == 'fs':
        return 'fs'
    if len(aa_three_str) % 3 != 0:
        log.warning('AA notation not a multiple of 3: %s', aa_three_str)
        return None

    aa_one_str = ''
    for i in range(0, len(aa_three_str), 3):
        aa = aa_three_str[i:i + 3]
        if aa not in AA_CODES:
            log.warning('AA notation %s not found within string %s', aa, aa_three_str)
            return None
        aa_one_str += AA_CODES[aa]
    return aa_one_str


def _parse_info_field(info_str):
    """Parse a VCF INFO string into a dict, safely skipping flag fields."""
    result = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            result[key] = value
    return result


def extract_info_df(df):
    """Extract selected INFO fields from a VCF DataFrame."""
    info_df = pd.DataFrame(list(df['INFO'].apply(_parse_info_field)))
    return info_df[INFO_FIELDS]


def extract_ann_df(info_df):
    """Extract and parse snpEff ANN sub-fields from the INFO DataFrame."""
    # Explode multi-allele ANN entries and split on pipe
    ann_df = info_df['ANN'].str.split(",").explode().str.split(r"\|", expand=True)
    ann_df.columns = SNPEFF_ANN_FIELDS

    aa_muts = ann_df['hgvs_p'].str.extract(r'p\.([^\.]+)').squeeze()

    # Mask for multi-base deletions (e.g. Ala100_Gly105del)
    mask_del_multi = aa_muts.str.match(r"[A-Za-z\d]+_[A-Za-z\d]+del$", na=False)

    # Split amino acid mutations into ref / pos / alt columns
    aa_muts_split = aa_muts.str.extract(r'^([A-z]+)(\d+)([A-z]+)$').fillna("")
    aa_muts_split[0] = aa_muts_split[0].apply(convert_aa_three_one)
    aa_muts_split[2] = aa_muts_split[2].apply(convert_aa_three_one)

    # Reformat multi-base deletions into "startPos-endPos" / "del"
    aa_muts_del = aa_muts.str.extract(r'[A-Za-z]+(\d+)_[A-Za-z]+(\d+)(del)$').fillna('')
    aa_muts_del[1] = aa_muts_del[0] + "-" + aa_muts_del[1]
    aa_muts_del[0] = ''

    aa_muts_split.loc[mask_del_multi] = aa_muts_del[mask_del_multi]
    aa_muts_split.columns = ['REF_AA', 'POS_AA', 'ALT_AA']

    return pd.concat([ann_df, aa_muts_split], axis=1)


def split_aa_muts(df):
    """Expand compound amino-acid mutations into semicolon-separated lists."""
    def _split_mutations(ref, pos, alt):
        if pos == '':
            return np.nan
        ref, alt = str(ref), str(alt)

        if alt == 'del':
            return pos + alt
        if alt == 'fs':
            return ref + pos + alt

        results = [ref[n] + str(int(pos) + n) + alt[n] for n in range(len(ref))]
        return ";".join(results)

    split_aa_df = (
        df[['REF_AA', 'POS_AA', 'ALT_AA']]
        .apply(lambda x: _split_mutations(x['REF_AA'], x['POS_AA'], x['ALT_AA']), axis=1)
        .to_frame(name='MUT_LIST_AA')
    )
    return split_aa_df


def reformat_vcf(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract snpEff annotations from INFO field and create separate columns.

    Args:
        df: VCF DataFrame with INFO column containing ANN field.

    Returns:
        DataFrame with parsed annotations.
    """
    info_df = extract_info_df(df)
    ann_df = extract_ann_df(info_df)

    # Merge annotation columns back into the original VCF rows
    info_df_merge = info_df.merge(ann_df, how='right', left_index=True, right_index=True)
    df_merge = df.merge(info_df_merge, how='right', left_index=True, right_index=True)

    split_aa_df = split_aa_muts(df_merge)
    df_nt_final = pd.concat([df_merge, split_aa_df], axis=1)

    df_nt_final['gene_name'] = df_nt_final['gene_name'].str.upper()

    # Move amino-acid columns to a consistent position
    aa_cols = ['REF_AA', 'POS_AA', 'ALT_AA', 'MUT_LIST_AA']
    for col_name in reversed(aa_cols):
        df_nt_final.insert(6, col_name, df_nt_final.pop(col_name))

    return df_nt_final


def filter_vcf(df):
    """Drop helper columns and rows that are not protein-coding mutations."""
    df = df.reset_index(drop=True)
    df = df.drop(columns=DROP_COLS)
    df = df[~df['gene_name'].str.contains("CHR")]
    df = df[df['hgvs_p'] != '']
    return df


def validate_outpath(outpath):
    if not outpath.endswith('.tsv') and not outpath.endswith('.csv'):
        raise ValueError("ERROR: Output path must end with either .tsv or .csv")
    return outpath


def get_args():
    parser = argparse.ArgumentParser(description="Reformat snpEff-annotated VCF into a mutation summary table.")
    parser.add_argument("-i", "--input", required=True, help="Path to the annotated VCF file")
    parser.add_argument("-o", "--output", required=True, type=validate_outpath, help="Path to the output TSV/CSV file")
    return parser.parse_args()


def main():
    args = get_args()

    log.info("Reading VCF ...")
    vcf_df, _ = read_vcf(args.input)

    log.info("Reformatting VCF ...")
    vcf_df = reformat_vcf(vcf_df)

    log.info("Filtering for relevant rows/columns ...")
    vcf_df_final = filter_vcf(vcf_df)

    vcf_df_final.to_csv(args.output, sep='\t', index=False)
    log.info("Complete.")


if __name__ == '__main__':
    main()
