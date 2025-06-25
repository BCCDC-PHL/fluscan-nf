#!/usr/bin/env python3
import os , sys
import pandas as pd 
import argparse
from typing import Tuple, List
from functools import partial
import gzip 
from subprocess import run, CalledProcessError
from tools import read_vcf, write_vcf

SEGMENTS = 'PB2 PB1 PA HA NP NA M NS'.split()
SEGMENTS = '13 14 15 16 17 18 19 20'.split()
def rename_vcf(df):
    df['CHROM'] = df['CHROM'].str.split("|").str[-1]
    return df

def run_snpeff(input_path, genome, output_path, config_path):
    try:
        command = ['snpEff', '-config', config_path, genome, input_path]
        process = run(command, capture_output=True, check=True, timeout=600)

    except FileNotFoundError as e:
        raise FileNotFoundError(f"File not found error: {str(e)}")
    
    except CalledProcessError as e:
        raise CalledProcessError(e.returncode, e.cmd,
            f"snpEff failed with return code {e.returncode}. Error: {e.stderr.decode('utf-8')}"
        )
    
    except OSError as e:
        raise OSError(f"Operating system error: {str(e)}")
    
    except UnicodeDecodeError as e:
        raise UnicodeDecodeError(f"Error decoding snpEff output: {str(e)}")
    
    except Exception as e:
        raise Exception(f"Unexpected error running snpEff: {str(e)}")
    
    with open(output_path,'w') as outfile:
	    outfile.write(process.stdout.decode('utf-8'))

def concat_vcf(vcf_path_dict):
    df_list = []

    for segment in SEGMENTS:
        if segment in vcf_path_dict:
            df, _ = read_vcf(vcf_path_dict[segment])
            df_list.append(df)

    return pd.concat(df_list)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', required=True, help="Input VCF")
    parser.add_argument('-c','--config', required=True, help="Path to the snpEff.config for the custom database")
    parser.add_argument('-o','--outname', required=True, help="Output VCF")
    return parser.parse_args()


def main():
    args = get_args()

    vcf_df, header_lines = read_vcf(args.input)

    annotated_vcfs = {}

    # split the input VCF by flu segment
    for ref_name, seg_df in vcf_df.groupby('CHROM'):
        print(ref_name)

        # parse out segment and subtype information from the header 
        #ref_name_fields = ref_name.split("|")
        segment = ref_name[-2:]
        subtype = 'h5n1'

        # rename the #CHROM column entries because these have to match with the snpEff database entry
        rename_vcf_df = rename_vcf(seg_df)
        
        # write a temp VCF file containing only the entries for one segment 
        write_vcf(rename_vcf_df, header_lines, '.tmp.vcf')

        # run snpEff on the temp VCF file to produce an annotated VCF for one segment only 
        output_vcf = f'{args.outname}_{segment}_anno.vcf'
        run_snpeff('.tmp.vcf', f'influenza_{subtype}', output_vcf, args.config)
        
        # save the output annotated VCF filepath for later concatentation 
        annotated_vcfs[segment] = output_vcf
    
    # extract header lines from last temp VCF file 
    _ , ann_header_lines = read_vcf('.tmp.vcf')
    
    print(annotated_vcfs.keys())
    # concatenate all the segment-based annotated VCFs together 
    final_vcf_df = concat_vcf(annotated_vcfs)

    # write the final concatenated VCF 
    final_vcf_path = f"{args.outname}_all_anno.vcf"
    print(f'Writing combined VCF to {final_vcf_path}')
    write_vcf(final_vcf_df, ann_header_lines, final_vcf_path)


if __name__ == '__main__':
    main()    
