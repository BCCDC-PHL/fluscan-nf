#%%
import numpy as np
import pandas as pd
from typing import Tuple, List
import gzip

BASE_DICT = dict(zip(list("ACGT"), range(4)))


SEGMENTS = 'PB2 PB1 PA HA NP NA M NS'.split()

def read_vcf(file_path: str, max_header_lines: int = 100_000) -> Tuple[pd.DataFrame, List[str]]:
    """
    Parse a Variant Call Format (VCF) file into a pandas DataFrame and header list.

    Transparently handles gzip-compressed (.gz) input. Every metadata ('##...') line and
    the '#CHROM' column line are returned verbatim in ``header_lines``. Data columns are
    read as strings with NA-detection disabled so values round-trip byte-for-byte: VCF '.'
    and 'NA'-like tokens and empty fields are preserved (not coerced to NaN), and POS stays
    a string — cast it (``.astype(int)``) when you need integers.

    Args:
        file_path: Path to the VCF file (optionally .gz).
        max_header_lines: Safety cap on the header scan before giving up.

    Returns:
        Tuple of (DataFrame of variant rows, list of header lines including '#CHROM').

    Raises:
        ValueError: if no '#CHROM' header line is found (missing/truncated/empty file).
    """
    opener = gzip.open if file_path.endswith('.gz') else open
    header_lines: List[str] = []

    with opener(file_path, 'rt') as infile:
        # Use readline() rather than iterating the handle: iteration uses a read-ahead
        # buffer that would consume part of the body, leaving nothing for pandas to read.
        line = infile.readline()
        while True:
            if line == '':
                raise ValueError(f"Reached end of {file_path} without a '#CHROM' header line")
            stripped = line.rstrip('\r\n')
            header_lines.append(stripped)
            if stripped.lstrip('#').startswith("CHROM"):
                break
            if len(header_lines) >= max_header_lines:
                raise ValueError(
                    f"No '#CHROM' header found within {max_header_lines} lines of {file_path}"
                )
            line = infile.readline()

        column_names = header_lines[-1].lstrip("#").split("\t")
        df = pd.read_csv(infile, sep='\t', names=column_names, dtype=str, na_filter=False)

    return df, header_lines


def write_vcf(df: pd.DataFrame, header_lines: List[str], output_path: str) -> None:
    """
    Write a VCF file from a DataFrame of variant rows and a list of header lines.

    The inverse of :func:`read_vcf`: emits each header line, then the DataFrame rows as
    tab-separated values (no index, no pandas header). Gzip-compresses the output when
    ``output_path`` ends in .gz.

    Args:
        df: Variant rows whose columns are already in VCF order.
        header_lines: Header lines to write verbatim (including the '#CHROM' line).
        output_path: Destination path (optionally .gz).
    """
    opener = gzip.open if str(output_path).endswith('.gz') else open
    with opener(output_path, 'wt') as outfile:
        outfile.write('\n'.join(header_lines) + '\n')
        df.to_csv(outfile, sep='\t', index=False, header=False)