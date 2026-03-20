#!/usr/bin/env python3

import argparse
import logging
import math
import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    from cyvcf2 import VCF
except ImportError as exc:
    raise SystemExit(
        "ERROR: cyvcf2 is required for this script. "
        "Install it in the active environment before running."
    ) from exc


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)


SNPEFF_ANN_FIELDS = [
    "allele",
    "annotation",
    "impact",
    "gene_name",
    "gene_id",
    "feature_type",
    "feature_id",
    "transcript_biotype",
    "rank",
    "hgvs_c",
    "hgvs_p",
    "cdna_position_cdna_len",
    "cds_position_cds_len",
    "protein_position_protein_len",
    "distance_to_feature",
    "ann_notes",
]


AA_CODES = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Asx": "B",
    "Cys": "C",
    "Glu": "E",
    "Gln": "Q",
    "Glx": "Z",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
    "del": "del",
    "fs": "fs",
}


GT_TYPE_LABELS = {
    0: "HOM_REF",
    1: "HET",
    2: "HOM_ALT",
    3: "UNKNOWN",
}


ALLELE_INFO_FIELDS = [
    "AF",
    "AO",
    "QA",
    "PAO",
    "PQA",
    "SAF",
    "SAR",
    "RPL",
    "RPR",
    "EPP",
    "SAP",
    "RPP",
    "MQM",
    "PAIRED",
    "MEANALT",
    "LEN",
    "TYPE",
    "CIGAR",
    "AB",
    "ABP",
    "RUN",
    "DPRA",
]


OUTPUT_COLUMNS = [
    "chrom",
    "pos",
    "id",
    "ref",
    "alt",
    "qual",
    "filter",
    "var_type",
    "is_snp",
    "is_indel",
    "is_mnp",
    "alt_index",
    "call_rate",
    "num_called",
    "num_het",
    "num_hom_ref",
    "num_hom_alt",
    "num_unknown",
    "sample_count",
    "samples",
    "sample_gt_bases",
    "sample_gt_types",
    "sample_gt_depths",
    "sample_gt_quals",
    "sample_ref_depths",
    "sample_alt_depths",
    "sample_alt_freqs",
    "dp",
    "ro",
    "qr",
    "srf",
    "srr",
    "mqmr",
    "pairedr",
    "af",
    "ao",
    "qa",
    "pao",
    "pqa",
    "saf",
    "sar",
    "rpl",
    "rpr",
    "epp",
    "sap",
    "rpp",
    "mqm",
    "paired",
    "meanalt",
    "allele_len",
    "caller_type",
    "cigar",
    "ab",
    "abp",
    "run",
    "dpra",
    "annotation",
    "impact",
    "gene_name",
    "gene_id",
    "feature_type",
    "feature_id",
    "transcript_biotype",
    "rank",
    "hgvs_c",
    "hgvs_p",
    "cdna_pos",
    "cdna_len",
    "cds_pos",
    "cds_len",
    "protein_pos",
    "protein_len",
    "aa_ref",
    "aa_pos",
    "aa_alt",
    "aa_change",
    "distance_to_feature",
    "ann_notes",
    "lof",
    "nmd",
]


def validate_outpath(outpath: str) -> str:
    if not outpath.endswith(".tsv"):
        raise ValueError("ERROR: Output path must end with .tsv")
    return outpath


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize a SnpEff-annotated VCF into a flattened TSV using cyvcf2 "
            "for VCF parsing and pandas for tabular output."
        )
    )
    parser.add_argument("-i", "--input", required=True, help="Path to the input VCF/BCF")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=validate_outpath,
        help="Path to the output TSV",
    )
    return parser.parse_args()


def normalize_value(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, bytes):
        value = value.decode()
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and math.isnan(value):
        return None
    if value in ("", "."):
        return None
    return value


def as_list(value: Any) -> List[Any]:
    value = normalize_value(value)
    if value is None:
        return []
    if isinstance(value, np.ndarray):
        return [normalize_value(x) for x in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [normalize_value(x) for x in value]
    if isinstance(value, str) and "," in value:
        return [normalize_value(x) for x in value.split(",")]
    return [value]


def collapse_values(values: Optional[Sequence[Any]]) -> Optional[str]:
    if values is None:
        return None
    if isinstance(values, np.ndarray):
        values = values.tolist()
    if not isinstance(values, (list, tuple)):
        value = normalize_value(values)
        return None if value is None else str(value)

    cleaned = [normalize_value(v) for v in values]
    cleaned = [v for v in cleaned if v is not None]
    if not cleaned:
        return None
    if len(cleaned) == 1:
        return str(cleaned[0])
    return "|".join(str(v) for v in cleaned)


def collapse_gt_types(values: Optional[Sequence[Any]]) -> Optional[str]:
    values = as_list(values)
    if not values:
        return None
    labels = [GT_TYPE_LABELS.get(int(value), str(value)) for value in values if value is not None]
    if not labels:
        return None
    if len(labels) == 1:
        return labels[0]
    return "|".join(labels)


def allele_info_value(variant, key: str, alt_index: int) -> Any:
    values = as_list(variant.INFO.get(key))
    if not values:
        return None
    if len(values) == 1:
        return values[0]
    if alt_index < len(values):
        return values[alt_index]
    return None


def split_position_length(value: Any) -> Tuple[Any, Any]:
    value = normalize_value(value)
    if value is None:
        return None, None
    parts = str(value).split("/", maxsplit=1)
    if len(parts) == 1:
        return normalize_value(parts[0]), None
    return normalize_value(parts[0]), normalize_value(parts[1])


def convert_aa_three_one(aa_three_str: str) -> Optional[str]:
    aa_three_str = str(aa_three_str)
    if aa_three_str in ("", "nan", "None"):
        return None
    if aa_three_str == "fs":
        return "fs"
    if len(aa_three_str) % 3 != 0:
        return None

    aa_one_str = ""
    for i in range(0, len(aa_three_str), 3):
        aa = aa_three_str[i : i + 3]
        if aa not in AA_CODES:
            return None
        aa_one_str += AA_CODES[aa]
    return aa_one_str


def parse_hgvs_p(hgvs_p: Any) -> Dict[str, Any]:
    hgvs_p = normalize_value(hgvs_p)
    parsed = {
        "aa_ref": None,
        "aa_pos": None,
        "aa_alt": None,
        "aa_change": None,
    }
    if hgvs_p is None:
        return parsed

    protein_change = str(hgvs_p).removeprefix("p.")
    protein_change = protein_change or None
    parsed["aa_change"] = protein_change

    if protein_change is None:
        return parsed

    simple_sub = re.fullmatch(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|del|fs)", protein_change)
    if simple_sub:
        ref_aa, aa_pos, alt_aa = simple_sub.groups()
        parsed["aa_ref"] = convert_aa_three_one(ref_aa)
        parsed["aa_pos"] = aa_pos
        parsed["aa_alt"] = convert_aa_three_one(alt_aa)
        if parsed["aa_ref"] and parsed["aa_pos"] and parsed["aa_alt"]:
            parsed["aa_change"] = f"{parsed['aa_ref']}{parsed['aa_pos']}{parsed['aa_alt']}"
        return parsed

    multi_del = re.fullmatch(
        r"([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)del",
        protein_change,
    )
    if multi_del:
        _, start_pos, _, end_pos = multi_del.groups()
        parsed["aa_pos"] = f"{start_pos}-{end_pos}"
        parsed["aa_alt"] = "del"
        return parsed

    pos_only = re.fullmatch(r"(\d+)", protein_change)
    if pos_only:
        parsed["aa_pos"] = pos_only.group(1)
        return parsed

    return parsed


def parse_ann_field(raw_ann: Any) -> List[Dict[str, Any]]:
    raw_ann = normalize_value(raw_ann)
    if raw_ann is None:
        return []

    annotations = []
    for entry in str(raw_ann).split(","):
        parts = entry.split("|")
        if len(parts) < len(SNPEFF_ANN_FIELDS):
            parts.extend([""] * (len(SNPEFF_ANN_FIELDS) - len(parts)))
        ann = dict(zip(SNPEFF_ANN_FIELDS, parts[: len(SNPEFF_ANN_FIELDS)]))
        ann = {key: normalize_value(value) for key, value in ann.items()}

        ann["gene_name"] = (
            str(ann["gene_name"]).upper() if ann.get("gene_name") is not None else None
        )
        ann["cdna_pos"], ann["cdna_len"] = split_position_length(ann.pop("cdna_position_cdna_len"))
        ann["cds_pos"], ann["cds_len"] = split_position_length(ann.pop("cds_position_cds_len"))
        ann["protein_pos"], ann["protein_len"] = split_position_length(
            ann.pop("protein_position_protein_len")
        )
        ann.update(parse_hgvs_p(ann.get("hgvs_p")))
        annotations.append(ann)
    return annotations


def normalize_effect_blob(value: Any) -> Optional[str]:
    values = as_list(value)
    if not values:
        return None

    cleaned = []
    for item in values:
        item = normalize_value(item)
        if item is None:
            continue
        item = str(item)
        if item.startswith("(") and item.endswith(")"):
            item = item[1:-1]
        cleaned.append(item)
    return None if not cleaned else ";".join(cleaned)


def annotations_for_alt(
    annotations: List[Dict[str, Any]],
    alt: str,
    total_alts: int,
) -> List[Dict[str, Any]]:
    matching = [ann for ann in annotations if ann.get("allele") == alt]
    if matching:
        return matching
    if total_alts == 1:
        return annotations
    return []


def summarize_vcf(vcf_path: str) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    vcf = VCF(vcf_path, gts012=True, strict_gt=True)

    try:
        samples = list(vcf.samples)

        for variant in vcf:
            alts = [alt for alt in variant.ALT if alt is not None]
            annotations = parse_ann_field(variant.INFO.get("ANN"))

            base_row = {
                "chrom": variant.CHROM,
                "pos": variant.POS,
                "id": normalize_value(variant.ID),
                "ref": variant.REF,
                "qual": normalize_value(variant.QUAL),
                "filter": "PASS" if variant.FILTER is None else collapse_values(variant.FILTERS),
                "var_type": normalize_value(variant.var_type),
                "is_snp": bool(variant.is_snp),
                "is_indel": bool(variant.is_indel),
                "is_mnp": bool(getattr(variant, "is_mnp", False)),
                "call_rate": normalize_value(variant.call_rate),
                "num_called": normalize_value(variant.num_called),
                "num_het": normalize_value(variant.num_het),
                "num_hom_ref": normalize_value(variant.num_hom_ref),
                "num_hom_alt": normalize_value(variant.num_hom_alt),
                "num_unknown": normalize_value(variant.num_unknown),
                "sample_count": len(samples),
                "samples": collapse_values(samples),
                "sample_gt_bases": collapse_values(getattr(variant, "gt_bases", None)),
                "sample_gt_types": collapse_gt_types(getattr(variant, "gt_types", None)),
                "sample_gt_depths": collapse_values(getattr(variant, "gt_depths", None)),
                "sample_gt_quals": collapse_values(getattr(variant, "gt_quals", None)),
                "sample_ref_depths": collapse_values(getattr(variant, "gt_ref_depths", None)),
                "sample_alt_depths": collapse_values(getattr(variant, "gt_alt_depths", None)),
                "sample_alt_freqs": collapse_values(getattr(variant, "gt_alt_freqs", None)),
                "dp": normalize_value(variant.INFO.get("DP")),
                "ro": normalize_value(variant.INFO.get("RO")),
                "qr": normalize_value(variant.INFO.get("QR")),
                "srf": normalize_value(variant.INFO.get("SRF")),
                "srr": normalize_value(variant.INFO.get("SRR")),
                "mqmr": normalize_value(variant.INFO.get("MQMR")),
                "pairedr": normalize_value(variant.INFO.get("PAIREDR")),
                "lof": normalize_effect_blob(variant.INFO.get("LOF")),
                "nmd": normalize_effect_blob(variant.INFO.get("NMD")),
            }

            for alt_index, alt in enumerate(alts):
                alt_row = dict(base_row)
                alt_row["alt"] = alt
                alt_row["alt_index"] = alt_index + 1

                for field in ALLELE_INFO_FIELDS:
                    out_key = field.lower()
                    if out_key == "len":
                        out_key = "allele_len"
                    elif out_key == "type":
                        out_key = "caller_type"
                    alt_row[out_key] = allele_info_value(variant, field, alt_index)

                alt_annotations = annotations_for_alt(annotations, alt, len(alts))
                if not alt_annotations:
                    rows.append(alt_row)
                    continue

                for ann in alt_annotations:
                    row = dict(alt_row)
                    row.update(ann)
                    rows.append(row)

        df = pd.DataFrame(rows)
    finally:
        vcf.close()

    for column in OUTPUT_COLUMNS:
        if column not in df.columns:
            df[column] = None

    df = df[OUTPUT_COLUMNS]
    if not df.empty:
        df = df.sort_values(["chrom", "pos", "alt_index", "gene_name", "annotation"], na_position="last")
    return df


def main() -> None:
    args = get_args()
    log.info("Reading and summarizing VCF ...")
    summary_df = summarize_vcf(args.input)

    log.info("Writing %s rows to %s", len(summary_df), args.output)
    summary_df.to_csv(args.output, sep="\t", index=False)
    log.info("Complete.")


if __name__ == "__main__":
    main()
