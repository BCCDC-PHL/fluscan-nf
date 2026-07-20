#!/usr/bin/env python3
"""Translation-based sanity check of Nextclade vs paftools AA calls.

This lifts the latest session prototype into a reusable CLI that reads a
Nextclade dataset directory directly rather than relying on hard-coded session
paths or contig names.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    import cyvcf2
except ImportError:
    cyvcf2 = None

try:
    import pandas as pd
except ImportError:
    pd = None

try:
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner, substitution_matrices
except ImportError:
    SeqIO = None
    PairwiseAligner = Any
    substitution_matrices = None

try:
    import plotly.express as px
    import plotly.graph_objects as go
except ImportError:  # pragma: no cover - optional dependency in this repo
    px = None
    go = None


AA_EFFECT_TERMS = {
    "missense_variant", "synonymous_variant", "stop_gained",
    "stop_lost", "start_lost", "stop_retained_variant",
    "initiator_codon_variant",
    "conservative_inframe_insertion", "conservative_inframe_deletion",
    "disruptive_inframe_insertion", "disruptive_inframe_deletion",
    "frameshift_variant", "protein_altering_variant",
}

THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*", "Xaa": "X", "Xle": "J",
    "Asx": "B", "Glx": "Z",
}

CLASS_MAP = {
    "missense_variant":               "missense",
    "synonymous_variant":             "synonymous",
    "stop_retained_variant":          "synonymous",
    "stop_gained":                    "stop_gained",
    "stop_lost":                      "stop_lost",
    "start_lost":                     "start_lost",
    "initiator_codon_variant":        "start_lost",
    "conservative_inframe_insertion": "inframe_indel",
    "conservative_inframe_deletion":  "inframe_indel",
    "disruptive_inframe_insertion":   "inframe_indel",
    "disruptive_inframe_deletion":    "inframe_indel",
    "frameshift_variant":             "frameshift",
    "protein_altering_variant":       "frameshift",
}

SNV_CLASSES = {"missense", "synonymous", "stop_gained", "stop_lost", "start_lost"}
ALL_CLASSES = set(CLASS_MAP.values())


@dataclass(frozen=True)
class SegmentInfo:
    label: str
    directory: str
    nextclade_contig: str
    fasta_tag: str


@dataclass(frozen=True)
class AAKey:
    gene: str
    hgvs_p: str
    klass: str


@dataclass
class CDSDef:
    gene: str
    exons: list[tuple[int, int]]


@dataclass
class PafHeader:
    status: str
    fraction: float
    range_start: int
    range_end: int


@dataclass
class BioCall:
    gene: str
    hgvs_p: str
    klass: str
    ref_aa: str
    alt_aa: str
    aa_pos: int
    codon_ref_nt: str
    codon_alt_nt: str
    codon_has_n: bool
    skipped: bool


@dataclass
class RunConfig:
    consensus_dir: Path
    reference_datasets_dir: Path
    nextclade_vcf_dir: Path
    paftools_vcf_dir: Path
    out_dir: Path
    sample_filter: set[str]
    skip_plots: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-datasets-dir",
        required=True,
        type=Path,
        help="Directory containing per-segment Nextclade datasets and config.json.",
    )
    parser.add_argument(
        "--consensus-dir",
        required=True,
        type=Path,
        help="Directory of *.consensus.fasta multi-FASTA files.",
    )
    parser.add_argument(
        "--nextclade-vcf-dir",
        required=True,
        type=Path,
        help="Directory containing *_nc_final.vcf files.",
    )
    parser.add_argument(
        "--paftools-vcf-dir",
        required=True,
        type=Path,
        help="Directory containing *_paf_final.vcf files.",
    )
    parser.add_argument(
        "--out-dir",
        default=Path("biopython_sanity"),
        type=Path,
        help="Output directory for CSV, report, and optional plots.",
    )
    parser.add_argument(
        "--sample",
        action="append",
        default=[],
        help="Limit analysis to one or more sample IDs. May be provided multiple times.",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip HTML plot generation even if plotly is installed.",
    )
    return parser.parse_args()


def build_run_config(args: argparse.Namespace) -> RunConfig:
    sample_filter = {sample.strip() for sample in args.sample if sample.strip()}
    return RunConfig(
        consensus_dir=args.consensus_dir.resolve(),
        reference_datasets_dir=args.reference_datasets_dir.resolve(),
        nextclade_vcf_dir=args.nextclade_vcf_dir.resolve(),
        paftools_vcf_dir=args.paftools_vcf_dir.resolve(),
        out_dir=args.out_dir.resolve(),
        sample_filter=sample_filter,
        skip_plots=args.skip_plots,
    )


def _parse_attrs(attrs_str: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for part in attrs_str.split(";"):
        if "=" in part:
            key, _, value = part.partition("=")
            result[key.strip()] = value.strip()
    return result


def load_segment_info(dataset_root: Path) -> list[SegmentInfo]:
    config_path = dataset_root / "config.json"
    if not config_path.exists():
        raise SystemExit(f"config.json not found in {dataset_root}")

    with open(config_path) as handle:
        config = json.load(handle)

    segments: list[SegmentInfo] = []
    for segment_cfg in config.get("segments", []):
        seg_dir = dataset_root / segment_cfg["dir"]
        pathogen_path = seg_dir / "pathogen.json"
        if not pathogen_path.exists():
            raise SystemExit(f"pathogen.json not found for segment {segment_cfg['name']}: {pathogen_path}")
        with open(pathogen_path) as handle:
            pathogen = json.load(handle)

        nextclade_contig = pathogen["attributes"]["reference accession"]
        label = segment_cfg["name"].upper()
        fasta_tag = f"segment{segment_cfg['segment_num']}_{label}"
        segments.append(
            SegmentInfo(
                label=label,
                directory=segment_cfg["dir"],
                nextclade_contig=nextclade_contig,
                fasta_tag=fasta_tag,
            )
        )

    if not segments:
        raise SystemExit(f"No segment definitions found in {config_path}")
    return segments


def parse_gff3(dataset_root: Path, segment: SegmentInfo) -> list[CDSDef]:
    gff_path = dataset_root / segment.directory / "genome_annotation.gff3"

    exons_by_protein: dict[str, list[tuple[int, int]]] = defaultdict(list)
    gene_by_protein: dict[str, str] = {}

    with open(gff_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            start = int(parts[3])
            end = int(parts[4])
            attrs = _parse_attrs(parts[8])
            protein_id = attrs.get("Name", "")
            gene_name = attrs.get("gene", "")
            if protein_id:
                exons_by_protein[protein_id].append((start, end))
                if gene_name:
                    gene_by_protein[protein_id] = gene_name

    cds_defs = [
        CDSDef(gene=gene_by_protein.get(protein_id, protein_id), exons=sorted(exons))
        for protein_id, exons in exons_by_protein.items()
    ]
    cds_defs.sort(key=lambda cds: cds.gene)
    return cds_defs


def load_reference(dataset_root: Path, segment: SegmentInfo) -> str:
    reference_path = dataset_root / segment.directory / "reference.fasta"
    record = next(SeqIO.parse(str(reference_path), "fasta"))
    return str(record.seq).upper()


def load_all_consensus(consensus_dir: Path, sample_filter: set[str]) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    fasta_paths = sorted(consensus_dir.glob("*.consensus.fasta"))
    if not fasta_paths:
        raise SystemExit(f"No *.consensus.fasta files found in {consensus_dir}")

    for fasta_path in fasta_paths:
        sample = fasta_path.stem.replace(".consensus", "")
        if sample_filter and sample not in sample_filter:
            continue
        sequences: dict[str, str] = {}
        for record in SeqIO.parse(str(fasta_path), "fasta"):
            tag = record.id[len(sample) + 1:]
            sequences[tag] = str(record.seq).upper()
        result[sample] = sequences

    if sample_filter:
        missing = sorted(sample_filter - set(result))
        if missing:
            raise SystemExit(f"Requested sample(s) not found in {consensus_dir}: {', '.join(missing)}")
    return result


def _make_dna_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    return aligner


def _make_protein_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner

def align_dna(ref_seq: str, alt_seq: str) -> tuple[str, str]:
    alignment = next(iter(dna_aligner().align(ref_seq, alt_seq)))
    return str(alignment[0]), str(alignment[1])


def align_proteins(ref_protein: str, alt_protein: str) -> tuple[str, str]:
    if not ref_protein or not alt_protein:
        return ref_protein, alt_protein
    alignment = next(iter(protein_aligner().align(ref_protein, alt_protein)))
    return str(alignment[0]), str(alignment[1])


def extract_cds_via_alignment(
    ref_aln: str,
    alt_aln: str,
    exons: list[tuple[int, int]],
) -> tuple[str, str]:
    exon_positions: set[int] = set()
    for start, end in exons:
        exon_positions.update(range(start, end + 1))

    ref_cds: list[str] = []
    alt_cds: list[str] = []
    ref_pos = 0

    for ref_char, alt_char in zip(ref_aln, alt_aln):
        if ref_char != "-":
            ref_pos += 1
            if ref_pos in exon_positions:
                ref_cds.append(ref_char)
                alt_cds.append(alt_char)
        elif ref_pos in exon_positions:
            ref_cds.append("-")
            alt_cds.append(alt_char)

    return "".join(ref_cds), "".join(alt_cds)


def _build_codon_table() -> dict[str, str]:
    from Bio.Data.CodonTable import standard_dna_table

    table: dict[str, str] = {}
    for codon, aa in standard_dna_table.forward_table.items():
        table[codon] = aa
    for stop_codon in standard_dna_table.stop_codons:
        table[stop_codon] = "*"
    return table


_CODON_TABLE: dict[str, str] = {}
_DNA_ALIGNER: PairwiseAligner | None = None
_PROTEIN_ALIGNER: PairwiseAligner | None = None


def ensure_runtime_dependencies() -> None:
    missing = []
    if cyvcf2 is None:
        missing.append("cyvcf2")
    if pd is None:
        missing.append("pandas")
    if SeqIO is None or substitution_matrices is None:
        missing.append("biopython")
    if missing:
        raise SystemExit(
            "Missing runtime dependencies: "
            + ", ".join(missing)
            + ". Activate the fluscan environment before running this script."
        )


def dna_aligner() -> PairwiseAligner:
    global _DNA_ALIGNER
    if _DNA_ALIGNER is None:
        _DNA_ALIGNER = _make_dna_aligner()
    return _DNA_ALIGNER


def protein_aligner() -> PairwiseAligner:
    global _PROTEIN_ALIGNER
    if _PROTEIN_ALIGNER is None:
        _PROTEIN_ALIGNER = _make_protein_aligner()
    return _PROTEIN_ALIGNER


def _codon_table() -> dict[str, str]:
    global _CODON_TABLE
    if not _CODON_TABLE:
        _CODON_TABLE = _build_codon_table()
    return _CODON_TABLE


def translate_cds(nt: str) -> tuple[str, list[bool]]:
    table = _codon_table()
    protein: list[str] = []
    has_n: list[bool] = []
    for index in range(0, len(nt) - 2, 3):
        codon = nt[index:index + 3].upper()
        if any(base not in "ACGT" for base in codon):
            protein.append("X")
            has_n.append(True)
        else:
            protein.append(table.get(codon, "X"))
            has_n.append(False)
    return "".join(protein), has_n


def cds_nt_to_genomic(cds_nt_pos: int, exons: list[tuple[int, int]]) -> int:
    pos = cds_nt_pos
    for start, end in exons:
        exon_len = end - start + 1
        if pos < exon_len:
            return start + pos
        pos -= exon_len
    return -1


def aa_genomic_start(aa_pos_1based: int, exons: list[tuple[int, int]]) -> int:
    return cds_nt_to_genomic((aa_pos_1based - 1) * 3, exons)


def call_aa_mutations(
    gene: str,
    ref_cds_gapped: str,
    alt_cds_gapped: str,
) -> tuple[list[BioCall], int]:
    ref_nt = ref_cds_gapped.replace("-", "")
    alt_nt = alt_cds_gapped.replace("-", "")

    ref_protein, _ = translate_cds(ref_nt)
    alt_protein, _ = translate_cds(alt_nt)

    calls: list[BioCall] = []
    n_skipped = 0

    indel_len = len(alt_nt) - len(ref_nt)
    is_frameshift = indel_len != 0 and indel_len % 3 != 0

    if is_frameshift:
        min_len = min(len(ref_protein), len(alt_protein))
        first_diff = next(
            (index for index in range(min_len) if ref_protein[index] != alt_protein[index]),
            min_len,
        )
        fs_pos = first_diff + 1
        ref_aa_at_fs = ref_protein[first_diff] if first_diff < len(ref_protein) else "?"
        downstream = alt_protein[first_diff:]
        stop_offset = downstream.find("*")
        hgvs_p = f"{ref_aa_at_fs}{fs_pos}fs*{stop_offset + 1}" if stop_offset >= 0 else f"{ref_aa_at_fs}{fs_pos}fs?"
        ref_codon = ref_nt[first_diff * 3:first_diff * 3 + 3] if (first_diff + 1) * 3 <= len(ref_nt) else "???"
        alt_codon = alt_nt[first_diff * 3:first_diff * 3 + 3] if (first_diff + 1) * 3 <= len(alt_nt) else "???"
        calls.append(BioCall(
            gene=gene,
            hgvs_p=hgvs_p,
            klass="frameshift",
            ref_aa=ref_aa_at_fs,
            alt_aa=alt_protein[first_diff] if first_diff < len(alt_protein) else "?",
            aa_pos=fs_pos,
            codon_ref_nt=ref_codon,
            codon_alt_nt=alt_codon,
            codon_has_n=False,
            skipped=False,
        ))
        return calls, n_skipped

    ref_stop_pos = ref_protein.find("*")
    if ref_stop_pos < 0:
        ref_stop_pos = len(ref_protein)

    alt_early_stop = next(
        (index for index, aa in enumerate(alt_protein) if aa == "*" and index < ref_stop_pos),
        -1,
    )

    ref_p_no_stop = ref_protein[:ref_stop_pos]
    if alt_early_stop >= 0:
        alt_p_for_aln = alt_protein[:alt_early_stop]
    else:
        alt_stop_pos = alt_protein.find("*")
        alt_p_for_aln = alt_protein[:alt_stop_pos] if alt_stop_pos >= 0 else alt_protein

    ref_p_aln, alt_p_aln = align_proteins(ref_p_no_stop, alt_p_for_aln)

    ref_idx = 0
    alt_idx = 0
    idx = 0
    while idx < len(ref_p_aln):
        ref_char = ref_p_aln[idx]
        alt_char = alt_p_aln[idx]

        if ref_char == "-" or alt_char == "-":
            gap_end = idx
            ref_gap = 0
            alt_gap = 0
            while gap_end < len(ref_p_aln) and (ref_p_aln[gap_end] == "-" or alt_p_aln[gap_end] == "-"):
                if ref_p_aln[gap_end] == "-":
                    ref_gap += 1
                if alt_p_aln[gap_end] == "-":
                    alt_gap += 1
                gap_end += 1

            if alt_gap > 0:
                start_pos = ref_idx + 1
                end_pos = ref_idx + alt_gap
                deleted = ref_p_no_stop[ref_idx:ref_idx + alt_gap]
                hgvs_p = f"{deleted[0]}{start_pos}del" if alt_gap == 1 else f"{deleted[0]}{start_pos}_{deleted[-1]}{end_pos}del"
                calls.append(BioCall(
                    gene=gene,
                    hgvs_p=hgvs_p,
                    klass="inframe_indel",
                    ref_aa=deleted,
                    alt_aa="-" * alt_gap,
                    aa_pos=start_pos,
                    codon_ref_nt="",
                    codon_alt_nt="",
                    codon_has_n=False,
                    skipped=False,
                ))

            if ref_gap > 0:
                inserted = alt_p_for_aln[alt_idx:alt_idx + ref_gap]
                anchor_pos = ref_idx
                anchor_aa = ref_p_no_stop[ref_idx - 1] if ref_idx > 0 else "?"
                next_aa = ref_p_no_stop[ref_idx] if ref_idx < len(ref_p_no_stop) else "?"
                hgvs_p = f"{anchor_aa}{anchor_pos}_{next_aa}{ref_idx + 1}ins{inserted}"
                calls.append(BioCall(
                    gene=gene,
                    hgvs_p=hgvs_p,
                    klass="inframe_indel",
                    ref_aa="",
                    alt_aa=inserted,
                    aa_pos=anchor_pos,
                    codon_ref_nt="",
                    codon_alt_nt="",
                    codon_has_n=False,
                    skipped=False,
                ))

            ref_idx += alt_gap
            alt_idx += ref_gap
            idx = gap_end
            continue

        aa_pos = ref_idx + 1
        ref_codon = ref_nt[ref_idx * 3:(ref_idx + 1) * 3] if (ref_idx + 1) * 3 <= len(ref_nt) else ""
        alt_codon = alt_nt[alt_idx * 3:(alt_idx + 1) * 3] if (alt_idx + 1) * 3 <= len(alt_nt) else ""

        if alt_codon and any(base not in "ACGT" for base in alt_codon):
            n_skipped += 1
            ref_idx += 1
            alt_idx += 1
            idx += 1
            continue

        if aa_pos == 1 and ref_char == "M" and alt_char not in ("M", "*"):
            calls.append(BioCall(
                gene=gene,
                hgvs_p=f"M1{alt_char}",
                klass="start_lost",
                ref_aa="M",
                alt_aa=alt_char,
                aa_pos=1,
                codon_ref_nt=ref_codon,
                codon_alt_nt=alt_codon,
                codon_has_n=False,
                skipped=False,
            ))
        elif ref_char != alt_char:
            calls.append(BioCall(
                gene=gene,
                hgvs_p=f"{ref_char}{aa_pos}{alt_char}",
                klass="missense",
                ref_aa=ref_char,
                alt_aa=alt_char,
                aa_pos=aa_pos,
                codon_ref_nt=ref_codon,
                codon_alt_nt=alt_codon,
                codon_has_n=False,
                skipped=False,
            ))
        elif ref_codon and alt_codon and ref_codon != alt_codon:
            calls.append(BioCall(
                gene=gene,
                hgvs_p=f"{ref_char}{aa_pos}{alt_char}",
                klass="synonymous",
                ref_aa=ref_char,
                alt_aa=alt_char,
                aa_pos=aa_pos,
                codon_ref_nt=ref_codon,
                codon_alt_nt=alt_codon,
                codon_has_n=False,
                skipped=False,
            ))

        ref_idx += 1
        alt_idx += 1
        idx += 1

    if alt_early_stop >= 0:
        pos = alt_early_stop + 1
        ref_aa_here = ref_p_no_stop[alt_early_stop] if alt_early_stop < len(ref_p_no_stop) else "?"
        ref_codon = ref_nt[alt_early_stop * 3:(alt_early_stop + 1) * 3]
        alt_codon = alt_nt[alt_early_stop * 3:(alt_early_stop + 1) * 3]
        calls.append(BioCall(
            gene=gene,
            hgvs_p=f"{ref_aa_here}{pos}*",
            klass="stop_gained",
            ref_aa=ref_aa_here,
            alt_aa="*",
            aa_pos=pos,
            codon_ref_nt=ref_codon,
            codon_alt_nt=alt_codon,
            codon_has_n=False,
            skipped=False,
        ))

    ref_has_stop = ref_stop_pos < len(ref_protein) and ref_protein[ref_stop_pos] == "*"
    alt_stop_idx = alt_protein.find("*")
    if ref_has_stop and (alt_stop_idx < 0 or alt_stop_idx >= ref_stop_pos + 3):
        stop_aa_pos = ref_stop_pos + 1
        alt_aa_here = alt_protein[ref_stop_pos] if ref_stop_pos < len(alt_protein) else "?"
        calls.append(BioCall(
            gene=gene,
            hgvs_p=f"*{stop_aa_pos}{alt_aa_here}",
            klass="stop_lost",
            ref_aa="*",
            alt_aa=alt_aa_here,
            aa_pos=stop_aa_pos,
            codon_ref_nt="",
            codon_alt_nt="",
            codon_has_n=False,
            skipped=False,
        ))

    if ref_has_stop and alt_stop_idx == ref_stop_pos and alt_early_stop < 0:
        stop_codon_start = ref_stop_pos * 3
        ref_stop_codon = ref_nt[stop_codon_start:stop_codon_start + 3]
        alt_stop_codon = alt_nt[stop_codon_start:stop_codon_start + 3]
        if ref_stop_codon and alt_stop_codon and ref_stop_codon != alt_stop_codon:
            stop_aa_pos = ref_stop_pos + 1
            calls.append(BioCall(
                gene=gene,
                hgvs_p=f"*{stop_aa_pos}*",
                klass="synonymous",
                ref_aa="*",
                alt_aa="*",
                aa_pos=stop_aa_pos,
                codon_ref_nt=ref_stop_codon,
                codon_alt_nt=alt_stop_codon,
                codon_has_n=False,
                skipped=False,
            ))

    return calls, n_skipped


def normalize_hgvs_p(value: str) -> str:
    value = value.strip()
    if value.startswith("p."):
        value = value[2:]
    return re.sub(r"[A-Z][a-z]{2}", lambda match: THREE_TO_ONE.get(match.group(0), match.group(0)), value)


def load_vcf_aa_sets(vcf_path: Path) -> dict[str, set[AAKey]]:
    result: dict[str, set[AAKey]] = defaultdict(set)
    if not vcf_path.exists():
        return dict(result)
    try:
        vcf = cyvcf2.VCF(str(vcf_path))
    except Exception as exc:
        print(f"WARNING: cannot open {vcf_path}: {exc}", file=sys.stderr)
        return dict(result)

    for record in vcf:
        if record.INFO.get("MNV_UNRESOLVED", False):
            continue
        ann_str = record.INFO.get("ANN")
        if not ann_str:
            continue
        for entry in ann_str.split(","):
            fields = entry.split("|")
            if len(fields) < 16:
                continue
            annotation = fields[1]
            gene_id = fields[4]
            gene_name = fields[3] or gene_id
            feature_type = fields[5]
            feature_id = fields[6]
            biotype = fields[7]
            hgvs_p_raw = fields[10]
            if feature_type != "transcript":
                continue
            if biotype != "protein_coding":
                continue
            if feature_id != gene_id:
                continue
            if annotation not in AA_EFFECT_TERMS:
                continue
            if not hgvs_p_raw.strip():
                continue
            klass = CLASS_MAP.get(annotation)
            if klass is None:
                continue
            hgvs_p = normalize_hgvs_p(hgvs_p_raw)
            result[gene_name].add(AAKey(gene=gene_name, hgvs_p=hgvs_p, klass=klass))
    vcf.close()
    return dict(result)


def parse_paf_header(vcf_path: Path) -> PafHeader:
    if not vcf_path.exists():
        return PafHeader("MISSING", 0.0, 0, 0)
    try:
        with open(vcf_path) as handle:
            for line in handle:
                if not line.startswith("##"):
                    break
                if line.startswith("##fluscan_alignment_coverage="):
                    content = line.strip().split("=", 1)[1]
                    parts: dict[str, str] = {}
                    for chunk in content.split(";"):
                        key, _, value = chunk.partition("=")
                        parts[key] = value
                    range_value = parts.get("range", "0-0")
                    range_start, _, range_end = range_value.partition("-")
                    return PafHeader(
                        status=parts.get("STATUS", "PASS"),
                        fraction=float(parts.get("fraction", "0")),
                        range_start=int(range_start),
                        range_end=int(range_end),
                    )
    except Exception as exc:
        print(f"WARNING: paf header parse error {vcf_path}: {exc}", file=sys.stderr)
    return PafHeader("PASS", 1.0, 0, 999_999)


def find_workflow_vcf(vcf_dir: Path, sample: str, segment_label: str, workflow_tag: str) -> Path | None:
    pattern = f"{sample}__*__*__{segment_label}_{workflow_tag}_final.vcf"
    matches = sorted(vcf_dir.glob(pattern))
    if not matches:
        return None
    if len(matches) > 1:
        print(
            f"WARNING: multiple matches for {sample}/{segment_label}/{workflow_tag}; using {matches[0].name}",
            file=sys.stderr,
        )
    return matches[0]


def _tpfpfn(workflow_set: set[AAKey], bio_set: set[AAKey]) -> tuple[int, int, int]:
    tp = len(workflow_set & bio_set)
    fp = len(workflow_set - bio_set)
    fn = len(bio_set - workflow_set)
    return tp, fp, fn


def _scoped(keys: set[AAKey], klasses: set[str]) -> set[AAKey]:
    return {key for key in keys if key.klass in klasses}


def compare_one(
    sample: str,
    segment: SegmentInfo,
    consensus_sequences: dict[str, str],
    ref_seq: str,
    cds_defs: list[CDSDef],
    config: RunConfig,
    skipped_pairs: list[str],
) -> tuple[list[dict], list[dict], list[dict]]:
    alt_seq = consensus_sequences.get(segment.fasta_tag)
    if alt_seq is None:
        print(f"WARNING: no consensus for {sample} {segment.label}", file=sys.stderr)
        return [], [], []

    nextclade_vcf = find_workflow_vcf(config.nextclade_vcf_dir, sample, segment.label, "nc")
    paftools_vcf = find_workflow_vcf(config.paftools_vcf_dir, sample, segment.label, "paf")

    if nextclade_vcf is None and paftools_vcf is None:
        skipped_pairs.append(f"{sample}/{segment.label}: both VCFs missing")
        return [], [], []
    if nextclade_vcf is None:
        skipped_pairs.append(f"{sample}/{segment.label}: NC VCF missing")
    if paftools_vcf is None:
        skipped_pairs.append(f"{sample}/{segment.label}: PAF VCF missing")

    nextclade_sets = load_vcf_aa_sets(nextclade_vcf) if nextclade_vcf else {}
    paftools_sets = load_vcf_aa_sets(paftools_vcf) if paftools_vcf else {}
    paftools_header = parse_paf_header(paftools_vcf) if paftools_vcf else PafHeader("MISSING", 0.0, 0, 0)

    ref_aln, alt_aln = align_dna(ref_seq, alt_seq)

    per_pair_rows: list[dict] = []
    disagreement_rows: list[dict] = []
    bio_calls_rows: list[dict] = []

    for cds in cds_defs:
        ref_cds_gapped, alt_cds_gapped = extract_cds_via_alignment(ref_aln, alt_aln, cds.exons)
        bio_calls, n_skipped = call_aa_mutations(cds.gene, ref_cds_gapped, alt_cds_gapped)

        bio_set = {
            AAKey(gene=call.gene, hgvs_p=call.hgvs_p, klass=call.klass)
            for call in bio_calls
            if not call.skipped
        }
        nextclade_set = nextclade_sets.get(cds.gene, set())
        paftools_set = paftools_sets.get(cds.gene, set())

        bio_outside_paf: set[AAKey] = set()
        if paftools_header.status == "WARN":
            for key in bio_set:
                call = next((item for item in bio_calls if item.hgvs_p == key.hgvs_p and item.gene == cds.gene), None)
                if call is None:
                    continue
                genomic_pos = aa_genomic_start(call.aa_pos, cds.exons)
                if genomic_pos < 0 or genomic_pos < paftools_header.range_start or genomic_pos > paftools_header.range_end:
                    bio_outside_paf.add(key)
        bio_for_paf = bio_set - bio_outside_paf

        scopes = {
            "miss": {"missense"},
            "snv": SNV_CLASSES,
            "all": ALL_CLASSES,
        }
        scope_metrics: dict[str, dict[str, int]] = {}
        for scope_key, klasses in scopes.items():
            bio_scoped = _scoped(bio_set, klasses)
            nc_scoped = _scoped(nextclade_set, klasses)
            paf_scoped = _scoped(paftools_set, klasses)
            paf_bio_scoped = _scoped(bio_for_paf, klasses)
            tp_nc, fp_nc, fn_nc = _tpfpfn(nc_scoped, bio_scoped)
            tp_paf, fp_paf, fn_paf = _tpfpfn(paf_scoped, paf_bio_scoped)
            scope_metrics[scope_key] = {
                f"n_bio_{scope_key}": len(bio_scoped),
                f"n_nc_{scope_key}": len(nc_scoped),
                f"n_paf_{scope_key}": len(paf_scoped),
                f"tp_nc_{scope_key}": tp_nc,
                f"fp_nc_{scope_key}": fp_nc,
                f"fn_nc_{scope_key}": fn_nc,
                f"tp_paf_{scope_key}": tp_paf,
                f"fp_paf_{scope_key}": fp_paf,
                f"fn_paf_{scope_key}": fn_paf,
            }

        row = {
            "sample": sample,
            "segment": segment.label,
            "contig": segment.nextclade_contig,
            "gene": cds.gene,
            "paf_status": paftools_header.status,
            "paf_fraction": paftools_header.fraction,
            "paf_range_start": paftools_header.range_start,
            "paf_range_end": paftools_header.range_end,
            "n_bio_outside_paf": len(bio_outside_paf),
            "n_codons_skipped": n_skipped,
        }
        for metrics in scope_metrics.values():
            row.update(metrics)
        per_pair_rows.append(row)

        for call in bio_calls:
            bio_calls_rows.append({
                "sample": sample,
                "segment": segment.label,
                "contig": segment.nextclade_contig,
                "gene": call.gene,
                "hgvs_p": call.hgvs_p,
                "annotation_class": call.klass,
                "ref_aa": call.ref_aa,
                "alt_aa": call.alt_aa,
                "aa_pos": call.aa_pos,
                "codon_ref_nt": call.codon_ref_nt,
                "codon_alt_nt": call.codon_alt_nt,
                "codon_has_n": call.codon_has_n,
                "skipped": call.skipped,
            })

        all_keys = bio_set | nextclade_set | paftools_set
        for key in all_keys:
            in_bio = key in bio_set
            in_nextclade = key in nextclade_set
            in_paftools = key in paftools_set
            in_bio_for_paf = key in bio_for_paf
            for workflow, in_workflow, bio_reference in [
                ("NC", in_nextclade, in_bio),
                ("PF", in_paftools, in_bio_for_paf),
            ]:
                if in_workflow and not bio_reference:
                    disagreement_rows.append({
                        "workflow": workflow,
                        "sample": sample,
                        "segment": segment.label,
                        "gene": key.gene,
                        "hgvs_p": key.hgvs_p,
                        "annotation_class": key.klass,
                        "kind": "FP",
                        "paf_status": paftools_header.status,
                    })
                elif not in_workflow and bio_reference:
                    disagreement_rows.append({
                        "workflow": workflow,
                        "sample": sample,
                        "segment": segment.label,
                        "gene": key.gene,
                        "hgvs_p": key.hgvs_p,
                        "annotation_class": key.klass,
                        "kind": "FN",
                        "paf_status": paftools_header.status,
                    })

    return per_pair_rows, disagreement_rows, bio_calls_rows


def _prf1(tp: int, fp: int, fn: int) -> tuple[float, float, float]:
    precision = tp / (tp + fp) if tp + fp > 0 else float("nan")
    recall = tp / (tp + fn) if tp + fn > 0 else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if precision + recall > 0 else float("nan")
    return precision, recall, f1


def write_headline(pp_df: pd.DataFrame, out_path: Path) -> None:
    rows = []
    for scope_key, scope_label in [("miss", "missense"), ("snv", "all_snv"), ("all", "all")]:
        for workflow, tp_col, fp_col, fn_col in [
            ("NC", f"tp_nc_{scope_key}", f"fp_nc_{scope_key}", f"fn_nc_{scope_key}"),
            ("PF_strict", f"tp_paf_{scope_key}", f"fp_paf_{scope_key}", f"fn_paf_{scope_key}"),
        ]:
            tp = int(pp_df[tp_col].sum())
            fp = int(pp_df[fp_col].sum())
            fn = int(pp_df[fn_col].sum())
            precision, recall, f1 = _prf1(tp, fp, fn)
            rows.append({
                "workflow": workflow,
                "scope": scope_label,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "precision": round(precision, 4) if precision == precision else None,
                "recall": round(recall, 4) if recall == recall else None,
                "f1": round(f1, 4) if f1 == f1 else None,
                "n_pairs": len(pp_df),
            })

        aligned = pp_df[pp_df["paf_status"] != "FAIL"]
        tp = int(aligned[f"tp_paf_{scope_key}"].sum())
        fp = int(aligned[f"fp_paf_{scope_key}"].sum())
        fn = int(aligned[f"fn_paf_{scope_key}"].sum())
        precision, recall, f1 = _prf1(tp, fp, fn)
        rows.append({
            "workflow": "PF_aligned",
            "scope": scope_label,
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": round(precision, 4) if precision == precision else None,
            "recall": round(recall, 4) if recall == recall else None,
            "f1": round(f1, 4) if f1 == f1 else None,
            "n_pairs": len(aligned),
        })

    pd.DataFrame(rows).to_csv(out_path, index=False)


def write_by_segment(pp_df: pd.DataFrame, out_path: Path, segment_labels: list[str]) -> None:
    rows = []
    for segment_label in segment_labels:
        subset = pp_df[pp_df["segment"] == segment_label]
        for scope_key, scope_label in [("miss", "missense"), ("snv", "all_snv"), ("all", "all")]:
            for workflow, tp_col, fp_col, fn_col in [
                ("NC", f"tp_nc_{scope_key}", f"fp_nc_{scope_key}", f"fn_nc_{scope_key}"),
                ("PF", f"tp_paf_{scope_key}", f"fp_paf_{scope_key}", f"fn_paf_{scope_key}"),
            ]:
                tp = int(subset[tp_col].sum())
                fp = int(subset[fp_col].sum())
                fn = int(subset[fn_col].sum())
                precision, recall, f1 = _prf1(tp, fp, fn)
                rows.append({
                    "segment": segment_label,
                    "workflow": workflow,
                    "scope": scope_label,
                    "tp": tp,
                    "fp": fp,
                    "fn": fn,
                    "precision": round(precision, 4) if precision == precision else None,
                    "recall": round(recall, 4) if recall == recall else None,
                    "f1": round(f1, 4) if f1 == f1 else None,
                    "n_pairs": len(subset),
                })
    pd.DataFrame(rows).to_csv(out_path, index=False)


def write_by_class(dis_df: pd.DataFrame, out_path: Path) -> None:
    rows = []
    for klass in sorted(ALL_CLASSES):
        for workflow in ["NC", "PF"]:
            subset = dis_df[(dis_df["annotation_class"] == klass) & (dis_df["workflow"] == workflow)]
            rows.append({
                "annotation_class": klass,
                "workflow": workflow,
                "fp": int((subset["kind"] == "FP").sum()),
                "fn": int((subset["kind"] == "FN").sum()),
            })
    pd.DataFrame(rows).to_csv(out_path, index=False)


def write_top_fp_fn(dis_df: pd.DataFrame, fp_path: Path, fn_path: Path) -> None:
    for kind, out_path in [("FP", fp_path), ("FN", fn_path)]:
        subset = dis_df[dis_df["kind"] == kind]
        if subset.empty:
            pd.DataFrame(columns=["workflow", "gene", "hgvs_p", "annotation_class", "sample_count"]).to_csv(out_path, index=False)
            continue
        grouped = (
            subset.groupby(["workflow", "gene", "hgvs_p", "annotation_class"])
            .agg(sample_count=("sample", "nunique"))
            .reset_index()
            .sort_values("sample_count", ascending=False)
        )
        grouped.to_csv(out_path, index=False)


def plot_pr_by_class(out_dir: Path) -> None:
    if px is None:
        return
    headline_path = out_dir / "headline.csv"
    if not headline_path.exists():
        return
    headline = pd.read_csv(headline_path, keep_default_na=False, na_values=[""])
    fig = px.bar(
        headline,
        x="scope",
        y="f1",
        color="workflow",
        barmode="group",
        title="F1 by Scope and Workflow",
        labels={"f1": "F1", "scope": "Scope"},
        color_discrete_map={"NC": "#1976D2", "PF_strict": "#D32F2F", "PF_aligned": "#E57373"},
    )
    plots_dir = out_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(plots_dir / "precision_recall_by_class.html"))


def plot_f1_by_segment(out_dir: Path, segment_labels: list[str]) -> None:
    if px is None:
        return
    segment_path = out_dir / "by_segment.csv"
    if not segment_path.exists():
        return
    segments = pd.read_csv(segment_path, keep_default_na=False, na_values=[""])
    missense = segments[(segments["scope"] == "missense") & (segments["workflow"].isin(["NC", "PF"]))]
    fig = px.bar(
        missense,
        x="segment",
        y="f1",
        color="workflow",
        barmode="group",
        category_orders={"segment": segment_labels},
        title="Missense F1 by Segment",
        labels={"f1": "F1 (missense)", "segment": "Segment"},
        color_discrete_map={"NC": "#1976D2", "PF": "#D32F2F"},
    )
    plots_dir = out_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(plots_dir / "f1_by_segment.html"))


def plot_fp_fn_heatmap(pp_df: pd.DataFrame, out_dir: Path, segment_labels: list[str]) -> None:
    if go is None:
        return
    samples = sorted(pp_df["sample"].unique())

    def build_matrix(fp_col: str, fn_col: str) -> list[list[int | None]]:
        matrix: list[list[int | None]] = []
        for sample in samples:
            row: list[int | None] = []
            for segment_label in segment_labels:
                subset = pp_df[(pp_df["sample"] == sample) & (pp_df["segment"] == segment_label)]
                if subset.empty:
                    row.append(None)
                else:
                    row.append(int(subset[fp_col].sum() + subset[fn_col].sum()))
            matrix.append(row)
        return matrix

    plots_dir = out_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    nc_matrix = build_matrix("fp_nc_miss", "fn_nc_miss")
    paf_matrix = build_matrix("fp_paf_miss", "fn_paf_miss")

    nc_fig = go.Figure(go.Heatmap(
        z=nc_matrix,
        x=segment_labels,
        y=samples,
        colorscale="Reds",
        colorbar=dict(title="FP+FN (missense)"),
    ))
    nc_fig.update_layout(
        title="Sample x Segment Missense (FP+FN) - Nextclade",
        height=max(400, len(samples) * 12),
        xaxis_title="Segment",
        yaxis_title="Sample",
    )
    nc_fig.write_html(str(plots_dir / "fp_fn_heatmap.html"))

    paf_fig = go.Figure(go.Heatmap(
        z=paf_matrix,
        x=segment_labels,
        y=samples,
        colorscale="Reds",
        colorbar=dict(title="FP+FN (missense)"),
    ))
    paf_fig.update_layout(
        title="Sample x Segment Missense (FP+FN) - paftools (strict)",
        height=max(400, len(samples) * 12),
        xaxis_title="Segment",
        yaxis_title="Sample",
    )
    paf_fig.write_html(str(plots_dir / "fp_fn_heatmap_paf.html"))


def write_report(pp_df: pd.DataFrame, skipped_pairs: list[str], out_dir: Path, segment_labels: list[str]) -> None:
    headline_path = out_dir / "headline.csv"
    headline = pd.read_csv(headline_path, keep_default_na=False, na_values=[""]) if headline_path.exists() else pd.DataFrame()

    lines: list[str] = [
        "# Biopython Sanity Check Report",
        "",
        f"Samples: {pp_df['sample'].nunique()}  |  Protein triples: {len(pp_df)}  |  Skipped pairs: {len(skipped_pairs)}",
        "",
    ]

    if not headline.empty:
        lines += ["## Headline Metrics", ""]
        lines.append("| Workflow | Scope | TP | FP | FN | Precision | Recall | F1 | N pairs |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        for _, row in headline.iterrows():
            lines.append(
                f"| {row['workflow']} | {row['scope']} | {row['tp']} | {row['fp']} | {row['fn']} "
                f"| {row['precision']} | {row['recall']} | {row['f1']} | {row['n_pairs']} |"
            )
        lines.append("")

        missense = headline[headline["scope"] == "missense"]
        nc_row = missense[missense["workflow"] == "NC"]
        pf_row = missense[missense["workflow"] == "PF_strict"]
        if not nc_row.empty and not pf_row.empty:
            nc = nc_row.iloc[0]
            pf = pf_row.iloc[0]
            winner = "Nextclade" if (nc["f1"] or 0) >= (pf["f1"] or 0) else "paftools"
            lines += [
                "## Final Recommendation",
                "",
                f"**{winner}** has higher missense F1 against the Biopython ground truth.",
                f"- Nextclade: F1={nc['f1']:.4f}, Recall={nc['recall']:.4f}",
                f"- paftools (strict): F1={pf['f1']:.4f}, Recall={pf['recall']:.4f}",
                "",
            ]

    by_segment_path = out_dir / "by_segment.csv"
    if by_segment_path.exists():
        by_segment = pd.read_csv(by_segment_path, keep_default_na=False, na_values=[""])
        missense = by_segment[(by_segment["scope"] == "missense") & (by_segment["workflow"].isin(["NC", "PF"]))]
        if not missense.empty:
            lines += ["## Missense F1 by Segment", ""]
            lines.append("| Segment | NC F1 | PF F1 | Better |")
            lines.append("|---|---|---|---|")
            for segment_label in segment_labels:
                nc_subset = missense[(missense["segment"] == segment_label) & (missense["workflow"] == "NC")]
                pf_subset = missense[(missense["segment"] == segment_label) & (missense["workflow"] == "PF")]
                nc_f1 = round(nc_subset.iloc[0]["f1"], 4) if not nc_subset.empty and nc_subset.iloc[0]["f1"] == nc_subset.iloc[0]["f1"] else None
                pf_f1 = round(pf_subset.iloc[0]["f1"], 4) if not pf_subset.empty and pf_subset.iloc[0]["f1"] == pf_subset.iloc[0]["f1"] else None
                better = "NC" if (nc_f1 or 0) >= (pf_f1 or 0) else "PF"
                lines.append(f"| {segment_label} | {nc_f1} | {pf_f1} | {better} |")
            lines.append("")

    for kind, heading in [("fp", "Top Systematic Over-calls"), ("fn", "Top Systematic Misses")]:
        csv_path = out_dir / f"top_{kind}.csv"
        if not csv_path.exists():
            continue
        table = pd.read_csv(csv_path, keep_default_na=False, na_values=[""])
        if table.empty:
            continue
        lines += [f"## {heading} (top 10, any workflow)", ""]
        for _, row in table.head(10).iterrows():
            lines.append(f"- [{row['workflow']}] {row['gene']} {row['hgvs_p']} ({row['annotation_class']}): {row['sample_count']} samples")
        lines.append("")

    status_counts = pp_df.groupby(["segment", "paf_status"]).size().unstack(fill_value=0)
    lines += ["## paftools Alignment Status by Segment", "", str(status_counts), ""]

    if skipped_pairs:
        lines += ["## Skipped Pairs", ""]
        lines.extend(f"- {pair}" for pair in sorted(skipped_pairs))
        lines.append("")

    lines += [
        "## Caveats",
        "- Biopython is an independent check, not absolute truth.",
        "- Indel and frameshift HGVS representation may differ between tools.",
        "- N-containing codons in consensus are skipped; see n_codons_skipped in per_pair.csv.",
        "- PF_strict counts FAIL segments as all-FN; PF_aligned excludes them.",
        "- WARN: Biopython calls outside the paftools coverage range are excluded from paftools FN.",
    ]

    (out_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    ensure_runtime_dependencies()
    config = build_run_config(args)

    segment_infos = load_segment_info(config.reference_datasets_dir)
    segment_labels = [segment.label for segment in segment_infos]

    print("Loading references and GFF3 ...", flush=True)
    refs = {segment.label: load_reference(config.reference_datasets_dir, segment) for segment in segment_infos}
    cdss = {segment.label: parse_gff3(config.reference_datasets_dir, segment) for segment in segment_infos}

    print("Loading consensus FASTAs ...", flush=True)
    all_consensus = load_all_consensus(config.consensus_dir, config.sample_filter)
    samples = sorted(all_consensus)
    print(f"  {len(samples)} samples", flush=True)

    per_pair_rows: list[dict] = []
    disagreement_rows: list[dict] = []
    bio_calls_rows: list[dict] = []
    skipped_pairs: list[str] = []

    total = len(samples) * len(segment_infos)
    done = 0
    for sample in samples:
        for segment in segment_infos:
            done += 1
            if done % 200 == 0 or done == total:
                print(f"  {done}/{total} ...", flush=True)
            pp_rows, dis_rows, bio_rows = compare_one(
                sample=sample,
                segment=segment,
                consensus_sequences=all_consensus[sample],
                ref_seq=refs[segment.label],
                cds_defs=cdss[segment.label],
                config=config,
                skipped_pairs=skipped_pairs,
            )
            per_pair_rows.extend(pp_rows)
            disagreement_rows.extend(dis_rows)
            bio_calls_rows.extend(bio_rows)

    print(f"Writing CSVs to {config.out_dir} ...", flush=True)
    config.out_dir.mkdir(parents=True, exist_ok=True)
    (config.out_dir / "plots").mkdir(exist_ok=True)

    pp_df = pd.DataFrame(per_pair_rows)
    dis_df = pd.DataFrame(disagreement_rows)
    bio_df = pd.DataFrame(bio_calls_rows)

    pp_df.to_csv(config.out_dir / "per_pair.csv", index=False)
    dis_df.to_csv(config.out_dir / "disagreements.csv", index=False)
    bio_df.to_csv(config.out_dir / "biopython_calls.csv", index=False)

    write_headline(pp_df, config.out_dir / "headline.csv")
    write_by_segment(pp_df, config.out_dir / "by_segment.csv", segment_labels)
    write_by_class(dis_df, config.out_dir / "by_class.csv")
    write_top_fp_fn(dis_df, config.out_dir / "top_fp.csv", config.out_dir / "top_fn.csv")
    write_report(pp_df, skipped_pairs, config.out_dir, segment_labels)

    if config.skip_plots:
        print("Skipping plots (--skip-plots).", flush=True)
    elif px is None or go is None:
        print("Skipping plots because plotly is not installed in the active environment.", flush=True)
    else:
        print("Writing plots ...", flush=True)
        plot_pr_by_class(config.out_dir)
        plot_f1_by_segment(config.out_dir, segment_labels)
        plot_fp_fn_heatmap(pp_df, config.out_dir, segment_labels)

    if not pp_df.empty:
        n_paf_fail = int((pp_df["paf_status"] == "FAIL").sum())
        n_paf_partial = int((pp_df["paf_status"] == "WARN").sum())
        n_paf_pass = int((pp_df["paf_status"] == "PASS").sum())
        print("\n=== Summary ===")
        print(f"Protein triples : {len(pp_df)}")
        print(f"  paf PASS      : {n_paf_pass}")
        print(f"  paf PARTIAL   : {n_paf_partial}")
        print(f"  paf NO_ALIGN  : {n_paf_fail}")
        print(f"Skipped pairs   : {len(skipped_pairs)}")
        for scope_key, scope_label in [("miss", "missense"), ("snv", "all_snv"), ("all", "all")]:
            tp_nc = int(pp_df[f"tp_nc_{scope_key}"].sum())
            fp_nc = int(pp_df[f"fp_nc_{scope_key}"].sum())
            fn_nc = int(pp_df[f"fn_nc_{scope_key}"].sum())
            tp_paf = int(pp_df[f"tp_paf_{scope_key}"].sum())
            fp_paf = int(pp_df[f"fp_paf_{scope_key}"].sum())
            fn_paf = int(pp_df[f"fn_paf_{scope_key}"].sum())
            _, recall_nc, f1_nc = _prf1(tp_nc, fp_nc, fn_nc)
            _, recall_paf, f1_paf = _prf1(tp_paf, fp_paf, fn_paf)
            print(f"\n  [{scope_label}]")
            print(f"    NC  : TP={tp_nc}  FP={fp_nc}  FN={fn_nc}  recall={recall_nc:.4f}  F1={f1_nc:.4f}")
            print(f"    PAF : TP={tp_paf} FP={fp_paf} FN={fn_paf} recall={recall_paf:.4f} F1={f1_paf:.4f}")
    print(f"\nOutputs in {config.out_dir}")


if __name__ == "__main__":
    main()