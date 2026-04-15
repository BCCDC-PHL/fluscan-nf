#!/usr/bin/env python3
"""fix_mnv_codons.py

Post-process a snpEff-annotated LoFreq VCF to correct amino acid annotations
where two or more adjacent SNVs fall within the same codon (MNV problem).

Background
----------
snpEff annotates each VCF record independently against the reference.  When
two SNVs at adjacent positions occupy the same codon, snpEff evaluates each
in isolation and reports incorrect amino acid changes.  For haploid influenza
consensus calls (AF >= --min-af), adjacent high-frequency SNVs in the same
codon are mathematically guaranteed to co-occur on the same viral genome
(pigeonhole principle: if f1 + f2 > 1, minimum co-occurrence >= f1+f2-1).

This script reconstructs the true combined codon and overwrites the HGVS.p,
annotation type (missense/synonymous), and impact (MODERATE/LOW) in the INFO
ANN field accordingly.

Algorithm
---------
1. First pass  — scan all PASS SNVs with AF >= min_af.  For each coding
   ANN entry, compute the codon number and offset from the CDS position.
   Group SNVs by (CHROM, snpEff feature_id, codon_number).

2. For each group with >= 2 SNVs:
     a. Reconstruct the 3-base reference codon from the VCF REF nucleotides
        already present in the ANN hgvs_c field.  For any codon position not
        covered by a variant, infer the reference nucleotide by elimination
        (try all 4 bases, pick the one consistent with the annotated ref AA).
        If inference is ambiguous, skip the codon and log a warning.
     b. Apply all ALT nucleotides to obtain the mutant codon.
     c. Translate both codons to determine the true amino acid change.

3. Second pass — rewrite VCF records that belong to a corrected codon group:
     • Update ANN field: HGVS.p, annotation type, impact for each affected
       transcript entry.
     • Set INFO flag   MNV_CORRECTED = True
     • Set INFO string MNV_ORIGINAL_HGVS_P = original HGVS.p values

4. Optionally write a flat TSV with one row per coding ANN entry.

Usage
-----
    fix_mnv_codons.py -i sample_anno.vcf \\
        --output-vcf sample_anno_mnv_corrected.vcf \\
        --output-tsv sample_anno_mnv_corrected.tsv \\
        --min-af 0.5

Pipeline integration (after snpEff, before reformat_vcf):
    run_snpeff.py -> fix_mnv_codons.py -> reformat_vcf_cyvcf2.py
"""

import argparse
import csv
import logging
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

try:
    from cyvcf2 import VCF, Writer
except ImportError as exc:
    raise SystemExit(
        "ERROR: cyvcf2 is required. Install it with: conda install -c bioconda cyvcf2"
    ) from exc

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# snpEff ANN field indices (0-based, pipe-delimited within each comma entry)
ANN_ALLELE_IDX = 0
ANN_ANNOTATION_IDX = 1
ANN_IMPACT_IDX = 2
ANN_GENE_NAME_IDX = 3
ANN_GENE_ID_IDX = 4
ANN_FEATURE_TYPE_IDX = 5
ANN_FEATURE_ID_IDX = 6
ANN_BIOTYPE_IDX = 7
ANN_RANK_IDX = 8
ANN_HGVS_C_IDX = 9
ANN_HGVS_P_IDX = 10
ANN_CDNA_POS_IDX = 11
ANN_CDS_POS_IDX = 12
ANN_PROT_POS_IDX = 13
ANN_DISTANCE_IDX = 14
ANN_NOTES_IDX = 15
ANN_FIELD_COUNT = 16

# Amino acid 3-letter → 1-letter
AA_3_TO_1: Dict[str, str] = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}

# Amino acid 1-letter → 3-letter (snpEff style)
AA_1_TO_3: Dict[str, str] = {v: k for k, v in AA_3_TO_1.items()}
AA_1_TO_3["*"] = "Ter"

# Standard DNA codon table (stop = "*")
_CODON_TABLE: Dict[str, str] = {
    "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
    "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
    "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
    "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F",
}

# (annotation_type, impact) for each outcome
_ANN_MISSENSE = ("missense_variant", "MODERATE")
_ANN_SYNONYMOUS = ("synonymous_variant", "LOW")
_ANN_STOP_GAINED = ("stop_gained", "HIGH")
_ANN_STOP_LOST = ("stop_lost", "HIGH")

# Regex for snpEff hgvs_c: c.263A>G or c.-10A>C (negative = UTR, skip those)
_HGVS_C_RE = re.compile(r"c\.(\d+)([ACGTacgt])>([ACGTacgt])", re.ASCII)

# Regex for snpEff hgvs_p —  two styles depending on snpEff flags:
#   3-letter: p.Asn88Ser   (default snpEff output)
#   1-letter: p.N88S        (when -hgvs1LetterAa is used — fluscan-nf default)
_HGVS_P_3L_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)", re.ASCII)
_HGVS_P_1L_RE = re.compile(r"p\.([A-Z\*])(\d+)", re.ASCII)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class CodonSNV:
    """One SNV contributing to a codon group for a specific transcript."""
    chrom: str
    vcf_pos: int           # 1-based VCF POS
    ref_nuc: str           # Reference nucleotide (from VCF REF, single base)
    alt_nuc: str           # Alternate nucleotide (from VCF ALT, single base)
    codon_offset: int      # 0=1st base, 1=2nd base, 2=3rd base of codon
    cds_pos: int           # 1-based CDS position from snpEff ANN
    protein_pos: int       # 1-based protein (AA) position
    ref_aa_1: Optional[str]       # 1-letter reference AA parsed from hgvs_p
    hgvs_c_ref_nuc: Optional[str] # Reference nucleotide from hgvs_c (should match REF)


@dataclass
class AnnCorrection:
    """Correction to apply to a specific (chrom, vcf_pos, feature_id) triplet."""
    feature_id: str
    protein_pos: int       # 1-based protein position
    ref_aa_1: str          # 1-letter reference AA
    true_alt_aa_1: str     # 1-letter corrected alt AA (combined codon result)
    codon_group_id: str    # descriptive group key for logging/INFO tags


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def translate_codon(codon: str) -> Optional[str]:
    """Translate a 3-base uppercase DNA codon to 1-letter AA, or None."""
    return _CODON_TABLE.get(codon.upper())


def parse_hgvs_c(hgvs_c: str) -> Tuple[Optional[int], Optional[str], Optional[str]]:
    """
    Parse snpEff hgvs_c like 'c.263A>G'.
    Returns (cds_pos, ref_nuc, alt_nuc).  Negative positions (UTR) return Nones.
    """
    if not hgvs_c:
        return None, None, None
    m = _HGVS_C_RE.match(hgvs_c)
    if not m:
        return None, None, None
    cds_pos = int(m.group(1))
    return cds_pos, m.group(2).upper(), m.group(3).upper()


def parse_hgvs_p_ref_aa_1(hgvs_p: str) -> Optional[str]:
    """
    Extract 1-letter reference AA from snpEff hgvs_p.
    Handles both 3-letter ('p.Asn88Ser') and 1-letter ('p.N88S') formats.
    Returns the 1-letter code (e.g. 'N'), or None.
    """
    if not hgvs_p:
        return None
    # Try 3-letter first (e.g. p.Asn88Ser)
    m3 = _HGVS_P_3L_RE.match(hgvs_p)
    if m3:
        return AA_3_TO_1.get(m3.group(1))
    # Try 1-letter (e.g. p.N88S) — fluscan-nf default with -hgvs1LetterAa
    m1 = _HGVS_P_1L_RE.match(hgvs_p)
    if m1:
        ref = m1.group(1)
        return ref if ref in AA_1_TO_3 or ref == "*" else None
    return None


def parse_pos_len(field: str) -> Tuple[Optional[int], Optional[str]]:
    """Parse '263/1707' into (263, '1707'), or ('', None) if invalid."""
    if not field or field in ("", "."):
        return None, None
    parts = str(field).split("/", 1)
    try:
        pos = int(parts[0])
        length = parts[1] if len(parts) > 1 else None
        return (pos if pos > 0 else None), length
    except (ValueError, IndexError):
        return None, None


def parse_ann_entries(raw_ann: str) -> List[List[str]]:
    """
    Split the ANN INFO value into entries, each a list of ANN_FIELD_COUNT fields.
    Returns [] for empty/None input.
    """
    if not raw_ann:
        return []
    entries = []
    for entry in str(raw_ann).split(","):
        parts = entry.split("|")
        # Pad to full field count
        parts.extend([""] * max(0, ANN_FIELD_COUNT - len(parts)))
        entries.append(parts[:ANN_FIELD_COUNT])
    return entries


def reassemble_ann(entries: List[List[str]]) -> str:
    """Reassemble list-of-field-lists back into a VCF ANN string."""
    return ",".join("|".join(e) for e in entries)


def _is_1letter_hgvs_p(hgvs_p: str) -> bool:
    """Return True if hgvs_p uses 1-letter AA codes (e.g. 'p.N88S'), False for 3-letter."""
    if not hgvs_p:
        return True  # default to 1-letter (fluscan-nf standard)
    return bool(_HGVS_P_1L_RE.match(hgvs_p) and not _HGVS_P_3L_RE.match(hgvs_p))


def build_hgvs_p(ref_aa_1: str, protein_pos: int, alt_aa_1: str, one_letter: bool = True) -> str:
    """
    Build snpEff HGVS.p string.  Matches the format used in the original annotation:
      one_letter=True  → 'p.N88R'     (fluscan-nf default with -hgvs1LetterAa)
      one_letter=False → 'p.Asn88Arg' (standard snpEff 3-letter default)
    """
    if one_letter:
        return f"p.{ref_aa_1}{protein_pos}{alt_aa_1}"
    ref_3 = AA_1_TO_3.get(ref_aa_1, ref_aa_1)
    alt_3 = AA_1_TO_3.get(alt_aa_1, alt_aa_1)
    return f"p.{ref_3}{protein_pos}{alt_3}"


def determine_annotation(ref_aa: str, alt_aa: str) -> Tuple[str, str]:
    """Return (annotation_type, impact) for an amino acid change."""
    if alt_aa == "*" and ref_aa != "*":
        return _ANN_STOP_GAINED
    if ref_aa == "*" and alt_aa != "*":
        return _ANN_STOP_LOST
    if ref_aa == alt_aa:
        return _ANN_SYNONYMOUS
    return _ANN_MISSENSE


# ---------------------------------------------------------------------------
# Codon reconstruction
# ---------------------------------------------------------------------------

def resolve_mnv_codon(
    snvs: List[CodonSNV],
) -> Optional[Tuple[str, str, str, str]]:
    """
    Resolve the reference and mutant codons for a group of SNVs in the same codon.
    Returns (ref_codon, mut_codon, ref_aa_1, true_mut_aa_1) or None.

    Strategy:
      - For codon positions covered by variants: use VCF REF/ALT nucleotides.
      - For any uncovered position: enumerate all 4 nucleotides, keep those
        consistent with the annotated reference amino acid.
      - If only one reference nucleotide is consistent -> unambiguous, proceed.
      - If multiple reference nucleotides are consistent but all produce the SAME
        mutant amino acid -> correction is still unambiguous, proceed.
        (E.g. Lys codon AAA/AAG + two other alts -> both produce Arg -> K->R is clear.)
      - Otherwise log a warning and return None.
    """
    known_ref: Dict[int, str] = {}
    known_alt: Dict[int, str] = {}

    for snv in snvs:
        offset = snv.codon_offset
        ref_nuc = snv.ref_nuc.upper()
        if snv.hgvs_c_ref_nuc and snv.hgvs_c_ref_nuc != ref_nuc:
            log.warning(
                "hgvs_c ref (%s) != VCF REF (%s) at CDS %d / POS %d — using VCF REF",
                snv.hgvs_c_ref_nuc, ref_nuc, snv.cds_pos, snv.vcf_pos,
            )
        if offset in known_ref and known_ref[offset] != ref_nuc:
            log.warning(
                "Conflicting REF bases at codon offset %d: %s vs %s — skipping",
                offset, known_ref[offset], ref_nuc,
            )
            return None
        known_ref[offset] = ref_nuc
        known_alt[offset] = snv.alt_nuc.upper()

    missing = [i for i in range(3) if i not in known_ref]

    if not missing:
        ref_codon = "".join(known_ref[i] for i in range(3))
        mut_codon = "".join(known_alt.get(i, known_ref[i]) for i in range(3))
        return ref_codon, mut_codon, translate_codon(ref_codon), translate_codon(mut_codon)

    if len(missing) > 1:
        log.debug("Too many uncovered codon positions (%d) — skipping", len(missing))
        return None

    ref_aa_1 = next((s.ref_aa_1 for s in snvs if s.ref_aa_1), None)
    if ref_aa_1 is None:
        log.debug("No reference AA available for codon with 1 uncovered position — skipping")
        return None

    pos = missing[0]
    valid_nucs = []
    for nuc in "ACGT":
        rc = [known_ref.get(i, nuc if i == pos else "N") for i in range(3)]
        if translate_codon("".join(rc)) == ref_aa_1:
            valid_nucs.append(nuc)

    if not valid_nucs:
        log.warning(
            "No nucleotide at codon offset %d translates to ref AA '%s' — skipping",
            pos, ref_aa_1,
        )
        return None

    mut_aa_set: set = set()
    per_nuc: Dict[str, Tuple[str, str]] = {}
    for nuc in valid_nucs:
        rc = [known_ref.get(i, nuc if i == pos else "N") for i in range(3)]
        mc = [known_alt.get(i, rc[i]) for i in range(3)]
        mut_aa = translate_codon("".join(mc))
        if mut_aa:
            mut_aa_set.add(mut_aa)
            per_nuc[nuc] = ("".join(rc), "".join(mc))

    if len(mut_aa_set) == 1:
        if len(valid_nucs) > 1:
            log.debug(
                "Reference codon offset %d is ambiguous (%s -> ref_aa '%s') but "
                "mutant AA is unambiguous ('%s'). Proceeding.",
                pos, "/".join(valid_nucs), ref_aa_1, next(iter(mut_aa_set)),
            )
        nuc = valid_nucs[0]
        ref_codon, mut_codon = per_nuc[nuc]
        return ref_codon, mut_codon, ref_aa_1, mut_aa_set.pop()

    log.warning(
        "Ambiguous MNV: offset %d could be %s (ref_aa '%s') and mutant AA would be %s. "
        "Provide --reference FASTA — skipping.",
        pos, "/".join(valid_nucs), ref_aa_1, "/".join(sorted(mut_aa_set)),
    )
    return None


# ---------------------------------------------------------------------------
# Pass 1 — collect SNVs and group by codon
# ---------------------------------------------------------------------------

def collect_codon_groups(
    vcf_path: str,
    min_af: float,
) -> Dict[Tuple[str, str, int], List[CodonSNV]]:
    """
    First pass: read all PASS SNVs with AF >= min_af; group by
    (chrom, feature_id, codon_number) using CDS positions from the ANN field.

    Returns the codon_groups dict.
    """
    codon_groups: Dict[Tuple[str, str, int], List[CodonSNV]] = defaultdict(list)

    vcf = VCF(vcf_path)
    n_total = n_skipped = 0
    try:
        for variant in vcf:
            n_total += 1

            # Only PASS records (cyvcf2: FILTER is None → PASS or ".")
            if variant.FILTER is not None:
                n_skipped += 1
                continue

            # Only SNVs; lofreq calls indels separately
            if not variant.is_snp:
                n_skipped += 1
                continue

            # Get AF
            af_raw = variant.INFO.get("AF")
            if af_raw is None:
                n_skipped += 1
                continue
            try:
                af = float(af_raw) if not isinstance(af_raw, (list, tuple)) else float(af_raw[0])
            except (TypeError, ValueError):
                n_skipped += 1
                continue
            if af < min_af:
                n_skipped += 1
                continue

            raw_ann = variant.INFO.get("ANN")
            if raw_ann is None:
                continue

            chrom = variant.CHROM
            pos = variant.POS
            ref_nuc = variant.REF.upper()

            # lofreq always emits one ALT per record for SNVs
            alts = [a for a in variant.ALT if a is not None]
            if len(alts) != 1:
                continue
            alt_nuc = alts[0].upper()

            ann_entries = parse_ann_entries(str(raw_ann))

            for entry in ann_entries:
                feature_id = entry[ANN_FEATURE_ID_IDX]
                hgvs_c = entry[ANN_HGVS_C_IDX]
                hgvs_p = entry[ANN_HGVS_P_IDX]
                cds_raw = entry[ANN_CDS_POS_IDX]
                prot_raw = entry[ANN_PROT_POS_IDX]

                cds_pos, _ = parse_pos_len(cds_raw)
                if cds_pos is None:
                    continue  # non-coding annotation

                protein_pos, _ = parse_pos_len(prot_raw)
                if protein_pos is None:
                    protein_pos = math.ceil(cds_pos / 3)

                codon_offset = (cds_pos - 1) % 3
                codon_number = math.ceil(cds_pos / 3)

                _, hgvs_c_ref, _ = parse_hgvs_c(hgvs_c)
                ref_aa_1_snv = parse_hgvs_p_ref_aa_1(hgvs_p)

                snv = CodonSNV(
                    chrom=chrom,
                    vcf_pos=pos,
                    ref_nuc=ref_nuc,
                    alt_nuc=alt_nuc,
                    codon_offset=codon_offset,
                    cds_pos=cds_pos,
                    protein_pos=protein_pos,
                    ref_aa_1=ref_aa_1_snv,
                    hgvs_c_ref_nuc=hgvs_c_ref,
                )
                codon_groups[(chrom, feature_id, codon_number)].append(snv)
    finally:
        vcf.close()

    log.debug("Pass 1: %d records read, %d skipped (non-PASS / non-SNV / low-AF)", n_total, n_skipped)
    return dict(codon_groups)


# ---------------------------------------------------------------------------
# Compute corrections from codon groups
# ---------------------------------------------------------------------------

def compute_corrections(
    codon_groups: Dict[Tuple[str, str, int], List[CodonSNV]],
) -> Dict[Tuple[str, int], List[AnnCorrection]]:
    """
    For each codon group with >= 2 SNVs, resolve the true combined amino acid
    change and build AnnCorrection objects keyed by (chrom, vcf_pos).
    """
    corrections: Dict[Tuple[str, int], List[AnnCorrection]] = defaultdict(list)
    n_groups = n_corrected = n_skipped = 0

    for (chrom, feature_id, codon_number), snvs in codon_groups.items():
        if len(snvs) < 2:
            continue
        n_groups += 1

        # Guard: all SNVs must map to the same codon number
        codon_nums = {math.ceil(s.cds_pos / 3) for s in snvs}
        if len(codon_nums) > 1:
            log.warning(
                "Mixed codon numbers %s in group (%s, %s, %d) — skipping",
                codon_nums, chrom, feature_id, codon_number,
            )
            n_skipped += 1
            continue

        # Guard: no two variants at the same codon offset
        offsets = [s.codon_offset for s in snvs]
        if len(set(offsets)) < len(offsets):
            log.warning(
                "Duplicate codon offsets %s in group (%s, %s, %d) — skipping",
                offsets, chrom, feature_id, codon_number,
            )
            n_skipped += 1
            continue

        # Resolve codons (handles ambiguous reference with unambiguous mutant)
        resolved = resolve_mnv_codon(snvs)
        if resolved is None:
            n_skipped += 1
            continue

        ref_codon, mut_codon, ref_aa_1, true_alt_aa_1 = resolved

        if true_alt_aa_1 is None:
            log.warning(
                "Cannot translate mutant codon '%s' for group (%s, %s, %d) — skipping",
                mut_codon, chrom, feature_id, codon_number,
            )
            n_skipped += 1
            continue

        protein_pos = snvs[0].protein_pos
        codon_group_id = f"{chrom}:{feature_id}:codon{codon_number}"
        vcf_positions = sorted(s.vcf_pos for s in snvs)

        log.info(
            "MNV  %s | %s%d → %s  |  ref_codon=%s  mut_codon=%s  |  "
            "VCF positions: %s",
            codon_group_id,
            ref_aa_1, protein_pos, true_alt_aa_1,
            ref_codon, mut_codon,
            vcf_positions,
        )

        correction = AnnCorrection(
            feature_id=feature_id,
            protein_pos=protein_pos,
            ref_aa_1=ref_aa_1,
            true_alt_aa_1=true_alt_aa_1,
            codon_group_id=codon_group_id,
        )
        for snv in snvs:
            corrections[(chrom, snv.vcf_pos)].append(correction)

        n_corrected += 1

    log.info(
        "Codon groups: %d with >=2 SNVs found, %d correctable, %d skipped",
        n_groups, n_corrected, n_skipped,
    )
    return dict(corrections)


# ---------------------------------------------------------------------------
# Apply corrections to ANN field
# ---------------------------------------------------------------------------

def apply_corrections_to_ann(
    raw_ann: str,
    corrections: List[AnnCorrection],
) -> str:
    """
    Update HGVS.p, annotation type, and impact in each ANN entry whose
    feature_id matches a correction.  Returns the modified ANN string.
    """
    entries = parse_ann_entries(raw_ann)
    by_feature = {c.feature_id: c for c in corrections}

    for entry in entries:
        fid = entry[ANN_FEATURE_ID_IDX]
        corr = by_feature.get(fid)
        if corr is None:
            continue

        # Preserve original HGVS.p format (1-letter vs 3-letter)
        one_letter = _is_1letter_hgvs_p(entry[ANN_HGVS_P_IDX])
        new_hgvs_p = build_hgvs_p(corr.ref_aa_1, corr.protein_pos, corr.true_alt_aa_1, one_letter)
        new_annotation, new_impact = determine_annotation(corr.ref_aa_1, corr.true_alt_aa_1)

        entry[ANN_HGVS_P_IDX] = new_hgvs_p
        entry[ANN_ANNOTATION_IDX] = new_annotation
        entry[ANN_IMPACT_IDX] = new_impact

    return reassemble_ann(entries)


# ---------------------------------------------------------------------------
# Pass 2 — write corrected VCF
# ---------------------------------------------------------------------------

def write_corrected_vcf(
    vcf_path: str,
    output_vcf_path: str,
    corrections: Dict[Tuple[str, int], List[AnnCorrection]],
) -> int:
    """
    Second pass: emit all records to a new VCF.  For records belonging to a
    corrected codon group, update the ANN field and add MNV_CORRECTED / 
    MNV_ORIGINAL_HGVS_P INFO tags.  Returns the count of corrected records.
    """
    vcf = VCF(vcf_path)

    vcf.add_info_to_header({
        "ID": "MNV_CORRECTED",
        "Number": "0",
        "Type": "Flag",
        "Description": (
            "Variant is part of an MNV codon group. The snpEff ANN HGVS.p, "
            "annotation type, and impact have been corrected to reflect the "
            "true combined amino acid effect of all co-occurring SNVs in the "
            "same codon (fix_mnv_codons.py)."
        ),
    })
    vcf.add_info_to_header({
        "ID": "MNV_ORIGINAL_HGVS_P",
        "Number": "1",
        "Type": "String",
        "Description": (
            "Original snpEff HGVS.p annotations before MNV correction. "
            "Format: featureID:hgvs_p pairs separated by pipe (|)."
        ),
    })

    writer = Writer(output_vcf_path, vcf)
    n_corrected = 0

    try:
        for variant in vcf:
            key = (variant.CHROM, variant.POS)
            correction_list = corrections.get(key)

            if correction_list:
                raw_ann = variant.INFO.get("ANN")
                if raw_ann is not None:
                    original_ann = str(raw_ann)

                    # Capture original HGVS.p values for the features we are about to change
                    orig_entries = parse_ann_entries(original_ann)
                    affected_features = {c.feature_id for c in correction_list}
                    orig_hgvs_parts = []
                    for orig_entry in orig_entries:
                        fid = orig_entry[ANN_FEATURE_ID_IDX]
                        hp = orig_entry[ANN_HGVS_P_IDX]
                        if fid in affected_features and hp:
                            orig_hgvs_parts.append(f"{fid}:{hp}")
                    orig_hgvs_p_str = "|".join(orig_hgvs_parts)

                    corrected_ann = apply_corrections_to_ann(original_ann, correction_list)

                    if corrected_ann != original_ann:
                        variant.INFO["ANN"] = corrected_ann
                        variant.INFO["MNV_CORRECTED"] = True
                        if orig_hgvs_p_str:
                            variant.INFO["MNV_ORIGINAL_HGVS_P"] = orig_hgvs_p_str
                        n_corrected += 1

            writer.write_record(variant)
    finally:
        writer.close()
        vcf.close()

    return n_corrected


# ---------------------------------------------------------------------------
# Optional TSV output
# ---------------------------------------------------------------------------

TSV_COLUMNS = [
    "chrom", "pos", "ref", "alt", "qual", "filter",
    "dp", "af", "sb", "dp4",
    "annotation", "impact", "gene_name", "feature_id",
    "hgvs_c", "hgvs_p", "cds_pos", "cds_len",
    "protein_pos", "protein_len",
    "aa_ref", "aa_pos", "aa_alt", "aa_change",
    "mnv_corrected", "mnv_original_hgvs_p",
]


def _parse_aa_parts(
    hgvs_p: str,
) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    """
    Return (aa_ref_1, aa_pos, aa_alt_1, aa_change_1letter) from
    snpEff hgvs_p like 'p.Asn88Arg'.  Returns Nones for unparseable input.
    """
    if not hgvs_p:
        return None, None, None, None
    protein_change = hgvs_p.removeprefix("p.")
    simple = re.fullmatch(
        r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|del|fs)", protein_change
    )
    if simple:
        ref3, pos, alt3 = simple.groups()
        ref1 = AA_3_TO_1.get(ref3, ref3)
        alt1 = AA_3_TO_1.get(alt3, alt3)
        return ref1, pos, alt1, f"{ref1}{pos}{alt1}"
    return None, None, None, protein_change


def write_corrected_tsv(corrected_vcf_path: str, output_tsv_path: str) -> None:
    """
    Write a flat TSV from the corrected VCF: one row per coding ANN entry
    (protein_coding biotype with a valid CDS position).
    """
    vcf = VCF(corrected_vcf_path)
    rows = []

    try:
        for variant in vcf:
            raw_ann = variant.INFO.get("ANN")
            if raw_ann is None:
                continue

            mnv_corrected = bool(variant.INFO.get("MNV_CORRECTED", False))
            mnv_orig_hp = str(variant.INFO.get("MNV_ORIGINAL_HGVS_P") or "")

            # DP4 serialisation (Number=4 Integer → list/array)
            dp4_raw = variant.INFO.get("DP4")
            dp4_str = (
                ",".join(str(int(x)) for x in dp4_raw)
                if dp4_raw is not None else None
            )

            # AF (Number=1 Float → scalar)
            af_raw = variant.INFO.get("AF")
            af_val = float(af_raw) if af_raw is not None else None

            filt_str = "PASS" if variant.FILTER is None else "|".join(variant.FILTERS or [])
            alt_str = ",".join(a for a in variant.ALT if a is not None)

            ann_entries = parse_ann_entries(str(raw_ann))

            for entry in ann_entries:
                biotype = entry[ANN_BIOTYPE_IDX]
                if biotype and biotype != "protein_coding":
                    continue  # skip non-protein-coding transcripts

                cds_pos_val, cds_len = parse_pos_len(entry[ANN_CDS_POS_IDX])
                if cds_pos_val is None:
                    continue  # skip non-coding entries (upstream, intergenic, etc.)

                prot_pos_val, prot_len = parse_pos_len(entry[ANN_PROT_POS_IDX])
                hgvs_p = entry[ANN_HGVS_P_IDX]
                feature_id = entry[ANN_FEATURE_ID_IDX]
                aa_ref, aa_pos, aa_alt, aa_change = _parse_aa_parts(hgvs_p)

                # mnv_original_hgvs_p — per-feature value from the stored string
                orig_hp_for_feature = ""
                if mnv_corrected and mnv_orig_hp:
                    for pair in mnv_orig_hp.split("|"):
                        if ":" in pair and pair.split(":", 1)[0] == feature_id:
                            orig_hp_for_feature = pair.split(":", 1)[1]
                            break

                rows.append({
                    "chrom": variant.CHROM,
                    "pos": variant.POS,
                    "ref": variant.REF,
                    "alt": alt_str,
                    "qual": variant.QUAL,
                    "filter": filt_str,
                    "dp": variant.INFO.get("DP"),
                    "af": af_val,
                    "sb": variant.INFO.get("SB"),
                    "dp4": dp4_str,
                    "annotation": entry[ANN_ANNOTATION_IDX],
                    "impact": entry[ANN_IMPACT_IDX],
                    "gene_name": entry[ANN_GENE_NAME_IDX],
                    "feature_id": feature_id,
                    "hgvs_c": entry[ANN_HGVS_C_IDX],
                    "hgvs_p": hgvs_p,
                    "cds_pos": cds_pos_val,
                    "cds_len": cds_len,
                    "protein_pos": prot_pos_val,
                    "protein_len": prot_len,
                    "aa_ref": aa_ref,
                    "aa_pos": aa_pos,
                    "aa_alt": aa_alt,
                    "aa_change": aa_change,
                    "mnv_corrected": mnv_corrected,
                    "mnv_original_hgvs_p": orig_hp_for_feature,
                })
    finally:
        vcf.close()

    with open(output_tsv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    log.info("Wrote %d coding TSV rows to %s", len(rows), output_tsv_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        metavar="VCF",
        help="Path to snpEff-annotated LoFreq VCF (must contain ANN INFO field).",
    )
    parser.add_argument(
        "--output-vcf",
        required=True,
        metavar="VCF",
        help="Path for the MNV-corrected output VCF.",
    )
    parser.add_argument(
        "--output-tsv",
        metavar="TSV",
        help=(
            "Optional path for a flat TSV of coding annotations from the corrected "
            "VCF.  If omitted, no TSV is produced."
        ),
    )
    parser.add_argument(
        "--min-af",
        type=float,
        default=0.5,
        metavar="FLOAT",
        help=(
            "Minimum allele frequency for two SNVs to be considered in-cis and "
            "eligible for MNV codon correction (default: 0.5). "
            "Both SNVs in a codon group must individually meet this threshold."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable DEBUG-level logging.",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    args = get_args()
    if args.verbose:
        log.setLevel(logging.DEBUG)

    log.info(
        "Pass 1: scanning %s for MNV codon groups (min_af=%.3f) ...",
        args.input, args.min_af,
    )
    codon_groups = collect_codon_groups(args.input, args.min_af)

    mnv_candidates = {k: v for k, v in codon_groups.items() if len(v) >= 2}
    log.info(
        "%d total codon groups (per-transcript), %d with >= 2 high-AF SNVs",
        len(codon_groups), len(mnv_candidates),
    )

    corrections = compute_corrections(mnv_candidates)
    log.info(
        "Computed corrections for %d variant positions across %d unique codon groups",
        len(corrections),
        len({c.codon_group_id for clist in corrections.values() for c in clist}),
    )

    log.info("Pass 2: writing corrected VCF → %s ...", args.output_vcf)
    n_records_corrected = write_corrected_vcf(args.input, args.output_vcf, corrections)
    log.info("%d VCF records updated with MNV corrections", n_records_corrected)

    if args.output_tsv:
        log.info("Writing corrected TSV → %s ...", args.output_tsv)
        write_corrected_tsv(args.output_vcf, args.output_tsv)

    if n_records_corrected == 0 and corrections:
        log.debug(
            "MNV codon groups were resolved but 0 VCF records needed annotation updates "
            "(combined codon changes may already be correctly annotated, or all groups "
            "were synonymous). Check WARNING messages above if this is unexpected."
        )
    elif n_records_corrected == 0 and mnv_candidates:
        log.warning(
            "MNV candidate groups were found but all were skipped. "
            "Check for ambiguous codon reconstruction warnings above."
        )

    log.info("Done.")


if __name__ == "__main__":
    main()
