#!/usr/bin/env python3
"""fix_mnv_codons.py

Post-process a snpEff-annotated VCF to correct amino acid annotations where
two or more adjacent SNVs fall within the same codon (MNV problem).

Supports any snpEff-annotated VCF containing an ANN INFO field — LoFreq-style
VCFs (with INFO/AF, DP, DP4, etc.) and alignment-based consensus VCFs (with
GT=1/1 and no AF) are both handled.  When INFO/AF is absent the variant is
treated as AF=1.0 (fixed consensus call).

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
     - Update ANN field: HGVS.p, annotation type, impact for each affected
       transcript entry.
     - Set INFO flag   MNV_CORRECTED = True
     - Set INFO string MNV_ORIGINAL_HGVS_P = original HGVS.p values

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
from dataclasses import dataclass
from enum import IntEnum
from typing import Dict, List, Optional, Tuple

try:
    from cyvcf2 import VCF, Writer
except ImportError as exc:
    raise SystemExit("cyvcf2 required: conda install -c bioconda cyvcf2") from exc

try:
    from Bio.Data.CodonTable import unambiguous_dna_by_name as _codon_tables
    from Bio.Data.IUPACData import protein_letters_3to1 as _aa_3to1
except ImportError as exc:
    raise SystemExit("biopython required: conda install -c conda-forge biopython") from exc

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

class ANN(IntEnum):
    """snpEff ANN field indices (0-based, pipe-delimited within each comma entry)."""
    ALLELE = 0; ANNOTATION = 1; IMPACT = 2; GENE_NAME = 3; GENE_ID = 4
    FEATURE_TYPE = 5; FEATURE_ID = 6; BIOTYPE = 7; RANK = 8
    HGVS_C = 9; HGVS_P = 10; CDNA_POS = 11; CDS_POS = 12
    PROT_POS = 13; DISTANCE = 14; NOTES = 15

ANN_FIELD_COUNT = len(ANN)

# Amino acid mappings — sourced from Biopython; 'Ter' added for snpEff stop notation
AA_3_TO_1 = {**_aa_3to1, "Ter": "*"}
AA_1_TO_3 = {v: k for k, v in AA_3_TO_1.items()}

# Standard DNA codon table (stop = "*") — sourced from Biopython
_t = _codon_tables["Standard"]
_CODON_TABLE = {**_t.forward_table, **{c: "*" for c in _t.stop_codons}}

# Regex for snpEff hgvs_c: c.263A>G (negative positions = UTR, skipped)
_HGVS_C_RE = re.compile(r"c\.(\d+)([ACGTacgt])>([ACGTacgt])", re.ASCII)
# Regex for snpEff hgvs_p — 3-letter (p.Asn88Ser) and 1-letter (p.N88S)
_HGVS_P_3L_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)", re.ASCII)
_HGVS_P_1L_RE = re.compile(r"p\.([A-Z\*])(\d+)", re.ASCII)

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class CodonSNV:
    """One SNV contributing to a codon group for a specific transcript."""
    chrom: str
    vcf_pos: int
    ref_nuc: str
    alt_nuc: str
    codon_offset: int  # 0-2
    cds_pos: int
    protein_pos: int
    ref_aa_1: Optional[str]
    hgvs_c_ref_nuc: Optional[str]

@dataclass(slots=True)
class AnnCorrection:
    """Correction to apply to a specific (chrom, vcf_pos, feature_id) triplet."""
    feature_id: str
    protein_pos: int
    ref_aa_1: str
    true_alt_aa_1: str
    codon_group_id: str

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def translate_codon(codon: str) -> Optional[str]:
    return _CODON_TABLE.get(codon.upper())

def parse_hgvs_c(hgvs_c: str) -> Tuple[Optional[int], Optional[str], Optional[str]]:
    """Parse 'c.263A>G' -> (cds_pos, ref_nuc, alt_nuc). Returns Nones for UTR/empty."""
    if not hgvs_c:
        return None, None, None
    m = _HGVS_C_RE.match(hgvs_c)
    if not m:
        return None, None, None
    return int(m.group(1)), m.group(2).upper(), m.group(3).upper()

def parse_hgvs_p_ref_aa(hgvs_p: str) -> Optional[str]:
    """Extract 1-letter reference AA from snpEff hgvs_p (3-letter or 1-letter format)."""
    if not hgvs_p:
        return None
    m3 = _HGVS_P_3L_RE.match(hgvs_p)
    if m3:
        return AA_3_TO_1.get(m3.group(1))
    m1 = _HGVS_P_1L_RE.match(hgvs_p)
    if m1:
        ref = m1.group(1)
        return ref if ref in AA_1_TO_3 else None
    return None

def parse_pos_len(field: str) -> Tuple[Optional[int], Optional[str]]:
    """Parse '263/1707' -> (263, '1707'). Returns (None, None) if invalid."""
    if not field or field in ("", "."):
        return None, None
    parts = str(field).split("/", 1)
    try:
        pos = int(parts[0])
    except (ValueError, IndexError):
        return None, None
    return (pos if pos > 0 else None), (parts[1] if len(parts) > 1 else None)

def parse_ann_entries(raw_ann: str) -> List[List[str]]:
    """Split ANN INFO value into list of field-lists (padded to ANN_FIELD_COUNT)."""
    if not raw_ann:
        return []
    entries = []
    for entry in str(raw_ann).split(","):
        parts = entry.split("|")
        parts.extend([""] * max(0, ANN_FIELD_COUNT - len(parts)))
        entries.append(parts[:ANN_FIELD_COUNT])
    return entries

def reassemble_ann(entries: List[List[str]]) -> str:
    return ",".join("|".join(e) for e in entries)

def _is_1letter_hgvs_p(hgvs_p: str) -> bool:
    if not hgvs_p:
        return True  # default to 1-letter (fluscan-nf standard)
    return bool(_HGVS_P_1L_RE.match(hgvs_p) and not _HGVS_P_3L_RE.match(hgvs_p))

def build_hgvs_p(ref_aa_1: str, protein_pos: int, alt_aa_1: str, one_letter: bool = True) -> str:
    """Build HGVS.p string in 1-letter ('p.N88R') or 3-letter ('p.Asn88Arg') format."""
    if one_letter:
        return f"p.{ref_aa_1}{protein_pos}{alt_aa_1}"
    return f"p.{AA_1_TO_3.get(ref_aa_1, ref_aa_1)}{protein_pos}{AA_1_TO_3.get(alt_aa_1, alt_aa_1)}"

def determine_annotation(ref_aa: str, alt_aa: str) -> Tuple[str, str]:
    """Return (annotation_type, impact) for an amino acid change."""
    if alt_aa == "*" and ref_aa != "*":
        return "stop_gained", "HIGH"
    if ref_aa == "*" and alt_aa != "*":
        return "stop_lost", "HIGH"
    if ref_aa == alt_aa:
        return "synonymous_variant", "LOW"
    return "missense_variant", "MODERATE"

# ---------------------------------------------------------------------------
# Codon reconstruction
# ---------------------------------------------------------------------------

def resolve_mnv_codon(snvs: List[CodonSNV]) -> Optional[Tuple[str, str, str, str]]:
    """Resolve ref/mut codons for a group of SNVs in the same codon.
    Returns (ref_codon, mut_codon, ref_aa_1, true_mut_aa_1) or None."""
    known_ref: Dict[int, str] = {}
    known_alt: Dict[int, str] = {}

    for snv in snvs:
        offset, ref_nuc = snv.codon_offset, snv.ref_nuc.upper()
        if snv.hgvs_c_ref_nuc and snv.hgvs_c_ref_nuc != ref_nuc:
            log.warning("hgvs_c ref (%s) != VCF REF (%s) at CDS %d / POS %d",
                        snv.hgvs_c_ref_nuc, ref_nuc, snv.cds_pos, snv.vcf_pos)
        if offset in known_ref and known_ref[offset] != ref_nuc:
            log.warning("Conflicting REF at offset %d: %s vs %s", offset, known_ref[offset], ref_nuc)
            return None
        known_ref[offset] = ref_nuc
        known_alt[offset] = snv.alt_nuc.upper()

    missing = [i for i in range(3) if i not in known_ref]

    # All three codon positions covered
    if not missing:
        ref_codon = "".join(known_ref[i] for i in range(3))
        mut_codon = "".join(known_alt.get(i, known_ref[i]) for i in range(3))
        return ref_codon, mut_codon, translate_codon(ref_codon), translate_codon(mut_codon)

    if len(missing) > 1:
        log.debug("Too many uncovered codon positions (%d)", len(missing))
        return None

    # One uncovered position — infer from reference AA
    ref_aa_1 = next((s.ref_aa_1 for s in snvs if s.ref_aa_1), None)
    if ref_aa_1 is None:
        log.debug("No reference AA for codon with 1 uncovered position")
        return None

    pos = missing[0]
    valid_nucs = []
    for nuc in "ACGT":
        rc = [known_ref.get(i, nuc if i == pos else "N") for i in range(3)]
        if translate_codon("".join(rc)) == ref_aa_1:
            valid_nucs.append(nuc)

    if not valid_nucs:
        log.warning("No nucleotide at offset %d gives ref AA '%s'", pos, ref_aa_1)
        return None

    # Check if all valid nucs produce the same mutant AA
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
            log.debug("Ambiguous ref offset %d (%s->'%s') but unambiguous mut AA '%s'",
                      pos, "/".join(valid_nucs), ref_aa_1, next(iter(mut_aa_set)))
        ref_codon, mut_codon = per_nuc[valid_nucs[0]]
        return ref_codon, mut_codon, ref_aa_1, mut_aa_set.pop()

    log.warning("Ambiguous MNV: offset %d=%s (ref_aa '%s'), mut AA=%s — skipping",
                pos, "/".join(valid_nucs), ref_aa_1, "/".join(sorted(mut_aa_set)))
    return None

# ---------------------------------------------------------------------------
# Pass 1 — collect SNVs and group by codon
# ---------------------------------------------------------------------------

def collect_codon_groups(vcf_path: str, min_af: float) -> Dict[Tuple[str, str, int], List[CodonSNV]]:
    """Read all PASS SNVs with AF >= min_af; group by (chrom, feature_id, codon_number)."""
    codon_groups: Dict[Tuple[str, str, int], List[CodonSNV]] = defaultdict(list)
    vcf = VCF(vcf_path)
    n_total = n_skipped = 0

    try:
        for variant in vcf:
            n_total += 1
            if variant.FILTER is not None or not variant.is_snp:
                n_skipped += 1
                continue

            af_raw = variant.INFO.get("AF")
            if af_raw is None:
                # Consensus-only callers (e.g. paftools.js call with GT=1/1) omit AF.
                af = 1.0
            else:
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

            chrom, pos, ref_nuc = variant.CHROM, variant.POS, variant.REF.upper()
            alts = [a for a in variant.ALT if a is not None]
            if len(alts) != 1:
                continue
            alt_nuc = alts[0].upper()

            for entry in parse_ann_entries(str(raw_ann)):
                cds_pos, _ = parse_pos_len(entry[ANN.CDS_POS])
                if cds_pos is None:
                    continue
                protein_pos, _ = parse_pos_len(entry[ANN.PROT_POS])
                if protein_pos is None:
                    protein_pos = math.ceil(cds_pos / 3)

                codon_number = math.ceil(cds_pos / 3)
                _, hgvs_c_ref, _ = parse_hgvs_c(entry[ANN.HGVS_C])

                codon_groups[(chrom, entry[ANN.FEATURE_ID], codon_number)].append(
                    CodonSNV(chrom=chrom, vcf_pos=pos, ref_nuc=ref_nuc, alt_nuc=alt_nuc,
                             codon_offset=(cds_pos - 1) % 3, cds_pos=cds_pos,
                             protein_pos=protein_pos,
                             ref_aa_1=parse_hgvs_p_ref_aa(entry[ANN.HGVS_P]),
                             hgvs_c_ref_nuc=hgvs_c_ref)
                )
    finally:
        vcf.close()

    log.debug("Pass 1: %d records, %d skipped", n_total, n_skipped)
    return dict(codon_groups)

# ---------------------------------------------------------------------------
# Compute corrections from codon groups
# ---------------------------------------------------------------------------

def compute_corrections(
    codon_groups: Dict[Tuple[str, str, int], List[CodonSNV]],
) -> Dict[Tuple[str, int], List[AnnCorrection]]:
    """For each codon group with >=2 SNVs, resolve the true combined AA change."""
    corrections: Dict[Tuple[str, int], List[AnnCorrection]] = defaultdict(list)
    n_groups = n_corrected = n_skipped = 0

    for (chrom, feature_id, codon_number), snvs in codon_groups.items():
        if len(snvs) < 2:
            continue
        n_groups += 1

        # Guard: consistent codon numbers and no duplicate offsets
        if len({math.ceil(s.cds_pos / 3) for s in snvs}) > 1:
            log.warning("Mixed codon numbers in group (%s, %s, %d)", chrom, feature_id, codon_number)
            n_skipped += 1
            continue
        offsets = [s.codon_offset for s in snvs]
        if len(set(offsets)) < len(offsets):
            log.warning("Duplicate offsets %s in group (%s, %s, %d)", offsets, chrom, feature_id, codon_number)
            n_skipped += 1
            continue

        resolved = resolve_mnv_codon(snvs)
        if resolved is None:
            n_skipped += 1
            continue

        ref_codon, mut_codon, ref_aa_1, true_alt_aa_1 = resolved
        if true_alt_aa_1 is None:
            log.warning("Cannot translate mut codon '%s' for (%s, %s, %d)",
                        mut_codon, chrom, feature_id, codon_number)
            n_skipped += 1
            continue

        protein_pos = snvs[0].protein_pos
        codon_group_id = f"{chrom}:{feature_id}:codon{codon_number}"
        log.info("MNV %s | %s%d->%s | %s->%s | pos=%s",
                 codon_group_id, ref_aa_1, protein_pos, true_alt_aa_1,
                 ref_codon, mut_codon, sorted(s.vcf_pos for s in snvs))

        correction = AnnCorrection(feature_id=feature_id, protein_pos=protein_pos,
                                   ref_aa_1=ref_aa_1, true_alt_aa_1=true_alt_aa_1,
                                   codon_group_id=codon_group_id)
        for snv in snvs:
            corrections[(chrom, snv.vcf_pos)].append(correction)
        n_corrected += 1

    log.info("Codon groups: %d found, %d correctable, %d skipped", n_groups, n_corrected, n_skipped)
    return dict(corrections)

# ---------------------------------------------------------------------------
# Apply corrections to ANN field
# ---------------------------------------------------------------------------

def apply_corrections_to_ann(raw_ann: str, corrections: List[AnnCorrection]) -> str:
    """Update HGVS.p, annotation type, and impact in matching ANN entries."""
    entries = parse_ann_entries(raw_ann)
    by_feature = {c.feature_id: c for c in corrections}

    for entry in entries:
        corr = by_feature.get(entry[ANN.FEATURE_ID])
        if corr is None:
            continue
        one_letter = _is_1letter_hgvs_p(entry[ANN.HGVS_P])
        entry[ANN.HGVS_P] = build_hgvs_p(corr.ref_aa_1, corr.protein_pos, corr.true_alt_aa_1, one_letter)
        entry[ANN.ANNOTATION], entry[ANN.IMPACT] = determine_annotation(corr.ref_aa_1, corr.true_alt_aa_1)

    return reassemble_ann(entries)

# ---------------------------------------------------------------------------
# Pass 2 — write corrected VCF
# ---------------------------------------------------------------------------

def write_corrected_vcf(
    vcf_path: str, output_vcf_path: str,
    corrections: Dict[Tuple[str, int], List[AnnCorrection]],
) -> int:
    """Emit all records to a new VCF, correcting ANN fields for MNV groups.
    Returns the count of corrected records."""
    vcf = VCF(vcf_path)
    vcf.add_info_to_header({
        "ID": "MNV_CORRECTED", "Number": "0", "Type": "Flag",
        "Description": ("Variant is part of an MNV codon group. ANN HGVS.p, annotation type, "
                        "and impact corrected for true combined amino acid effect (fix_mnv_codons.py)."),
    })
    vcf.add_info_to_header({
        "ID": "MNV_ORIGINAL_HGVS_P", "Number": "1", "Type": "String",
        "Description": "Original snpEff HGVS.p before MNV correction (featureID:hgvs_p pipe-separated).",
    })

    writer = Writer(output_vcf_path, vcf)
    n_corrected = 0

    try:
        for variant in vcf:
            correction_list = corrections.get((variant.CHROM, variant.POS))
            if correction_list:
                raw_ann = variant.INFO.get("ANN")
                if raw_ann is not None:
                    original_ann = str(raw_ann)
                    affected = {c.feature_id for c in correction_list}
                    orig_hgvs = "|".join(
                        f"{e[ANN.FEATURE_ID]}:{e[ANN.HGVS_P]}"
                        for e in parse_ann_entries(original_ann)
                        if e[ANN.FEATURE_ID] in affected and e[ANN.HGVS_P]
                    )
                    corrected_ann = apply_corrections_to_ann(original_ann, correction_list)
                    if corrected_ann != original_ann:
                        variant.INFO["ANN"] = corrected_ann
                        variant.INFO["MNV_CORRECTED"] = True
                        if orig_hgvs:
                            variant.INFO["MNV_ORIGINAL_HGVS_P"] = orig_hgvs
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

def _parse_aa_parts(hgvs_p: str) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    """Parse 'p.Asn88Arg' -> (ref_1, pos, alt_1, change_str). Returns Nones if unparseable."""
    if not hgvs_p:
        return None, None, None, None
    protein_change = hgvs_p.removeprefix("p.")
    m = re.fullmatch(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|del|fs)", protein_change)
    if m:
        ref3, pos, alt3 = m.groups()
        ref1, alt1 = AA_3_TO_1.get(ref3, ref3), AA_3_TO_1.get(alt3, alt3)
        return ref1, pos, alt1, f"{ref1}{pos}{alt1}"
    return None, None, None, protein_change

def write_corrected_tsv(corrected_vcf_path: str, output_tsv_path: str) -> None:
    """Write flat TSV from corrected VCF: one row per coding ANN entry."""
    vcf = VCF(corrected_vcf_path)
    rows = []

    try:
        for variant in vcf:
            raw_ann = variant.INFO.get("ANN")
            if raw_ann is None:
                continue

            mnv_corrected = bool(variant.INFO.get("MNV_CORRECTED", False))
            mnv_orig_hp = str(variant.INFO.get("MNV_ORIGINAL_HGVS_P") or "")
            dp4_raw = variant.INFO.get("DP4")
            dp4_str = ",".join(str(int(x)) for x in dp4_raw) if dp4_raw is not None else None
            af_raw = variant.INFO.get("AF")
            af_val = float(af_raw) if af_raw is not None else None
            filt_str = "PASS" if variant.FILTER is None else "|".join(variant.FILTERS or [])
            alt_str = ",".join(a for a in variant.ALT if a is not None)

            for entry in parse_ann_entries(str(raw_ann)):
                biotype = entry[ANN.BIOTYPE]
                if biotype and biotype != "protein_coding":
                    continue
                cds_pos_val, cds_len = parse_pos_len(entry[ANN.CDS_POS])
                if cds_pos_val is None:
                    continue
                prot_pos_val, prot_len = parse_pos_len(entry[ANN.PROT_POS])
                hgvs_p = entry[ANN.HGVS_P]
                feature_id = entry[ANN.FEATURE_ID]
                aa_ref, aa_pos, aa_alt, aa_change = _parse_aa_parts(hgvs_p)

                # Extract per-feature original hgvs_p
                orig_hp_for_feature = ""
                if mnv_corrected and mnv_orig_hp:
                    for pair in mnv_orig_hp.split("|"):
                        if ":" in pair and pair.split(":", 1)[0] == feature_id:
                            orig_hp_for_feature = pair.split(":", 1)[1]
                            break

                rows.append({
                    "chrom": variant.CHROM, "pos": variant.POS,
                    "ref": variant.REF, "alt": alt_str,
                    "qual": variant.QUAL, "filter": filt_str,
                    "dp": variant.INFO.get("DP"), "af": af_val,
                    "sb": variant.INFO.get("SB"), "dp4": dp4_str,
                    "annotation": entry[ANN.ANNOTATION], "impact": entry[ANN.IMPACT],
                    "gene_name": entry[ANN.GENE_NAME], "feature_id": feature_id,
                    "hgvs_c": entry[ANN.HGVS_C], "hgvs_p": hgvs_p,
                    "cds_pos": cds_pos_val, "cds_len": cds_len,
                    "protein_pos": prot_pos_val, "protein_len": prot_len,
                    "aa_ref": aa_ref, "aa_pos": aa_pos,
                    "aa_alt": aa_alt, "aa_change": aa_change,
                    "mnv_corrected": mnv_corrected,
                    "mnv_original_hgvs_p": orig_hp_for_feature,
                })
    finally:
        vcf.close()

    with open(output_tsv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote %d coding TSV rows -> %s", len(rows), output_tsv_path)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, metavar="VCF",
                   help="snpEff-annotated LoFreq VCF (must contain ANN INFO field).")
    p.add_argument("--output-vcf", required=True, metavar="VCF",
                   help="Path for MNV-corrected output VCF.")
    p.add_argument("--output-tsv", metavar="TSV",
                   help="Optional flat TSV of coding annotations from corrected VCF.")
    p.add_argument("--min-af", type=float, default=0.5, metavar="FLOAT",
                   help=("Minimum AF for MNV codon correction eligibility (default: 0.5). "
                         "Variants with no INFO/AF are treated as AF=1.0."))
    p.add_argument("--verbose", action="store_true", help="Enable DEBUG logging.")
    return p.parse_args()

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    args = get_args()
    if args.verbose:
        log.setLevel(logging.DEBUG)

    log.info("Pass 1: scanning %s (min_af=%.3f)", args.input, args.min_af)
    codon_groups = collect_codon_groups(args.input, args.min_af)
    mnv_candidates = {k: v for k, v in codon_groups.items() if len(v) >= 2}
    log.info("%d codon groups, %d with >=2 high-AF SNVs", len(codon_groups), len(mnv_candidates))

    corrections = compute_corrections(mnv_candidates)
    log.info("%d positions to correct across %d codon groups",
             len(corrections), len({c.codon_group_id for cl in corrections.values() for c in cl}))

    log.info("Pass 2: writing corrected VCF -> %s", args.output_vcf)
    n = write_corrected_vcf(args.input, args.output_vcf, corrections)
    log.info("%d VCF records updated", n)

    if args.output_tsv:
        write_corrected_tsv(args.output_vcf, args.output_tsv)

    if n == 0 and mnv_candidates:
        log.warning("MNV candidates found but 0 records updated — check warnings above")

if __name__ == "__main__":
    main()
