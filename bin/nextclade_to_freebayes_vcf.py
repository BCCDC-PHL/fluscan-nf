#!/usr/bin/env python3
"""Convert nextclade JSON or TSV output to a freebayes-format VCF.

Synthesises depth-1 (consensus-quality) freebayes INFO/FORMAT fields so that
reformat_vcf_cyvcf2.py can parse the output unchanged.

For JSON input, emits snpEff-compatible ANN= annotations directly in the INFO
column, sourced from Nextclade's per-CDS amino-acid data.  This bypasses the
need for snpEff and correctly handles Influenza's overlapping reading frames
(M1/M2, NS1/NEP, PB1/PB1-F2).

TSV input does not provide per-CDS resolution of overlapping ORFs and therefore
does not emit ANN annotations.

Requires the reference FASTA used during nextclade run to resolve the contig ID
and anchor bases for insertion/deletion records.

Pipeline position: after nextclade run, before reformat_vcf_cyvcf2.py.

Critical evaluation — pandas / cyvcf2:
  pandas:  csv.DictReader is sufficient for row-by-row TSV streaming; pandas
           adds dependency weight with no benefit here.
  cyvcf2:  Writer requires a parsed template header; creating Variant objects
           from scratch is not supported. Plain text output is cleaner and has
           no extra dependencies.
"""

import argparse
import csv
import json
import logging
import re
import sys
from datetime import date
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Freebayes header definitions (verbatim from paftools_to_freebayes_vcf.py)
# ---------------------------------------------------------------------------

_FREEBAYES_INFO = [
    ("NS",      "1", "Integer", "Number of samples with data"),
    ("DP",      "1", "Integer", "Total read depth at the locus"),
    ("DPB",     "1", "Float",   "Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype"),
    ("AC",      "A", "Integer", "Total number of alternate alleles in called genotypes"),
    ("AN",      "1", "Integer", "Total number of alleles in called genotypes"),
    ("AF",      "A", "Float",   "Estimated allele frequency in the range (0,1]"),
    ("RO",      "1", "Integer", "Count of full observations of the reference haplotype."),
    ("AO",      "A", "Integer", "Count of full observations of this alternate haplotype."),
    ("PRO",     "1", "Float",   "Reference allele observation count, with partial observations recorded fractionally"),
    ("PAO",     "A", "Float",   "Alternate allele observations, with partial observations recorded fractionally"),
    ("QR",      "1", "Integer", "Reference allele quality sum in phred"),
    ("QA",      "A", "Integer", "Alternate allele quality sum in phred"),
    ("PQR",     "1", "Float",   "Reference allele quality sum in phred for partial observations"),
    ("PQA",     "A", "Float",   "Alternate allele quality sum in phred for partial observations"),
    ("SRF",     "1", "Integer", "Number of reference observations on the forward strand"),
    ("SRR",     "1", "Integer", "Number of reference observations on the reverse strand"),
    ("SAF",     "A", "Integer", "Number of alternate observations on the forward strand"),
    ("SAR",     "A", "Integer", "Number of alternate observations on the reverse strand"),
    ("SRP",     "1", "Float",   "Strand balance probability for the reference allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SRF and SRR given E(SRF/SRR) ~ 0.5, derived using Hoeffding's inequality"),
    ("SAP",     "A", "Float",   "Strand balance probability for the alternate allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SAF and SAR given E(SAF/SAR) ~ 0.5, derived using Hoeffding's inequality"),
    ("AB",      "A", "Float",   "Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous"),
    ("ABP",     "A", "Float",   "Allele balance probability at heterozygous sites: Phred-scaled upper-bounds estimate of the probability of observing the deviation between ABR and ABA given E(ABR/ABA) ~ 0.5, derived using Hoeffding's inequality"),
    ("RUN",     "A", "Integer", "Run length: the number of consecutive repeats of the alternate allele in the reference genome"),
    ("RPP",     "A", "Float",   "Read Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality"),
    ("RPPR",    "1", "Float",   "Read Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality"),
    ("RPL",     "A", "Float",   "Reads Placed Left: number of reads supporting the alternate balanced to the left (5') of the alternate allele"),
    ("RPR",     "A", "Float",   "Reads Placed Right: number of reads supporting the alternate balanced to the right (3') of the alternate allele"),
    ("EPP",     "A", "Float",   "End Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality"),
    ("EPPR",    "1", "Float",   "End Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality"),
    ("DPRA",    "A", "Float",   "Alternate allele depth ratio.  Ratio between depth in samples with each called alternate allele and those without."),
    ("ODDS",    "1", "Float",   "The log odds ratio of the best genotype combination to the second-best."),
    ("GTI",     "1", "Integer", "Number of genotyping iterations required to reach convergence or bailout."),
    ("TYPE",    "A", "String",  "The type of allele, either snp, mnp, ins, del, or complex."),
    ("CIGAR",   "A", "String",  "The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR."),
    ("NUMALT",  "1", "Integer", "Number of unique non-reference alleles in called genotypes at this position."),
    ("MEANALT", "A", "Float",   "Mean number of unique non-reference allele observations per sample with the corresponding alternate alleles."),
    ("LEN",     "A", "Integer", "allele length"),
    ("MQM",     "A", "Float",   "Mean mapping quality of observed alternate alleles"),
    ("MQMR",    "1", "Float",   "Mean mapping quality of observed reference alleles"),
    ("PAIRED",  "A", "Float",   "Proportion of observed alternate alleles which are supported by properly paired read fragments"),
    ("PAIREDR", "1", "Float",   "Proportion of observed reference alleles which are supported by properly paired read fragments"),
]

_FREEBAYES_FORMAT = [
    ("GT",  "1", "String",  "Genotype"),
    ("GQ",  "1", "Float",   "Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype"),
    ("GL",  "G", "Float",   "Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy"),
    ("DP",  "1", "Integer", "Read Depth"),
    ("AD",  "R", "Integer", "Number of observation for each allele"),
    ("RO",  "1", "Integer", "Reference allele observation count"),
    ("QR",  "1", "Integer", "Sum of quality of the reference observations"),
    ("AO",  "A", "Integer", "Alternate allele observation count"),
    ("QA",  "A", "Integer", "Sum of quality of the alternate observations"),
]

# ---------------------------------------------------------------------------
# Amino acid 1-letter to 3-letter code map
# ---------------------------------------------------------------------------

_AA1_TO_AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter", "X": "Xaa", "-": "del",
}

# ---------------------------------------------------------------------------
# Variant classification (verbatim from paftools_to_freebayes_vcf.py)
# ---------------------------------------------------------------------------

def classify(ref: str, alt: str) -> tuple[str, str, int]:
    """Return (type_str, cigar_str, len_val) for a REF/ALT pair."""
    if len(ref) == 1 and len(alt) == 1:
        return "snp", "1X", 1
    if len(ref) < len(alt) and alt.startswith(ref):
        ins = len(alt) - len(ref)
        return "ins", f"{len(ref)}M{ins}I", ins
    if len(ref) > len(alt) and ref.startswith(alt):
        dele = len(ref) - len(alt)
        return "del", f"{len(alt)}M{dele}D", dele
    if len(ref) == len(alt):
        return "mnp", f"{len(ref)}X", len(ref)
    diff = abs(len(alt) - len(ref))
    direction = "I" if len(alt) > len(ref) else "D"
    return "complex", f"{min(len(ref), len(alt))}M{diff}{direction}", max(len(ref), len(alt))

# ---------------------------------------------------------------------------
# Reference FASTA
# ---------------------------------------------------------------------------

def load_reference(ref_path: str) -> tuple[str, str]:
    """Return (chrom_id, sequence_uppercase) from a single-record FASTA."""
    rec = SeqIO.read(ref_path, "fasta")
    return rec.id, str(rec.seq).upper()

# ---------------------------------------------------------------------------
# VCF header
# ---------------------------------------------------------------------------

def build_header(chrom: str, ref_len: int, ref_path: str, sample: str) -> str:
    lines = [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        f"##fileDate={date.today().strftime('%Y%m%d')}",
        "##source=nextclade_to_freebayes_vcf.py",
        f"##reference={ref_path}",
        f"##contig=<ID={chrom},length={ref_len}>",
        "##phasing=none",
    ]
    for fid, number, typ, desc in _FREEBAYES_INFO:
        lines.append(f'##INFO=<ID={fid},Number={number},Type={typ},Description="{desc}">')
    lines.append(
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: '
        "'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type"
        " | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos /"
        " cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance |"
        " ERRORS / WARNINGS / INFO'\">"
    )
    for fid, number, typ, desc in _FREEBAYES_FORMAT:
        lines.append(f'##FORMAT=<ID={fid},Number={number},Type={typ},Description="{desc}">')
    lines.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}")
    return "\n".join(lines)

# ---------------------------------------------------------------------------
# Single variant row
# ---------------------------------------------------------------------------

def _info(type_str: str, cigar: str, length: int, qual: int) -> str:
    return (
        f"NS=1;DP=1;DPB=1.0;AC=1;AN=1;AF=1.0;RO=0;AO=1;PRO=0.0;PAO=0.0;"
        f"QR=0;QA={qual};PQR=0.0;PQA=0.0;SRF=0;SRR=0;SAF=0;SAR=0;SRP=0.0;SAP=0.0;"
        f"AB=0.0;ABP=0.0;RUN=1;RPP=0.0;RPPR=0.0;RPL=0.0;RPR=0.0;EPP=0.0;EPPR=0.0;"
        f"DPRA=0.0;ODDS=0.0;GTI=0;TYPE={type_str};CIGAR={cigar};NUMALT=1;MEANALT=1.0;"
        f"LEN={length};MQM=60.0;MQMR=60.0;PAIRED=0.0;PAIREDR=0.0"
    )


def _sample(qual: int) -> str:
    gl = f"{-qual / 10:.3f},0"
    return f"1:1:0,1:0:0:1:{qual}:{gl}"


def variant_row(chrom: str, pos: int, ref: str, alt: str, qual: int, ann: str = "") -> str:
    type_str, cigar, length = classify(ref, alt)
    info = _info(type_str, cigar, length, qual)
    if ann:
        info = f"{info};{ann}"
    return (
        f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\tPASS"
        f"\t{info}"
        f"\tGT:DP:AD:RO:QR:AO:QA:GL\t{_sample(qual)}"
    )

# ---------------------------------------------------------------------------
# Indel anchor helpers
# ---------------------------------------------------------------------------

def _deletion_alleles(ref_seq: str, begin: int, end: int) -> tuple[int, str, str]:
    """Return (vcf_pos_1based, ref_allele, alt_allele) for a deletion."""
    if begin > 0:
        anchor = begin - 1
        return anchor + 1, ref_seq[anchor] + ref_seq[begin:end], ref_seq[anchor]
    # deletion at position 0: use trailing anchor
    return 1, ref_seq[begin:end] + ref_seq[end], ref_seq[end]


def _insertion_alleles(ref_seq: str, pos: int, ins_seq: str) -> tuple[int, str, str]:
    """Return (vcf_pos_1based, ref_allele, alt_allele) for an insertion."""
    anchor = pos  # 0-based position of anchor base
    return anchor + 1, ref_seq[anchor], ref_seq[anchor] + ins_seq

# ---------------------------------------------------------------------------
# CDS index helpers
# ---------------------------------------------------------------------------

def _sample_to_ref(s: int, del_ranges: list[tuple[int, int]]) -> int:
    """Convert a 0-based sample position to the corresponding 0-based reference position.

    annotation.genes[].range is in the sample's coordinate space (after deletions).
    Nucleotide positions in substitutions/nucToAaMuts are in the reference coordinate
    space.  This function bridges the gap so the CDS index uses reference coordinates.

    del_ranges must be sorted by start position.
    """
    r = s
    for _ in range(200):
        deleted = 0
        for b, e in del_ranges:
            if b >= r:
                break
            deleted += min(e, r) - b
        new_r = s + deleted
        if new_r == r:
            break
        r = new_r
    # Advance past any deleted range that starts at r (deleted positions have
    # no sample equivalent; the next non-deleted position is the correct answer).
    changed = True
    while changed:
        changed = False
        for b, e in del_ranges:
            if b <= r < e:
                r = e
                changed = True
                break
    return r


def _build_cds_index(annotation: dict, deletions: list[dict], ref_seq: str) -> list[dict]:
    """Parse annotation.genes into a flat list of CDS segment records.

    annotation.genes[].range is in sample coordinates; this function converts
    them to reference coordinates using the deletions list so that CDS positions
    are directly comparable with substitutions.pos / nucToAaMuts keys.

    Nextclade CDS annotations sometimes include 5' UTR before the ATG.  The
    'atg_begin' field stores the first ATG at or after 'begin' (up to 60 nt
    ahead) and is used for HGVS nucleotide position calculations so that HGVS.c
    numbering is consistent with the AA positions reported by nucToAaMuts.
    """
    del_ranges = sorted(
        [(d["range"]["begin"], d["range"]["end"]) for d in deletions]
    )
    index = []
    for gene in annotation.get("genes", []):
        for cds in gene.get("cdses", []):
            for seg in cds.get("segments", []):
                ref_begin = _sample_to_ref(seg["range"]["begin"], del_ranges)
                ref_end   = _sample_to_ref(seg["range"]["end"],   del_ranges)
                strand    = seg.get("strand", "+")
                length    = ref_end - ref_begin
                aa_length = length // 3

                # Find the true ATG start for HGVS.c numbering.
                # For + strand, search forward from ref_begin up to 60 nt.
                # For - strand, fall back to ref_begin (rare in flu).
                atg_begin = ref_begin
                if strand == "+":
                    for i in range(ref_begin, min(ref_begin + 60, ref_end - 2)):
                        if ref_seq[i:i + 3].upper() == "ATG":
                            atg_begin = i
                            break

                index.append({
                    "gene_name": gene["name"],
                    "gene_id":   gene["id"],
                    "cds_id":    cds["id"],
                    "cds_name":  cds["name"],
                    "begin":     ref_begin,
                    "end":       ref_end,
                    "atg_begin": atg_begin,
                    "strand":    strand,
                    "length":    length,
                    "aa_length": aa_length,
                })
    return index


def _locate_cds(cds_index: list[dict], pos_0: int) -> list[dict]:
    """Return every CDS record whose range contains 0-based position pos_0."""
    return [c for c in cds_index if c["begin"] <= pos_0 < c["end"]]


def _cds_local_pos(cds: dict, pos_0: int, use_atg: bool = False) -> int:
    """Return 1-based CDS-local nucleotide position for genomic 0-based pos_0.

    When use_atg=True, numbering is relative to the ATG start codon (correct
    for HGVS.c).  When False, numbering is relative to the annotation start
    (used for CDS boundary arithmetic).
    """
    if cds["strand"] == "-":
        return cds["end"] - pos_0
    start = cds["atg_begin"] if use_atg else cds["begin"]
    return pos_0 - start + 1


def _translate_codon(codon: str) -> str:
    try:
        return str(Seq(codon).translate())
    except Exception:
        return "X"


def _synonymous_aa_change(
    ref_seq: str, cds: dict, pos_0: int, ref_nt: str, alt_nt: str
) -> tuple[str, str, int]:
    """Return (refAa, altAa, aa_pos_1based) for a substitution at pos_0.

    Computes codon translation from the reference sequence.  If the result
    turns out to be non-synonymous (refAa != altAa), the caller should emit
    a missense_variant annotation rather than synonymous_variant.
    """
    if cds["strand"] == "+":
        atg = cds["atg_begin"]
        cds_local_0 = pos_0 - atg
        codon_idx   = cds_local_0 // 3
        codon_start = atg + codon_idx * 3
        if codon_start + 3 > len(ref_seq):
            return "X", "X", codon_idx + 1
        ref_codon = ref_seq[codon_start:codon_start + 3]
        offset    = cds_local_0 % 3
        alt_codon = ref_codon[:offset] + alt_nt + ref_codon[offset + 1:]
        return _translate_codon(ref_codon), _translate_codon(alt_codon), codon_idx + 1
    else:
        # - strand: codon is the reverse complement of the genomic triplet
        cds_local_0 = cds["end"] - 1 - pos_0
        codon_idx   = cds_local_0 // 3
        codon_start = cds["end"] - (codon_idx + 1) * 3
        if codon_start < 0:
            return "X", "X", codon_idx + 1
        ref_codon_fwd = ref_seq[codon_start:codon_start + 3]
        rc_ref = str(Seq(ref_codon_fwd).reverse_complement())
        offset_fwd    = pos_0 - codon_start
        alt_codon_fwd = ref_codon_fwd[:offset_fwd] + alt_nt + ref_codon_fwd[offset_fwd + 1:]
        rc_alt = str(Seq(alt_codon_fwd).reverse_complement())
        return _translate_codon(rc_ref), _translate_codon(rc_alt), codon_idx + 1


def _aa_pos_to_nuc_range(cds: dict, aa_pos_0: int) -> tuple[int, int]:
    """Return (begin_0, end_0) genomic nuc range for a 0-based AA position in CDS."""
    if cds["strand"] == "+":
        begin = cds["begin"] + aa_pos_0 * 3
        return begin, begin + 3
    else:
        end = cds["end"] - aa_pos_0 * 3
        return end - 3, end

# ---------------------------------------------------------------------------
# HGVS formatters
# ---------------------------------------------------------------------------

def _format_hgvs_c_snp(cds_local_1: int, ref_nt: str, alt_nt: str) -> str:
    return f"c.{cds_local_1}{ref_nt}>{alt_nt}"


def _format_hgvs_c_del(start_1: int, end_1: int) -> str:
    if start_1 == end_1:
        return f"c.{start_1}del"
    return f"c.{start_1}_{end_1}del"


def _format_hgvs_c_ins(anchor_1: int, ins_seq: str) -> str:
    return f"c.{anchor_1}_{anchor_1 + 1}ins{ins_seq}"


def _format_hgvs_p_sub(ref_aa: str, aa_pos_1: int, alt_aa: str) -> str:
    ref3 = _AA1_TO_AA3.get(ref_aa, ref_aa)
    alt3 = _AA1_TO_AA3.get(alt_aa, alt_aa)
    return f"p.{ref3}{aa_pos_1}{alt3}"


def _format_hgvs_p_syn(ref_aa: str, aa_pos_1: int) -> str:
    ref3 = _AA1_TO_AA3.get(ref_aa, ref_aa)
    return f"p.{ref3}{aa_pos_1}="


def _format_hgvs_p_del(aa_changes: list[tuple[str, int]]) -> str:
    """Return HGVS.p for an in-frame deletion given list of (refAa, aa_pos_1)."""
    if not aa_changes:
        return ""
    sorted_ch = sorted(aa_changes, key=lambda x: x[1])
    first = sorted_ch[0]
    last  = sorted_ch[-1]
    first3 = _AA1_TO_AA3.get(first[0], first[0])
    if len(sorted_ch) == 1:
        return f"p.{first3}{first[1]}del"
    last3 = _AA1_TO_AA3.get(last[0], last[0])
    return f"p.{first3}{first[1]}_{last3}{last[1]}del"


def _format_hgvs_p_ins(aa_pos_1: int, ins_aa_seq: str) -> str:
    ins3 = "".join(_AA1_TO_AA3.get(aa, aa) for aa in ins_aa_seq)
    return f"p.{aa_pos_1}_{aa_pos_1 + 1}ins{ins3}"


def _format_hgvs_p_fs(ref_aa: str, aa_pos_1: int) -> str:
    ref3 = _AA1_TO_AA3.get(ref_aa, ref_aa)
    return f"p.{ref3}{aa_pos_1}fs"

# ---------------------------------------------------------------------------
# Consequence classification
# ---------------------------------------------------------------------------

def _classify_aa_consequence(ref_aa: str, alt_aa: str, aa_pos_1: int) -> tuple[str, str]:
    """Return (SO_term, impact) for a substitution-derived AA change."""
    if ref_aa == alt_aa:
        return "synonymous_variant", "LOW"
    if alt_aa == "*":
        return "stop_gained", "HIGH"
    if ref_aa == "*":
        return "stop_lost", "HIGH"
    if aa_pos_1 == 1 and ref_aa == "M":
        return "start_lost", "HIGH"
    return "missense_variant", "MODERATE"

# ---------------------------------------------------------------------------
# ANN entry builder
# ---------------------------------------------------------------------------

def _build_ann_entry(
    allele:            str,
    annotation:        str,
    impact:            str,
    gene_name:         str,
    gene_id:           str,
    feature_id:        str,
    transcript_biotype: str,
    hgvs_c:            str,
    hgvs_p:            str,
    cds_pos:           int,
    cds_len:           int,
    aa_pos:            int,
    aa_len:            int,
    distance:          str = "",
    feature_type:      str = "transcript",
) -> str:
    cds_str  = f"{cds_pos}/{cds_len}" if cds_len > 0 else ""
    aa_str   = f"{aa_pos}/{aa_len}"   if aa_len  > 0 else ""
    rank_str = "1/1" if feature_type == "transcript" else ""
    return "|".join([
        allele,
        annotation,
        impact,
        gene_name,
        gene_id,
        feature_type,
        feature_id,
        transcript_biotype,
        rank_str,
        hgvs_c,
        hgvs_p,
        cds_str,   # cDNA.pos/cDNA.length
        cds_str,   # CDS.pos/CDS.length
        aa_str,
        distance,
        "",        # ERRORS/WARNINGS/INFO
    ])

# ---------------------------------------------------------------------------
# Per-variant ANN field builders
# ---------------------------------------------------------------------------

def _ann_field_for_sub(
    pos_0:      int,
    ref_nt:     str,
    alt_nt:     str,
    nuc_to_aa:  dict,
    cds_index:  list[dict],
    ref_seq:    str,
) -> str:
    """Return 'ANN=...' INFO sub-field for a nucleotide substitution, or ''."""
    entries: list[str] = []

    aa_changes = nuc_to_aa.get(str(pos_0), [])
    if aa_changes:
        for change in aa_changes:
            cds_name = change["cdsName"]
            aa_pos_0 = change["pos"]
            ref_aa   = change["refAa"]
            qry_aa   = change["qryAa"]
            aa_pos_1 = aa_pos_0 + 1

            cds_list = [c for c in cds_index if c["cds_name"] == cds_name]
            if not cds_list:
                continue
            cds = cds_list[0]

            # use_atg=True: HGVS.c numbered from ATG, consistent with AA pos
            cds_local_1            = _cds_local_pos(cds, pos_0, use_atg=True)
            annotation, impact     = _classify_aa_consequence(ref_aa, qry_aa, aa_pos_1)
            hgvs_c                 = _format_hgvs_c_snp(cds_local_1, ref_nt, alt_nt)
            if annotation == "synonymous_variant":
                hgvs_p = _format_hgvs_p_syn(ref_aa, aa_pos_1)
            else:
                hgvs_p = _format_hgvs_p_sub(ref_aa, aa_pos_1, qry_aa)

            entries.append(_build_ann_entry(
                allele=alt_nt,
                annotation=annotation,
                impact=impact,
                gene_name=cds["gene_name"],
                gene_id=cds["gene_id"],
                feature_id=cds["cds_id"],
                transcript_biotype="protein_coding",
                hgvs_c=hgvs_c,
                hgvs_p=hgvs_p,
                cds_pos=cds_local_1,
                cds_len=cds["length"],
                aa_pos=aa_pos_1,
                aa_len=cds["aa_length"],
            ))
    else:
        cds_hits = _locate_cds(cds_index, pos_0)
        if cds_hits:
            for cds in cds_hits:
                # Positions in [begin, atg_begin) are 5' UTR within the gene annotation.
                if pos_0 < cds["atg_begin"]:
                    entries.append(_build_ann_entry(
                        allele=alt_nt,
                        annotation="5_prime_UTR_variant",
                        impact="MODIFIER",
                        gene_name=cds["gene_name"],
                        gene_id=cds["gene_id"],
                        feature_id=cds["cds_id"],
                        transcript_biotype="protein_coding",
                        hgvs_c=f"c.-{cds['atg_begin'] - pos_0}{ref_nt}>{alt_nt}",
                        hgvs_p="",
                        cds_pos=0,
                        cds_len=cds["length"],
                        aa_pos=0,
                        aa_len=cds["aa_length"],
                    ))
                    continue
                ref_aa, alt_aa, aa_pos_1 = _synonymous_aa_change(
                    ref_seq, cds, pos_0, ref_nt, alt_nt
                )
                cds_local_1        = _cds_local_pos(cds, pos_0, use_atg=True)
                annotation, impact = _classify_aa_consequence(ref_aa, alt_aa, aa_pos_1)
                hgvs_c             = _format_hgvs_c_snp(cds_local_1, ref_nt, alt_nt)
                if annotation == "synonymous_variant":
                    hgvs_p = _format_hgvs_p_syn(ref_aa, aa_pos_1)
                else:
                    hgvs_p = _format_hgvs_p_sub(ref_aa, aa_pos_1, alt_aa)

                entries.append(_build_ann_entry(
                    allele=alt_nt,
                    annotation=annotation,
                    impact=impact,
                    gene_name=cds["gene_name"],
                    gene_id=cds["gene_id"],
                    feature_id=cds["cds_id"],
                    transcript_biotype="protein_coding",
                    hgvs_c=hgvs_c,
                    hgvs_p=hgvs_p,
                    cds_pos=cds_local_1,
                    cds_len=cds["length"],
                    aa_pos=aa_pos_1,
                    aa_len=cds["aa_length"],
                ))
        else:
            entries.append(_build_ann_entry(
                allele=alt_nt,
                annotation="intergenic_region",
                impact="MODIFIER",
                gene_name="",
                gene_id="",
                feature_id="",
                transcript_biotype="",
                hgvs_c=f"n.{pos_0 + 1}{ref_nt}>{alt_nt}",
                hgvs_p="",
                cds_pos=0,
                cds_len=0,
                aa_pos=0,
                aa_len=0,
                feature_type="intergenic_region",
            ))

    if not entries:
        return ""
    return "ANN=" + ",".join(entries)


def _ann_field_for_del(
    begin_0:     int,
    end_0:       int,
    ref_allele:  str,
    alt_allele:  str,
    aa_dels:     list,
    frame_shifts: list,
    cds_index:   list[dict],
) -> str:
    """Return 'ANN=...' INFO sub-field for a nucleotide deletion, or ''."""
    entries: list[str] = []
    del_len = end_0 - begin_0

    overlapping = [
        c for c in cds_index if c["begin"] < end_0 and c["end"] > begin_0
    ]

    if not overlapping:
        entries.append(_build_ann_entry(
            allele=alt_allele,
            annotation="intergenic_region",
            impact="MODIFIER",
            gene_name="", gene_id="", feature_id="",
            transcript_biotype="",
            hgvs_c=f"n.{begin_0 + 1}_{end_0}del",
            hgvs_p="", cds_pos=0, cds_len=0, aa_pos=0, aa_len=0,
            feature_type="intergenic_region",
        ))
        return "ANN=" + ",".join(entries)

    for cds in overlapping:
        is_frameshift = (del_len % 3 != 0) or any(
            fs.get("cdsName") == cds["cds_name"] for fs in frame_shifts
        )

        cds_start_1 = _cds_local_pos(cds, max(begin_0, cds["begin"]), use_atg=True)
        cds_end_1   = _cds_local_pos(cds, min(end_0 - 1, cds["end"] - 1), use_atg=True)
        hgvs_c      = _format_hgvs_c_del(cds_start_1, cds_end_1)
        aa_pos_1    = (max(begin_0, cds["atg_begin"]) - cds["atg_begin"]) // 3 + 1

        if is_frameshift:
            matching = [d for d in aa_dels if d["cdsName"] == cds["cds_name"]]
            if matching:
                first = min(matching, key=lambda d: d["pos"])
                ref_aa   = first.get("refAa", "X")
                aa_pos_1 = first["pos"] + 1
            else:
                ref_aa = "X"
            hgvs_p = _format_hgvs_p_fs(ref_aa, aa_pos_1)
            entries.append(_build_ann_entry(
                allele=alt_allele,
                annotation="frameshift_variant",
                impact="HIGH",
                gene_name=cds["gene_name"],
                gene_id=cds["gene_id"],
                feature_id=cds["cds_id"],
                transcript_biotype="protein_coding",
                hgvs_c=hgvs_c,
                hgvs_p=hgvs_p,
                cds_pos=cds_start_1,
                cds_len=cds["length"],
                aa_pos=aa_pos_1,
                aa_len=cds["aa_length"],
            ))
        else:
            matching = [
                d for d in aa_dels
                if d["cdsName"] == cds["cds_name"]
                and cds["atg_begin"] + d["pos"] * 3 >= begin_0
                and cds["atg_begin"] + d["pos"] * 3 < end_0
            ]
            if matching:
                aa_changes = [(d.get("refAa", "X"), d["pos"] + 1) for d in matching]
                hgvs_p = _format_hgvs_p_del(aa_changes)
            else:
                hgvs_p = f"p.{aa_pos_1}del"

            entries.append(_build_ann_entry(
                allele=alt_allele,
                annotation="conservative_inframe_deletion",
                impact="MODERATE",
                gene_name=cds["gene_name"],
                gene_id=cds["gene_id"],
                feature_id=cds["cds_id"],
                transcript_biotype="protein_coding",
                hgvs_c=hgvs_c,
                hgvs_p=hgvs_p,
                cds_pos=cds_start_1,
                cds_len=cds["length"],
                aa_pos=aa_pos_1,
                aa_len=cds["aa_length"],
            ))

    return "ANN=" + ",".join(entries)


def _ann_field_for_ins(
    pos_0:       int,
    ins_seq:     str,
    alt_allele:  str,
    aa_ins:      list,
    frame_shifts: list,
    cds_index:   list[dict],
) -> str:
    """Return 'ANN=...' INFO sub-field for a nucleotide insertion, or ''."""
    entries: list[str] = []
    ins_len = len(ins_seq)

    overlapping = _locate_cds(cds_index, pos_0)

    if not overlapping:
        entries.append(_build_ann_entry(
            allele=alt_allele,
            annotation="intergenic_region",
            impact="MODIFIER",
            gene_name="", gene_id="", feature_id="",
            transcript_biotype="",
            hgvs_c=f"n.{pos_0}_{pos_0 + 1}ins{ins_seq}",
            hgvs_p="", cds_pos=0, cds_len=0, aa_pos=0, aa_len=0,
            feature_type="intergenic_region",
        ))
        return "ANN=" + ",".join(entries)

    for cds in overlapping:
        is_frameshift = (ins_len % 3 != 0) or any(
            fs.get("cdsName") == cds["cds_name"] for fs in frame_shifts
        )

        cds_anchor_1 = _cds_local_pos(cds, pos_0, use_atg=True)
        aa_pos_1     = (pos_0 - cds["atg_begin"]) // 3 + 1
        hgvs_c       = _format_hgvs_c_ins(cds_anchor_1, ins_seq)

        if is_frameshift:
            matching = [i for i in aa_ins if i["cdsName"] == cds["cds_name"]]
            if matching:
                first    = min(matching, key=lambda i: i["pos"])
                aa_pos_1 = first["pos"] + 1
                ref_aa   = "X"
            else:
                ref_aa = "X"
            hgvs_p = _format_hgvs_p_fs(ref_aa, aa_pos_1)
            entries.append(_build_ann_entry(
                allele=alt_allele,
                annotation="frameshift_variant",
                impact="HIGH",
                gene_name=cds["gene_name"],
                gene_id=cds["gene_id"],
                feature_id=cds["cds_id"],
                transcript_biotype="protein_coding",
                hgvs_c=hgvs_c,
                hgvs_p=hgvs_p,
                cds_pos=cds_anchor_1,
                cds_len=cds["length"],
                aa_pos=aa_pos_1,
                aa_len=cds["aa_length"],
            ))
        else:
            matching = [
                i for i in aa_ins
                if i["cdsName"] == cds["cds_name"]
                and cds["atg_begin"] + i["pos"] * 3 == pos_0 + 1
            ]
            if matching:
                first    = matching[0]
                aa_pos_1 = first["pos"] + 1
                hgvs_p   = _format_hgvs_p_ins(aa_pos_1, first.get("ins", ""))
            else:
                hgvs_p = f"p.{aa_pos_1}_{aa_pos_1 + 1}ins"

            entries.append(_build_ann_entry(
                allele=alt_allele,
                annotation="conservative_inframe_insertion",
                impact="MODERATE",
                gene_name=cds["gene_name"],
                gene_id=cds["gene_id"],
                feature_id=cds["cds_id"],
                transcript_biotype="protein_coding",
                hgvs_c=hgvs_c,
                hgvs_p=hgvs_p,
                cds_pos=cds_anchor_1,
                cds_len=cds["length"],
                aa_pos=aa_pos_1,
                aa_len=cds["aa_length"],
            ))

    return "ANN=" + ",".join(entries)

# ---------------------------------------------------------------------------
# Mutation extraction: JSON
# ---------------------------------------------------------------------------

def variants_from_json_result(
    result: dict, ref_seq: str, chrom: str, qual: int
) -> list[tuple[int, str]]:
    """Return sorted list of (vcf_pos, vcf_row) for one nextclade result."""
    rows: list[tuple[int, str]] = []

    cds_index    = _build_cds_index(result.get("annotation", {}), result.get("deletions", []), ref_seq)
    nuc_to_aa    = result.get("nucToAaMuts", {})
    aa_dels      = result.get("aaDeletions", [])
    aa_ins       = result.get("aaInsertions", [])
    frame_shifts = result.get("frameShifts", [])

    for sub in result.get("substitutions", []):
        vcf_pos = sub["pos"] + 1  # nextclade pos is 0-based
        ann     = _ann_field_for_sub(
            sub["pos"], sub["refNuc"], sub["qryNuc"],
            nuc_to_aa, cds_index, ref_seq,
        )
        rows.append((vcf_pos, variant_row(chrom, vcf_pos, sub["refNuc"], sub["qryNuc"], qual, ann=ann)))

    for dele in result.get("deletions", []):
        begin = dele["range"]["begin"]  # 0-based inclusive
        end   = dele["range"]["end"]    # 0-based exclusive
        vcf_pos, ref_a, alt_a = _deletion_alleles(ref_seq, begin, end)
        ann = _ann_field_for_del(begin, end, ref_a, alt_a, aa_dels, frame_shifts, cds_index)
        rows.append((vcf_pos, variant_row(chrom, vcf_pos, ref_a, alt_a, qual, ann=ann)))

    for ins in result.get("insertions", []):
        # pos: 0-based position of the last ref base before the insertion
        vcf_pos, ref_a, alt_a = _insertion_alleles(ref_seq, ins["pos"], ins["ins"])
        ann = _ann_field_for_ins(ins["pos"], ins["ins"], alt_a, aa_ins, frame_shifts, cds_index)
        rows.append((vcf_pos, variant_row(chrom, vcf_pos, ref_a, alt_a, qual, ann=ann)))

    rows.sort(key=lambda x: x[0])
    return rows

# ---------------------------------------------------------------------------
# Mutation extraction: TSV
# ---------------------------------------------------------------------------

_SUB_RE = re.compile(r"([A-Z])(\d+)([A-Z])")
_DEL_RE = re.compile(r"(\d+)-(\d+)")
_INS_RE = re.compile(r"(\d+):([A-Za-z]+)")


def variants_from_tsv_row(
    row: dict, ref_seq: str, chrom: str, qual: int
) -> list[tuple[int, str]]:
    """Return sorted list of (vcf_pos, vcf_row) for one TSV row.

    NOTE: TSV input does not provide per-CDS amino-acid resolution and therefore
    does not emit ANN= annotations.  Use JSON input for annotated output.
    """
    rows: list[tuple[int, str]] = []

    subs_str = (row.get("substitutions") or "").strip()
    for m in _SUB_RE.finditer(subs_str):
        vcf_pos = int(m.group(2))  # TSV positions are 1-based
        rows.append((vcf_pos, variant_row(chrom, vcf_pos, m.group(1), m.group(3), qual)))

    dels_str = (row.get("deletions") or "").strip()
    for m in _DEL_RE.finditer(dels_str):
        begin = int(m.group(1)) - 1  # TSV is 1-based inclusive → 0-based
        end   = int(m.group(2))      # TSV end is 1-based inclusive → 0-based exclusive
        vcf_pos, ref_a, alt_a = _deletion_alleles(ref_seq, begin, end)
        rows.append((vcf_pos, variant_row(chrom, vcf_pos, ref_a, alt_a, qual)))

    ins_str = (row.get("insertions") or "").strip()
    for m in _INS_RE.finditer(ins_str):
        # TSV format: "{1-based-last-ref-pos}:{inserted-sequence}"
        anchor_0 = int(m.group(1)) - 1  # convert to 0-based
        vcf_pos, ref_a, alt_a = _insertion_alleles(ref_seq, anchor_0, m.group(2))
        rows.append((vcf_pos, variant_row(chrom, vcf_pos, ref_a, alt_a, qual)))

    rows.sort(key=lambda x: x[0])
    return rows

# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

def write_vcf(
    output_path: str,
    chrom: str,
    ref_len: int,
    ref_path: str,
    sample: str,
    rows: list[tuple[int, str]],
) -> int:
    header = build_header(chrom, ref_len, ref_path, sample)
    with open(output_path, "w") as fh:
        fh.write(header)
        fh.write("\n")
        for _, row in rows:
            fh.write(row)
            fh.write("\n")
    return len(rows)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Convert nextclade JSON or TSV output to a freebayes-format VCF. "
            "Synthesises depth-1 freebayes INFO/FORMAT fields for use with "
            "reformat_vcf_cyvcf2.py.  JSON input emits snpEff-compatible ANN= "
            "annotations; TSV input does not."
        )
    )
    p.add_argument("-i", "--input",     required=True,
                   help="nextclade JSON (default) or TSV (use --tsv) input file")
    p.add_argument("-r", "--reference", required=True,
                   help="Reference FASTA used during nextclade run (single sequence)")
    p.add_argument("-o", "--output",    required=True,
                   help="Output VCF file (used as prefix when input has multiple samples)")
    p.add_argument("--tsv",    action="store_true",
                   help="Treat input as TSV instead of JSON")
    p.add_argument("--sample", default=None,
                   help="Override sample name (single-sample inputs only)")
    p.add_argument("--qual",   type=int, default=60,
                   help="QUAL score for all records (default: 60)")
    return p.parse_args()


def main() -> None:
    args = get_args()
    chrom, ref_seq = load_reference(args.reference)
    ref_len = len(ref_seq)
    log.info("Loaded reference %s (length %d) from %s", chrom, ref_len, args.reference)

    results: list[tuple[str, list[tuple[int, str]]]] = []

    if args.tsv:
        with open(args.input) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("errors"):
                    log.warning("Skipping %s: errors=%s", row.get("seqName", "?"), row["errors"])
                    continue
                sample = args.sample or row.get("seqName", "SAMPLE")
                results.append((sample, variants_from_tsv_row(row, ref_seq, chrom, args.qual)))
    else:
        with open(args.input) as fh:
            data = json.load(fh)
        for result in data.get("results", []):
            if result.get("errors"):
                log.warning("Skipping %s: errors=%s", result.get("seqName", "?"), result["errors"])
                continue
            sample = args.sample or result.get("seqName", "SAMPLE")
            results.append((sample, variants_from_json_result(result, ref_seq, chrom, args.qual)))

    if not results:
        log.error("No results found in %s", args.input)
        sys.exit(1)

    output_path = Path(args.output)

    if len(results) == 1:
        sample, rows = results[0]
        n = write_vcf(str(output_path), chrom, ref_len, args.reference, sample, rows)
        log.info("Wrote %d records for '%s' -> %s", n, sample, output_path)
    else:
        if args.sample:
            log.warning("--sample ignored for multi-sample input")
        for sample, rows in results:
            safe = re.sub(r"[^\w.-]", "_", sample)
            out = output_path.parent / f"{output_path.stem}_{safe}.vcf"
            n = write_vcf(str(out), chrom, ref_len, args.reference, sample, rows)
            log.info("Wrote %d records for '%s' -> %s", n, sample, out)


if __name__ == "__main__":
    main()
