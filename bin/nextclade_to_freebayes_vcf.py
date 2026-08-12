#!/usr/bin/env python3
"""Convert nextclade JSON output to a freebayes-format VCF with snpEff-compatible ANN= annotations.

Reads a per-segment GFF3 annotation (the same file Nextclade consumed) to resolve CDS
structure, then annotates each variant with per-transcript HGVS.c / HGVS.p entries for
main CDS transcripts, mature-peptide transcripts, intronic positions, upstream regions,
and intergenic positions — matching the annotation style produced by snpEff in the
paftools pipeline.

Handles overlapping influenza ORFs (M1/M2, NS1/NEP, PB1/PB1-F2, PA/PA-X) and
spliced CDSes (M2, NEP, PA-X) natively from GFF3 segment data.

Requires:
  -i    Nextclade v3 JSON (per-segment)
  -d    Nextclade reference dataset directory (must contain reference.fasta and genome_annotation.gff3)

--output-no-mature writes a second VCF alongside the first, carrying the same records but
no mature-peptide ANN entries: HA is annotated as HA0 only, with no HA1/HA2 entries.

Pipeline position: after nextclade run, before fix_mnv_codons.py.
"""

import argparse
import json
import logging
import re
import sys
from datetime import date
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Freebayes header definitions
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
    ("DP",  "1", "Integer", "Read Depth"),
    ("AD",  "R", "Integer", "Number of observation for each allele"),
    ("RO",  "1", "Integer", "Reference allele observation count"),
    ("QR",  "1", "Integer", "Sum of quality of the reference observations"),
    ("AO",  "A", "Integer", "Alternate allele observation count"),
    ("QA",  "A", "Integer", "Sum of quality of the alternate observations"),
    ("GL",  "G", "Float",   "Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy"),
    ("GQ",  "1", "Float",   "Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype"),
]

# ---------------------------------------------------------------------------
# Variant classification
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
# GFF3 parsing
# ---------------------------------------------------------------------------

def _parse_gff3_attrs(attr_str: str) -> dict:
    attrs = {}
    for part in attr_str.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def load_gff3(gff_path: str) -> tuple[list, list]:
    """Parse GFF3; return (genes, transcripts).

    genes: [{gene_name, gene_id, begin_0, end_0}]
    transcripts: [{transcript_kind, gene_name, gene_id, feature_id, protein_id,
                   strand, segments, is_mature, warnings, cds_length, aa_length}]

    segments is a sorted list of (begin_0, end_0) half-open intervals.
    protein_id matches Nextclade's cdsName key in nucToAaMuts/aaSubstitutions.
    For main CDS transcripts, feature_id == gene_id (locus_tag).
    For mature peptide transcripts, feature_id == protein_id.
    aa_length excludes the stop codon (cds_length // 3 - 1).

    Tolerates custom references that lack locus_tag and protein_id by falling
    back to gene name and the GFF3 ID attribute respectively.
    """
    genes: list = []
    gene_by_locus_tag: dict = {}     # locus_tag → (gene_name, gene_id)
    cds_segments: dict = {}          # cds_id (GFF3 ID) → [(begin_0, end_0)]
    cds_meta: dict = {}              # cds_id → {gene_name, gene_id, strand, protein_id}
    mature_list: list = []           # raw mature peptide rows

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype  = parts[2]
            b0     = int(parts[3]) - 1   # GFF3 1-based inclusive → 0-based
            e0     = int(parts[4])        # 1-based inclusive → 0-based half-open
            strand = parts[6]
            attrs  = _parse_gff3_attrs(parts[8])

            if ftype == "gene":
                gname = attrs.get("gene", attrs.get("Name", ""))
                # Fall back to gene name when locus_tag is absent (custom references).
                ltag  = attrs.get("locus_tag") or gname
                genes.append({"gene_name": gname, "gene_id": ltag,
                               "begin_0": b0, "end_0": e0})
                gene_by_locus_tag[ltag] = (gname, ltag)

            elif ftype == "CDS":
                pid_attr = attrs.get("protein_id", "")
                cds_id   = attrs.get("ID", "")
                parent   = attrs.get("Parent", "")
                ltag_p   = parent[len("gene-"):] if parent.startswith("gene-") else parent
                gname, gid = gene_by_locus_tag.get(ltag_p, ("", ltag_p))
                # Use the full GFF3 ID as the storage key so the Parent attribute of
                # mature_protein_region_of_CDS rows matches verbatim (e.g. "cds-HA"
                # or "cds-YP_308669.1"), avoiding the key mismatch when protein_id is absent.
                storage_key = cds_id or pid_attr or attrs.get("Name", "")
                # Nextclade uses protein_id as cdsName when present; falls back to the
                # GFF3 ID when protein_id is absent (verified against chicken_2021 JSON).
                nextclade_match = pid_attr or cds_id or attrs.get("Name", "")
                if storage_key not in cds_segments:
                    cds_segments[storage_key] = []
                    cds_meta[storage_key] = {
                        "gene_name":  gname,
                        "gene_id":    gid,
                        "strand":     strand,
                        "protein_id": nextclade_match,
                    }
                cds_segments[storage_key].append((b0, e0))

            elif ftype == "mature_protein_region_of_CDS":
                pid    = attrs.get("protein_id", attrs.get("product", attrs.get("Name", "")))
                parent = attrs.get("Parent", "")
                # parent is already the full GFF3 ID (e.g. "cds-HA") — use it as-is
                # so the lookup into cds_meta succeeds without stripping "cds-".
                mature_list.append({
                    "protein_id":     pid,
                    "parent_cds_key": parent,
                    "begin_0":        b0,
                    "end_0":          e0,
                    "strand":         strand,
                })

    transcripts: list = []

    for storage_key, segs in cds_segments.items():
        meta    = cds_meta[storage_key]
        segs_s  = sorted(segs)
        cds_len = sum(e - b for b, e in segs_s)
        transcripts.append({
            "transcript_kind": "cds",
            "gene_name":   meta["gene_name"],
            "gene_id":     meta["gene_id"],
            "feature_id":  meta["gene_id"],    # Feature_ID = locus_tag, matching NCBI VCF shape
            "protein_id":  meta["protein_id"], # matches Nextclade cdsName
            "strand":      meta["strand"],
            "segments":    segs_s,
            "is_mature":   False,
            "warnings":    "",
            "cds_length":  cds_len,
            "aa_length":   max(cds_len // 3 - 1, 0),
        })

    for me in mature_list:
        meta    = cds_meta.get(me["parent_cds_key"], {"gene_name": "", "gene_id": "", "protein_id": ""})
        seg_len = me["end_0"] - me["begin_0"]
        transcripts.append({
            "transcript_kind": "mature_peptide",
            "gene_name":   meta["gene_name"],
            "gene_id":     meta["gene_id"],
            "feature_id":  me["protein_id"],   # mature peptide uses its own protein_id
            "protein_id":  me["protein_id"],
            "strand":      me["strand"],
            "segments":    [(me["begin_0"], me["end_0"])],
            "is_mature":   True,
            "warnings":    "WARNING_TRANSCRIPT_NO_START_CODON",
            "cds_length":  seg_len,
            "aa_length":   max(seg_len // 3 - 1, 0),
        })

    return genes, transcripts

# ---------------------------------------------------------------------------
# CDS position helpers
# ---------------------------------------------------------------------------

def _cds_pos_from_genomic(tx: dict, pos_0: int) -> tuple:
    """Return (cds_local_1, exon_rank, exon_count) if pos_0 is inside a CDS segment, else (None,None,None)."""
    segments = tx["segments"]
    running  = 0
    for i, (b, e) in enumerate(segments):
        if b <= pos_0 < e:
            if tx["strand"] == "+":
                local_1 = running + (pos_0 - b) + 1
            else:
                local_1 = running + (e - 1 - pos_0) + 1
            return local_1, i + 1, len(segments)
        running += (e - b)
    return None, None, None


def _intron_info(tx: dict, pos_0: int) -> Optional[tuple]:
    """For a spliced transcript, return (cds_end_1, offset, intron_rank, n_introns) if pos_0 is in an intron."""
    segs     = tx["segments"]
    n_introns = len(segs) - 1
    if n_introns == 0:
        return None
    running = 0
    for i, (b, e) in enumerate(segs[:-1]):
        next_b = segs[i + 1][0]
        running += (e - b)
        if e <= pos_0 < next_b:
            offset = pos_0 - e + 1   # 1-based distance into intron
            return running, offset, i + 1, n_introns
    return None


def _is_in_gene(genes: list, pos_0: int) -> bool:
    return any(g["begin_0"] <= pos_0 < g["end_0"] for g in genes)


def _intergenic_name(genes: list, pos_0: int) -> tuple[str, str, str]:
    """Return (gene_name_str, gene_id_str, feature_id_str) for an intergenic entry."""
    sorted_g    = sorted(genes, key=lambda g: g["begin_0"])
    upstream    = [g for g in sorted_g if g["end_0"] <= pos_0]
    downstream  = [g for g in sorted_g if g["begin_0"] > pos_0]
    prev_name   = upstream[-1]["gene_name"] if upstream else "CHR_START"
    prev_id     = upstream[-1]["gene_id"]   if upstream else "CHR_START"
    next_name   = downstream[0]["gene_name"] if downstream else "CHR_END"
    next_id     = downstream[0]["gene_id"]   if downstream else "CHR_END"
    gene_name   = f"{prev_name}-{next_name}"
    gene_id     = f"{prev_id}-{next_id}"
    return gene_name, gene_id, gene_id

# ---------------------------------------------------------------------------
# Codon helpers
# ---------------------------------------------------------------------------

def _translate_codon(codon: str) -> str:
    try:
        return str(Seq(codon).translate())
    except Exception:
        return "X"


def _cds_to_genomic(cds_pos_0: int, segments: list, strand: str) -> Optional[int]:
    """Map 0-based CDS-local position to 0-based genomic position."""
    running = 0
    for b, e in segments:
        seg_len = e - b
        if running + seg_len > cds_pos_0:
            if strand == "+":
                return b + (cds_pos_0 - running)
            else:
                return e - 1 - (cds_pos_0 - running)
        running += seg_len
    return None


def _synonymous_aa_change(ref_seq: str, tx: dict, pos_0: int, ref_nt: str, alt_nt: str) -> tuple[str, str, int]:
    """Return (refAa, altAa, aa_pos_1) by direct codon translation from the reference.

    Falls back to ('X','X',0) when the codon cannot be determined (e.g. minus strand,
    segment-boundary ambiguity, or out-of-range position).
    """
    segs   = tx["segments"]
    strand = tx["strand"]

    # CDS-local 0-based position
    running = 0
    cds_local_0 = None
    for b, e in segs:
        if b <= pos_0 < e:
            cds_local_0 = (running + (pos_0 - b)) if strand == "+" else (running + (e - 1 - pos_0))
            break
        running += (e - b)
    if cds_local_0 is None:
        return "X", "X", 0

    codon_idx    = cds_local_0 // 3
    codon_offset = cds_local_0 % 3
    aa_pos_1     = codon_idx + 1

    # Genomic positions for all 3 codon bases
    gpos = [_cds_to_genomic(codon_idx * 3 + i, segs, strand) for i in range(3)]
    if None in gpos:
        return "X", "X", aa_pos_1

    if strand == "+":
        ref_codon = "".join(ref_seq[g] for g in gpos)
        alt_codon = ref_codon[:codon_offset] + alt_nt + ref_codon[codon_offset + 1:]
    else:
        # minus strand: reading direction is reversed; just fall back
        return "X", "X", aa_pos_1

    return _translate_codon(ref_codon), _translate_codon(alt_codon), aa_pos_1

# ---------------------------------------------------------------------------
# HGVS formatters  (1-letter amino-acid codes, matching paftools/snpEff output)
# ---------------------------------------------------------------------------

def _format_hgvs_c_snp(cds_local_1: int, ref_nt: str, alt_nt: str) -> str:
    return f"c.{cds_local_1}{ref_nt}>{alt_nt}"


def _format_hgvs_c_del(start_1: int, end_1: int) -> str:
    if start_1 == end_1:
        return f"c.{start_1}del"
    return f"c.{start_1}_{end_1}del"


def _format_hgvs_c_ins(anchor_1: int, ins_seq: str) -> str:
    return f"c.{anchor_1}_{anchor_1 + 1}ins{ins_seq}"


def _format_hgvs_c_intron(cds_end_1: int, offset: int, ref_nt: str, alt_nt: str) -> str:
    return f"c.{cds_end_1}+{offset}{ref_nt}>{alt_nt}"


def _format_hgvs_p_sub(ref_aa: str, aa_pos_1: int, alt_aa: str) -> str:
    return f"p.{ref_aa}{aa_pos_1}{alt_aa}"


def _format_hgvs_p_syn(ref_aa: str, aa_pos_1: int) -> str:
    return f"p.{ref_aa}{aa_pos_1}{ref_aa}"


def _format_hgvs_p_del(aa_changes: list) -> str:
    """Return HGVS.p for an in-frame deletion given list of (refAa, aa_pos_1)."""
    if not aa_changes:
        return ""
    sc = sorted(aa_changes, key=lambda x: x[1])
    if len(sc) == 1:
        return f"p.{sc[0][0]}{sc[0][1]}del"
    return f"p.{sc[0][0]}{sc[0][1]}_{sc[-1][0]}{sc[-1][1]}del"


def _format_hgvs_p_ins(aa_pos_1: int, ins_aa_seq: str) -> str:
    return f"p.{aa_pos_1}_{aa_pos_1 + 1}ins{ins_aa_seq}"


def _format_hgvs_p_fs(ref_aa: str, aa_pos_1: int) -> str:
    return f"p.{ref_aa}{aa_pos_1}fs"

# ---------------------------------------------------------------------------
# Consequence classification
# ---------------------------------------------------------------------------

def _classify_aa_consequence(ref_aa: str, alt_aa: str, aa_pos_1: int) -> tuple[str, str]:
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
    allele:             str,
    annotation:         str,
    impact:             str,
    gene_name:          str,
    gene_id:            str,
    feature_id:         str,
    transcript_biotype: str,
    hgvs_c:             str,
    hgvs_p:             str,
    cds_pos:            int,
    cds_len:            int,
    aa_pos:             int,
    aa_len:             int,
    distance:           str = "",
    feature_type:       str = "transcript",
    rank_str:           str = "",
    warnings:           str = "",
) -> str:
    cds_str = f"{cds_pos}/{cds_len}" if cds_len > 0 else ""
    aa_str  = f"{aa_pos}/{aa_len}"   if aa_len  > 0 else ""
    if not rank_str and feature_type == "transcript" and cds_len > 0:
        rank_str = "1/1"
    return "|".join([
        allele, annotation, impact,
        gene_name, gene_id,
        feature_type, feature_id, transcript_biotype,
        rank_str, hgvs_c, hgvs_p,
        cds_str, cds_str,  # cDNA.pos/length repeated for CDS.pos/length
        aa_str, distance, warnings,
    ])

# ---------------------------------------------------------------------------
# Indel anchor helpers
# ---------------------------------------------------------------------------

def _deletion_alleles(ref_seq: str, begin: int, end: int) -> tuple[int, str, str]:
    if begin > 0:
        anchor = begin - 1
        return anchor + 1, ref_seq[anchor] + ref_seq[begin:end], ref_seq[anchor]
    return 1, ref_seq[begin:end] + ref_seq[end], ref_seq[end]


def _insertion_alleles(ref_seq: str, pos: int, ins_seq: str) -> tuple[int, str, str]:
    if pos < 0:
        # Defensive fallback: terminal-overhang insertions (nextclade pos == -1) are
        # normally filtered out upstream in variants_from_json_result, so this branch is
        # not reached in the standard pipeline. Kept to guarantee valid VCF if ever called.
        # Right-anchor to POS 1: ALT = ins_seq + ref[0], REF = ref[0].
        # _deletion_alleles already handles its analogous begin==0 edge case;
        # without this guard, anchor = -1 causes POS = 0 (invalid VCF) and
        # ref_seq[-1] silently returns the *last* reference base.
        return 1, ref_seq[0], ins_seq + ref_seq[0]
    anchor = pos
    return anchor + 1, ref_seq[anchor], ref_seq[anchor] + ins_seq

# ---------------------------------------------------------------------------
# VCF row helpers
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
    return f"1:1:0,1:0:0:1:{qual}:{gl}:{qual}"


def variant_row(chrom: str, pos: int, ref: str, alt: str, qual: int, ann: str = "", filter: str = "PASS") -> str:
    type_str, cigar, length = classify(ref, alt)
    info = _info(type_str, cigar, length, qual)
    if ann:
        info = f"{info};{ann}"
    return (
        f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter}"
        f"\t{info}"
        f"\tGT:DP:AD:RO:QR:AO:QA:GL:GQ\t{_sample(qual)}"
    )

# ---------------------------------------------------------------------------
# VCF header
# ---------------------------------------------------------------------------

def build_header(chrom: str, ref_len: int, ref_path: str, sample: str) -> str:
    lines = [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=AmbiguousInsert,Description="Insertion contains ambiguous (N) bases in the inserted sequence">',
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
# Per-variant ANN field builders
# ---------------------------------------------------------------------------

def _ann_field_for_sub(
    pos_0:      int,
    ref_nt:     str,
    alt_nt:     str,
    nuc_to_aa:  dict,
    genes:      list,
    transcripts: list,
    ref_seq:    str,
) -> str:
    """Return 'ANN=...' INFO sub-field for a nucleotide substitution."""
    entries: list[str] = []

    # --- 1. Main CDS transcripts: in-CDS or intronic ---
    for tx in transcripts:
        if tx["is_mature"]:
            continue
        cds_local_1, exon_rank, exon_count = _cds_pos_from_genomic(tx, pos_0)
        if cds_local_1 is not None:
            rank_str = f"{exon_rank}/{exon_count}"
            aa_changes = [c for c in nuc_to_aa.get(str(pos_0), [])
                          if c["cdsName"] == tx["protein_id"]]
            if aa_changes:
                ch      = aa_changes[0]
                ref_aa  = ch["refAa"]
                alt_aa  = ch["qryAa"]
                aa_pos_1 = ch["pos"] + 1
            else:
                ref_aa, alt_aa, aa_pos_1 = _synonymous_aa_change(ref_seq, tx, pos_0, ref_nt, alt_nt)

            annotation, impact = _classify_aa_consequence(ref_aa, alt_aa, aa_pos_1)
            hgvs_c = _format_hgvs_c_snp(cds_local_1, ref_nt, alt_nt)
            hgvs_p = (_format_hgvs_p_syn(ref_aa, aa_pos_1)
                      if annotation == "synonymous_variant"
                      else _format_hgvs_p_sub(ref_aa, aa_pos_1, alt_aa))

            entries.append(_build_ann_entry(
                allele=alt_nt, annotation=annotation, impact=impact,
                gene_name=tx["gene_name"], gene_id=tx["gene_id"],
                feature_id=tx["feature_id"],
                transcript_biotype="protein_coding",
                hgvs_c=hgvs_c, hgvs_p=hgvs_p,
                cds_pos=cds_local_1, cds_len=tx["cds_length"],
                aa_pos=aa_pos_1, aa_len=tx["aa_length"],
                rank_str=rank_str,
            ))
            continue

        # Check for intron position in spliced transcripts
        intron = _intron_info(tx, pos_0)
        if intron is not None:
            cds_end_1, offset, intron_rank, n_introns = intron
            rank_str = f"{intron_rank}/{n_introns}"
            hgvs_c = _format_hgvs_c_intron(cds_end_1, offset, ref_nt, alt_nt)
            entries.append(_build_ann_entry(
                allele=alt_nt, annotation="intron_variant", impact="MODIFIER",
                gene_name=tx["gene_name"], gene_id=tx["gene_id"],
                feature_id=tx["feature_id"],
                transcript_biotype="protein_coding",
                hgvs_c=hgvs_c, hgvs_p="",
                cds_pos=0, cds_len=0, aa_pos=0, aa_len=0,
                rank_str=rank_str,
            ))

    # --- 2. Mature peptide transcripts: in-region or upstream ---
    for tx in transcripts:
        if not tx["is_mature"]:
            continue
        seg_begin = tx["segments"][0][0]
        seg_end   = tx["segments"][-1][1]

        if seg_begin <= pos_0 < seg_end:
            cds_local_1, exon_rank, exon_count = _cds_pos_from_genomic(tx, pos_0)
            if cds_local_1 is None:
                continue
            rank_str = f"{exon_rank}/{exon_count}"
            ref_aa, alt_aa, aa_pos_1 = _synonymous_aa_change(ref_seq, tx, pos_0, ref_nt, alt_nt)
            annotation, impact = _classify_aa_consequence(ref_aa, alt_aa, aa_pos_1)
            hgvs_c = _format_hgvs_c_snp(cds_local_1, ref_nt, alt_nt)
            hgvs_p = (_format_hgvs_p_syn(ref_aa, aa_pos_1)
                      if annotation == "synonymous_variant"
                      else _format_hgvs_p_sub(ref_aa, aa_pos_1, alt_aa))
            entries.append(_build_ann_entry(
                allele=alt_nt, annotation=annotation, impact=impact,
                gene_name=tx["gene_name"], gene_id=tx["gene_id"],
                feature_id=tx["feature_id"],
                transcript_biotype="protein_coding",
                hgvs_c=hgvs_c, hgvs_p=hgvs_p,
                cds_pos=cds_local_1, cds_len=tx["cds_length"],
                aa_pos=aa_pos_1, aa_len=tx["aa_length"],
                rank_str=rank_str,
                warnings=tx["warnings"],
            ))
        elif pos_0 < seg_begin:
            distance = seg_begin - pos_0
            hgvs_c   = f"c.-{distance}{ref_nt}>{alt_nt}"
            entries.append(_build_ann_entry(
                allele=alt_nt, annotation="upstream_gene_variant", impact="MODIFIER",
                gene_name=tx["gene_name"], gene_id=tx["gene_id"],
                feature_id=tx["feature_id"],
                transcript_biotype="protein_coding",
                hgvs_c=hgvs_c, hgvs_p="",
                cds_pos=0, cds_len=0, aa_pos=0, aa_len=0,
                distance=str(distance), warnings=tx["warnings"],
            ))

    # --- 3. Upstream entries for main CDS transcripts ---
    for tx in transcripts:
        if tx["is_mature"]:
            continue
        tx_start = tx["segments"][0][0]
        if pos_0 < tx_start:
            distance = tx_start - pos_0
            hgvs_c   = f"c.-{distance}{ref_nt}>{alt_nt}"
            entries.append(_build_ann_entry(
                allele=alt_nt, annotation="upstream_gene_variant", impact="MODIFIER",
                gene_name=tx["gene_name"], gene_id=tx["gene_id"],
                feature_id=tx["feature_id"],
                transcript_biotype="protein_coding",
                hgvs_c=hgvs_c, hgvs_p="",
                cds_pos=0, cds_len=0, aa_pos=0, aa_len=0,
                distance=str(distance),
            ))

    # --- 4. Intergenic (when outside all gene boundaries) ---
    if not _is_in_gene(genes, pos_0):
        gname, gid, fid = _intergenic_name(genes, pos_0)
        entries.append(_build_ann_entry(
            allele=alt_nt, annotation="intergenic_region", impact="MODIFIER",
            gene_name=gname, gene_id=gid, feature_id=fid,
            transcript_biotype="",
            hgvs_c=f"n.{pos_0 + 1}{ref_nt}>{alt_nt}", hgvs_p="",
            cds_pos=0, cds_len=0, aa_pos=0, aa_len=0,
            feature_type="intergenic_region",
        ))

    if not entries:
        return ""
    return "ANN=" + ",".join(entries)


def _ann_field_for_del(
    begin_0:      int,
    end_0:        int,
    ref_allele:   str,
    alt_allele:   str,
    aa_dels:      list,
    frame_shifts: list,
    genes:        list,
    transcripts:  list,
) -> str:
    """Return 'ANN=...' INFO sub-field for a nucleotide deletion."""
    entries: list[str] = []
    del_len = end_0 - begin_0

    for tx in transcripts:
        if tx["is_mature"]:
            continue
        # Check if deletion overlaps any segment of this transcript
        overlap = any(b < end_0 and begin_0 < e for b, e in tx["segments"])
        if not overlap:
            continue

        is_frameshift = (del_len % 3 != 0) or any(
            fs.get("cdsName") == tx["protein_id"] for fs in frame_shifts
        )

        # CDS-local 1-based start of deletion
        del_in_seg_start = max(begin_0, tx["segments"][0][0])
        cds_start_1, _, _ = _cds_pos_from_genomic(tx, del_in_seg_start)
        if cds_start_1 is None:
            cds_start_1 = 1  # fallback

        del_in_seg_end = min(end_0 - 1, tx["segments"][-1][1] - 1)
        cds_end_1, _, _ = _cds_pos_from_genomic(tx, del_in_seg_end)
        if cds_end_1 is None:
            cds_end_1 = cds_start_1

        hgvs_c  = _format_hgvs_c_del(cds_start_1, cds_end_1)
        aa_pos_1 = (cds_start_1 - 1) // 3 + 1

        if is_frameshift:
            matching = [d for d in aa_dels if d["cdsName"] == tx["protein_id"]]
            if matching:
                first    = min(matching, key=lambda d: d["pos"])
                ref_aa   = first.get("refAa", "X")
                aa_pos_1 = first["pos"] + 1
            else:
                ref_aa = "X"
            hgvs_p = _format_hgvs_p_fs(ref_aa, aa_pos_1)
            entries.append(_build_ann_entry(
                allele=alt_allele, annotation="frameshift_variant", impact="HIGH",
                gene_name=tx["gene_name"], gene_id=tx["gene_id"],
                feature_id=tx["feature_id"],
                transcript_biotype="protein_coding",
                hgvs_c=hgvs_c, hgvs_p=hgvs_p,
                cds_pos=cds_start_1, cds_len=tx["cds_length"],
                aa_pos=aa_pos_1, aa_len=tx["aa_length"],
            ))
        else:
            matching = [
                d for d in aa_dels
                if d["cdsName"] == tx["protein_id"]
                and tx["segments"][0][0] + d["pos"] * 3 >= begin_0
                and tx["segments"][0][0] + d["pos"] * 3 < end_0
            ]
            if matching:
                aa_changes = [(d.get("refAa", "X"), d["pos"] + 1) for d in matching]
                hgvs_p = _format_hgvs_p_del(aa_changes)
            else:
                hgvs_p = f"p.{aa_pos_1}del"
            entries.append(_build_ann_entry(
                allele=alt_allele, annotation="conservative_inframe_deletion", impact="MODERATE",
                gene_name=tx["gene_name"], gene_id=tx["gene_id"],
                feature_id=tx["feature_id"],
                transcript_biotype="protein_coding",
                hgvs_c=hgvs_c, hgvs_p=hgvs_p,
                cds_pos=cds_start_1, cds_len=tx["cds_length"],
                aa_pos=aa_pos_1, aa_len=tx["aa_length"],
            ))

    if not entries:
        if not _is_in_gene(genes, begin_0):
            gname, gid, fid = _intergenic_name(genes, begin_0)
            entries.append(_build_ann_entry(
                allele=alt_allele, annotation="intergenic_region", impact="MODIFIER",
                gene_name=gname, gene_id=gid, feature_id=fid,
                transcript_biotype="",
                hgvs_c=f"n.{begin_0 + 1}_{end_0}del", hgvs_p="",
                cds_pos=0, cds_len=0, aa_pos=0, aa_len=0,
                feature_type="intergenic_region",
            ))

    return ("ANN=" + ",".join(entries)) if entries else ""


def _ann_field_for_ins(
    pos_0:        int,
    ins_seq:      str,
    alt_allele:   str,
    aa_ins:       list,
    frame_shifts: list,
    genes:        list,
    transcripts:  list,
) -> str:
    """Return 'ANN=...' INFO sub-field for a nucleotide insertion."""
    entries: list[str] = []
    ins_len = len(ins_seq)

    for tx in transcripts:
        if tx["is_mature"]:
            continue
        cds_anchor_1, exon_rank, exon_count = _cds_pos_from_genomic(tx, pos_0)
        if cds_anchor_1 is None:
            continue

        rank_str      = f"{exon_rank}/{exon_count}"
        aa_pos_1      = (pos_0 - tx["segments"][0][0]) // 3 + 1
        hgvs_c        = _format_hgvs_c_ins(cds_anchor_1, ins_seq)
        is_frameshift = (ins_len % 3 != 0) or any(
            fs.get("cdsName") == tx["protein_id"] for fs in frame_shifts
        )

        if is_frameshift:
            matching = [i for i in aa_ins if i["cdsName"] == tx["protein_id"]]
            if matching:
                first    = min(matching, key=lambda i: i["pos"])
                aa_pos_1 = first["pos"] + 1
                ref_aa   = "X"
            else:
                ref_aa = "X"
            hgvs_p = _format_hgvs_p_fs(ref_aa, aa_pos_1)
            entries.append(_build_ann_entry(
                allele=alt_allele, annotation="frameshift_variant", impact="HIGH",
                gene_name=tx["gene_name"], gene_id=tx["gene_id"],
                feature_id=tx["feature_id"],
                transcript_biotype="protein_coding",
                hgvs_c=hgvs_c, hgvs_p=hgvs_p,
                cds_pos=cds_anchor_1, cds_len=tx["cds_length"],
                aa_pos=aa_pos_1, aa_len=tx["aa_length"],
                rank_str=rank_str,
            ))
        else:
            matching = [
                i for i in aa_ins
                if i["cdsName"] == tx["protein_id"]
                and tx["segments"][0][0] + i["pos"] * 3 == pos_0 + 1
            ]
            if matching:
                first    = matching[0]
                aa_pos_1 = first["pos"] + 1
                hgvs_p   = _format_hgvs_p_ins(aa_pos_1, first.get("ins", ""))
            else:
                hgvs_p = f"p.{aa_pos_1}_{aa_pos_1 + 1}ins"
            entries.append(_build_ann_entry(
                allele=alt_allele, annotation="conservative_inframe_insertion", impact="MODERATE",
                gene_name=tx["gene_name"], gene_id=tx["gene_id"],
                feature_id=tx["feature_id"],
                transcript_biotype="protein_coding",
                hgvs_c=hgvs_c, hgvs_p=hgvs_p,
                cds_pos=cds_anchor_1, cds_len=tx["cds_length"],
                aa_pos=aa_pos_1, aa_len=tx["aa_length"],
                rank_str=rank_str,
            ))

    if not entries:
        if not _is_in_gene(genes, pos_0):
            gname, gid, fid = _intergenic_name(genes, pos_0)
            entries.append(_build_ann_entry(
                allele=alt_allele, annotation="intergenic_region", impact="MODIFIER",
                gene_name=gname, gene_id=gid, feature_id=fid,
                transcript_biotype="",
                hgvs_c=f"n.{pos_0}_{pos_0 + 1}ins{ins_seq}", hgvs_p="",
                cds_pos=0, cds_len=0, aa_pos=0, aa_len=0,
                feature_type="intergenic_region",
            ))

    return ("ANN=" + ",".join(entries)) if entries else ""

# ---------------------------------------------------------------------------
# Mutation extraction: JSON
# ---------------------------------------------------------------------------

def variants_from_json_result(
    result:      dict,
    ref_seq:     str,
    chrom:       str,
    genes:       list,
    transcripts: list,
    qual:        int,
) -> list[tuple[int, str]]:
    """Return sorted list of (vcf_pos, vcf_row) for one nextclade result."""
    rows: list[tuple[int, str]] = []

    nuc_to_aa    = result.get("nucToAaMuts", {})
    aa_dels      = result.get("aaDeletions", [])
    aa_ins       = [
        {**i, "cdsName": i.get("cdsName") or i.get("cds", "")}
        for i in result.get("aaInsertions", [])
    ]
    frame_shifts = result.get("frameShifts", [])

    for sub in result.get("substitutions", []):
        vcf_pos = sub["pos"] + 1
        ann     = _ann_field_for_sub(
            sub["pos"], sub["refNuc"], sub["qryNuc"],
            nuc_to_aa, genes, transcripts, ref_seq,
        )
        rows.append((vcf_pos, variant_row(chrom, vcf_pos, sub["refNuc"], sub["qryNuc"], qual, ann=ann)))

    for dele in result.get("deletions", []):
        begin   = dele["range"]["begin"]
        end     = dele["range"]["end"]
        vcf_pos, ref_a, alt_a = _deletion_alleles(ref_seq, begin, end)
        ann = _ann_field_for_del(begin, end, ref_a, alt_a, aa_dels, frame_shifts, genes, transcripts)
        rows.append((vcf_pos, variant_row(chrom, vcf_pos, ref_a, alt_a, qual, ann=ann)))

    for ins in result.get("insertions", []):
        raw_pos = ins["pos"]
        ins_seq = ins["ins"]
        # Skip terminal-overhang insertions: 5' inserts before the first reference base
        # (raw_pos < 0) and 3' inserts after the last reference base (raw_pos >= ref_len-1).
        # These encode query overhang / N-padding outside the reference span; nucmer/paftools
        # align only within the reference span and never emit them, so keeping them produces
        # terminal false positives. Genuine internal insertions (raw_pos in 0 .. ref_len-2)
        # are unaffected.
        if raw_pos < 0 or raw_pos >= len(ref_seq) - 1:
            continue
        vcf_pos, ref_a, alt_a = _insertion_alleles(ref_seq, raw_pos, ins_seq)
        ann_pos = raw_pos if raw_pos >= 0 else 0
        ann = _ann_field_for_ins(ann_pos, ins_seq, alt_a, aa_ins, frame_shifts, genes, transcripts)
        filter_val = "AmbiguousInsert" if "N" in ins_seq.upper() else "PASS"
        rows.append((vcf_pos, variant_row(chrom, vcf_pos, ref_a, alt_a, qual, ann=ann, filter=filter_val)))

    rows.sort(key=lambda x: x[0])
    return rows

# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

def write_vcf(
    output_path: str,
    chrom:       str,
    ref_len:     int,
    ref_path:    str,
    sample:      str,
    rows:        list[tuple[int, str]],
) -> int:
    header = build_header(chrom, ref_len, ref_path, sample)
    with open(output_path, "w") as fh:
        fh.write(header)
        fh.write("\n")
        for _, row in rows:
            fh.write(row)
            fh.write("\n")
    return len(rows)


def write_results(
    output_path: Path,
    results:     list[tuple[str, list[tuple[int, str]]]],
    chrom:       str,
    ref_len:     int,
    ref_path:    str,
    named_sample: bool,
) -> None:
    """Write one VCF per result; multi-sample input appends the sample name to the stem."""
    if len(results) == 1:
        sample, rows = results[0]
        n = write_vcf(str(output_path), chrom, ref_len, ref_path, sample, rows)
        log.info("Wrote %d records for '%s' -> %s", n, sample, output_path)
        return

    if named_sample:
        log.warning("--sample ignored for multi-sample input")
    for sample, rows in results:
        safe = re.sub(r"[^\w.-]", "_", sample)
        out  = output_path.parent / f"{output_path.stem}_{safe}.vcf"
        n    = write_vcf(str(out), chrom, ref_len, ref_path, sample, rows)
        log.info("Wrote %d records for '%s' -> %s", n, sample, out)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Convert nextclade JSON output to a freebayes-format VCF with "
            "snpEff-compatible ANN= annotations sourced from a GFF3 annotation file."
        )
    )
    p.add_argument("-i", "--input",   required=True,
                   help="Nextclade v3 JSON input file")
    p.add_argument("-d", "--dataset", required=True,
                   help="Nextclade reference dataset directory (must contain reference.fasta and genome_annotation.gff3)")
    p.add_argument("-o", "--output",  required=True,
                   help="Output VCF file (used as prefix when input has multiple samples)")
    p.add_argument("--output-no-mature", default=None, dest="output_no_mature",
                   help="Optional second output VCF with the same records but no "
                        "mature-peptide (HA1/HA2) ANN entries")
    p.add_argument("--sample", default=None,
                   help="Override sample name (single-sample inputs only)")
    p.add_argument("--qual",   type=int, default=60,
                   help="QUAL score for all records (default: 60)")
    return p.parse_args()


_SEED_COVER_RE = re.compile(r"seed alignment covers (\d+(?:\.\d+)?%)")


def _extract_seed_cover(error_messages: list[str]) -> Optional[str]:
    for msg in error_messages:
        m = _SEED_COVER_RE.search(msg)
        if m:
            return m.group(1)
    return None


def main() -> None:
    args = get_args()

    dataset_dir = Path(args.dataset)
    ref_path    = dataset_dir / "reference.fasta"
    gff_path    = dataset_dir / "genome_annotation.gff3"

    for p in (ref_path, gff_path):
        if not p.exists():
            log.error("Expected file not found in dataset directory: %s", p)
            sys.exit(1)

    chrom, ref_seq = load_reference(str(ref_path))
    ref_len = len(ref_seq)
    log.info("Loaded reference %s (length %d) from %s", chrom, ref_len, ref_path)

    genes, transcripts = load_gff3(str(gff_path))
    log.info(
        "Loaded GFF3 %s: %d gene(s), %d transcript(s)",
        gff_path, len(genes), len(transcripts),
    )

    # The no-mature output is the same annotation pass run against a transcript list with
    # the mature peptides removed, so it matches what a GFF3 without them would produce.
    transcripts_no_mature = [tx for tx in transcripts if not tx["is_mature"]]
    if args.output_no_mature:
        log.info(
            "Also writing a no-mature VCF -> %s (%d mature-peptide transcript(s) omitted)",
            args.output_no_mature, len(transcripts) - len(transcripts_no_mature),
        )

    with open(args.input) as fh:
        data = json.load(fh)

    results:           list[tuple[str, list[tuple[int, str]]]] = []
    results_no_mature: list[tuple[str, list[tuple[int, str]]]] = []

    for result in data.get("results", []):
        if result.get("errors"):
            log.warning("Skipping %s: errors=%s", result.get("seqName", "?"), result["errors"])
            continue
        sample = args.sample or result.get("seqName", "SAMPLE")
        results.append((
            sample,
            variants_from_json_result(result, ref_seq, chrom, genes, transcripts, args.qual),
        ))
        if args.output_no_mature:
            results_no_mature.append((
                sample,
                variants_from_json_result(
                    result, ref_seq, chrom, genes, transcripts_no_mature, args.qual,
                ),
            ))

    for error_entry in data.get("errors", []):
        seq_name = error_entry.get("seqName", "UNKNOWN")
        sample   = args.sample or seq_name
        msgs     = error_entry.get("errors", [])
        if any("Unable to align" in m for m in msgs):
            seed_cover = _extract_seed_cover(msgs)
            if seed_cover:
                log.warning("Unable to align '%s' (seed alignment covers %s); writing empty VCF",
                            sample, seed_cover)
            else:
                log.warning("Unable to align '%s'; writing empty VCF", sample)
        else:
            log.warning("Nextclade errors for '%s': %s; writing empty VCF",
                        sample, "; ".join(str(m) for m in msgs))
        results.append((sample, []))
        if args.output_no_mature:
            results_no_mature.append((sample, []))

    if not results:
        log.error("No results found in %s", args.input)
        sys.exit(1)

    write_results(Path(args.output), results, chrom, ref_len, str(ref_path), bool(args.sample))
    if args.output_no_mature:
        # named_sample=False: the first call has already warned about --sample if needed.
        write_results(Path(args.output_no_mature), results_no_mature,
                      chrom, ref_len, str(ref_path), False)


if __name__ == "__main__":
    main()
