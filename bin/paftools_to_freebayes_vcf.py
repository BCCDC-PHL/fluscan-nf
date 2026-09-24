#!/usr/bin/env python3
"""Convert a paftools.js-format VCF to a freebayes-format VCF.

Paftools VCFs are produced from consensus FASTA input (depth=1) and carry only
QNAME/QSTART/QSTRAND plus snpEff annotations.  The downstream summarizer
(summarize_vcf_cyvcf2.py) expects the richer freebayes INFO/FORMAT schema.

This script synthesises all missing freebayes fields with plausible depth-1
defaults (DP=1, AO=1, RO=0, AF=1, MQM=60, …) while preserving every
paftools-specific field (QNAME, QSTART, QSTRAND, ANN, LOF, NMD, MNV_*).

Pipeline position: after paftools.js + snpEff (+ optional fix_mnv_codons.py),
before summarize_vcf_cyvcf2.py.
"""

import argparse
import logging

import numpy as np

try:
    from cyvcf2 import VCF, Writer
except ImportError as exc:
    raise SystemExit("cyvcf2 is required: conda install -c bioconda cyvcf2") from exc

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Freebayes header definitions to inject
# Each tuple: (ID, Number, Type, Description)
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
    ("GQ",  "1", "Integer", "Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype"),
    ("GL",  "G", "Float",   "Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy"),
    ("DP",  "1", "Integer", "Read Depth"),
    ("AD",  "R", "Integer", "Number of observation for each allele"),
    ("RO",  "1", "Integer", "Reference allele observation count"),
    ("QR",  "1", "Integer", "Sum of quality of the reference observations"),
    ("AO",  "A", "Integer", "Alternate allele observation count"),
    ("QA",  "A", "Integer", "Sum of quality of the alternate observations"),
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
# Main conversion
# ---------------------------------------------------------------------------

def _existing_header_ids(vcf: VCF) -> set[str]:
    ids = set()
    for h in vcf.header_iter():
        info = h.info()
        if info.get("HeaderType") in ("INFO", "FORMAT"):
            ids.add(info.get("ID", ""))
    return ids


def convert(input_path: str, output_path: str) -> int:
    vcf = VCF(input_path)
    existing = _existing_header_ids(vcf)

    for field_id, number, typ, desc in _FREEBAYES_INFO:
        if field_id not in existing:
            vcf.add_info_to_header({"ID": field_id, "Number": number, "Type": typ, "Description": desc})

    for field_id, number, typ, desc in _FREEBAYES_FORMAT:
        if field_id not in existing:
            vcf.add_format_to_header({"ID": field_id, "Number": number, "Type": typ, "Description": desc})

    vcf.add_to_header("##source=paftools_to_freebayes_vcf.py")

    writer = Writer(output_path, vcf)
    n_written = 0

    try:
        for variant in vcf:
            alts = [a for a in variant.ALT if a is not None]
            if not alts:
                writer.write_record(variant)
                n_written += 1
                continue

            if len(alts) > 1:
                log.warning(
                    "%s:%d has %d ALT alleles; freebayes-format fields will reflect only the first",
                    variant.CHROM, variant.POS, len(alts),
                )

            alt0 = alts[0]
            qa = int(variant.QUAL or 60)
            type_str, cigar_str, len_val = classify(variant.REF, alt0)

            # --- INFO: scalar values (cyvcf2 accepts scalars for Number=A when one alt) ---
            variant.INFO["NS"]      = 1
            variant.INFO["DP"]      = 1
            variant.INFO["DPB"]     = 1.0
            variant.INFO["AC"]      = 1
            variant.INFO["AN"]      = 1
            variant.INFO["AF"]      = 1.0
            variant.INFO["RO"]      = 0
            variant.INFO["AO"]      = 1
            variant.INFO["PRO"]     = 0.0
            variant.INFO["PAO"]     = 0.0
            variant.INFO["QR"]      = 0
            variant.INFO["QA"]      = qa
            variant.INFO["PQR"]     = 0.0
            variant.INFO["PQA"]     = 0.0
            variant.INFO["SRF"]     = 0
            variant.INFO["SRR"]     = 0
            variant.INFO["SAF"]     = 0
            variant.INFO["SAR"]     = 0
            variant.INFO["SRP"]     = 0.0
            variant.INFO["SAP"]     = 0.0
            variant.INFO["AB"]      = 0.0
            variant.INFO["ABP"]     = 0.0
            variant.INFO["RUN"]     = 1
            variant.INFO["RPP"]     = 0.0
            variant.INFO["RPPR"]    = 0.0
            variant.INFO["RPL"]     = 0.0
            variant.INFO["RPR"]     = 0.0
            variant.INFO["EPP"]     = 0.0
            variant.INFO["EPPR"]    = 0.0
            variant.INFO["DPRA"]    = 0.0
            variant.INFO["ODDS"]    = 0.0
            variant.INFO["GTI"]     = 0
            variant.INFO["TYPE"]    = type_str
            variant.INFO["CIGAR"]   = cigar_str
            variant.INFO["NUMALT"]  = 1
            variant.INFO["MEANALT"] = 1.0
            variant.INFO["LEN"]     = len_val
            variant.INFO["MQM"]     = 60.0
            variant.INFO["MQMR"]    = 60.0
            variant.INFO["PAIRED"]  = 0.0
            variant.INFO["PAIREDR"] = 0.0

            # --- FORMAT: numpy arrays, one row per sample ---
            # AD: Number=R → shape (n_samples, n_alleles) = (1, 2) for ref + one alt
            ad = np.array([[0, 1]], dtype=np.int32)
            variant.set_format("DP", np.array([[1]], dtype=np.int32))
            variant.set_format("AD", ad)
            variant.set_format("RO", np.array([[0]], dtype=np.int32))
            variant.set_format("QR", np.array([[0]], dtype=np.int32))
            variant.set_format("AO", np.array([[1]], dtype=np.int32))
            variant.set_format("QA", np.array([[qa]], dtype=np.int32))
            # GL: Number=G → for haploid with 1 alt, 2 genotype likelihoods [P(REF), P(ALT)]
            variant.set_format("GL", np.array([[-float(qa) / 10.0, 0.0]], dtype=np.float32))
            variant.set_format("GQ", np.array([[qa]], dtype=np.int32))

            # --- GT: rewrite paftools 1/1 → freebayes haploid 1 ---
            variant.genotypes = [[1, False] for _ in vcf.samples]

            writer.write_record(variant)
            n_written += 1
    finally:
        writer.close()
        vcf.close()

    return n_written

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Convert a paftools.js+snpEff VCF to a freebayes-format VCF. "
            "Synthesises depth-1 freebayes INFO/FORMAT fields so that "
            "summarize_vcf_cyvcf2.py can parse the output unchanged."
        )
    )
    p.add_argument("-i", "--input",  required=True, help="Input paftools VCF")
    p.add_argument("-o", "--output", required=True, help="Output freebayes-format VCF")
    return p.parse_args()


def main() -> None:
    args = get_args()
    log.info("Converting %s -> %s", args.input, args.output)
    n = convert(args.input, args.output)
    log.info("Wrote %d records to %s", n, args.output)


if __name__ == "__main__":
    main()
