"""Unit tests for bin/fix_mnv_codons.py: corrections land only on the SNVs that make the codon.

Run with:  pytest tests/test_fix_mnv_codons.py

Reference ATG AAA CGT TTT GGG CCC TAA on chr1, one CDS over all of it. Codon 2 (AAA, pos 4-6)
carries A4C + A5T, resolvable to p.K2L. Codon 3 (CGT, pos 7-9) carries C7A + G8A, which is
AA? and ambiguous (K or N) without --reference, so it is MNV_UNRESOLVED. A deletion anchored
on pos 4 and on pos 7, and a low-AF second allele at pos 5, share a POS with a grouped SNV.
"""

import pytest

pytest.importorskip("cyvcf2", reason="cyvcf2 is required by the script under test")

import fix_mnv_codons as fmc  # noqa: E402  (needs cyvcf2, checked above)

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=21>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=ANN,Number=.,Type=String,Description="snpEff annotation">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def _ann(alt, effect, impact, hgvs_c, hgvs_p, cds_pos, aa_pos):
    return (f"{alt}|{effect}|{impact}|G|G|transcript|TX1|protein_coding|1/1|{hgvs_c}|{hgvs_p}|"
            f"{cds_pos}/21|{cds_pos}/21|{aa_pos}/6||")


RECORDS = [
    (4, "A",    "C", 1.0, _ann("C", "missense_variant", "MODERATE", "c.4A>C", "p.K2Q", 4, 2)),
    (4, "AAAC", "A", 1.0, _ann("A", "frameshift_variant", "HIGH", "c.5_7delAAC", "p.K2fs", 5, 2)),
    (5, "A",    "T", 1.0, _ann("T", "missense_variant", "MODERATE", "c.5A>T", "p.K2I", 5, 2)),
    (5, "A",    "G", 0.2, _ann("G", "missense_variant", "MODERATE", "c.5A>G", "p.K2R", 5, 2)),
    (7, "C",    "A", 1.0, _ann("A", "missense_variant", "MODERATE", "c.7C>A", "p.R3S", 7, 3)),
    (7, "CGTT", "C", 1.0, _ann("C", "frameshift_variant", "HIGH", "c.8_10delGTT", "p.R3fs", 8, 3)),
    (8, "G",    "A", 1.0, _ann("A", "missense_variant", "MODERATE", "c.8G>A", "p.R3H", 8, 3)),
]


@pytest.fixture
def written(tmp_path):
    """Run both passes and return {(pos, ref, alt): INFO} from the output VCF."""
    vcf_in, vcf_out = tmp_path / "in.vcf", tmp_path / "out.vcf"
    vcf_in.write_text(HEADER + "".join(
        f"chr1\t{pos}\t.\t{ref}\t{alt}\t60\t.\tAF={af};ANN={ann}\n" for pos, ref, alt, af, ann in RECORDS))
    corrections, unresolved = fmc.compute_corrections(fmc.collect_codon_groups(str(vcf_in), 0.5))
    fmc.write_corrected_vcf(str(vcf_in), str(vcf_out), corrections, unresolved)
    rows = [line.split("\t") for line in vcf_out.read_text().splitlines() if not line.startswith("#")]
    return {(int(r[1]), r[3], r[4]): r[7] for r in rows}


def test_snvs_in_one_codon_get_the_combined_change(written):
    for key in [(4, "A", "C"), (5, "A", "T")]:
        assert "MNV_CORRECTED" in written[key]
        assert "|p.K2L|" in written[key]


@pytest.mark.parametrize("key", [(4, "AAAC", "A"), (5, "A", "G")])
def test_other_records_at_a_corrected_position_are_untouched(written, key):
    original = next(ann for pos, ref, alt, _, ann in RECORDS if (pos, ref, alt) == key)
    assert written[key].endswith(f"ANN={original}")
    assert "MNV_" not in written[key]


def test_unresolved_flag_only_on_the_grouped_snvs(written):
    assert "MNV_UNRESOLVED" in written[(7, "C", "A")]
    assert "MNV_UNRESOLVED" in written[(8, "G", "A")]
    assert "MNV_" not in written[(7, "CGTT", "C")]
