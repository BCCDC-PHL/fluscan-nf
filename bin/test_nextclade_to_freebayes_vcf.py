"""Unit tests for nextclade_to_freebayes_vcf.py helper functions."""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from nextclade_to_freebayes_vcf import _insertion_alleles, _deletion_alleles


REF = "GCATGCAGGCAATCCT"  # 16-nt toy reference; first base G, last base T (distinct)


class TestInsertionAlleles:
    def test_interior(self):
        """Standard in-body insertion: anchor is the base at pos, POS is 1-based."""
        vcf_pos, ref_a, alt_a = _insertion_alleles(REF, 4, "AAA")
        assert vcf_pos == 5
        assert ref_a == REF[4]          # 'G'
        assert alt_a == REF[4] + "AAA"  # 'GAAA'

    def test_leading_pos_minus_one(self):
        """5' insertion: nextclade pos=-1 must not give POS=0 or use ref[-1] as anchor."""
        vcf_pos, ref_a, alt_a = _insertion_alleles(REF, -1, "NNNAGCAGG")
        # POS must be valid (>= 1)
        assert vcf_pos == 1, f"expected POS=1, got {vcf_pos}"
        # REF must be the actual first base, not the last
        assert ref_a == REF[0], f"expected ref[0]={REF[0]!r}, got {ref_a!r}"
        assert ref_a != REF[-1], "REF must not be the last base (Python -1 wraparound)"
        # ALT = ins_seq right-anchored before the first base
        assert alt_a == "NNNAGCAGG" + REF[0]

    def test_last_position(self):
        """Insertion after the last reference base: standard left-anchor."""
        last = len(REF) - 1
        vcf_pos, ref_a, alt_a = _insertion_alleles(REF, last, "GG")
        assert vcf_pos == len(REF)
        assert ref_a == REF[last]
        assert alt_a == REF[last] + "GG"


class TestDeletionAlleles:
    def test_interior(self):
        vcf_pos, ref_a, alt_a = _deletion_alleles(REF, 3, 6)
        assert vcf_pos == 3           # anchor = begin-1 = 2, POS = 3
        assert ref_a == REF[2] + REF[3:6]
        assert alt_a == REF[2]

    def test_leading(self):
        """Deletion at the start (begin=0): already guarded in the existing code."""
        vcf_pos, ref_a, alt_a = _deletion_alleles(REF, 0, 3)
        assert vcf_pos >= 1


if __name__ == "__main__":
    import traceback
    suites = [TestInsertionAlleles, TestDeletionAlleles]
    passed = failed = 0
    for cls in suites:
        obj = cls()
        for name in [n for n in dir(cls) if n.startswith("test_")]:
            try:
                getattr(obj, name)()
                print(f"  PASS  {cls.__name__}.{name}")
                passed += 1
            except Exception:
                print(f"  FAIL  {cls.__name__}.{name}")
                traceback.print_exc()
                failed += 1
    print(f"\n{passed} passed, {failed} failed")
    sys.exit(failed)
