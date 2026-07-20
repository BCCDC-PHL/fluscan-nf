"""Unit tests for bin/consensus_identity_qc.py (per-segment identity QC).

Covers the QC scenarios required by the consensus-identity QC plan: PASS, WARN,
FAIL, a missing/empty consensus segment, a missing/empty reference segment,
ambiguous bases, single-record FASTA loading, and the end-to-end CLI (TSV
schema + single data row).

Identity in the controlled tests is driven by replacing reference bases with
ambiguous ``N`` rather than by substituting other bases. ``N`` never
exact-matches a base, and the surrounding bases stay matched in place, so the
ungapped global alignment is optimal and identity is exactly ``(L - k) / L``.
Substituting real bases instead lets the global aligner recover spurious
matches via gaps, which makes identity non-deterministic.

Run with:  pytest tests/test_consensus_identity_qc.py
"""

import csv
import importlib.util
import random
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "bin" / "consensus_identity_qc.py"

PASS_PCT = 90.0
WARN_PCT = 70.0


def _load_module():
    spec = importlib.util.spec_from_file_location("consensus_identity_qc", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


qc = _load_module()

# A deterministic, non-repetitive 300 bp reference used across identity tests.
REF = "".join(random.Random(42).choice("ACGT") for _ in range(300))


def with_ns(seq, k):
    """Return seq with the first k bases replaced by ambiguous N.

    Yields identity exactly (len(seq) - k) / len(seq) against the original seq.
    """
    return "N" * k + seq[k:]


def write_fasta(path, records):
    """records: list of (header, sequence). Empty list -> empty file."""
    with open(path, "w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n{seq}\n")
    return path


# --------------------------------------------------------------------------
# alignment_identity_pct: core scoring
# --------------------------------------------------------------------------

def test_identity_identical_is_100():
    assert qc.alignment_identity_pct(REF, REF) == pytest.approx(100.0)


def test_identity_80pct():
    # 60 ambiguous bases over 300 columns -> exactly 80% identity.
    assert qc.alignment_identity_pct(with_ns(REF, 60), REF) == pytest.approx(80.0)


def test_identity_50pct():
    assert qc.alignment_identity_pct(with_ns(REF, 150), REF) == pytest.approx(50.0)


def test_identity_empty_sequence_is_zero():
    assert qc.alignment_identity_pct("", REF) == 0.0
    assert qc.alignment_identity_pct(REF, "") == 0.0


def test_ambiguous_bases_count_as_mismatch():
    # 75 N substitutions -> 75% identity (N never exact-matches a base) -> WARN.
    identity = qc.alignment_identity_pct(with_ns(REF, 75), REF)
    assert identity == pytest.approx(75.0)
    assert qc.segment_status(identity, PASS_PCT, WARN_PCT) == "WARN"


# --------------------------------------------------------------------------
# segment_status thresholds
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "identity,expected",
    [(100.0, "PASS"), (90.0, "PASS"), (89.999, "WARN"), (70.0, "WARN"), (69.999, "FAIL"), (0.0, "FAIL")],
)
def test_segment_status_thresholds(identity, expected):
    assert qc.segment_status(identity, PASS_PCT, WARN_PCT) == expected


# --------------------------------------------------------------------------
# load_single_sequence
# --------------------------------------------------------------------------

def test_load_single_sequence_reads_first_record(tmp_path):
    path = write_fasta(tmp_path / "c.fasta", [("SAMPLE_HA", REF)])
    rec_id, seq = qc.load_single_sequence(path)
    assert rec_id == "SAMPLE_HA"
    assert seq == REF


def test_load_single_sequence_empty_file(tmp_path):
    path = write_fasta(tmp_path / "empty.fasta", [])
    rec_id, seq = qc.load_single_sequence(path)
    assert rec_id is None
    assert seq == ""


def test_load_single_sequence_missing_file(tmp_path):
    rec_id, seq = qc.load_single_sequence(tmp_path / "does_not_exist.fasta")
    assert rec_id is None
    assert seq == ""


# --------------------------------------------------------------------------
# build_qc_row: PASS / WARN / FAIL / missing
# --------------------------------------------------------------------------

def build_row(consensus, reference, segment="HA"):
    return qc.build_qc_row(
        sample_id="SAMPLE",
        segment=segment,
        consensus=consensus,
        reference=reference,
        pass_pct=PASS_PCT,
        warn_pct=WARN_PCT,
    )


def test_row_pass():
    row = build_row(("SAMPLE_HA", REF), ("REF_HA", REF))
    assert row["status"] == "PASS"
    assert row["identity_pct"] == pytest.approx(100.0)
    assert row["segment"] == "HA"
    assert row["consensus_length"] == 300
    assert row["reference_length"] == 300
    assert row["reference_id"] == "REF_HA"


def test_row_warn():
    row = build_row(("SAMPLE_HA", with_ns(REF, 60)), ("REF_HA", REF))
    assert row["status"] == "WARN"
    assert row["identity_pct"] == pytest.approx(80.0)


def test_row_fail_low_identity():
    row = build_row(("SAMPLE_HA", with_ns(REF, 150)), ("REF_HA", REF))
    assert row["status"] == "FAIL"
    assert row["identity_pct"] == pytest.approx(50.0)


def test_row_missing_consensus():
    row = build_row((None, ""), ("REF_HA", REF))
    assert row["status"] == "FAIL"
    assert row["identity_pct"] is None
    assert row["consensus_length"] == 0
    assert row["reference_length"] == 300


def test_row_missing_reference():
    row = build_row(("SAMPLE_HA", REF), (None, ""))
    assert row["status"] == "FAIL"
    assert row["identity_pct"] is None
    assert row["consensus_length"] == 300
    assert row["reference_length"] == 0
    assert row["reference_id"] is None


# --------------------------------------------------------------------------
# End-to-end CLI: schema + single data row via subprocess
# --------------------------------------------------------------------------

def run_cli(tmp_path, consensus_records, reference_records, segment="HA"):
    consensus = write_fasta(tmp_path / "consensus.fasta", consensus_records)
    reference = write_fasta(tmp_path / "reference.fasta", reference_records)
    out_tsv = tmp_path / "out.tsv"
    subprocess.run(
        [
            sys.executable, str(SCRIPT_PATH),
            "--consensus", str(consensus),
            "--reference", str(reference),
            "--sample-id", "SAMPLE",
            "--segment", segment,
            "--pass-pct", str(PASS_PCT),
            "--warn-pct", str(WARN_PCT),
            "--out-tsv", str(out_tsv),
        ],
        check=True,
    )
    with open(out_tsv) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows


def test_cli_schema_and_single_row(tmp_path):
    rows = run_cli(tmp_path, [("SAMPLE_HA", REF)], [("REF_HA", REF)])
    assert len(rows) == 1
    row = rows[0]
    assert list(row.keys()) == list(qc.OUTPUT_COLUMNS)
    assert row["sample_id"] == "SAMPLE"
    assert row["segment"] == "HA"
    assert row["status"] == "PASS"
    assert row["identity_pct"] == "100.000"
    assert row["reference_id"] == "REF_HA"


@pytest.mark.parametrize(
    "n_ns,expected_status,expected_pct",
    [(0, "PASS", "100.000"), (60, "WARN", "80.000"), (150, "FAIL", "50.000")],
)
def test_cli_pass_warn_fail(tmp_path, n_ns, expected_status, expected_pct):
    rows = run_cli(tmp_path, [("SAMPLE_NA", with_ns(REF, n_ns))], [("REF_NA", REF)], segment="NA")
    assert rows[0]["status"] == expected_status
    assert rows[0]["identity_pct"] == expected_pct


def test_cli_missing_consensus_segment(tmp_path):
    # Empty consensus FASTA -> still writes a row, FAIL/NA (report-only).
    rows = run_cli(tmp_path, [], [("REF_PB2", REF)], segment="PB2")
    assert rows[0]["status"] == "FAIL"
    assert rows[0]["identity_pct"] == "NA"
    assert rows[0]["consensus_length"] == "0"


def test_cli_missing_reference_segment(tmp_path):
    rows = run_cli(tmp_path, [("SAMPLE_NS", REF)], [], segment="NS")
    assert rows[0]["status"] == "FAIL"
    assert rows[0]["identity_pct"] == "NA"
    assert rows[0]["reference_length"] == "0"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
