"""Unit tests for bin/prepare_consensus.py (consensus FASTA preparation).

Covers segment-label parsing across the accepted header forms, sequence cleaning
and its counts, the PASS/WARN/FAIL/MISSING QC statuses, and the CLI: output paths
are the ones the caller names, and headers are rewritten to the supplied
--sample-id no matter what the input headers claim.

Run with:  pytest tests/test_prepare_consensus.py
"""

import csv
import sys

import pytest
from Bio import SeqIO

# conftest.py puts bin/ on sys.path.
import prepare_consensus as prep

# Lengths inside each segment's accepted range, so a test sequence passes on length.
IN_RANGE = {"PB2": 2000, "PB1": 2000, "PA": 2200, "HA": 1700,
            "NP": 1500, "NA": 1400, "M": 1000, "NS": 900}


def write_fasta(path, records):
    """records: list of (header, sequence)."""
    with open(path, "w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n{seq}\n")
    return path


def run_cli(tmp_path, records, sample_id="SAMPLE-001", extra=()):
    """Run the CLI over records, returning (exit_code, fasta_path, qc_path)."""
    fasta_in = write_fasta(tmp_path / "in.fasta", records)
    fasta_out = tmp_path / f"{sample_id}.fasta"
    qc_out = tmp_path / f"{sample_id}_qc_report.tsv"
    code = prep.main([
        "--input", str(fasta_in),
        "--sample-id", sample_id,
        "--output-fasta", str(fasta_out),
        "--output-qc", str(qc_out),
        *extra,
    ])
    return code, fasta_out, qc_out


def qc_by_segment(qc_path):
    with open(qc_path) as handle:
        return {row["segment"]: row for row in csv.DictReader(handle, delimiter="\t")}


# --------------------------------------------------------------------------- #
# Segment labels
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("header,expected", [
    ("WIN-AH-2024-FAV-0415-1-OS_segment1_PB2", "PB2"),
    ("WIN-AH-2024-FAV-0415-1-OS_PB2", "PB2"),
    ("WIN-AH-2022-FAV-1687-1_OS_2_segment4_HA", "HA"),
    ("EPI_ISL_20399898|A_/_H5N1|NA|6|EPI5261179|2.3.4.4b_NA", "NA"),
    ("EPI_ISL_20399898|A_/_H5N1|HA|4|EPI5261179|2.3.4.4b", "HA"),
    ("sample_mp", "M"),
    ("sample_MP", "M"),
])
def test_segment_of_recognised(header, expected):
    assert prep.segment_of(header) == expected


@pytest.mark.parametrize("header", ["sample_UNKNOWN", "no-delimiters-here", ""])
def test_segment_of_unrecognised(header):
    assert prep.segment_of(header) is None


# --------------------------------------------------------------------------- #
# Sample stems (the one-file-one-sample check)
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("header,expected", [
    ("WIN-AH-2024-FAV-0415-1-OS_segment1_PB2", "WIN-AH-2024-FAV-0415-1-OS"),
    ("WIN-AH-2024-FAV-0415-1-OS_PB2", "WIN-AH-2024-FAV-0415-1-OS"),
    ("WIN-AH-2022-FAV-1687-1_OS_2_segment4_HA", "WIN-AH-2022-FAV-1687-1"),
    ("EPI_ISL_20399898|A_/_H5N1|NA|6|EPI5261179|2.3.4.4b_NA", "EPI_ISL_20399898"),
])
def test_sample_stem_is_constant_across_a_sample(header, expected):
    assert prep.sample_stem(header) == expected


# --------------------------------------------------------------------------- #
# Cleaning
# --------------------------------------------------------------------------- #

def clean(raw, *, mask_degenerate=True):
    record, qc = prep.process_segment("S", "HA", raw, mask_degenerate=mask_degenerate,
                                      max_n_frac=1.0)
    return str(record.seq), qc


def test_sequence_is_uppercased():
    sequence, _ = clean("acgtacgt")
    assert sequence == "ACGTACGT"


def test_degenerate_bases_masked():
    sequence, qc = clean("ACGT" + "".join(sorted(prep.DEGENERATE_BASES)))
    assert sequence == "ACGT" + "N" * len(prep.DEGENERATE_BASES)
    assert qc.n_replaced == len(prep.DEGENERATE_BASES)


def test_degenerate_bases_kept_when_not_masking():
    sequence, qc = clean("ACGTRY", mask_degenerate=False)
    assert sequence == "ACGTRY"
    assert qc.n_replaced == 0


def test_gaps_are_stripped():
    sequence, qc = clean("ACG-TAC.GT")
    assert sequence == "ACGTACGT"
    assert qc.n_gaps_removed == 2
    assert qc.n_invalid == 0


def test_unknown_characters_become_n():
    sequence, qc = clean("ACGT123X")
    assert sequence == "ACGTNNNN"
    assert qc.n_invalid == 4


def test_n_proportions_before_and_after():
    _, qc = clean("ACGTRYNN")
    assert qc.n_prop_before == pytest.approx(2 / 8)
    assert qc.n_prop_after == pytest.approx(4 / 8)


# --------------------------------------------------------------------------- #
# QC statuses
# --------------------------------------------------------------------------- #

def qc_for(segment, raw, **kwargs):
    kwargs.setdefault("mask_degenerate", True)
    kwargs.setdefault("max_n_frac", 1.0)
    return prep.process_segment("SAMPLE-001", segment, raw, **kwargs)[1]


def test_status_pass_in_range():
    qc = qc_for("HA", "A" * IN_RANGE["HA"])
    assert (qc.status, qc.length_ok, qc.notes) == (prep.STATUS_PASS, True, "")


def test_status_fail_outside_length_range():
    qc = qc_for("HA", "A" * 100)
    assert qc.status == prep.STATUS_FAIL
    assert qc.length_ok is False
    assert "outside expected range" in qc.notes


def test_status_fail_all_n():
    qc = qc_for("HA", "N" * IN_RANGE["HA"])
    assert qc.status == prep.STATUS_FAIL
    assert "all-N" in qc.notes


def test_status_warn_above_max_n_frac():
    raw = "N" * 170 + "A" * (IN_RANGE["HA"] - 170)          # 10% N
    qc = qc_for("HA", raw, max_n_frac=0.05)
    assert qc.status == prep.STATUS_WARN
    assert "exceeds" in qc.notes


def test_status_warn_on_stripped_gaps():
    qc = qc_for("HA", "-" * 5 + "A" * IN_RANGE["HA"])
    assert qc.status == prep.STATUS_WARN
    assert qc.n_gaps_removed == 5


def test_seq_id_uses_sample_id():
    record, qc = prep.process_segment("SAMPLE-001", "HA", "A" * IN_RANGE["HA"],
                                      mask_degenerate=True, max_n_frac=1.0)
    assert record.id == "SAMPLE-001_HA"
    assert record.description == ""
    assert qc.seq_id == "SAMPLE-001_HA"


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def test_cli_writes_the_named_outputs(tmp_path):
    code, fasta_out, qc_out = run_cli(tmp_path, [("anything_HA", "A" * IN_RANGE["HA"])])
    assert code == prep.EXIT_OK
    assert fasta_out.exists() and qc_out.exists()


def test_cli_headers_come_from_sample_id_not_the_input(tmp_path):
    """The input headers name a different isolate; the output must ignore that."""
    records = [
        ("EPI_ISL_20399898|A_/_H5N1|PB2|1|EPI5261178|2.3.4.4b_PB2", "A" * IN_RANGE["PB2"]),
        ("EPI_ISL_20399898|A_/_H5N1|NA|6|EPI5261179|2.3.4.4b_NA", "C" * IN_RANGE["NA"]),
    ]
    _, fasta_out, _ = run_cli(tmp_path, records, sample_id="SAMPLE-ABC")
    assert [rec.id for rec in SeqIO.parse(fasta_out, "fasta")] == \
        ["SAMPLE-ABC_PB2", "SAMPLE-ABC_NA"]


def test_cli_reports_every_expected_segment(tmp_path):
    _, _, qc_out = run_cli(tmp_path, [("x_HA", "A" * IN_RANGE["HA"])])
    rows = qc_by_segment(qc_out)
    assert list(rows) == prep.EXPECTED_SEGMENTS
    assert rows["HA"]["status"] == prep.STATUS_PASS
    assert rows["PB2"]["status"] == prep.STATUS_MISSING
    assert rows["PB2"]["notes"] == "segment not found in input"


def test_cli_qc_proportions_are_rounded(tmp_path):
    raw = "N" * 7 + "A" * (IN_RANGE["HA"] - 7)
    _, _, qc_out = run_cli(tmp_path, [("x_HA", raw)])
    assert qc_by_segment(qc_out)["HA"]["n_prop_after"] == "0.004118"


def test_cli_duplicate_segment_is_fatal(tmp_path):
    """Two records for one segment means two samples in one file."""
    records = [("A-1_HA", "A" * IN_RANGE["HA"]), ("B-2_HA", "C" * IN_RANGE["HA"])]
    code, fasta_out, qc_out = run_cli(tmp_path, records)
    assert code == prep.EXIT_INPUT_ERROR
    assert not fasta_out.exists() and not qc_out.exists()


def test_cli_disagreeing_headers_are_fatal(tmp_path):
    """Eight distinct segments from two samples must not be relabelled as one."""
    records = [("A-1_PB2", "A" * IN_RANGE["PB2"]), ("B-2_HA", "C" * IN_RANGE["HA"])]
    code, fasta_out, qc_out = run_cli(tmp_path, records)
    assert code == prep.EXIT_INPUT_ERROR
    assert not fasta_out.exists() and not qc_out.exists()


def test_cli_no_recognisable_segments_is_fatal(tmp_path):
    code, fasta_out, _ = run_cli(tmp_path, [("mystery-header", "ACGT" * 400)])
    assert code == prep.EXIT_INPUT_ERROR
    assert not fasta_out.exists()


def test_cli_no_mask_degenerate(tmp_path):
    raw = "R" * 10 + "A" * (IN_RANGE["HA"] - 10)
    _, fasta_out, qc_out = run_cli(tmp_path, [("x_HA", raw)], extra=["--no-mask-degenerate"])
    assert str(next(SeqIO.parse(fasta_out, "fasta")).seq).startswith("R" * 10)
    assert qc_by_segment(qc_out)["HA"]["n_replaced"] == "0"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
