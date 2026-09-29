"""Unit tests for bin/nextclade_post_process.py: per-sample splits and identity QC."""

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

import nextclade_post_process as npp

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "bin" / "nextclade_post_process.py"

TSV_FIELDS = ["index", "seqName", "alignmentStart", "alignmentEnd", "totalSubstitutions",
              "totalDeletions", "totalInsertions", "totalMissing", "totalNonACGTNs", "errors"]


def _row(name, start=1, end=100, subs=0, dels=0, ins=0, missing=0, non_acgtn=0, errors=""):
    return {"index": 0, "seqName": name, "alignmentStart": start, "alignmentEnd": end,
            "totalSubstitutions": subs, "totalDeletions": dels, "totalInsertions": ins,
            "totalMissing": missing, "totalNonACGTNs": non_acgtn, "errors": errors}


def _as_tsv_row(row):
    return {k: str(v) for k, v in row.items()}


# ---------------------------------------------------------------------------
# identity_pct
# ---------------------------------------------------------------------------

def test_full_length_identical_is_100():
    assert npp.identity_pct(_as_tsv_row(_row("S")), 100) == 100.0


def test_substitutions_missing_and_ambiguous_bases_are_mismatches():
    row = _as_tsv_row(_row("S", subs=5, missing=3, non_acgtn=2))
    assert npp.identity_pct(row, 100) == 90.0


def test_insertions_widen_the_alignment():
    row = _as_tsv_row(_row("S", ins=25))
    assert npp.identity_pct(row, 100) == 80.0


def test_deletions_and_unaligned_ends_are_mismatches():
    # 10 reference bases before the query starts and 10 deleted inside it.
    row = _as_tsv_row(_row("S", start=11, end=100, dels=10))
    assert npp.identity_pct(row, 100) == 80.0


def test_blank_counts_read_as_zero():
    row = _as_tsv_row(_row("S"))
    row["totalMissing"] = ""
    assert npp.identity_pct(row, 100) == 100.0


# ---------------------------------------------------------------------------
# write_identity_qc
# ---------------------------------------------------------------------------

@pytest.fixture
def batch(tmp_path):
    rows = [
        _row("PASSER"),
        _row("WARNER", subs=20),
        _row("FAILER", subs=40),
        _row("BOUNDARY", subs=10),
        {**_row("UNALIGNED"), "alignmentStart": "", "alignmentEnd": "",
         "errors": "Unable to align"},
    ]
    tsv = tmp_path / "batch.tsv"
    with open(tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TSV_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    fasta = tmp_path / "batch.fasta"
    fasta.write_text("".join(f">{r['seqName']}\n{'A' * 60}\n{'A' * 38}\n" for r in rows))
    reference = tmp_path / "reference.fasta"
    reference.write_text(">REF1 some description\n" + "A" * 100 + "\n")
    return tsv, fasta, reference


def _identity_rows(tmp_path, monkeypatch, batch):
    monkeypatch.chdir(tmp_path)
    npp.write_identity_qc(*batch, "HA", pass_pct=90, warn_pct=70)
    return {p.name.split("_")[0]: next(csv.DictReader(open(p), delimiter="\t"))
            for p in tmp_path.glob("*_HA_identity_qc.tsv")}


def test_statuses_follow_the_thresholds(tmp_path, monkeypatch, batch):
    rows = _identity_rows(tmp_path, monkeypatch, batch)
    assert {name: row["status"] for name, row in rows.items()} == {
        "PASSER": "PASS", "WARNER": "WARN", "FAILER": "FAIL", "BOUNDARY": "PASS", "UNALIGNED": "FAIL",
    }


def test_row_schema_and_values(tmp_path, monkeypatch, batch):
    row = _identity_rows(tmp_path, monkeypatch, batch)["WARNER"]
    assert list(row) == npp.IDENTITY_COLUMNS
    assert row == {"sample_id": "WARNER", "segment": "HA", "status": "WARN",
                   "identity_pct": "80.000", "consensus_length": "98",
                   "reference_length": "100", "reference_id": "REF1"}


def test_unaligned_sample_has_no_identity(tmp_path, monkeypatch, batch):
    assert _identity_rows(tmp_path, monkeypatch, batch)["UNALIGNED"]["identity_pct"] == "NA"


# ---------------------------------------------------------------------------
# CLI: splits, status summary and identity QC together
# ---------------------------------------------------------------------------

def test_cli_writes_every_per_sample_output(tmp_path, batch):
    tsv, fasta, reference = batch
    results = [{"seqName": name} for name in ("PASSER", "WARNER", "FAILER", "BOUNDARY")]
    (tmp_path / "batch.json").write_text(json.dumps({
        "results": results, "errors": [{"seqName": "UNALIGNED", "errors": ["Unable to align"]}],
    }))
    out = tmp_path / "out"
    out.mkdir()
    subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--json", str(tmp_path / "batch.json"), "--tsv", str(tsv),
         "--segment", "HA", "--fasta", str(fasta), "--reference", str(reference),
         "--pass-pct", "90", "--warn-pct", "70", "--status-summary", "HA_nextclade_status.tsv"],
        check=True, capture_output=True, text=True, cwd=out,
    )
    names = {p.name for p in out.iterdir()}
    assert "UNALIGNED_HA_nextclade.json" not in names
    assert "UNALIGNED_HA_nextclade.tsv" in names
    assert "UNALIGNED_HA_identity_qc.tsv" in names
    assert {f"{s}_HA_nextclade.json" for s in ("PASSER", "WARNER", "FAILER", "BOUNDARY")} <= names
    status = {r["sample"]: r["status"] for r in csv.DictReader(open(out / "HA_nextclade_status.tsv"), delimiter="\t")}
    assert status["UNALIGNED"] == "fail" and status["PASSER"] == "pass"
