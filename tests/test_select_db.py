"""Unit tests for bin/select_db.py (reference bundle selection from HA/NA BLAST hits).

Run with:  pytest tests/test_select_db.py
"""

import csv

import pytest

# conftest.py puts bin/ on sys.path.
import select_db

SAMPLE = "SAMPLE-001"


def hit(segment, db, subtype, pident=98.0, bitscore=3000, qlen=1700, aligned=1700, sample=SAMPLE):
    """One BLAST row in the select_db -outfmt column order."""
    accession = f"{db}-{segment}"
    return [f"{sample}_{segment}", f"{db}|{accession}|{subtype}|{segment}", pident, aligned,
            qlen, qlen, 1, aligned, 1, aligned, 0.0, bitscore, bitscore]


def run(tmp_path, rows, segments=("PB2", "HA", "NA", "M"), min_identity=80, min_coverage=0.5):
    blast = tmp_path / "blast.tsv"
    blast.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
    query = tmp_path / "query.fasta"
    query.write_text("".join(f">{SAMPLE}_{s}\nACGT\n" for s in segments))
    out_tsv, out_db = tmp_path / "selection.tsv", tmp_path / "db.txt"
    assert select_db.main([
        str(blast), "--query", str(query), "--sample-id", SAMPLE,
        "--min-identity", str(min_identity), "--min-coverage", str(min_coverage),
        "--out-tsv", str(out_tsv), "--out-db", str(out_db),
    ]) == 0
    with out_tsv.open() as handle:
        (row,) = list(csv.DictReader(handle, delimiter="\t"))
    assert out_db.read_text() == row["db"]
    return row


H1N1 = [hit("HA", "wisconsin588", "H1N1"), hit("NA", "wisconsin588", "H1N1", bitscore=2500),
        hit("HA", "chicken2021", "H5N1", pident=60, bitscore=400),
        hit("NA", "chicken2021", "H5N1", pident=85, bitscore=1800)]


def test_assigns_bundle_best_for_both_segments(tmp_path):
    row = run(tmp_path, H1N1)
    assert (row["status"], row["db"], row["subtype"]) == ("assigned", "wisconsin588", "H1N1")
    assert row["HA_ref"] == "wisconsin588-HA"
    assert (row["runner_up_db"], row["runner_up_bitscore"]) == ("chicken2021", "2200")


def test_other_segments_do_not_vote(tmp_path):
    row = run(tmp_path, H1N1 + [hit("PB2", "chicken2021", "H5N1", bitscore=99999)])
    assert row["db"] == "wisconsin588"


def test_cross_segment_hits_ignored(tmp_path):
    # An HA query hitting an NA subject must not count as an NA hit.
    rows = [hit("HA", "wisconsin588", "H1N1"),
            [f"{SAMPLE}_HA", "chicken2021|x|H5N1|NA", 99, 1700, 1700, 1700, 1, 1700, 1, 1700, 0, 9999, 9999]]
    assert run(tmp_path, rows)["status"] == "no_hit"


def test_disagreement_is_unassigned(tmp_path):
    # An H5N5: HA closest to the H5N1 bundle, NA closest to something else.
    rows = [hit("HA", "chicken2021", "H5N1"), hit("NA", "chicken2021", "H5N1", pident=50, bitscore=100),
            hit("HA", "wisconsin588", "H1N1", pident=60, bitscore=300), hit("NA", "wisconsin588", "H1N1", bitscore=200)]
    row = run(tmp_path, rows)
    assert (row["status"], row["db"]) == ("disagree", "")
    assert (row["HA_db"], row["NA_db"]) == ("chicken2021", "wisconsin588")


@pytest.mark.parametrize("kwargs", [dict(pident=79.9), dict(aligned=800)])
def test_below_threshold_is_unassigned(tmp_path, kwargs):
    rows = [hit("HA", "wisconsin588", "H1N1", **kwargs), hit("NA", "wisconsin588", "H1N1")]
    assert run(tmp_path, rows)["status"] == "below_threshold"


def test_missing_segment(tmp_path):
    row = run(tmp_path, [hit("HA", "wisconsin588", "H1N1")], segments=("HA", "PB2"))
    assert (row["status"], row["db"]) == ("missing_segment", "")


def test_no_hits(tmp_path):
    row = run(tmp_path, [])
    assert (row["status"], row["db"], row["runner_up_db"]) == ("no_hit", "", "")


def test_tie_breaks_by_bundle_name(tmp_path):
    rows = [hit("HA", "b", "H1N1"), hit("NA", "b", "H1N1"), hit("HA", "a", "H1N1"), hit("NA", "a", "H1N1")]
    assert run(tmp_path, rows)["db"] == "a"


def test_best_hsp_per_bundle_used(tmp_path):
    rows = [hit("HA", "wisconsin588", "H1N1", bitscore=100), hit("HA", "wisconsin588", "H1N1", bitscore=3000),
            hit("NA", "wisconsin588", "H1N1"), hit("HA", "darwin6", "H3N2", bitscore=2000),
            hit("NA", "darwin6", "H3N2", bitscore=10)]
    row = run(tmp_path, rows)
    assert (row["db"], row["HA_bitscore"]) == ("wisconsin588", "3000")
