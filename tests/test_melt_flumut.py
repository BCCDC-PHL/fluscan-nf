"""CLI tests for bin/melt_flumut.py."""

import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "bin" / "melt_flumut.py"


@pytest.mark.parametrize("wide, expected", [
    ("Sample\nS1\n", "sample\tgene\tmutation\tallele\n"),
    ("Sample\tHA1:K156E\tPB2:E627K\nS1\tE\t-\n",
     "sample\tgene\tmutation\tallele\nS1\tHA1\tK156E\tE\nS1\tPB2\tE627K\t-\n"),
])
def test_melt(tmp_path, wide, expected):
    (tmp_path / "wide.tsv").write_text(wide)
    subprocess.run([sys.executable, SCRIPT_PATH, "-i", tmp_path / "wide.tsv", "-o", tmp_path / "long.tsv"], check=True)
    assert (tmp_path / "long.tsv").read_text() == expected
