"""Unit tests for scripts/biopython_sanity_check.py dataset discovery.

The script compares two workflows' amino-acid calls against an independent Biopython
translation, so it has to find every segment's reference and annotation. It reads the
bundle layout directly: uppercase per-segment directories in the canonical order, each
carrying the `pathogen.json` Nextclade itself reads. A dataset `config.json` is
`scripts/build_dataset.py`'s input manifest and must not be required here.

Only the discovery path is exercised; the comparison needs cyvcf2, pandas and real VCFs.
"""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "biopython_sanity_check.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("biopython_sanity_check", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    # Registered before execution because @dataclass resolves annotations through
    # sys.modules[cls.__module__]; without this the class body raises on import.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


checker = _load_module()

# Canonical order; the segment number in a consensus FASTA tag is the position in it.
SEGMENT_NUMBERS = {
    "PB2": 1, "PB1": 2, "PA": 3, "HA": 4, "NP": 5, "NA": 6, "M": 7, "NS": 8,
}


def write_datasets(root: Path, segments=checker.SEGMENTS) -> Path:
    for segment in segments:
        directory = root / segment
        directory.mkdir(parents=True)
        (directory / "pathogen.json").write_text(
            json.dumps({"attributes": {"reference accession": f"NC_ref_{segment}"}})
        )
    return root


def test_all_eight_segments_are_discovered_without_a_config(tmp_path):
    root = write_datasets(tmp_path / "nextclade")

    segments = checker.load_segment_info(root)

    assert [segment.label for segment in segments] == list(checker.SEGMENTS)
    assert not (root / "config.json").exists()


def test_a_config_json_is_ignored_rather_than_believed(tmp_path):
    """The layout is the interface. A stale or lowercase manifest must not redirect the
    lookup -- macOS would resolve `ha/` against `HA/` and hide the mistake until Linux."""
    root = write_datasets(tmp_path / "nextclade")
    (root / "config.json").write_text(
        json.dumps({"segments": [{"name": "ha", "dir": "ha", "segment_num": 4}]})
    )

    segments = checker.load_segment_info(root)

    assert [segment.directory for segment in segments] == list(checker.SEGMENTS)


def test_segment_numbers_follow_the_canonical_order(tmp_path):
    """The fasta_tag is matched against consensus headers, so a wrong number silently
    finds no sequence for that segment rather than failing."""
    root = write_datasets(tmp_path / "nextclade")

    for segment in checker.load_segment_info(root):
        expected = SEGMENT_NUMBERS[segment.label]
        assert segment.fasta_tag == f"segment{expected}_{segment.label}"


def test_reference_accession_comes_from_each_pathogen_json(tmp_path):
    root = write_datasets(tmp_path / "nextclade")

    segments = checker.load_segment_info(root)

    assert {s.label: s.nextclade_contig for s in segments} == {
        segment: f"NC_ref_{segment}" for segment in checker.SEGMENTS
    }


def test_a_missing_segment_directory_is_named(tmp_path):
    root = write_datasets(tmp_path / "nextclade")
    (root / "NA" / "pathogen.json").unlink()
    (root / "NA").rmdir()

    with pytest.raises(SystemExit, match="Segment dataset directory not found"):
        checker.load_segment_info(root)


def test_a_dataset_without_pathogen_json_is_named(tmp_path):
    root = write_datasets(tmp_path / "nextclade")
    (root / "M" / "pathogen.json").unlink()

    with pytest.raises(SystemExit, match="pathogen.json not found for segment M"):
        checker.load_segment_info(root)
