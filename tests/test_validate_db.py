"""Unit tests for scripts/validate_db.py.

The bundle validator exists to catch reference databases that resolve cleanly and then
annotate nothing, so most assertions here are that a *silent* defect is reported: a
snpEff genome with no built database, a BLAST header that does not parse, a Nextclade
dataset directory that is missing a segment. The counterpart assertions are that a
legitimately partial bundle -- the FluMut bundle carries only nextclade/ -- stays valid,
because a validator that fails on those would just be turned off.

Run with:  pytest tests/test_validate_db.py
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "validate_db.py"

SEGMENTS = ["pb2", "pb1", "pa", "ha", "np", "na", "m", "ns"]


def _load_module():
    spec = importlib.util.spec_from_file_location("validate_db", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    # Registered before exec: @dataclass resolves the `from __future__` string
    # annotations through sys.modules[cls.__module__], which is None otherwise.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


vdb = _load_module()


# ---------------------------------------------------------------------------
# Bundle fixtures
# ---------------------------------------------------------------------------

BLAST_FASTA = "".join(
    f">NC_00735{i}.1|H5N1|{seg.upper()}\nACGTACGTACGT\n"
    for i, seg in enumerate(SEGMENTS)
)

SNPEFF_CONFIG = """\
data.dir = ./data

influenza_h5n1.genome : Influenza_A_H5N1
influenza_h5n1.chromosomes : NC_007357.1
influenza_h5n1.NC_007357.1.codonTable : Standard
"""


def make_nextclade(root: Path, segments=SEGMENTS) -> None:
    for segment in segments:
        directory = root / "nextclade" / segment
        directory.mkdir(parents=True)
        (directory / "reference.fasta").write_text(f">NC_ref_{segment}\nACGT\n")
        (directory / "genome_annotation.gff3").write_text("##gff-version 3\n")
        (directory / "pathogen.json").write_text("{}\n")


def make_snpeff(root: Path, genomes=("influenza_h5n1",), build=True) -> None:
    directory = root / "snpeff"
    directory.mkdir(parents=True)
    config = "data.dir = ./data\n\n" + "".join(
        f"{genome}.genome : {genome}\n" for genome in genomes
    )
    (directory / "snpEff.config").write_text(config)
    for genome in genomes:
        data = directory / "data" / genome
        data.mkdir(parents=True)
        (data / "genes.gbk").write_text("LOCUS test\n//\n")
        if build:
            (data / "snpEffectPredictor.bin").write_bytes(b"\x00")


@pytest.fixture
def complete_bundle(tmp_path: Path) -> Path:
    root = tmp_path / "gsgd"
    root.mkdir()
    (root / "blastn_db.fa").write_text(BLAST_FASTA)
    (root / "db.yml").write_text(
        'name: gsgd\nversion: "2025-07-07"\nsnpeff_genomes: [influenza_h5n1]\n'
    )
    make_nextclade(root)
    make_snpeff(root)
    return root


@pytest.fixture
def nextclade_only_bundle(tmp_path: Path) -> Path:
    """Shaped like the real FluMut bundle: datasets only, no BLAST or SnpEff."""
    root = tmp_path / "flumut"
    root.mkdir()
    make_nextclade(root)
    return root


def run_cli(root: Path, *extra: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), str(root), *extra],
        capture_output=True, text=True,
    )


# ---------------------------------------------------------------------------
# Valid bundles
# ---------------------------------------------------------------------------

def test_complete_bundle_is_valid(complete_bundle):
    report = vdb.validate(complete_bundle, [])
    assert report.ok, report.problems
    assert report.components["blastn"].startswith("OK")
    assert report.components["nextclade"].startswith("OK")
    assert report.components["snpeff"].startswith("OK")


def test_partial_bundle_is_valid_and_reports_skips(nextclade_only_bundle):
    """A bundle only the Nextclade workflow consumes must not fail for what it omits."""
    report = vdb.validate(nextclade_only_bundle, [])
    assert report.ok, report.problems
    assert report.components["blastn"].startswith("SKIP")
    assert report.components["snpeff"].startswith("SKIP")
    assert report.components["nextclade"].startswith("OK")


def test_require_fails_on_absent_component(nextclade_only_bundle):
    report = vdb.validate(nextclade_only_bundle, ["snpeff"])
    assert not report.ok
    assert any("--require snpeff" in p for p in report.problems)


def test_missing_bundle_is_reported(tmp_path):
    report = vdb.validate(tmp_path / "nope", [])
    assert not report.ok
    assert any("does not exist" in p for p in report.problems)


# ---------------------------------------------------------------------------
# SnpEff: the silent-failure case this validator exists for
# ---------------------------------------------------------------------------

def test_declared_genome_without_built_database_fails(complete_bundle):
    """snpEff would emit empty ANN= rather than error, so nothing downstream notices."""
    (complete_bundle / "snpeff" / "data" / "influenza_h5n1" / "snpEffectPredictor.bin").unlink()
    report = vdb.validate(complete_bundle, [])
    assert not report.ok
    assert any("no built database" in p for p in report.problems)
    assert report.components["snpeff"] == "FAIL"


def test_second_genome_without_built_database_fails(tmp_path):
    root = tmp_path / "multi"
    root.mkdir()
    make_snpeff(root, genomes=("influenza_h5n1", "influenza_h1n1"))
    (root / "snpeff" / "data" / "influenza_h1n1" / "snpEffectPredictor.bin").unlink()
    report = vdb.validate(root, [])
    assert not report.ok
    assert any("influenza_h1n1" in p for p in report.problems)
    assert not any("influenza_h5n1," in p for p in report.problems)


def test_snpeff_dir_without_config_fails(complete_bundle):
    (complete_bundle / "snpeff" / "snpEff.config").unlink()
    report = vdb.validate(complete_bundle, [])
    assert not report.ok
    assert any("snpEff.config is missing" in p for p in report.problems)


def test_config_with_no_genome_entry_fails(complete_bundle):
    (complete_bundle / "snpeff" / "snpEff.config").write_text("data.dir = ./data\n")
    report = vdb.validate(complete_bundle, [])
    assert not report.ok
    assert any("declares no" in p for p in report.problems)


def test_commented_genome_entry_is_not_counted(complete_bundle):
    (complete_bundle / "snpeff" / "snpEff.config").write_text(
        "data.dir = ./data\n# influenza_h9n2.genome : old\ninfluenza_h5n1.genome : x\n"
    )
    report = vdb.validate(complete_bundle, [])
    assert report.ok, report.problems
    assert "influenza_h9n2" not in report.components["snpeff"]


def test_stale_absolute_data_dir_warns_but_does_not_fail(complete_bundle):
    """-dataDir overrides it at runtime, so this is a moved-bundle smell, not a failure."""
    (complete_bundle / "snpeff" / "snpEff.config").write_text(
        "data.dir = /nonexistent/snpeff/gsgd\ninfluenza_h5n1.genome : x\n"
    )
    report = vdb.validate(complete_bundle, [])
    assert report.ok, report.problems
    assert any("data.dir points outside the bundle" in w for w in report.warnings)


# ---------------------------------------------------------------------------
# BLAST database
# ---------------------------------------------------------------------------

def test_malformed_blast_header_fails(complete_bundle):
    with (complete_bundle / "blastn_db.fa").open("a") as handle:
        handle.write(">no_pipes_here\nACGT\n")
    report = vdb.validate(complete_bundle, [])
    assert not report.ok
    assert any("accession|subtype|segment" in p for p in report.problems)


def test_empty_blast_database_fails(complete_bundle):
    (complete_bundle / "blastn_db.fa").write_text("")
    report = vdb.validate(complete_bundle, [])
    assert not report.ok
    assert any("no sequences" in p for p in report.problems)


def test_missing_segment_in_blast_database_warns_only(complete_bundle):
    """A subtype-specific database need not carry all eight segments."""
    kept = [line for line in BLAST_FASTA.splitlines(keepends=True)]
    (complete_bundle / "blastn_db.fa").write_text("".join(kept[:-2]))
    report = vdb.validate(complete_bundle, [])
    assert report.ok, report.problems
    assert any("no sequence for segment" in w for w in report.warnings)


# ---------------------------------------------------------------------------
# Nextclade datasets
# ---------------------------------------------------------------------------

def test_missing_segment_directory_fails(complete_bundle):
    import shutil
    shutil.rmtree(complete_bundle / "nextclade" / "pb1")
    report = vdb.validate(complete_bundle, [])
    assert not report.ok
    assert any("PB1" in p for p in report.problems)


def test_segment_missing_annotation_fails(complete_bundle):
    (complete_bundle / "nextclade" / "ha" / "genome_annotation.gff3").unlink()
    report = vdb.validate(complete_bundle, [])
    assert not report.ok
    assert any("genome_annotation.gff3" in p for p in report.problems)


def test_reference_that_is_not_fasta_fails(complete_bundle):
    (complete_bundle / "nextclade" / "na" / "reference.fasta").write_text("ACGT\n")
    report = vdb.validate(complete_bundle, [])
    assert not report.ok
    assert any("not FASTA" in p for p in report.problems)


def test_missing_pathogen_json_warns_only(complete_bundle):
    (complete_bundle / "nextclade" / "m" / "pathogen.json").unlink()
    report = vdb.validate(complete_bundle, [])
    assert report.ok, report.problems
    assert any("pathogen.json" in w for w in report.warnings)


# ---------------------------------------------------------------------------
# db.yml
# ---------------------------------------------------------------------------

def test_absent_db_yml_is_valid_with_basename_fallback(nextclade_only_bundle):
    report = vdb.validate(nextclade_only_bundle, [])
    assert report.ok
    assert report.metadata == {}
    assert "no db.yml" in " ".join(report.warnings)


def test_malformed_db_yml_fails(complete_bundle):
    (complete_bundle / "db.yml").write_text("name: [unclosed\n")
    report = vdb.validate(complete_bundle, [])
    assert not report.ok
    assert any("not valid YAML" in p for p in report.problems)


def test_non_mapping_db_yml_fails(complete_bundle):
    (complete_bundle / "db.yml").write_text("- just\n- a list\n")
    report = vdb.validate(complete_bundle, [])
    assert not report.ok
    assert any("must contain a mapping" in p for p in report.problems)


def test_unquoted_date_version_warns(complete_bundle):
    """PyYAML and snakeyaml both turn 2025-07-07 into a date; quoting keeps it a string."""
    (complete_bundle / "db.yml").write_text("name: gsgd\nversion: 2025-07-07\n")
    report = vdb.validate(complete_bundle, [])
    assert report.ok, report.problems
    assert any("quote it" in w for w in report.warnings)


def test_unknown_db_yml_key_warns(complete_bundle):
    (complete_bundle / "db.yml").write_text("name: gsgd\nnextclade_dir: nope\n")
    report = vdb.validate(complete_bundle, [])
    assert report.ok, report.problems
    assert any("unrecognised key" in w for w in report.warnings)


def test_snpeff_genomes_not_in_config_warns(complete_bundle):
    (complete_bundle / "db.yml").write_text(
        "name: gsgd\nsnpeff_genomes: [influenza_h5n1, influenza_h7n9]\n"
    )
    report = vdb.validate(complete_bundle, [])
    assert report.ok, report.problems
    assert any("influenza_h7n9" in w for w in report.warnings)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def test_cli_exit_zero_on_valid_bundle(complete_bundle):
    result = run_cli(complete_bundle)
    assert result.returncode == 0, result.stdout
    assert "VALID" in result.stdout


def test_cli_exit_one_on_invalid_bundle(complete_bundle):
    (complete_bundle / "snpeff" / "data" / "influenza_h5n1" / "snpEffectPredictor.bin").unlink()
    result = run_cli(complete_bundle)
    assert result.returncode == 1
    assert "INVALID" in result.stdout


def test_cli_json_output(complete_bundle):
    import json
    result = run_cli(complete_bundle, "--json")
    assert result.returncode == 0
    payload = json.loads(result.stdout)
    assert payload["valid"] is True
    assert payload["name"] == "gsgd"
    assert payload["version"] == "2025-07-07"
    assert set(payload["components"]) == {"blastn", "nextclade", "snpeff"}
