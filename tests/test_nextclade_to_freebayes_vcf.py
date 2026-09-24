"""Unit tests for bin/nextclade_to_freebayes_vcf.py, covering the two-output split.

Run with:  pytest tests/test_nextclade_to_freebayes_vcf.py

The script writes the VCF the rest of the workflow consumes and, when
--output-no-mature is given, a second VCF annotated without mature peptides -- HA as HA0
only, no HA1/HA2 entries. The contract the downstream consumer relies on is that the two
files differ *only* in their ANN sub-field: same records, same order, same alleles. That is
what is asserted here, alongside the guarantee that the standard output is byte-identical
whether or not the second one was requested.

The dataset is synthesised rather than taken from a reference bundle so the suite needs
nothing outside the repo. An HA-shaped annotation is used: one gene, one CDS, and two
mature peptides standing in for HA1 and HA2.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT   = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "bin" / "nextclade_to_freebayes_vcf.py"

pytest.importorskip("Bio", reason="Biopython is required by the script under test")

CHROM = "TEST_REF"

# CDS spans 0-based [20, 620); HA1 stands in at [23, 320), HA2 at [320, 617).
CDS_BEGIN, CDS_END = 20, 620
MAT1_BEGIN, MAT1_END = 23, 320
MAT2_BEGIN, MAT2_END = 320, 617

GFF3_WITH_MATURE = f"""\
##gff-version 3
{CHROM}\tTest\tgene\t{CDS_BEGIN + 1}\t{CDS_END}\t.\t+\t.\tID=gene-TESTLOCUS;Name=HA;gbkey=Gene;gene=HA;locus_tag=TESTLOCUS
{CHROM}\tTest\tCDS\t{CDS_BEGIN + 1}\t{CDS_END}\t.\t+\t0\tID=cds-TESTHA;Parent=gene-TESTLOCUS;gene=HA;locus_tag=TESTLOCUS;protein_id=TESTHA
{CHROM}\tTest\tmature_protein_region_of_CDS\t{MAT1_BEGIN + 1}\t{MAT1_END}\t.\t+\t.\tID=matpep-1;Parent=cds-TESTHA;gbkey=Prot;product=HA1
{CHROM}\tTest\tmature_protein_region_of_CDS\t{MAT2_BEGIN + 1}\t{MAT2_END}\t.\t+\t.\tID=matpep-2;Parent=cds-TESTHA;gbkey=Prot;product=HA2
"""

# Same annotation with the mature peptides left out: the no-mature output must reproduce it.
GFF3_NO_MATURE = "\n".join(GFF3_WITH_MATURE.splitlines()[:3]) + "\n"

# Positions chosen to land in each region of interest: upstream of the CDS, inside HA1,
# and inside HA2. All are 0-based, as nextclade reports them.
SUB_POSITIONS = [5, 100, 400]


def _reference_sequence(length: int = 700) -> str:
    """A reference with an in-frame ATG at the CDS start and no internal stop codons."""
    seq = list("ACGT" * (length // 4 + 1))[:length]
    seq[CDS_BEGIN:CDS_BEGIN + 3] = list("ATG")
    return "".join(seq)


REF_SEQ = _reference_sequence()


def _write_dataset(directory: Path, gff3: str) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "reference.fasta").write_text(f">{CHROM}\n{REF_SEQ}\n")
    (directory / "genome_annotation.gff3").write_text(gff3)
    return directory


def _write_json(path: Path, positions=SUB_POSITIONS, sample: str = "SAMPLE01") -> Path:
    """A minimal nextclade v3 result carrying one substitution per position."""
    substitutions = [
        {"pos": pos, "refNuc": REF_SEQ[pos], "qryNuc": "A" if REF_SEQ[pos] != "A" else "C"}
        for pos in positions
    ]
    path.write_text(json.dumps({
        "results": [{
            "seqName":       sample,
            "substitutions": substitutions,
            "deletions":     [],
            "insertions":    [],
            "nucToAaMuts":   {},
            "aaDeletions":   [],
            "aaInsertions":  [],
            "frameShifts":   [],
        }],
    }))
    return path


def _run(dataset: Path, json_path: Path, output: Path, no_mature: Path | None = None) -> None:
    cmd = [sys.executable, str(SCRIPT_PATH), "-i", str(json_path),
           "-d", str(dataset), "-o", str(output), "--sample", "SAMPLE01"]
    if no_mature is not None:
        cmd += ["--output-no-mature", str(no_mature)]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def _records(vcf: Path) -> list[str]:
    return [ln for ln in vcf.read_text().splitlines() if not ln.startswith("#")]


def _ann_entries(record: str) -> list[list[str]]:
    info = record.split("\t")[7]
    for token in info.split(";"):
        if token.startswith("ANN="):
            return [entry.split("|") for entry in token[4:].split(",") if entry]
    return []


def _feature_ids(vcf: Path) -> set[str]:
    return {fields[6] for record in _records(vcf) for fields in _ann_entries(record)}


@pytest.fixture
def dataset(tmp_path):
    return _write_dataset(tmp_path / "dataset", GFF3_WITH_MATURE)


@pytest.fixture
def nextclade_json(tmp_path):
    return _write_json(tmp_path / "nextclade.json")


@pytest.fixture
def both_outputs(tmp_path, dataset, nextclade_json):
    standard  = tmp_path / "out.vcf"
    no_mature = tmp_path / "out_nomat.vcf"
    _run(dataset, nextclade_json, standard, no_mature)
    return standard, no_mature


# ---------------------------------------------------------------------------
# The split itself
# ---------------------------------------------------------------------------

def test_standard_output_keeps_mature_peptide_entries(both_outputs):
    standard, _ = both_outputs
    assert {"HA1", "HA2"} <= _feature_ids(standard)


def test_no_mature_output_drops_mature_peptide_entries(both_outputs):
    _, no_mature = both_outputs
    feature_ids = _feature_ids(no_mature)
    assert not ({"HA1", "HA2"} & feature_ids)
    # The parent CDS must survive: an HA0-only VCF, not an unannotated one.
    assert "TESTLOCUS" in feature_ids


def test_upstream_mature_entries_are_dropped_too(both_outputs):
    """Position 5 precedes the CDS, so HA1/HA2 appear there as upstream_gene_variant."""
    standard, no_mature = both_outputs
    upstream_std = [fields for fields in _ann_entries(_records(standard)[0])
                    if fields[1] == "upstream_gene_variant"]
    assert {"HA1", "HA2"} <= {fields[6] for fields in upstream_std}
    assert not ({"HA1", "HA2"} & {fields[6] for fields in _ann_entries(_records(no_mature)[0])})


def test_only_the_ann_field_differs(both_outputs):
    """POS/REF/ALT/QUAL/FILTER and the record set are the downstream contract."""
    standard, no_mature = both_outputs
    std_records, nomat_records = _records(standard), _records(no_mature)
    assert len(std_records) == len(nomat_records) == len(SUB_POSITIONS)
    for std, nomat in zip(std_records, nomat_records):
        assert std.split("\t")[:7] == nomat.split("\t")[:7]


def test_no_mature_entries_are_a_subset_of_the_standard_ones(both_outputs):
    """Nothing is rewritten on the way out; entries are only ever removed."""
    standard, no_mature = both_outputs
    for std, nomat in zip(_records(standard), _records(no_mature)):
        std_entries = ["|".join(fields) for fields in _ann_entries(std)]
        for entry in ("|".join(fields) for fields in _ann_entries(nomat)):
            assert entry in std_entries


def test_no_mature_output_matches_a_dataset_without_mature_peptides(tmp_path, nextclade_json):
    """The second output is the annotation a GFF3 with no mature peptides would produce."""
    with_mature = _write_dataset(tmp_path / "with", GFF3_WITH_MATURE)
    without     = _write_dataset(tmp_path / "without", GFF3_NO_MATURE)

    from_flag = tmp_path / "from_flag.vcf"
    _run(with_mature, nextclade_json, tmp_path / "ignored.vcf", from_flag)

    from_gff = tmp_path / "from_gff.vcf"
    _run(without, nextclade_json, from_gff)

    assert _records(from_flag) == _records(from_gff)


# ---------------------------------------------------------------------------
# The existing single-output behaviour
# ---------------------------------------------------------------------------

def test_standard_output_is_unchanged_by_the_new_argument(tmp_path, dataset, nextclade_json):
    with_flag = tmp_path / "with_flag.vcf"
    _run(dataset, nextclade_json, with_flag, tmp_path / "nomat.vcf")

    without_flag = tmp_path / "without_flag.vcf"
    _run(dataset, nextclade_json, without_flag)

    assert with_flag.read_text() == without_flag.read_text()


def test_no_second_file_without_the_argument(tmp_path, dataset, nextclade_json):
    standard = tmp_path / "out.vcf"
    _run(dataset, nextclade_json, standard)
    assert standard.exists()
    assert [p.name for p in tmp_path.glob("*.vcf")] == ["out.vcf"]


def test_multi_sample_input_splits_both_outputs(tmp_path, dataset):
    json_path = tmp_path / "multi.json"
    json_path.write_text(json.dumps({
        "results": [
            {"seqName": name, "substitutions": [
                {"pos": 100, "refNuc": REF_SEQ[100], "qryNuc": "A"}],
             "deletions": [], "insertions": [], "nucToAaMuts": {},
             "aaDeletions": [], "aaInsertions": [], "frameShifts": []}
            for name in ("SAMPLE01", "SAMPLE02")
        ],
    }))
    subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "-i", str(json_path), "-d", str(dataset),
         "-o", str(tmp_path / "out.vcf"),
         "--output-no-mature", str(tmp_path / "out_nomat.vcf")],
        check=True, capture_output=True, text=True,
    )
    for sample in ("SAMPLE01", "SAMPLE02"):
        assert (tmp_path / f"out_{sample}.vcf").exists()
        assert (tmp_path / f"out_nomat_{sample}.vcf").exists()
        assert not _feature_ids(tmp_path / f"out_nomat_{sample}.vcf") & {"HA1", "HA2"}


def test_unalignable_sample_writes_both_empty_vcfs(tmp_path, dataset):
    json_path = tmp_path / "errors.json"
    json_path.write_text(json.dumps({
        "results": [],
        "errors":  [{"seqName": "SAMPLE01",
                     "errors": ["Unable to align: seed alignment covers 3.0% of the query"]}],
    }))
    standard  = tmp_path / "out.vcf"
    no_mature = tmp_path / "out_nomat.vcf"
    _run(dataset, json_path, standard, no_mature)
    assert _records(standard) == []
    assert _records(no_mature) == []
