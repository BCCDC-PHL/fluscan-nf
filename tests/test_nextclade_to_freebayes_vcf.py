"""Unit tests for bin/nextclade_to_freebayes_vcf.py: the two-output split and shared codons.

Run with:  pytest tests/test_nextclade_to_freebayes_vcf.py

The script writes one VCF per sample, the file the rest of the workflow consumes, and, when
--no-mature-suffix is given, a second VCF annotated without mature peptides -- HA as HA0
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


def _write_json(path: Path, positions=SUB_POSITIONS, sample: str = "SAMPLE01",
                qry: dict[int, str] | None = None, nuc_to_aa: dict | None = None) -> Path:
    """A minimal nextclade v3 result carrying one substitution per position (qry sets the base)."""
    qry = qry or {}
    substitutions = [
        {"pos": pos, "refNuc": REF_SEQ[pos],
         "qryNuc": qry.get(pos, "A" if REF_SEQ[pos] != "A" else "C")}
        for pos in positions
    ]
    path.write_text(json.dumps({
        "results": [{
            "seqName":       sample,
            "substitutions": substitutions,
            "deletions":     [],
            "insertions":    [],
            "nucToAaMuts":   nuc_to_aa or {},
            "aaDeletions":   [],
            "aaInsertions":  [],
            "frameShifts":   [],
        }],
    }))
    return path


def _run(dataset: Path, json_path: Path, workdir: Path, no_mature: bool = False) -> tuple[Path, Path]:
    """Convert json_path inside workdir; return SAMPLE01's standard and no-mature VCF paths."""
    workdir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(SCRIPT_PATH), "-i", str(json_path),
           "-d", str(dataset), "--suffix", ".vcf"]
    if no_mature:
        cmd += ["--no-mature-suffix", "_nomat.vcf"]
    subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=workdir)
    return workdir / "SAMPLE01.vcf", workdir / "SAMPLE01_nomat.vcf"


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
    return _run(dataset, nextclade_json, tmp_path / "out", no_mature=True)


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


def test_gq_is_integer(both_outputs):
    for vcf in both_outputs:
        assert "##FORMAT=<ID=GQ,Number=1,Type=Integer," in vcf.read_text()
        for record in _records(vcf):
            fields = record.split("\t")
            assert fields[9].split(":")[fields[8].split(":").index("GQ")].isdigit()


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

    _, from_flag = _run(with_mature, nextclade_json, tmp_path / "flag", no_mature=True)
    from_gff, _ = _run(without, nextclade_json, tmp_path / "gff")

    assert _records(from_flag) == _records(from_gff)


# ---------------------------------------------------------------------------
# Substitutions sharing a codon
# ---------------------------------------------------------------------------

# 0-based 107-109 is a TAC (Tyr) codon: CDS codon 30, HA1 codon 29. T>G alone gives GAC (Asp),
# A>C alone gives TCC (Ser), together GCC (Ala).
MNV_POSITIONS = [107, 108, 400]
MNV_QRY = {107: "G", 108: "C"}


def _info(record: str) -> dict[str, str]:
    return dict(token.partition("=")[::2] for token in record.split("\t")[7].split(";"))


def _hgvs_p(record: str) -> dict[str, str]:
    return {fields[6]: fields[10] for fields in _ann_entries(record) if fields[10]}


@pytest.fixture
def mnv_records(tmp_path, dataset):
    json_path = _write_json(tmp_path / "mnv.json", positions=MNV_POSITIONS, qry=MNV_QRY)
    out, _ = _run(dataset, json_path, tmp_path / "mnv")
    return _records(out)


def test_shared_codon_is_annotated_with_the_combined_change(mnv_records):
    for record in mnv_records[:2]:
        assert _hgvs_p(record) == {"TESTLOCUS": "p.Y30A", "HA1": "p.Y29A"}
        assert "missense_variant" in record


def test_shared_codon_records_the_single_base_calls(mnv_records):
    first, second, _ = mnv_records
    assert "MNV_CORRECTED" in _info(first)
    assert _info(first)["MNV_ORIGINAL_HGVS_P"] == "TESTLOCUS:p.Y30D|HA1:p.Y29D"
    assert _info(second)["MNV_ORIGINAL_HGVS_P"] == "TESTLOCUS:p.Y30S|HA1:p.Y29S"


def test_lone_substitution_carries_no_mnv_fields(mnv_records):
    info = _info(mnv_records[2])
    assert "MNV_CORRECTED" not in info and "MNV_ORIGINAL_HGVS_P" not in info


def _calls(alt_aa: str) -> dict:
    """nucToAaMuts giving Nextclade's TESTHA call alt_aa for the shared codon at 107/108."""
    call = [{"cdsName": "TESTHA", "refAa": "Y", "qryAa": alt_aa, "pos": 29}]
    return {"107": call, "108": call}


def test_shared_codon_overrides_a_disagreeing_nextclade_call(tmp_path, dataset):
    """Nextclade's peptide alignment can gap a codon its nucleotide alignment keeps (HA K345)."""
    lone = {"400": [{"cdsName": "TESTHA", "refAa": "V", "qryAa": "-", "pos": 126}]}
    json_path = _write_json(tmp_path / "mnv.json", positions=MNV_POSITIONS, qry=MNV_QRY,
                            nuc_to_aa=_calls("-") | lone)
    first, second, third = _records(_run(dataset, json_path, tmp_path / "out")[0])
    for record in (first, second):
        assert _hgvs_p(record) == {"TESTLOCUS": "p.Y30A", "HA1": "p.Y29A"}
        assert "MNV_CORRECTED" in _info(record)
    assert _info(first)["MNV_ORIGINAL_HGVS_P"] == "TESTLOCUS:p.Y30-|HA1:p.Y29D"
    # A lone substitution keeps Nextclade's call.
    assert _hgvs_p(third)["TESTLOCUS"] == "p.V127-"


def test_shared_codon_matching_nextclade_is_not_flagged(tmp_path):
    """MNV_CORRECTED compares against Nextclade's call, not the single-base change."""
    dataset = _write_dataset(tmp_path / "dataset", GFF3_NO_MATURE)
    json_path = _write_json(tmp_path / "mnv.json", positions=MNV_POSITIONS[:2], qry=MNV_QRY,
                            nuc_to_aa=_calls("A"))
    for record in _records(_run(dataset, json_path, tmp_path / "out")[0]):
        assert _hgvs_p(record) == {"TESTLOCUS": "p.Y30A"}
        assert "MNV_CORRECTED" not in _info(record) and "MNV_ORIGINAL_HGVS_P" not in _info(record)


def test_mnv_fields_are_declared(tmp_path, dataset, nextclade_json):
    out, _ = _run(dataset, nextclade_json, tmp_path / "out")
    header = out.read_text()
    assert "##INFO=<ID=MNV_CORRECTED,Number=0,Type=Flag," in header
    assert "##INFO=<ID=MNV_ORIGINAL_HGVS_P,Number=1,Type=String," in header


# ---------------------------------------------------------------------------
# Batch input
# ---------------------------------------------------------------------------

def test_standard_output_is_unchanged_by_the_no_mature_suffix(tmp_path, dataset, nextclade_json):
    with_flag, _ = _run(dataset, nextclade_json, tmp_path / "with", no_mature=True)
    without_flag, _ = _run(dataset, nextclade_json, tmp_path / "without")
    assert with_flag.read_text() == without_flag.read_text()


def test_no_second_file_without_the_suffix(tmp_path, dataset, nextclade_json):
    standard, _ = _run(dataset, nextclade_json, tmp_path / "out")
    assert [p.name for p in standard.parent.glob("*.vcf")] == ["SAMPLE01.vcf"]


def _result(name: str) -> dict:
    return {"seqName": name, "substitutions": [{"pos": 100, "refNuc": REF_SEQ[100], "qryNuc": "A"}],
            "deletions": [], "insertions": [], "nucToAaMuts": {},
            "aaDeletions": [], "aaInsertions": [], "frameShifts": []}


def test_batch_input_writes_one_named_vcf_per_sample(tmp_path, dataset):
    json_path = tmp_path / "batch.json"
    json_path.write_text(json.dumps({"results": [_result("SAMPLE01"), _result("SAMPLE02")]}))
    subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "-i", str(json_path), "-d", str(dataset),
         "--suffix", "__REF__HA_nc_final.vcf", "--no-mature-suffix", "__REF__HA_nc_nomat_final.vcf"],
        check=True, capture_output=True, text=True, cwd=tmp_path,
    )
    for sample in ("SAMPLE01", "SAMPLE02"):
        standard = tmp_path / f"{sample}__REF__HA_nc_final.vcf"
        assert len(_records(standard)) == 1
        assert f"FORMAT\t{sample}" in standard.read_text()
        assert not _feature_ids(tmp_path / f"{sample}__REF__HA_nc_nomat_final.vcf") & {"HA1", "HA2"}


def test_unalignable_sample_gets_no_vcf(tmp_path, dataset):
    json_path = tmp_path / "batch.json"
    json_path.write_text(json.dumps({
        "results": [_result("SAMPLE02")],
        "errors":  [{"seqName": "SAMPLE01",
                     "errors": ["Unable to align: seed alignment covers 3.0% of the query"]}],
    }))
    _run(dataset, json_path, tmp_path / "out", no_mature=True)
    assert sorted(p.name for p in (tmp_path / "out").glob("*.vcf")) == ["SAMPLE02.vcf", "SAMPLE02_nomat.vcf"]


def test_batch_where_nothing_aligned_writes_nothing(tmp_path, dataset):
    json_path = tmp_path / "batch.json"
    json_path.write_text(json.dumps({"results": [], "errors": [{"seqName": "SAMPLE01", "errors": ["x"]}]}))
    _run(dataset, json_path, tmp_path / "out")
    assert list((tmp_path / "out").glob("*.vcf")) == []
