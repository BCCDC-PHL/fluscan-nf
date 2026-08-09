"""Unit tests for bin/build_reference_coordinate_map.py.

Three things are pinned here:

1. **HA subunits are read from either annotation shape.** FluMutDB numbers HA by its
   HA1/HA2 subunits. A target dataset may express that split as two real CDS features
   (chicken2021) or as `mature_protein_region_of_CDS` children of one HA CDS (the
   NCBI/RefSeq model, the GsGD dataset among them). The two must produce the same
   coordinate map -- the second one produced no HA map at all until the mature rows were
   read, putting all 76 of FluMutDB's HA mappings out of reach on the default dataset.

2. **QC reporting does not go silent on the worst input.** `safe_pct` yields NaN when
   nothing aligned, and every `nan < threshold` comparison is False, so a plain comparison
   stays quiet on exactly the case that most needs reporting. A questionable map is
   reported, never withheld.

3. **The published artifacts are an interface.** Filenames, JSON keys and TSV columns are
   consumed downstream and by hand; they are asserted here so a refactor cannot quietly
   drop one.

Only Biopython is required; the datasets are synthesised in tmp_path.
"""

import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "bin" / "build_reference_coordinate_map.py"

pytest.importorskip("Bio", reason="build_reference_coordinate_map requires Biopython")


def _load_module():
    name = "build_reference_coordinate_map"
    spec = importlib.util.spec_from_file_location(name, SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    # Registered before execution because @dataclass resolves annotations through
    # sys.modules[cls.__module__]; without this the class body raises on import.
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


builder = _load_module()

# HA is modelled on chicken2021: HA1 at 77..1063 (987 nt) and HA2 at 1064..1732 (669 nt),
# adjacent, inside a 1776 nt segment. Because they abut, one CDS spanning 77..1732 with two
# mature children over the same intervals is an exact re-encoding of the same protein
# structure -- which is what makes the equivalence assertion below meaningful.
HA_LENGTH = 1776
HA1_RANGE = (77, 1063)
HA2_RANGE = (1064, 1732)

# Every sense codon; no stops, so a translated model never truncates unexpectedly.
CODONS = [a + b + c
          for a in "ACGT" for b in "ACGT" for c in "ACGT"
          if a + b + c not in ("TAA", "TAG", "TGA")]


def synthetic_sequence(length: int, seed: int) -> str:
    """A deterministic, stop-free nucleotide sequence of the requested length."""
    out = []
    value = seed
    while len(out) * 3 < length:
        value = (value * 1103515245 + 12345) % 2147483648
        out.append(CODONS[value % len(CODONS)])
    return "".join(out)[:length]


def gff_line(chrom: str, feature: str, start: int, end: int, attrs: str) -> str:
    return f"{chrom}\t.\t{feature}\t{start}\t{end}\t.\t+\t0\t{attrs}"


def write_segment(seg_dir: Path, chrom: str, sequence: str, gff_lines: list[str]) -> None:
    seg_dir.mkdir(parents=True, exist_ok=True)
    (seg_dir / "reference.fasta").write_text(f">{chrom}\n{sequence}\n")
    (seg_dir / "genome_annotation.gff3").write_text("\n".join(["##gff-version 3", *gff_lines]) + "\n")


def ha_gff_split_cds(chrom: str) -> list[str]:
    """HA1/HA2 as two genuine CDS features (the chicken2021 / flumut shape)."""
    lines = []
    for gene, (start, end) in (("HA1", HA1_RANGE), ("HA2", HA2_RANGE)):
        lines.append(gff_line(chrom, "gene", start, end,
                              f"ID=gene-{gene};Name={gene};gene={gene}"))
        lines.append(gff_line(chrom, "CDS", start, end,
                              f"ID=cds-{gene};Parent=gene-{gene};Name={gene};gene={gene}"))
    return lines


def ha_gff_mature_peptides(chrom: str) -> list[str]:
    """One HA CDS with HA1/HA2 as mature peptide children (the NCBI/RefSeq / GsGD shape)."""
    start, end = HA1_RANGE[0], HA2_RANGE[1]
    lines = [
        gff_line(chrom, "gene", start, end, "ID=gene-HA;Name=HA;gene=HA"),
        gff_line(chrom, "CDS", start, end,
                 "ID=cds-HA;Parent=gene-HA;Name=HA;gene=HA;protein_id=SYN_HA"),
    ]
    for product, (m_start, m_end) in (("HA1", HA1_RANGE), ("HA2", HA2_RANGE)):
        lines.append(gff_line(chrom, "mature_protein_region_of_CDS", m_start, m_end,
                              f"ID=id-{product};Parent=cds-HA;product={product};"
                              f"protein_id=SYN_{product}"))
    return lines


def write_dataset(root: Path, prefix: str, ha_gff, ha_sequence: str) -> Path:
    """A minimal 8-segment reference dataset. Only HA carries a meaningful annotation; the
    other seven exist because build_map requires the full canonical segment set."""
    root.mkdir(parents=True, exist_ok=True)
    for segment in builder.SEGMENTS:
        chrom = f"{prefix}_{segment}"
        seg_dir = root / segment
        if segment == "HA":
            write_segment(seg_dir, chrom, ha_sequence, ha_gff(chrom))
            continue
        sequence = synthetic_sequence(300, seed=hash(segment) % 9973 + 1)
        write_segment(seg_dir, chrom, sequence, [
            gff_line(chrom, "gene", 1, 300, f"ID=gene-{segment};Name={segment};gene={segment}"),
            gff_line(chrom, "CDS", 1, 300, f"ID=cds-{segment};Parent=gene-{segment};gene={segment}"),
        ])
    return root


def build(tmp_path: Path, source: Path, target: Path, name: str, **kwargs) -> dict:
    out_dir = tmp_path / name
    args = builder.parse_args([
        "--source-ref-dir", str(source),
        "--target-ref-dir", str(target),
        "--out-dir", str(out_dir),
        "--source-name", "src",
        "--target-name", "tgt",
        *[a for a in kwargs.get("extra", [])],
    ])
    builder.build_map(args)
    return json.loads((out_dir / "coordinate_map.json").read_text())


@pytest.fixture
def datasets(tmp_path):
    """A source dataset plus two targets that differ only in how HA is annotated."""
    source_ha = synthetic_sequence(HA_LENGTH, seed=11)
    # A different strain: same structure, different residues, so the alignment is real work
    # rather than a self-map.
    target_ha = synthetic_sequence(HA_LENGTH, seed=11)[:900] + synthetic_sequence(HA_LENGTH, seed=29)[900:]

    source = write_dataset(tmp_path / "source", "SRC", ha_gff_split_cds, source_ha)
    split = write_dataset(tmp_path / "target_split", "TGT", ha_gff_split_cds, target_ha)
    mature = write_dataset(tmp_path / "target_mature", "TGT", ha_gff_mature_peptides, target_ha)
    return source, split, mature


# ---------------------------------------------------------------------------
# HA subunit recovery from mature peptides
# ---------------------------------------------------------------------------

def test_mature_peptides_yield_the_same_ha_map_as_split_cds(tmp_path, datasets):
    source, split, mature = datasets

    split_ha = build(tmp_path, source, split, "split")["segments"]["HA"]
    mature_ha = build(tmp_path, source, mature, "mature")["segments"]["HA"]

    assert sorted(split_ha["genes"]) == ["HA1", "HA2"]
    assert sorted(mature_ha["genes"]) == ["HA1", "HA2"]
    for gene in ("HA1", "HA2"):
        assert (mature_ha["genes"][gene]["aa_source_to_target"]
                == split_ha["genes"][gene]["aa_source_to_target"])
        assert (mature_ha["genes"][gene]["source_length"]
                == split_ha["genes"][gene]["source_length"])


def test_mature_peptide_target_is_not_flagged_stale(tmp_path, datasets):
    """The dataset is not out of date -- the NCBI mature-peptide model is the RefSeq
    standard. Only an HA with no subunit features of any kind is stale."""
    source, _, mature = datasets
    ha = build(tmp_path, source, mature, "mature")["segments"]["HA"]

    assert ha["stale_annotation"] is False


def test_mature_peptide_aliases_are_recorded_for_the_merge(tmp_path, datasets):
    """The merge keys on these: a mature peptide's ANN entry names the parent gene in
    Gene_Name, so only Feature_ID identifies the subunit."""
    source, _, mature = datasets
    ha = build(tmp_path, source, mature, "mature")["segments"]["HA"]

    assert "SYN_HA1" in ha["genes"]["HA1"]["aliases"]
    assert "SYN_HA2" in ha["genes"]["HA2"]["aliases"]


def test_subunit_parent_is_named_rather_than_reported_missing(tmp_path, datasets):
    """HA has no source counterpart once its subunits carry the map. Listing it as missing
    would warn, every run, about the arrangement the build just handled correctly."""
    source, _, mature = datasets
    ha = build(tmp_path, source, mature, "mature")["segments"]["HA"]

    assert ha["subunit_parent_genes"] == ["HA"]
    assert "HA" not in ha["missing_target_genes"]


def test_split_cds_target_records_no_subunit_metadata(tmp_path, datasets):
    """Datasets whose subunits are real CDS features are untouched by any of this."""
    source, split, _ = datasets
    ha = build(tmp_path, source, split, "split")["segments"]["HA"]

    assert ha["subunit_parent_genes"] == []
    assert all("aliases" not in gene_doc for gene_doc in ha["genes"].values())


def test_a_real_cds_wins_over_a_mature_peptide_of_the_same_name(tmp_path):
    """If both are present, the CDS defines the key -- a mature row must never displace it."""
    chrom = "X_HA"
    gff = ha_gff_split_cds(chrom) + [
        gff_line(chrom, "mature_protein_region_of_CDS", 100, 400,
                 "ID=id-HA1;Parent=cds-HA1;product=HA1;protein_id=DECOY"),
    ]
    path = tmp_path / "ha.gff3"
    path.write_text("\n".join(gff) + "\n")

    models = builder.read_cds_models(path, "HA")

    assert models["HA1"].parts == (HA1_RANGE,)
    assert models["HA1"].is_mature is False
    assert models["HA1"].aliases == ()


# ---------------------------------------------------------------------------
# Mature-peptide reading frame
# ---------------------------------------------------------------------------
#
# A mature peptide is translated from its own first base, so treating it as a CDS is only
# valid when it begins on a codon boundary of its parent. A peptide starting mid-codon still
# translates, still aligns, and still yields a position map -- of entirely different
# residues. Nothing errors, so the check is the only thing standing between a malformed
# annotation and confidently wrong markers.

def mature_ha_gff(chrom: str, ha1_start: int, parent: str = "cds-HA") -> list[str]:
    return [
        gff_line(chrom, "gene", 22, 1728, "ID=gene-HA;Name=HA;gene=HA"),
        gff_line(chrom, "CDS", 22, 1728, "ID=cds-HA;Parent=gene-HA;gene=HA;protein_id=P_HA"),
        gff_line(chrom, "mature_protein_region_of_CDS", ha1_start, 1059,
                 f"ID=id-HA1;Parent={parent};product=HA1;protein_id=P_HA1"),
    ]


@pytest.mark.parametrize("ha1_start,in_frame", [
    (70, True),    # 70 - 22 = 48 = exactly 16 codons, the real RefSeq arrangement
    (71, False),   # one base into a codon
    (72, False),   # two bases in
    (73, True),    # 51 = 17 codons
])
def test_mature_peptide_must_be_codon_aligned_to_its_parent(tmp_path, capsys, ha1_start, in_frame):
    path = tmp_path / f"ha_{ha1_start}.gff3"
    path.write_text("\n".join(mature_ha_gff("X_HA", ha1_start)) + "\n")

    models = builder.read_cds_models(path, "HA")
    err = capsys.readouterr().err

    assert ("HA1" in models) is in_frame
    if in_frame:
        assert models["HA1"].is_mature is True
        assert "mature peptide HA1" not in err
    else:
        assert "into a codon of its parent" in err and "skipping it" in err


def test_mature_peptide_with_an_unresolvable_parent_is_refused(tmp_path, capsys):
    """Without the parent CDS there is no way to know where its codons begin, so the frame
    cannot be checked at all. Accepting the one case we cannot verify would defeat the check
    -- the peptide is skipped, and HA then falls to the stale-annotation path, which reports
    the loss loudly rather than emitting markers of unknown validity."""
    path = tmp_path / "ha.gff3"
    path.write_text("\n".join(mature_ha_gff("X_HA", 70, parent="cds-TYPO")) + "\n")

    models = builder.read_cds_models(path, "HA")
    err = capsys.readouterr().err

    assert "HA1" not in models
    assert "no resolvable parent CDS" in err and "cannot be verified" in err
    # HA1/HA2 absent => the existing stale check fires, so the operator is told.
    assert builder.detect_stale_ha_annotation("HA", {"HA1", "HA2"}, set(models)) is True


def test_frame_offset_follows_a_spliced_parent(tmp_path):
    """Offset is measured along the spliced CDS, not by genomic subtraction: for a parent
    with an intron, the two disagree. No influenza reference exercises this today, which is
    exactly why it is pinned here."""
    parent = builder.CdsModel(
        segment="M", gene="M2", contig="X", strand="+",
        parts=((1, 26), (715, 1000)),
    )
    # Spliced coding sequence is 1..26 then 715..1000, so 715 is the 27th coding base.
    assert builder.mature_peptide_frame_offset(parent, 715, 1000) == 26
    assert builder.mature_peptide_frame_offset(parent, 716, 1000) == 27
    # Genomic subtraction would have said 714 and 715 -- a different frame verdict.
    assert (715 - 1) % 3 != 26 % 3
    # A position inside the intron belongs to no codon of the parent.
    assert builder.mature_peptide_frame_offset(parent, 300, 1000) is None


def test_frame_offset_counts_from_the_end_on_the_minus_strand(tmp_path):
    parent = builder.CdsModel(
        segment="X", gene="G", contig="X", strand="-",
        parts=((100, 200),),
    )
    # Translation starts at 200 and runs down, so the peptide ending at 200 is at offset 0.
    assert builder.mature_peptide_frame_offset(parent, 150, 200) == 0
    assert builder.mature_peptide_frame_offset(parent, 150, 197) == 3


def test_v1_target_cds_payload_is_gone(tmp_path, datasets):
    """The block existed so v1 could place synthetic rows; v1 is retired and v2 never read
    it. Its removal is what drops the map builder's dependency on the ANN producer."""
    source, split, _ = datasets
    doc = build(tmp_path, source, split, "split")

    assert doc["schema_version"] == "2.2"
    for seg_doc in doc["segments"].values():
        for gene_doc in seg_doc["genes"].values():
            assert "target_cds" not in gene_doc


# ---------------------------------------------------------------------------
# QC reporting
# ---------------------------------------------------------------------------

def alignment(identities: int, aligned: int, source_length: int) -> "builder.AlignmentResult":
    """An AlignmentResult carrying only the counts the quality properties read."""
    return builder.AlignmentResult(
        source_seq="X" * source_length,
        target_seq="X" * source_length,
        source_aligned="",
        target_aligned="",
        source_to_target={},
        identities=identities,
        aligned=aligned,
        gaps=0,
        co_optimal=1,
    )


def test_nothing_aligned_is_reported_rather_than_passing_silently():
    """safe_pct returns NaN here, and `nan < 80.0` is False -- so a plain threshold test
    stays silent on a map with zero aligned positions, the worst possible outcome."""
    nothing_aligned = alignment(identities=0, aligned=0, source_length=100)

    assert math.isnan(builder.safe_pct(0, 0))
    assert math.isnan(nothing_aligned.identity_pct)
    assert nothing_aligned.is_questionable is True
    assert nothing_aligned.status == "warning"


@pytest.mark.parametrize("identities,aligned,source_length,expected", [
    (100, 100, 100, False),   # perfect identity, full coverage
    (79, 100, 100, True),     # identity below 80%
    (94, 94, 100, True),      # coverage below 95%
    (0, 0, 100, True),        # nothing aligned: both percentages are NaN
])
def test_a_questionable_map_is_flagged(identities, aligned, source_length, expected):
    result = alignment(identities, aligned, source_length)
    assert result.is_questionable is expected
    assert result.status == ("warning" if expected else "ok")


def test_qc_tsv_leaves_an_unmeasurable_identity_blank():
    """`nan` as a literal string is awkward downstream, and inconsistent with the
    missing_in_source / extra_in_source rows, which already leave the column empty."""
    assert builder.format_pct(math.nan) == ""
    assert builder.format_pct(87.234) == "87.2340"


def test_qc_tsv_carries_the_co_optimal_alignment_count(tmp_path, datasets):
    """Co-optimal alignments make gap placement arbitrary; the count is what lets a reviewer
    see that the map they are reading was one choice among several."""
    source, split, _ = datasets
    out_dir = tmp_path / "qc"
    builder.build_map(builder.parse_args([
        "--source-ref-dir", str(source), "--target-ref-dir", str(split),
        "--out-dir", str(out_dir), "--source-name", "src", "--target-name", "tgt",
    ]))

    header, *rows = (out_dir / "coordinate_map.qc.tsv").read_text().splitlines()
    assert "co_optimal_alignments" in header.split("\t")
    index = header.split("\t").index("co_optimal_alignments")
    assert all(int(row.split("\t")[index]) >= 1 for row in rows if row.split("\t")[2] == "nt")


def test_a_questionable_map_warns_and_is_written_anyway(tmp_path, datasets, capsys):
    """Real reference pairs sit below these thresholds -- flumut -> gsgd reports NS/NS1 at
    66% identity -- so refusing to write one would reject a dataset the pipeline documents
    as supported. The warning and the QC row are the whole report."""
    source, split, _ = datasets
    out_dir = tmp_path / "questionable"

    builder.build_map(builder.parse_args([
        "--source-ref-dir", str(source), "--target-ref-dir", str(split),
        "--out-dir", str(out_dir), "--source-name", "src", "--target-name", "tgt",
    ]))
    err = capsys.readouterr().err

    assert "WARNING" in err and "map has identity=" in err

    header, *rows = (out_dir / "coordinate_map.qc.tsv").read_text().splitlines()
    status = header.split("\t").index("status")
    assert any(row.split("\t")[status] == "warning" for row in rows)
    assert (out_dir / "coordinate_map.json").exists()


def test_per_position_status_follows_the_mapped_coordinate(tmp_path, datasets):
    """`status` is a reading of the coordinate, not a second record kept alongside it: a
    source position with no target counterpart is exactly the one written as unmapped."""
    assert builder.map_status(42) == "mapped"
    assert builder.map_status(None) == "unmapped"

    source, split, _ = datasets
    out_dir = tmp_path / "status"
    builder.build_map(builder.parse_args([
        "--source-ref-dir", str(source), "--target-ref-dir", str(split),
        "--out-dir", str(out_dir), "--source-name", "src", "--target-name", "tgt",
    ]))

    header, *rows = (out_dir / "coordinate_map.nt.tsv").read_text().splitlines()
    columns = header.split("\t")
    target_pos = columns.index("target_pos")
    status = columns.index("status")
    for row in rows:
        fields = row.split("\t")
        assert fields[status] == ("mapped" if fields[target_pos] else "unmapped")


# ---------------------------------------------------------------------------
# Published artifacts
# ---------------------------------------------------------------------------

TSV_COLUMNS = {
    "coordinate_map.nt.tsv": [
        "source_reference", "segment", "source_chrom", "source_pos", "source_base",
        "target_chrom", "target_pos", "target_base", "status",
    ],
    "coordinate_map.aa.tsv": [
        "source_reference", "segment", "gene", "source_aa_pos", "source_residue",
        "target_aa_pos", "target_residue", "status",
    ],
    "coordinate_map.qc.tsv": [
        "source_reference", "segment", "layer", "gene", "source_length", "target_length",
        "aligned_positions", "gap_columns", "identity_pct", "co_optimal_alignments",
        "status",
    ],
}

SEGMENT_KEYS = {
    "source_chrom", "target_chrom", "source_length", "target_length", "nt_identity_pct",
    "nt_aligned_positions", "nt_gap_columns", "nt_co_optimal_alignments",
    "nt_source_to_target", "missing_target_genes", "extra_source_genes",
    "stale_annotation", "subunit_parent_genes", "genes",
}

GENE_KEYS = {
    "source_length", "target_length", "aa_identity_pct", "aa_aligned_positions",
    "aa_gap_columns", "aa_co_optimal_alignments", "aa_source_to_target",
}


@pytest.fixture
def built(tmp_path, datasets):
    """One real build, for the assertions that only read its outputs."""
    source, split, _ = datasets
    out_dir = tmp_path / "artifacts"
    builder.build_map(builder.parse_args([
        "--source-ref-dir", str(source), "--target-ref-dir", str(split),
        "--out-dir", str(out_dir), "--source-name", "src", "--target-name", "tgt",
    ]))
    return out_dir


def test_every_documented_artifact_is_written(built):
    expected = {
        "coordinate_map.json",
        *(f"coordinate_map.{segment}.json" for segment in builder.SEGMENTS),
        "coordinate_map.nt.tsv",
        "coordinate_map.aa.tsv",
        "coordinate_map.qc.tsv",
        "src_alignment.fasta",
        "src_aa_alignment.fasta",
    }

    assert expected <= {path.name for path in built.iterdir()}


def test_tsv_columns_are_unchanged(built):
    for name, columns in TSV_COLUMNS.items():
        header = (built / name).read_text().splitlines()[0]
        assert header.split("\t") == columns, name


def test_json_fields_are_unchanged(built):
    doc = json.loads((built / "coordinate_map.json").read_text())

    assert doc["schema_version"] == "2.2"
    assert set(doc) == {
        "schema_version", "script_version", "source_reference", "target_reference",
        "alignment_parameters", "segments",
    }
    assert sorted(doc["segments"]) == sorted(builder.SEGMENTS)
    for seg_doc in doc["segments"].values():
        assert set(seg_doc) == SEGMENT_KEYS
        for gene_doc in seg_doc["genes"].values():
            # `aliases` is present only for genes recovered from mature peptides.
            assert set(gene_doc) - {"aliases"} == GENE_KEYS


def test_shards_drop_only_the_nucleotide_position_map(built):
    """Shards are what the merge actually reads. They keep the same segments[SEG] nesting
    so a consumer works unchanged, minus the large nt map nothing downstream uses."""
    doc = json.loads((built / "coordinate_map.json").read_text())

    for segment in builder.SEGMENTS:
        shard = json.loads((built / f"coordinate_map.{segment}.json").read_text())
        assert set(shard["segments"]) == {segment}
        seg_doc = shard["segments"][segment]
        assert set(seg_doc) == SEGMENT_KEYS - {"nt_source_to_target"}
        assert seg_doc == {
            k: v for k, v in doc["segments"][segment].items() if k != "nt_source_to_target"
        }


# ---------------------------------------------------------------------------
# Dataset discovery
# ---------------------------------------------------------------------------

def test_config_json_is_not_consulted(tmp_path, datasets):
    """config.json is scripts/build_dataset.py's input manifest, not a runtime interface.
    A bundle need not carry one, and a broken one must not reach the builder at all."""
    source, split, _ = datasets
    (split / "config.json").write_text("{ this is not json")

    doc = build(tmp_path, source, split, "no_config")

    assert sorted(doc["segments"]) == sorted(builder.SEGMENTS)
    assert doc["source_reference"] == "src"
    assert doc["target_reference"] == "tgt"
