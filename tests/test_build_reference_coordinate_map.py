"""Unit tests for bin/build_reference_coordinate_map.py.

Two things are pinned here:

1. **HA subunits are read from either annotation shape.** FluMutDB numbers HA by its
   HA1/HA2 subunits. A target dataset may express that split as two real CDS features
   (chicken2021) or as `mature_protein_region_of_CDS` children of one HA CDS (the
   NCBI/RefSeq model, the GsGD dataset among them). The two must produce the same
   coordinate map -- the second one produced no HA map at all until the mature rows were
   read, putting all 76 of FluMutDB's HA mappings out of reach on the default dataset.

2. **QC reporting does not go silent on the worst input.** `safe_pct` yields NaN when
   nothing aligned, and every `nan < threshold` comparison is False, so a plain comparison
   stays quiet on exactly the case that most needs reporting.

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
        seg_dir = root / segment.lower()
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
        parts=((1, 26), (715, 1000)), cds_ids=(),
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
        parts=((100, 200),), cds_ids=(),
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

def test_nothing_aligned_is_reported_rather_than_passing_silently():
    """safe_pct returns NaN here, and `nan < 80.0` is False -- so a plain threshold test
    stays silent on a map with zero aligned positions, the worst possible outcome."""
    assert math.isnan(builder.safe_pct(0, 0))
    assert builder.below_threshold(builder.safe_pct(0, 0), 100.0) is True
    assert builder.compute_status(builder.safe_pct(0, 0), 100.0) == "warning"


@pytest.mark.parametrize("identity,coverage,expected", [
    (100.0, 100.0, False),
    (79.9, 100.0, True),
    (100.0, 94.9, True),
    (math.nan, 100.0, True),
    (100.0, math.nan, True),
])
def test_below_threshold(identity, coverage, expected):
    assert builder.below_threshold(identity, coverage) is expected


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


def test_low_identity_fails_only_when_asked(tmp_path, datasets):
    """Real reference pairs sit below these thresholds (flumut -> gsgd reports NS1 at 66%),
    so failing by default would reject a dataset the pipeline documents as supported."""
    source, split, _ = datasets
    common = ["--source-ref-dir", str(source), "--target-ref-dir", str(split),
              "--source-name", "src", "--target-name", "tgt"]

    builder.build_map(builder.parse_args([*common, "--out-dir", str(tmp_path / "warn_only")]))

    with pytest.raises(SystemExit, match="Low-quality coordinate map"):
        builder.build_map(builder.parse_args(
            [*common, "--out-dir", str(tmp_path / "strict"), "--fail-on-low-identity"]))


# ---------------------------------------------------------------------------
# Reference naming
# ---------------------------------------------------------------------------

def test_explicit_name_beats_config_json(tmp_path):
    """The caller's name wins so the process's stub block can reproduce the real output
    filenames; a dataset carrying an id_prefix used to make the two diverge."""
    ref_dir = tmp_path / "chicken2021" / "nextclade"
    ref_dir.mkdir(parents=True)
    (ref_dir / "config.json").write_text(json.dumps({"id_prefix": "CHICKEN2021"}))

    assert builder.resolve_reference_name(ref_dir, "chicken2021") == "chicken2021"
    assert builder.resolve_reference_name(ref_dir, None) == "CHICKEN2021"


def test_parent_directory_is_the_last_resort(tmp_path):
    ref_dir = tmp_path / "gsgd" / "nextclade"
    ref_dir.mkdir(parents=True)

    assert builder.resolve_reference_name(ref_dir, None) == "gsgd"
