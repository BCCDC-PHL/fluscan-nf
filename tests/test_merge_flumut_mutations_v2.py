"""Unit tests for bin/merge_flumut_mutations_v2.py (FluMutDB VCF annotation).

Covers the behaviours the v2 merge is specified around: matching on
(protein, AA position, alt AA), translation through the coordinate map's
aa_source_to_target block, the shape of the FLUMUT INFO field, and the guarantee that no
row is added, removed, reordered or otherwise altered.

The real FluMutDB shipped by the flumutdb conda package is used rather than a fixture, so
these assertions double as a check that the queries still fit the released schema. The
concrete markers asserted on (FluMutDB 6.6):

    PB2:E627K   reference PB2, protein PB2, position 627, alt K
    M1:N30D     reference MP,  protein M1,  position 30,  alt D
    M2:A30V/T/S reference MP,  protein M2,  position 30   (no alt D -- see the M1/M2 test)
    HA1-5:T156A reference HA,  protein HA1, position 156, alt A

Run with:  pytest tests/test_merge_flumut_mutations_v2.py

The database is located the same way the script does, via the installed flumutdb package.
Set FLUMUT_DB to a flumut_db.sqlite path to run the suite from an environment that has
pytest but not flumut installed; the whole module skips when neither is available.
"""

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "bin" / "merge_flumut_mutations_v2.py"

VCF_COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
BASE_INFO = "NS=1;DP=1;AF=1.0;TYPE=snp;LEN=1"
FORMAT_COLS = "GT:DP\t1:1"


def _load_module():
    spec = importlib.util.spec_from_file_location("merge_flumut_mutations_v2", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


merge_v2 = _load_module()


def _find_db() -> str | None:
    if os.environ.get("FLUMUT_DB"):
        return os.environ["FLUMUT_DB"]
    try:
        return merge_v2.locate_db(None)
    except SystemExit:
        return None


FLUMUT_DB = _find_db()

pytestmark = pytest.mark.skipif(
    FLUMUT_DB is None,
    reason="flumutdb is not installed and FLUMUT_DB is not set",
)


# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------

def ann(gene: str, hgvs_p: str, annotation: str = "missense_variant", alt: str = "A") -> str:
    """A 16-field snpEff-style ANN entry, matching what the converter emits."""
    return "|".join([
        alt, annotation, "MODERATE", gene, f"LOCUS_{gene}", "transcript",
        f"LOCUS_{gene}", "protein_coding", "1/1", "c.1A>C", hgvs_p,
        "1/999", "1/999", "1/333", "", "",
    ])


def vcf_row(pos: int, ann_entries: list[str], ref: str = "A", alt: str = "C") -> str:
    info = BASE_INFO
    if ann_entries:
        info += ";ANN=" + ",".join(ann_entries)
    return f"CHROM_X\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t{info}\t{FORMAT_COLS}"


def write_vcf(path: Path, rows: list[str]) -> Path:
    header = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=CHROM_X,length=2400>",
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations">',
        VCF_COLUMNS,
    ]
    path.write_text("\n".join(header + rows) + "\n")
    return path


def write_coord_map(path: Path, segment: str, genes: dict, stale: bool = False) -> Path:
    """A coordinate-map shard. `genes` maps gene key -> {target_pos: source_pos} pairs
    expressed the way the real map stores them: aa_source_to_target, keyed by source."""
    doc = {
        "schema_version": "2.1",
        "source_reference": "FLUMUT",
        "target_reference": "TEST",
        "segments": {
            segment: {
                "stale_annotation": stale,
                "genes": {
                    gene: {"aa_source_to_target": {str(src): tgt for src, tgt in mapping.items()}}
                    for gene, mapping in genes.items()
                },
            }
        },
    }
    path.write_text(json.dumps(doc, indent=2))
    return path


def identity_map(n: int = 800) -> dict:
    return {i: i for i in range(1, n + 1)}


def run_script(tmp_path: Path, vcf: Path, coord_map: Path, segment: str,
               extra: list[str] | None = None) -> tuple[subprocess.CompletedProcess, Path]:
    out = tmp_path / "out.vcf"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--vcf", str(vcf), "--coord-map", str(coord_map),
         "--segment", segment, "--sample", "SAMPLE_1", "-o", str(out),
         "--flumut-db", FLUMUT_DB] + (extra or []),
        capture_output=True, text=True,
    )
    return proc, out


def flumut_fields(vcf_path: Path) -> list[str | None]:
    """The FLUMUT INFO value of each data row (None when the row has no FLUMUT key)."""
    values: list[str | None] = []
    for line in vcf_path.read_text().splitlines():
        if line.startswith("#"):
            continue
        info = line.split("\t")[7]
        hit = [t[len("FLUMUT="):] for t in info.split(";") if t.startswith("FLUMUT=")]
        values.append(hit[0] if hit else None)
    return values


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------

def test_matching_marker_is_annotated_with_name_and_effects(tmp_path):
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(1881, [ann("PB2", "p.E627K")])])
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": identity_map()})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr

    (field,) = flumut_fields(out)
    protein, pos, alt, names, effects = field.split("|")
    assert (protein, pos, alt, names) == ("PB2", "627", "K", "PB2:E627K")
    assert "Increased_virulence_in_mice" in effects.split("&")
    assert " " not in field and "," not in field  # legal in a VCF INFO value


def test_marker_matches_only_its_own_protein(tmp_path):
    """FluMutDB has M1:N30D but no M2 marker with alt D at position 30.

    Both ANN entries carry position 30 / alt D on the same segment reference (MP), so a
    match keyed on (position, alt) alone would emit two entries. Only M1 may match.
    """
    row = vcf_row(120, [ann("M1", "p.N30D"), ann("M2", "p.A30D")])
    vcf = write_vcf(tmp_path / "in.vcf", [row])
    cmap = write_coord_map(tmp_path / "map.json", "M",
                           {"M1": identity_map(300), "M2": identity_map(100)})

    proc, out = run_script(tmp_path, vcf, cmap, "M")
    assert proc.returncode == 0, proc.stderr

    (field,) = flumut_fields(out)
    assert field.count("|") == 4, f"expected exactly one entry, got {field}"
    assert field.startswith("M1|30|D|M1:N30D|")


def test_ha1_marker_matches_subunit_numbering(tmp_path):
    """HA markers are numbered against the HA1/HA2 subunits, not the full HA CDS."""
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(543, [ann("HA1", "p.T156A")])])
    cmap = write_coord_map(tmp_path / "map.json", "HA",
                           {"HA1": identity_map(340), "HA2": identity_map(230)})

    proc, out = run_script(tmp_path, vcf, cmap, "HA")
    assert proc.returncode == 0, proc.stderr

    (field,) = flumut_fields(out)
    assert field.startswith("HA1|156|A|HA1-5:T156A|")


def test_positions_are_translated_through_the_coordinate_map(tmp_path):
    """A shifted map must move the match, not merely pass positions through.

    Models a two-residue insertion in the target: source 1..600 align 1:1, everything from
    601 on shifts by two. FluMutDB position 627 therefore lands on target 629, so the row
    at 629 matches and an otherwise identical row at 627 does not.
    """
    shifted = {src: (src if src <= 600 else src + 2) for src in range(1, 801)}
    rows = [vcf_row(10, [ann("PB2", "p.E629K")]), vcf_row(20, [ann("PB2", "p.E627K")])]
    vcf = write_vcf(tmp_path / "in.vcf", rows)
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": shifted})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr

    translated, untranslated = flumut_fields(out)
    assert translated.startswith("PB2|629|K|PB2:E627K|")
    assert untranslated is None


def test_no_effects_flag_omits_the_effect_field(tmp_path):
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(1881, [ann("PB2", "p.E627K")])])
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": identity_map()})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2", extra=["--no-effects"])
    assert proc.returncode == 0, proc.stderr

    (field,) = flumut_fields(out)
    assert field == "PB2|627|K|PB2:E627K"


# ---------------------------------------------------------------------------
# Non-matching rows
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("entries", [
    [ann("PB2", "p.E627E", annotation="synonymous_variant")],   # synonymous
    [ann("PB2", "", annotation="intron_variant")],              # no HGVS.p
    [ann("PB2", "p.E627R")],                                    # position curated, alt is not
    [ann("PB2", "p.A5T")],                                      # position not curated
], ids=["synonymous", "intronic", "wrong-alt", "uncurated-position"])
def test_non_matching_rows_are_left_alone(tmp_path, entries):
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(100, entries)])
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": identity_map()})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out) == [None]


def test_unmapped_gene_is_reported_rather_than_silently_skipped(tmp_path):
    """A VCF annotated against a single-HA-gene dataset cannot match FluMutDB's HA1/HA2."""
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(543, [ann("HA", "p.T172A")])])
    cmap = write_coord_map(tmp_path / "map.json", "HA",
                           {"HA1": identity_map(340), "HA2": identity_map(230)})

    proc, out = run_script(tmp_path, vcf, cmap, "HA")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out) == [None]
    assert "HA" in proc.stderr and "no FluMutDB counterpart" in proc.stderr


# ---------------------------------------------------------------------------
# Degenerate inputs
# ---------------------------------------------------------------------------

def test_header_only_vcf_passes_through(tmp_path):
    vcf = write_vcf(tmp_path / "in.vcf", [])
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": identity_map()})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out) == []
    assert any(line.startswith("##INFO=<ID=FLUMUT") for line in out.read_text().splitlines())


def test_segment_without_genes_passes_through(tmp_path):
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(1881, [ann("PB2", "p.E627K")])])
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out) == [None]


def test_segment_absent_from_coordinate_map_passes_through(tmp_path):
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(1881, [ann("PB2", "p.E627K")])])
    cmap = write_coord_map(tmp_path / "map.json", "NA", {"NA": identity_map()})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out) == [None]


def test_stale_annotation_flag_is_warned_about(tmp_path):
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(543, [ann("HA", "p.T172A")])])
    cmap = write_coord_map(tmp_path / "map.json", "HA", {}, stale=True)

    proc, _ = run_script(tmp_path, vcf, cmap, "HA")
    assert proc.returncode == 0, proc.stderr
    assert "stale" in proc.stderr


# ---------------------------------------------------------------------------
# The core guarantee: nothing but INFO changes
# ---------------------------------------------------------------------------

def test_rows_are_annotated_in_place_and_never_added_or_reordered(tmp_path):
    rows = [
        vcf_row(10, [ann("PB2", "p.A5T")]),
        vcf_row(1881, [ann("PB2", "p.E627K")]),
        vcf_row(30, [ann("PB2", "p.E627E", annotation="synonymous_variant")]),
        vcf_row(2000, []),
    ]
    vcf = write_vcf(tmp_path / "in.vcf", rows)
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": identity_map()})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr

    def data_lines(text):
        return [ln for ln in text.splitlines() if not ln.startswith("#")]

    before, after = data_lines(vcf.read_text()), data_lines(out.read_text())
    assert len(after) == len(before) == 4

    for old, new in zip(before, after):
        old_cols, new_cols = old.split("\t"), new.split("\t")
        # Every column except INFO is untouched (POS order included, hence zip).
        assert old_cols[:7] == new_cols[:7]
        assert old_cols[8:] == new_cols[8:]
        # INFO may only gain a ';FLUMUT=...' suffix, so ANN survives verbatim.
        assert new_cols[7] == old_cols[7] or new_cols[7].startswith(old_cols[7] + ";FLUMUT=")

    # Only the curated row gains a FLUMUT key; the others are returned untouched.
    fields = flumut_fields(out)
    assert [f is None for f in fields] == [True, False, True, True]
    assert fields[1].startswith("PB2|627|K|PB2:E627K|")


# ---------------------------------------------------------------------------
# Helper-level checks
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("protein,available,expected", [
    ("NS-1", {"NS1", "NEP"}, "NS1"),
    ("NS-2", {"NS1", "NEP"}, "NEP"),
    ("NS-2", {"NS1", "NS2"}, "NS2"),
    ("HA1", {"HA1", "HA2"}, "HA1"),
    ("HA1", {"HA"}, None),
    ("PB1-F2", {"PB1", "PB1-F2"}, "PB1-F2"),
])
def test_resolve_gene(protein, available, expected):
    assert merge_v2.resolve_gene(protein, available) == expected


def test_database_is_found_without_an_explicit_path():
    """The default path: locate_db() falls back to the installed flumutdb package."""
    pytest.importorskip("flumut")
    assert Path(merge_v2.locate_db(None)).name == "flumut_db.sqlite"


def test_sanitize_produces_a_legal_info_value():
    dirty = "Decreased virus binding to α2–3 (weak); n=2, maybe|not"
    clean = merge_v2.sanitize(dirty)
    assert not any(ch in clean for ch in " \t;,=|")


def test_aa_changes_from_info_skips_synonymous_and_deduplicates():
    info = BASE_INFO + ";ANN=" + ",".join([
        ann("PB2", "p.E627K"),
        ann("PB2", "p.E627K"),                                    # duplicate feature entry
        ann("PB2", "p.G100G", annotation="synonymous_variant"),
        ann("PB1", "p.N375S"),
    ])
    assert merge_v2.aa_changes_from_info(info) == [("PB2", 627, "K"), ("PB1", 375, "S")]
