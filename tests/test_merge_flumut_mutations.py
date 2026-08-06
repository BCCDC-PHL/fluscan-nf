"""Unit tests for bin/merge_flumut_mutations.py (FluMutDB VCF annotation).

Covers the behaviours the merge is specified around: matching on
(protein, AA position, alt AA), translation through the coordinate map's
aa_source_to_target block, canonicalisation of ANN gene names onto the coordinate map's
gene keys, the shape of the FLUMUT INFO field including the direct/associated effect
split, and the guarantee that no row is added, removed, reordered or otherwise altered.

The real FluMutDB shipped by the flumutdb conda package is used rather than a fixture, so
these assertions double as a check that the queries still fit the released schema. The
concrete markers asserted on (FluMutDB 6.6):

    PB2:E627K   reference PB2, protein PB2, position 627, alt K  (both effect classes)
    PB2:S715N   reference PB2, protein PB2, position 715, alt N  (direct only)
    NS-1:A53D   reference NS,  protein NS-1, position 53, alt D  (associated only)
    NS-2:T48A   reference NS,  protein NS-2, position 48, alt A  (associated only; the
                marker ISSUES used to reproduce the NS2/NEP gene-alias bug)
    M1:N30D     reference MP,  protein M1,  position 30,  alt D
    M2:A30V/T/S reference MP,  protein M2,  position 30   (no alt D -- see the M1/M2 test)
    HA1-5:T156A reference HA,  protein HA1, position 156, alt A

Run with:  pytest tests/test_merge_flumut_mutations.py

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
SCRIPT_PATH = REPO_ROOT / "bin" / "merge_flumut_mutations.py"

VCF_COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
BASE_INFO = "NS=1;DP=1;AF=1.0;TYPE=snp;LEN=1"
FORMAT_COLS = "GT:DP\t1:1"


def _load_module():
    spec = importlib.util.spec_from_file_location("merge_flumut_mutations", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


merge_mod = _load_module()


def _find_db() -> str | None:
    if os.environ.get("FLUMUT_DB"):
        return os.environ["FLUMUT_DB"]
    try:
        return merge_mod.locate_db(None)
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

def ann(gene: str, hgvs_p: str, annotation: str = "missense_variant", alt: str = "A",
        feature_id: str | None = None) -> str:
    """A 16-field snpEff-style ANN entry, matching what the converter emits.

    `feature_id` defaults to a locus tag derived from the gene, as a main CDS entry carries.
    Mature-peptide entries pass it explicitly: their Gene_Name is the *parent* gene, and only
    Feature_ID names the subunit.
    """
    return "|".join([
        alt, annotation, "MODERATE", gene, f"LOCUS_{gene}", "transcript",
        feature_id if feature_id is not None else f"LOCUS_{gene}",
        "protein_coding", "1/1", "c.1A>C", hgvs_p,
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


def write_coord_map(path: Path, segment: str, genes: dict, stale: bool = False,
                    aliases: dict | None = None,
                    subunit_parent_genes: list | None = None) -> Path:
    """A coordinate-map shard. `genes` maps gene key -> {target_pos: source_pos} pairs
    expressed the way the real map stores them: aa_source_to_target, keyed by source.

    `aliases` maps a gene key to the ANN Feature_ID values that also identify it (what the
    builder records for a mature peptide); `subunit_parent_genes` names target genes left
    unmapped because their subunits carry the map instead.
    """
    aliases = aliases or {}
    doc = {
        "schema_version": "2.2",
        "source_reference": "FLUMUT",
        "target_reference": "TEST",
        "segments": {
            segment: {
                "stale_annotation": stale,
                "subunit_parent_genes": subunit_parent_genes or [],
                "genes": {
                    gene: {
                        "aa_source_to_target": {str(src): tgt for src, tgt in mapping.items()},
                        **({"aliases": aliases[gene]} if gene in aliases else {}),
                    }
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
    protein, pos, alt, names, direct, associated = field.split("|")
    assert (protein, pos, alt, names) == ("PB2", "627", "K", "PB2:E627K")
    assert "H5N1:Increased_virulence_in_mice" in direct.split("&")
    assert " " not in field and "," not in field  # legal in a VCF INFO value


def test_effects_are_split_into_direct_and_associated(tmp_path):
    """PB2:E627K is carried by both single- and multi-mutation markers.

    Every value in either class must be subtype-qualified, and the two classes are
    populated independently -- an effect may be directly attributable in one subtype and
    only associated in another.
    """
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(1881, [ann("PB2", "p.E627K")])])
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": identity_map()})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr

    (field,) = flumut_fields(out)
    _, _, _, _, direct, associated = field.split("|")
    assert direct != "." and associated != "."
    for value in direct.split("&") + associated.split("&"):
        subtype, _, effect = value.partition(":")
        assert subtype and effect, f"{value} is not Subtype:Effect"
    # Each class is sorted and deduplicated.
    for values in (direct.split("&"), associated.split("&")):
        assert values == sorted(set(values))


def test_single_mutation_marker_populates_only_direct_effect(tmp_path):
    """PB2:S715N defines a marker on its own, so its effect is directly attributable."""
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(2159, [ann("PB2", "p.S715N")])])
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": identity_map()})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr

    (field,) = flumut_fields(out)
    assert field == "PB2|715|N|PB2:S715N|H5N1:Decreased_virulence_in_mice|."


def test_multi_mutation_marker_populates_only_associated_effect(tmp_path):
    """NS-1:A53D appears only in a marker requiring further mutations.

    The mutation is still detected and named -- 118 of FluMutDB 6.6's 333 mutations have no
    single-mutation marker, and suppressing them would discard their only evidence -- but
    the effect lands in AssociatedEffect, which does not claim the whole marker is present.
    """
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(200, [ann("NS1", "p.A53D")])])
    cmap = write_coord_map(tmp_path / "map.json", "NS",
                           {"NS1": identity_map(240), "NEP": identity_map(130)})

    proc, out = run_script(tmp_path, vcf, cmap, "NS")
    assert proc.returncode == 0, proc.stderr

    (field,) = flumut_fields(out)
    assert field == "NS-1|53|D|NS-1:A53D|.|H1N1:Increased_virulence_in_mice"


def test_mutation_without_any_effect_is_still_detected(tmp_path):
    """PB2:M535I is curated but carries no published effect in FluMutDB 6.6.

    Detection is keyed on the mutation, not on the effect join, so the marker is named and
    both effect classes are empty. A non-LEFT join anywhere in MARKER_QUERY would lose it.
    """
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(1605, [ann("PB2", "p.M535I")])])
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": identity_map()})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out) == ["PB2|535|I|PB2:M535I|.|."]


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
    assert field.count("|") == 5, f"expected exactly one entry, got {field}"
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


# ---------------------------------------------------------------------------
# Gene-name canonicalisation (ANN Gene_Name -> coordinate-map gene key)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("ann_gene", ["NS2", "NS-2", "NEP"])
def test_nep_annotates_under_every_spelling(tmp_path, ann_gene):
    """The map builder canonicalises NS2/NS-2 to NEP; the merge must do the same.

    Without matching canonicalisation a target dataset spelling the gene NS2 -- which the
    flumut source reference itself does -- produces a map keyed NEP and a VCF gene keyed
    NS2, and every NEP marker is dropped with no indication of why.
    """
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(300, [ann(ann_gene, "p.T48A")])])
    cmap = write_coord_map(tmp_path / "map.json", "NS",
                           {"NEP": identity_map(130), "NS1": identity_map(240)})

    proc, out = run_script(tmp_path, vcf, cmap, "NS")
    assert proc.returncode == 0, proc.stderr

    (field,) = flumut_fields(out)
    assert field == "NS-2|48|A|NS-2:T48A|.|H5N1:Decreased_antiviral_response_in_ferrets"
    assert "no FluMutDB counterpart" not in proc.stderr


@pytest.mark.parametrize("ann_gene", ["NS1", "NS-1"])
def test_ns1_annotates_under_every_spelling(tmp_path, ann_gene):
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(200, [ann(ann_gene, "p.A53D")])])
    cmap = write_coord_map(tmp_path / "map.json", "NS",
                           {"NS1": identity_map(240), "NEP": identity_map(130)})

    proc, out = run_script(tmp_path, vcf, cmap, "NS")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out) == ["NS-1|53|D|NS-1:A53D|.|H1N1:Increased_virulence_in_mice"]


def test_unmapped_gene_warning_leads_with_the_general_condition(tmp_path):
    """An unresolvable non-HA gene must not be explained in terms of HA.

    The HA hint is the single most likely cause but not the only one; leading with it sent
    the ISSUES reproduction of the NS2/NEP bug looking at the wrong segment.
    """
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(300, [ann("BOGUS", "p.T48A")])])
    cmap = write_coord_map(tmp_path / "map.json", "NS",
                           {"NEP": identity_map(130), "NS1": identity_map(240)})

    proc, out = run_script(tmp_path, vcf, cmap, "NS")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out) == [None]
    assert "BOGUS" in proc.stderr and "must match one of the map's gene keys" in proc.stderr
    assert "HA1/HA2" not in proc.stderr


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
# Mature-peptide subunits (ANN Feature_ID -> coordinate-map gene key)
# ---------------------------------------------------------------------------
#
# NCBI/RefSeq flu annotations -- the GsGD dataset among them -- model HA as one CDS with
# HA1/HA2 as mature_protein_region_of_CDS children. The converter emits those entries with
# the *parent* gene in Gene_Name and the peptide's own identifier in Feature_ID, with AA
# numbering already relative to the subunit, i.e. FluMutDB's own numbering. Keying on
# Gene_Name alone put all 76 HA mappings out of reach on that dataset.

def test_mature_peptide_resolves_through_its_protein_id(tmp_path):
    vcf = write_vcf(tmp_path / "in.vcf",
                    [vcf_row(600, [ann("HA", "p.T195I", feature_id="YP_529486.1")])])
    cmap = write_coord_map(tmp_path / "map.json", "HA",
                           {"HA1": identity_map(340), "HA2": identity_map(230)},
                           aliases={"HA1": ["YP_529486.1"]},
                           subunit_parent_genes=["HA"])

    proc, out = run_script(tmp_path, vcf, cmap, "HA")
    assert proc.returncode == 0, proc.stderr

    (field,) = flumut_fields(out)
    protein, pos, alt, names = field.split("|")[:4]
    assert (protein, pos, alt, names) == ("HA1", "195", "I", "HA1-5:T195I")


def test_mature_peptide_resolves_through_its_product_name(tmp_path):
    """A mature row with no protein_id: load_gff3 falls back to `product`, so Feature_ID is
    the subunit name itself and must resolve as a gene key in its own right."""
    vcf = write_vcf(tmp_path / "in.vcf",
                    [vcf_row(600, [ann("HA", "p.T195I", feature_id="HA1")])])
    cmap = write_coord_map(tmp_path / "map.json", "HA",
                           {"HA1": identity_map(340), "HA2": identity_map(230)},
                           subunit_parent_genes=["HA"])

    proc, out = run_script(tmp_path, vcf, cmap, "HA")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out)[0].startswith("HA1|195|I|HA1-5:T195I")


def test_parent_cds_entry_alongside_a_subunit_entry_is_not_warned_about(tmp_path):
    """Every in-HA1 row on a mature-peptide dataset also carries the parent HA CDS's entry,
    which cannot resolve by design. Warning about it would flag rows annotated correctly."""
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(600, [
        ann("HA", "p.T211I"),                                 # parent CDS, HA numbering
        ann("HA", "p.T195I", feature_id="YP_529486.1"),       # HA1 mature peptide
    ])])
    cmap = write_coord_map(tmp_path / "map.json", "HA",
                           {"HA1": identity_map(340), "HA2": identity_map(230)},
                           aliases={"HA1": ["YP_529486.1"]},
                           subunit_parent_genes=["HA"])

    proc, out = run_script(tmp_path, vcf, cmap, "HA")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out)[0].startswith("HA1|195|I|HA1-5:T195I")
    assert "no FluMutDB counterpart" not in proc.stderr


def test_unresolvable_gene_is_still_reported_when_nothing_on_the_row_resolves(tmp_path):
    """The per-row suppression above must not silence a genuinely unmapped gene."""
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(600, [ann("BOGUS", "p.T195I")])])
    cmap = write_coord_map(tmp_path / "map.json", "HA",
                           {"HA1": identity_map(340), "HA2": identity_map(230)},
                           subunit_parent_genes=["HA"])

    proc, out = run_script(tmp_path, vcf, cmap, "HA")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out) == [None]
    assert "BOGUS" in proc.stderr and "no FluMutDB counterpart" in proc.stderr


# ---------------------------------------------------------------------------
# HGVS.p parsing
# ---------------------------------------------------------------------------

def test_three_letter_hgvs_p_annotates_identically_to_one_letter(tmp_path):
    """snpEff's default spelling must not be silently unmatchable: the merge process is
    meant to be reusable from the snpEff-annotated branches."""
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": identity_map()})

    fields = []
    for hgvs_p in ("p.E627K", "p.Glu627Lys"):
        vcf = write_vcf(tmp_path / f"in_{hgvs_p}.vcf", [vcf_row(1881, [ann("PB2", hgvs_p)])])
        proc, out = run_script(tmp_path, vcf, cmap, "PB2")
        assert proc.returncode == 0, proc.stderr
        fields.append(flumut_fields(out))

    assert fields[0] == fields[1]
    assert fields[0][0].startswith("PB2|627|K|PB2:E627K")


def test_three_letter_synonymous_change_is_not_annotated(tmp_path):
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(1881, [ann("PB2", "p.Glu627Glu",
                                                             annotation="synonymous_variant")])])
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": identity_map()})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out) == [None]


def test_unparseable_hgvs_p_is_reported_rather_than_silently_ignored(tmp_path):
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(1881, [ann("PB2", "p.Xyzzy627Plugh")])])
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": identity_map()})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out) == [None]
    assert "cannot parse" in proc.stderr and "p.Xyzzy627Plugh" in proc.stderr


@pytest.mark.parametrize("hgvs_p", ["p.Lys5del", "p.Glu627=", "p.Ala10fs"])
def test_non_substitution_hgvs_p_is_not_reported_as_a_parse_failure(tmp_path, hgvs_p):
    """del/ins/fs and snpEff's '=' synonymous shorthand are legitimately skipped; counting
    them as failures would make the warning fire on ordinary VCFs."""
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(1881, [ann("PB2", hgvs_p)])])
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": identity_map()})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out) == [None]
    assert "cannot parse" not in proc.stderr


# ---------------------------------------------------------------------------
# MNV provenance
# ---------------------------------------------------------------------------

def test_mnv_unresolved_rows_are_skipped(tmp_path):
    """fix_mnv_codons.py sets this flag on rows whose HGVS.p it knows to be wrong but could
    not correct. Matching a curated marker against one manufactures a false positive."""
    row = vcf_row(1881, [ann("PB2", "p.E627K")]).replace(BASE_INFO,
                                                         BASE_INFO + ";MNV_UNRESOLVED")
    vcf = write_vcf(tmp_path / "in.vcf", [row])
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": identity_map()})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out) == [None]
    assert "MNV_UNRESOLVED" in proc.stderr


def test_mnv_corrected_rows_are_still_annotated(tmp_path):
    """Only the *unresolved* flag suppresses a row; a corrected MNV carries a HGVS.p the
    pipeline stands behind. The substring must not match the wrong key."""
    row = vcf_row(1881, [ann("PB2", "p.E627K")]).replace(BASE_INFO,
                                                         BASE_INFO + ";MNV_CORRECTED")
    vcf = write_vcf(tmp_path / "in.vcf", [row])
    cmap = write_coord_map(tmp_path / "map.json", "PB2", {"PB2": identity_map()})

    proc, out = run_script(tmp_path, vcf, cmap, "PB2")
    assert proc.returncode == 0, proc.stderr
    assert flumut_fields(out)[0].startswith("PB2|627|K|PB2:E627K")


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
    # The bare stem is tried last, so a specific match always wins over it. A dataset naming
    # its NS1 gene 'NS' offers that stem to NS-1 and NS-2 alike; resolving NS-2 to it would
    # query NEP markers through NS1's coordinate map -- a false positive, not a silent drop.
    ("NS-2", {"NS", "NEP"}, "NEP"),
    ("NS-1", {"NS", "NS1"}, "NS1"),
    ("HA1-5", {"HA1", "HA2"}, "HA1"),
    # Only when nothing more specific exists does the stem win.
    ("NS-1", {"NS"}, "NS"),
])
def test_resolve_gene(protein, available, expected):
    assert merge_mod.resolve_gene(protein, available) == expected


def test_colliding_gene_keys_are_refused_rather_than_silently_overwritten(tmp_path):
    """Two FluMutDB proteins resolving to one gene key must not overwrite each other.

    NS1 (230 aa) and NEP (121 aa) are unrelated in numbering, so annotating one through the
    other's map produces wrong positions rather than none -- the worst failure mode here.
    """
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(300, [ann("NS", "p.T48A")])])
    cmap = write_coord_map(tmp_path / "map.json", "NS", {"NS": identity_map(240)})

    proc, out = run_script(tmp_path, vcf, cmap, "NS")
    assert proc.returncode == 0, proc.stderr
    assert "both resolve to coordinate-map gene key NS" in proc.stderr
    # NS-1 wins (PROTEIN_QUERY is ORDER BY'd), so the NS-2 marker at 48 must not appear.
    assert "NS-2" not in (flumut_fields(out)[0] or "")


@pytest.mark.parametrize("raw,expected", [
    ("NS2", "NEP"),
    ("NS-2", "NEP"),
    ("NEP", "NEP"),
    ("NS-1", "NS1"),
    ("NS1", "NS1"),
    ("HA1", "HA1"),
    ("PB2", "PB2"),
])
def test_canonical_gene(raw, expected):
    assert merge_mod.canonical_gene(raw) == expected


def test_canonical_gene_comes_from_the_shared_module():
    """The merge must not carry its own copy of the alias table.

    A second copy is how the two sides drift: the map ends up keyed under a name the merge
    never looks up, which is exactly how every NEP marker went missing on NS2 datasets.
    bin/ is already on sys.path -- merge_flumut_mutations put it there when it imported
    this same module.
    """
    from lib import flu_gene_names

    assert merge_mod.canonical_gene is flu_gene_names.canonical_gene


def test_ann_parsing_comes_from_the_shared_module():
    """The merge must not carry its own copy of the ANN/HGVS.p mechanics either.

    bin/mutation_watcher.py asks the same questions of the same VCFs against a different
    curated source. Two copies are free to drift over which HGVS.p spellings parse and
    which non-substitution forms are skipped, and a disagreement there is silent -- the
    same failure mode as a second alias table.
    """
    from lib import aa_markers

    assert merge_mod.parse_hgvs_p is aa_markers.parse_hgvs_p
    assert merge_mod.ann_entries is aa_markers.ann_entries
    assert merge_mod.has_info_flag is aa_markers.has_info_flag
    assert merge_mod.aa_changes_from_info is aa_markers.aa_changes_from_info


def test_database_is_found_without_an_explicit_path():
    """The default path: locate_db() falls back to the installed flumutdb package."""
    pytest.importorskip("flumut")
    assert Path(merge_mod.locate_db(None)).name == "flumut_db.sqlite"


def test_sanitize_produces_a_legal_info_value():
    dirty = "Decreased virus binding to α2–3 (weak); n=2, maybe|not"
    clean = merge_mod.sanitize(dirty)
    assert not any(ch in clean for ch in " \t;,=|")


def test_aa_changes_from_info_skips_synonymous_and_deduplicates():
    info = BASE_INFO + ";ANN=" + ",".join([
        ann("PB2", "p.E627K"),
        ann("PB2", "p.E627K"),                                    # duplicate feature entry
        ann("PB2", "p.G100G", annotation="synonymous_variant"),
        ann("PB1", "p.N375S"),
    ])
    assert merge_mod.aa_changes_from_info(info) == [
        ("PB2", "LOCUS_PB2", 627, "K"),
        ("PB1", "LOCUS_PB1", 375, "S"),
    ]


def test_aa_changes_keeps_one_entry_per_feature_of_the_same_gene():
    """A mature peptide and its parent CDS share a Gene_Name but are different features, so
    de-duplication must key on the feature too or the subunit entry is dropped."""
    info = BASE_INFO + ";ANN=" + ",".join([
        ann("HA", "p.T211I"),
        ann("HA", "p.T195I", feature_id="YP_529486.1"),
    ])
    assert merge_mod.aa_changes_from_info(info) == [
        ("HA", "LOCUS_HA", 211, "I"),
        ("HA", "YP_529486.1", 195, "I"),
    ]


@pytest.mark.parametrize("hgvs_p,expected", [
    ("p.E627K", ("E", 627, "K")),
    ("p.Glu627Lys", ("E", 627, "K")),
    ("p.Ter5Gln", ("*", 5, "Q")),
    ("p.Lys5del", None),
    ("p.Glu627=", None),
    ("", None),
])
def test_parse_hgvs_p(hgvs_p, expected):
    assert merge_mod.parse_hgvs_p(hgvs_p) == expected


@pytest.mark.parametrize("info,key,expected", [
    ("NS=1;MNV_UNRESOLVED;DP=3", "MNV_UNRESOLVED", True),
    ("NS=1;MNV_UNRESOLVED=1", "MNV_UNRESOLVED", True),
    ("NS=1;MNV_CORRECTED;DP=3", "MNV_UNRESOLVED", False),
    # A substring test would match this; a token test must not.
    ("NS=1;ANN=A|x|y|MNV_UNRESOLVED", "MNV_UNRESOLVED", False),
])
def test_has_info_flag(info, key, expected):
    assert merge_mod.has_info_flag(info, key) is expected
