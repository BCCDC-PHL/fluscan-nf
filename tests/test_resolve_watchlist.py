"""Tests for scripts/resolve_watchlist.py.

The bundle here is synthetic: eight segment annotations built from known proteins, so a
position's residue is a fact the test states rather than one it looks up. HA carries its
subunits as mature peptides, which is the shape every real bundle uses.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


def _load_module():
    """Load the script by path; scripts/ is not a package and is not on sys.path.

    Registered in sys.modules before exec because @dataclass resolves __module__ through
    it while the class body is still running.
    """
    path = Path(__file__).resolve().parents[1] / "scripts" / "resolve_watchlist.py"
    spec = importlib.util.spec_from_file_location("resolve_watchlist", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


rw = _load_module()

# One codon per residue, standard table, so a protein round-trips through translate_model.
CODONS = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
    "I": "ATT", "K": "AAA", "L": "CTT", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA",
    "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}

# HA is 5 signal residues, then a 10-residue HA1, then an 8-residue HA2. HA1 position 3 is
# D while the parent HA position 3 is T, which is what separates the two numberings.
HA_SIGNAL, HA1, HA2 = "MKTII", "ACDEFGHIKL", "MNPQRSTV"

PROTEINS = {
    "PB2": "MERIKELRDLM",
    "PB1": "MDVNPTLLFLK",
    "PA": "MEDFVRQCFNP",
    "HA": HA_SIGNAL + HA1 + HA2,
    "NP": "MASQGTKRSYE",
    "NA": "MKTIIALSYIFCLVFAQKLP",
    "M": "MSLLTEVETYV",
    "NS": "MDPNTVSSFQV",
}

# The gene each segment's single CDS is named for; M and NS carry one product here rather
# than their real two, so M2 and NEP double as genes a bundle can lack.
GENES = {"PB2": "PB2", "PB1": "PB1", "PA": "PA", "HA": "HA",
         "NP": "NP", "NA": "NA", "M": "M1", "NS": "NS1"}


def back_translate(protein: str) -> str:
    return "".join(CODONS[residue] for residue in protein) + "TAA"


def write_segment(directory: Path, segment: str, protein: str) -> None:
    directory.mkdir(parents=True)
    sequence = back_translate(protein)
    contig = f"TEST_{segment}"
    (directory / "reference.fasta").write_text(f">{contig}\n{sequence}\n")

    gene = GENES[segment]
    rows = [
        f"{contig}\tTest\tCDS\t1\t{len(sequence)}\t.\t+\t0\t"
        f"ID=cds-{gene};gene={gene};product={gene}"
    ]
    if segment == "HA":
        # Mature peptides must start on a codon boundary within the parent to be promoted.
        ha1_start = len(HA_SIGNAL) * 3 + 1
        ha1_end = ha1_start + len(HA1) * 3 - 1
        ha2_end = ha1_end + len(HA2) * 3
        rows += [
            f"{contig}\tTest\tmature_protein_region_of_CDS\t{ha1_start}\t{ha1_end}\t.\t+\t.\t"
            f"ID=id-HA1;Parent=cds-HA;product=HA1",
            f"{contig}\tTest\tmature_protein_region_of_CDS\t{ha1_end + 1}\t{ha2_end}\t.\t+\t.\t"
            f"ID=id-HA2;Parent=cds-HA;product=HA2",
        ]
    (directory / "genome_annotation.gff3").write_text("##gff-version 3\n" + "\n".join(rows) + "\n")


@pytest.fixture
def bundle(tmp_path) -> Path:
    root = tmp_path / "testref"
    (root).mkdir()
    (root / "db.yml").write_text('name: testref\nversion: "2026-01-01"\n')
    for segment, protein in PROTEINS.items():
        write_segment(root / "nextclade" / segment, segment, protein)
    return root


@pytest.fixture
def index(bundle) -> dict:
    return rw.build_gene_index(bundle)


def resolve_one(index, gene, mutation, offsets=(0,), extras=None, row=None):
    row = row or {"gene": gene, "mutation": mutation}
    mapping = {"gene": "gene", "mutation": "mutation"}
    return rw.resolve(row, mapping, list(extras or []), index, list(offsets))


def write_source(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


# --- the gene index -------------------------------------------------------------------

def test_the_gene_index_is_discovered_from_every_segment(index):
    assert set(index) == {"PB2", "PB1", "PA", "HA", "HA1", "HA2", "NP", "NA", "M1", "NS1"}


def test_ha_subunits_resolve_as_their_own_proteins(index):
    """The HA trap: HA1 numbering starts after the signal peptide, so it is a third protein."""
    assert index["HA"][1] == PROTEINS["HA"]
    assert index["HA1"][1] == HA1
    assert index["HA2"][1] == HA2
    assert index["HA1"][0] == "HA"


def test_a_subunit_position_is_not_the_parent_position(index):
    """HA1 3 and HA 3 are different residues, which is the whole reason they are separate."""
    assert resolve_one(index, "HA1", "D3A").status == "ok"
    assert resolve_one(index, "HA", "T3A").status == "ok"
    # Keying an HA1 marker on the parent would land on T, not D.
    assert resolve_one(index, "HA", "D3A").status == "ref_differs"


def test_a_gene_the_bundle_lacks_is_reported_not_skipped(index):
    """wisconsin588 has no PB1-F2 annotation; reporting nothing silently would hide that."""
    assert resolve_one(index, "PB1-F2", "N66S").status == "gene_not_in_reference"


def test_a_gene_on_two_segments_is_an_error(bundle):
    extra = bundle / "nextclade" / "NP" / "genome_annotation.gff3"
    extra.write_text(
        extra.read_text() + "TEST_NP\tTest\tCDS\t1\t9\t.\t+\t0\tID=cds-NA;gene=NA;product=NA\n"
    )
    with pytest.raises(rw.SourceError, match="annotated on both"):
        rw.build_gene_index(bundle)


# --- statuses -------------------------------------------------------------------------

def test_an_agreeing_entry_resolves_unchanged(index):
    resolution = resolve_one(index, "NA", "I10V")
    assert (resolution.status, resolution.target_position, resolution.mutation) == ("ok", 10, "I10V")


def test_a_differing_wild_type_is_kept_and_rewritten(index):
    """What a sample ends up with is the point, so lineage drift is not a reason to drop it."""
    resolution = resolve_one(index, "NA", "A10V")
    assert resolution.status == "ref_differs"
    assert resolution.mutation == "I10V"
    assert "native=A10V" in resolution.note


def test_an_entry_the_reference_already_carries_is_excluded(index):
    """No variant can be called where the sample already matches, so it cannot be detected."""
    resolution = resolve_one(index, "NA", "A10I")
    assert resolution.status == "ref_is_alt"
    assert resolution.status not in rw.INCLUDED_STATUSES


def test_a_position_past_the_protein_is_reported(index):
    assert resolve_one(index, "NA", "I500V").status == "out_of_range"


def test_an_unparseable_or_unknown_row_is_reported_not_fatal(index):
    """A curation run reports every bad row at once rather than stopping at the first."""
    assert resolve_one(index, "NA", "Del 247-250").status == "unparseable"
    assert resolve_one(index, "NA", "I10I").status == "synonymous"
    assert resolve_one(index, "XYZ", "I10V").status == "unknown_gene"


# --- offsets --------------------------------------------------------------------------

def test_the_offset_that_matches_the_reference_residue_is_chosen(index):
    """NA 5 is I and NA 10 is I too, but only one candidate carries the claimed Y."""
    # NA: M K T I I A L S Y I F C L V F A Q K L P -- position 9 is Y, 4 is I.
    resolution = resolve_one(index, "NA", "Y4A", offsets=(0, 5))
    assert (resolution.status, resolution.offset, resolution.target_position) == ("ok", 5, 9)


def test_several_fitting_offsets_are_excluded_rather_than_guessed(index):
    """NA 4 and NA 5 are both I, so an I-at-4 entry with offset 0 or 1 is unresolvable."""
    resolution = resolve_one(index, "NA", "I4V", offsets=(0, 1))
    assert resolution.status == "ambiguous"


def test_no_fitting_offset_is_excluded(index):
    resolution = resolve_one(index, "NA", "W4V", offsets=(0, 5))
    assert resolution.status == "no_fit"


def test_a_single_declared_offset_is_trusted_rather_than_requiring_a_fit(index):
    """With nothing to disambiguate, the declared frame is taken at its word and reported."""
    resolution = resolve_one(index, "NA", "W4V", offsets=(5,))
    assert (resolution.status, resolution.target_position) == ("ref_differs", 9)


# --- deduplication --------------------------------------------------------------------

def test_duplicates_merge_instead_of_being_silently_dropped(index):
    """The watcher collapses repeated keys itself and loses the second note."""
    first = resolve_one(index, "NA", "I10V")
    second = rw.resolve(
        {"gene": "NA", "mutation": "A10V"}, {"gene": "gene", "mutation": "mutation"},
        [], index, [0],
    )
    merged = rw.deduplicate([first, second])
    assert len(merged) == 1
    assert "also_in=A10V" in merged[0].note


# --- the source list contract ---------------------------------------------------------

def test_columns_may_be_named_anything_and_extras_are_carried(tmp_path, index):
    source = write_source(
        tmp_path / "custom.csv",
        "# numbering: testref\n"
        "# columns: gene=Protein, mutation=Substitution\n"
        "Protein,Substitution,Drug,Citation\n"
        "NA,I10V,oseltamivir,Smith 2024\n",
    )
    directives, rows, extras = rw.load_source(source)
    assert extras == ["Drug", "Citation"]
    resolution = rw.resolve(rows[0], directives["_columns"], extras, index, [0])
    assert resolution.status == "ok"
    assert "Drug=oseltamivir" in resolution.note
    assert "Citation=Smith 2024" in resolution.note


def test_a_source_without_a_numbering_header_is_fatal(tmp_path):
    source = write_source(tmp_path / "s.csv", "gene,mutation\nNA,I10V\n")
    with pytest.raises(rw.SourceError, match="numbering"):
        rw.load_source(source)


def test_a_column_named_in_the_header_but_absent_is_fatal(tmp_path):
    source = write_source(
        tmp_path / "s.csv",
        "# numbering: testref\n# columns: gene=Protein, mutation=Sub\ngene,mutation\nNA,I10V\n",
    )
    with pytest.raises(rw.SourceError, match="not in the file"):
        rw.load_source(source)


def test_offsets_default_to_zero_and_parse_a_list(tmp_path):
    assert rw.parse_offsets("0, 20") == [0, 20]
    with pytest.raises(rw.SourceError, match="comma-separated integers"):
        rw.parse_offsets("twenty")


# --- the bundle contract --------------------------------------------------------------

def test_a_bundle_without_db_yml_is_refused(tmp_path):
    """Legacy bundles carry no db.yml, and a list cannot record a version it never saw."""
    root = tmp_path / "legacy"
    root.mkdir()
    with pytest.raises(rw.SourceError, match="no db.yml"):
        rw.read_db_metadata(root)


def test_db_yml_must_name_and_version_the_bundle(tmp_path):
    root = tmp_path / "partial"
    root.mkdir()
    (root / "db.yml").write_text("name: partial\n")
    with pytest.raises(rw.SourceError, match="version"):
        rw.read_db_metadata(root)


# --- end to end -----------------------------------------------------------------------

def test_a_generated_list_carries_the_header_the_watcher_requires(tmp_path, bundle, monkeypatch):
    source = write_source(
        tmp_path / "panel_source.csv",
        "# numbering: testref\n# source: a paper\ngene,mutation\nNA,I10V\nPB2,R3K\nHA1,D3A\n",
    )
    out = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", [
        "resolve_watchlist.py", "--source", str(source), "--db", str(bundle),
        "--out-dir", str(out),
    ])
    assert rw.main() == 0

    generated = (out / "panel_testref_watchlist.csv").read_text().splitlines()
    assert "# numbering: testref" in generated
    assert "# db-version: 2026-01-01" in generated
    assert "# source: a paper" in generated
    assert "gene,mutation,note" in generated
    assert any(line.startswith("HA1,D3A") for line in generated)
    assert (out / "panel_testref_reference_map.tsv").exists()
