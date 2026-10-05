"""Unit tests for bin/mutation_watcher.py.

Run with:  pytest tests/test_mutation_watcher.py

Most of what is asserted here fails *silently* when it breaks -- the tool reports zero
hits, which is also what a clean sample reports. That is why the gene-resolution and
MNV_UNRESOLVED cases are tested as thoroughly as the matching itself.
"""

import csv
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "bin" / "mutation_watcher.py"

VCF_COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
BASE_INFO = "NS=1;DP=1;AF=1.0;TYPE=snp;LEN=1"
FORMAT_COLS = "GT:DP\t1:1"


def _load_module():
    sys.path.insert(0, str(REPO_ROOT / "bin"))
    spec = importlib.util.spec_from_file_location("mutation_watcher", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    # Registered before exec_module because @dataclass resolves its class's __module__
    # through sys.modules while the class body is being processed.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


watcher = _load_module()


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def ann(gene: str, hgvs_p: str, feature_id: str | None = None,
        annotation: str = "missense_variant") -> str:
    """A 16-field snpEff-style ANN entry.

    `feature_id` defaults to a locus tag. Mature-peptide entries pass it explicitly: their
    Gene_Name is the *parent* gene, and only Feature_ID names the subunit.
    """
    return "|".join([
        "C", annotation, "MODERATE", gene, f"LOCUS_{gene}", "transcript",
        feature_id if feature_id is not None else f"LOCUS_{gene}",
        "protein_coding", "1/1", "c.1A>C", hgvs_p,
        "1/999", "1/999", "1/333", "", "",
    ])


def vcf_row(pos: int, ann_entries: list[str], extra_info: str = "") -> str:
    info = BASE_INFO
    if extra_info:
        info += f";{extra_info}"
    if ann_entries:
        info += ";ANN=" + ",".join(ann_entries)
    return f"CHROM_X\t{pos}\t.\tA\tC\t60\tPASS\t{info}\t{FORMAT_COLS}"


def write_vcf(path: Path, rows: list[str]) -> Path:
    header = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=CHROM_X,length=2400>",
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations">',
        VCF_COLUMNS,
    ]
    path.write_text("\n".join(header + rows) + "\n")
    return path


def write_watchlist(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def scan(tmp_path, rows: list[str], watchlist: str, sample="SAMPLE-1", segment="PB2"):
    """Run a scan. Returns (hit rows, observed reference residues, rows skipped)."""
    vcf = write_vcf(tmp_path / "in.vcf", rows)
    entries = watcher.load_watchlist(write_watchlist(tmp_path / "w.csv", watchlist))
    watched = {entry.gene for entry in entries.values()}
    changes, observed_ref, skipped = watcher.index_vcf(vcf, watched)
    return watcher.find_hits(entries, changes, sample, segment), observed_ref, skipped


def run_cli(*args) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *map(str, args)],
        capture_output=True, text=True,
    )


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------

def test_a_watched_substitution_is_a_hit(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")])],
                      "gene,mutation,note\nPB2,E627K,mammalian\n")
    assert len(hits) == 1
    assert hits[0]["gene"] == "PB2"
    assert hits[0]["mutation"] == "E627K"
    assert hits[0]["note"] == "mammalian"


def test_the_sample_id_and_segment_are_reported_verbatim(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")])],
                      "gene,mutation\nPB2,E627K\n", sample="WIN-AH-2024_1-OS", segment="PB2")
    assert hits[0]["sample_id"] == "WIN-AH-2024_1-OS"
    assert hits[0]["segment"] == "PB2"


def test_a_wrong_alt_residue_is_not_a_hit(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E627R")])],
                      "gene,mutation\nPB2,E627K\n")
    assert hits == []


def test_a_wrong_position_is_not_a_hit(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E628K")])],
                      "gene,mutation\nPB2,E627K\n")
    assert hits == []


def test_the_three_letter_spelling_gives_the_identical_hit(tmp_path):
    """snpEff emits p.Glu627Lys; the nextclade converter emits p.E627K. Both must match."""
    one, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")])],
                     "gene,mutation\nPB2,E627K\n")
    three, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.Glu627Lys")])],
                       "gene,mutation\nPB2,Glu627Lys\n")
    assert one[0]["gene"] == three[0]["gene"]
    assert one[0]["position"] == three[0]["position"]
    assert one[0]["alt_aa"] == three[0]["alt_aa"]


def test_a_synonymous_change_is_not_a_hit(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E627E")])],
                      "gene,mutation\nPB2,E627K\n")
    assert hits == []


@pytest.mark.parametrize("hgvs_p", ["p.Glu627del", "p.Glu627fs", "p.Glu627="])
def test_non_substitutions_are_never_hits(tmp_path, hgvs_p):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", hgvs_p)])],
                      "gene,mutation\nPB2,E627K\n")
    assert hits == []


def test_a_gene_the_watchlist_does_not_name_is_ignored(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("CHROM_X", "p.E627K")])],
                      "gene,mutation\nPB2,E627K\n")
    assert hits == []


# ---------------------------------------------------------------------------
# Gene resolution
# ---------------------------------------------------------------------------

def test_a_gene_alias_matches(tmp_path):
    """The watchlist's NS-2 and the VCF's NEP are the same protein."""
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("NEP", "p.M51I")])],
                      "gene,mutation\nNS-2,M51I\n", segment="NS")
    assert [h["gene"] for h in hits] == ["NEP"]


def test_the_reverse_alias_direction_also_matches(tmp_path):
    """Canonicalisation runs on the annotation too, so an NCBI-style NS2 resolves."""
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("NS2", "p.M51I")])],
                      "gene,mutation\nNEP,M51I\n", segment="NS")
    assert [h["gene"] for h in hits] == ["NEP"]


def test_a_mature_peptide_resolves_through_feature_id(tmp_path):
    """An NCBI-style HA annotation names the parent gene; only Feature_ID names the subunit."""
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("HA", "p.T195I", feature_id="HA1")])],
                      "gene,mutation\nHA1,T195I\n", segment="HA")
    assert [h["gene"] for h in hits] == ["HA1"]


def test_a_parent_gene_never_matches_a_subunit_watchlist_entry(tmp_path):
    """HA1 numbering starts after HA's 16-residue signal peptide, so this must not match."""
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("HA", "p.T195I")])],
                      "gene,mutation\nHA1,T195I\n", segment="HA")
    assert hits == []


def test_a_subunit_annotation_never_matches_a_parent_watchlist_entry(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("HA1", "p.T195I")])],
                      "gene,mutation\nHA,T195I\n", segment="HA")
    assert hits == []


def test_the_feature_id_wins_when_the_list_carries_both_parent_and_subunit(tmp_path):
    """Gene_Name=HA with Feature_ID=HA1 is an HA1-numbered change, so it resolves to HA1.

    Resolving it to bare HA instead would be wrong by the signal peptide at 93% of
    positions, and a list carrying both keys is the only case where the order is visible.
    """
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("HA", "p.Q226L", feature_id="HA1")])],
                      "gene,mutation\nHA,Q226L\nHA1,Q226L\n", segment="HA")
    assert [h["gene"] for h in hits] == ["HA1"]


def test_a_row_carrying_both_the_parent_and_the_subunit_hits_once(tmp_path):
    hits, _, _ = scan(
        tmp_path,
        [vcf_row(100, [ann("HA", "p.T195I"), ann("HA", "p.T195I", feature_id="HA1")])],
        "gene,mutation\nHA1,T195I\n", segment="HA",
    )
    assert len(hits) == 1


# ---------------------------------------------------------------------------
# The reference residue
# ---------------------------------------------------------------------------

def test_a_differing_reference_residue_is_still_a_hit(tmp_path):
    """Matching is on position and ALT: ending up with Y at 275 is what matters."""
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("NA", "p.N275Y")])],
                      "gene,mutation\nNA,H275Y\n", segment="NA")
    assert len(hits) == 1


def test_the_observed_reference_residue_reaches_the_output(tmp_path):
    """The sample's own reference is the fact; the watchlist's REF is a claim."""
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("NA", "p.N275Y")])],
                      "gene,mutation\nNA,H275Y\n", segment="NA")
    assert hits[0]["ref_aa"] == "N"
    assert hits[0]["mutation"] == "N275Y"


def test_the_hgvs_p_is_carried_through_to_the_output(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.Glu627Lys")])],
                      "gene,mutation\nPB2,E627K\n")
    assert hits[0]["hgvs_p"] == "p.Glu627Lys"
    assert hits[0]["chrom"] == "CHROM_X"
    assert hits[0]["pos"] == "100"


def test_a_reference_mismatch_warns_naming_both_residues(tmp_path, caplog):
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(100, [ann("NA", "p.N275Y")])])
    path = write_watchlist(tmp_path / "w.csv", "gene,mutation\nNA,H275Y\n")
    entries = watcher.load_watchlist(path)
    _, observed_ref, _ = watcher.index_vcf(vcf, {"NA"})
    with caplog.at_level("WARNING"):
        watcher.warn_ref_mismatches(entries, observed_ref, path)
    assert "reference N at NA:275, not H" in caplog.text


def test_a_reference_mismatch_warns_even_without_a_hit(tmp_path, caplog):
    """The check covers every annotated position, not only those that produced a hit."""
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(100, [ann("NA", "p.N275D")])])
    path = write_watchlist(tmp_path / "w.csv", "gene,mutation\nNA,H275Y\n")
    entries = watcher.load_watchlist(path)
    changes, observed_ref, _ = watcher.index_vcf(vcf, {"NA"})
    assert watcher.find_hits(entries, changes, "S", "NA") == []
    with caplog.at_level("WARNING"):
        watcher.warn_ref_mismatches(entries, observed_ref, path)
    assert "reference N at NA:275, not H" in caplog.text


def test_an_agreeing_reference_residue_is_silent(tmp_path, caplog):
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(100, [ann("NA", "p.H275D")])])
    path = write_watchlist(tmp_path / "w.csv", "gene,mutation\nNA,H275Y\n")
    entries = watcher.load_watchlist(path)
    _, observed_ref, _ = watcher.index_vcf(vcf, {"NA"})
    with caplog.at_level("WARNING"):
        watcher.warn_ref_mismatches(entries, observed_ref, path)
    assert caplog.text == ""


def test_an_unannotated_position_is_not_reported_as_a_mismatch(tmp_path, caplog):
    """A VCF carries variants, not every residue, so silence here is the only answer."""
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(100, [ann("NA", "p.V100I")])])
    path = write_watchlist(tmp_path / "w.csv", "gene,mutation\nNA,H275Y\n")
    entries = watcher.load_watchlist(path)
    _, observed_ref, _ = watcher.index_vcf(vcf, {"NA"})
    with caplog.at_level("WARNING"):
        watcher.warn_ref_mismatches(entries, observed_ref, path)
    assert caplog.text == ""


# ---------------------------------------------------------------------------
# MNV_UNRESOLVED
# ---------------------------------------------------------------------------

def test_an_mnv_unresolved_row_is_skipped(tmp_path):
    """fix_mnv_codons.py has already determined this row's HGVS.p is wrong."""
    hits, _, skipped = scan(
        tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")], extra_info="MNV_UNRESOLVED")],
        "gene,mutation\nPB2,E627K\n",
    )
    assert hits == []
    assert skipped == 1


def test_the_same_row_without_the_flag_is_a_hit(tmp_path):
    hits, _, skipped = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")])],
                            "gene,mutation\nPB2,E627K\n")
    assert len(hits) == 1
    assert skipped == 0


def test_the_flag_is_matched_as_a_token_not_a_substring(tmp_path):
    hits, _, skipped = scan(
        tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")], extra_info="NOT_MNV_UNRESOLVED_X=1")],
        "gene,mutation\nPB2,E627K\n",
    )
    assert len(hits) == 1
    assert skipped == 0


@pytest.mark.parametrize("info,key,expected", [
    ("NS=1;MNV_UNRESOLVED;DP=3", "MNV_UNRESOLVED", True),
    ("NS=1;MNV_UNRESOLVED=1", "MNV_UNRESOLVED", True),
    ("NS=1;MNV_CORRECTED;DP=3", "MNV_UNRESOLVED", False),
    # A substring test would match this; a token test must not.
    ("NS=1;ANN=A|x|y|MNV_UNRESOLVED", "MNV_UNRESOLVED", False),
])
def test_has_info_flag(info, key, expected):
    assert watcher.has_info_flag(info, key) is expected


# ---------------------------------------------------------------------------
# Loading the watchlist
# ---------------------------------------------------------------------------

def test_a_minimal_watchlist_loads(tmp_path):
    entries = watcher.load_watchlist(
        write_watchlist(tmp_path / "w.csv", "gene,mutation\nPB2,E627K\n")
    )
    entry = entries[("PB2", 627, "K")]
    assert (entry.gene, entry.ref_aa, entry.position, entry.alt_aa) == ("PB2", "E", 627, "K")


def test_the_note_column_is_optional(tmp_path):
    entries = watcher.load_watchlist(
        write_watchlist(tmp_path / "w.csv", "gene,mutation\nPB2,E627K\n")
    )
    assert entries[("PB2", 627, "K")].note == ""


def test_gene_names_are_canonicalised_on_load(tmp_path):
    entries = watcher.load_watchlist(
        write_watchlist(tmp_path / "w.csv", "gene,mutation\nNS-2,M51I\n")
    )
    assert list(entries) == [("NEP", 51, "I")]


def test_the_three_letter_spelling_is_accepted_on_load(tmp_path):
    entries = watcher.load_watchlist(
        write_watchlist(tmp_path / "w.csv", "gene,mutation\nPB2,Glu627Lys\n")
    )
    assert list(entries) == [("PB2", 627, "K")]


def test_comment_lines_are_ignored(tmp_path):
    entries = watcher.load_watchlist(write_watchlist(
        tmp_path / "w.csv", "# where this came from\ngene,mutation\n# a note\nPB2,E627K\n"
    ))
    assert list(entries) == [("PB2", 627, "K")]


def test_line_numbers_survive_comment_stripping(tmp_path):
    """An error must point at the line to edit, not a post-filter row index."""
    entries = watcher.load_watchlist(write_watchlist(
        tmp_path / "w.csv", "# one\n# two\ngene,mutation\n\nPB2,E627K\n"
    ))
    assert entries[("PB2", 627, "K")].line == 5


def test_duplicate_entries_collapse_to_one_hit(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")])],
                      "gene,mutation,note\nPB2,E627K,first\nPB2,Glu627Lys,second\n")
    assert len(hits) == 1
    assert hits[0]["note"] == "first"


@pytest.mark.parametrize("text,expected", [
    ("gene,note\nPB2,x\n", "missing required column"),
    ("gene,mutation\nPB2,E627\n", "cannot parse mutation"),
    ("gene,mutation\nPB2,E627E\n", "synonymous"),
    ("gene,mutation\nHA3,E627K\n", "unknown gene"),
    ("gene,mutation\n,E627K\n", "both 'gene' and 'mutation' are required"),
    ("# only a comment\n", "no data rows"),
    ("gene,mutation\n", "no data rows"),
])
def test_a_malformed_watchlist_is_rejected(tmp_path, text, expected):
    with pytest.raises(watcher.WatchlistError, match=expected):
        watcher.load_watchlist(write_watchlist(tmp_path / "w.csv", text))


def test_the_rejection_names_the_line_to_edit(tmp_path):
    with pytest.raises(watcher.WatchlistError, match=r"w\.csv:4:"):
        watcher.load_watchlist(write_watchlist(
            tmp_path / "w.csv", "# a comment\ngene,mutation\nPB2,E627K\nHA3,T195I\n"
        ))


# ---------------------------------------------------------------------------
# Output and CLI
# ---------------------------------------------------------------------------

def test_hits_are_reported_in_watchlist_order(tmp_path):
    """Watchlist order, not VCF order, so the output is predictable from the input list."""
    hits, _, _ = scan(
        tmp_path,
        [vcf_row(100, [ann("NEP", "p.M51I")]), vcf_row(200, [ann("NS1", "p.D92E")])],
        "gene,mutation\nNS1,D92E\nNS-2,M51I\n", segment="NS",
    )
    assert [h["gene"] for h in hits] == ["NS1", "NEP"]


def test_the_output_has_the_documented_columns(tmp_path):
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(100, [ann("PB2", "p.E627K")])])
    watchlist = write_watchlist(tmp_path / "w.csv", "gene,mutation,note\nPB2,E627K,x\n")
    out = tmp_path / "hits.tsv"
    result = run_cli("--vcf", vcf, "--watchlist", watchlist,
                     "--sample", "S1", "--segment", "PB2", "-o", out)
    assert result.returncode == 0, result.stderr
    with out.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert tuple(rows[0]) == watcher.COLUMNS


def test_no_hits_writes_a_header_only_file(tmp_path):
    """'Scanned, found nothing' must not read as 'did not run'."""
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(100, [ann("PB2", "p.A100V")])])
    watchlist = write_watchlist(tmp_path / "w.csv", "gene,mutation\nPB2,E627K\n")
    out = tmp_path / "hits.tsv"
    result = run_cli("--vcf", vcf, "--watchlist", watchlist,
                     "--sample", "S1", "--segment", "PB2", "-o", out)
    assert result.returncode == 0, result.stderr
    assert out.read_text().splitlines() == ["\t".join(watcher.COLUMNS)]


def test_an_empty_vcf_is_not_an_error(tmp_path):
    hits, _, skipped = scan(tmp_path, [], "gene,mutation\nPB2,E627K\n")
    assert hits == []
    assert skipped == 0


def test_a_malformed_watchlist_exits_nonzero(tmp_path):
    """Never a silent empty scan: that is indistinguishable from a clean sample."""
    vcf = write_vcf(tmp_path / "in.vcf", [vcf_row(100, [ann("PB2", "p.E627K")])])
    watchlist = write_watchlist(tmp_path / "w.csv", "gene,note\nPB2,x\n")
    result = run_cli("--vcf", vcf, "--watchlist", watchlist,
                     "--sample", "S1", "-o", tmp_path / "hits.tsv")
    assert result.returncode == 1
    assert "missing required column" in result.stderr


def test_the_cli_reports_a_version():
    """The provenance block shells out to this, so it must not drift."""
    result = run_cli("--version")
    assert result.returncode == 0
    assert watcher.SCRIPT_VERSION in result.stdout


def test_the_watcher_uses_the_shared_parsers():
    """A second ANN parser is how the two consumers of a VCF quietly disagree."""
    from lib import aa_markers

    assert watcher.annotated_changes is aa_markers.annotated_changes
    assert watcher.parse_hgvs_p is aa_markers.parse_hgvs_p
    assert watcher.has_info_flag is aa_markers.has_info_flag
