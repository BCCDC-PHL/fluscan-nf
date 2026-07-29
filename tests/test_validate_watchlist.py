"""Unit tests for bin/validate_watchlist.py.

The watchlist is hand-maintained, so almost every assertion here is about a defect being
*reported* rather than absorbed. The implementation this replaces caught every load error
with `except Exception: return dict()`, which turned a malformed file into "no mutations
of concern" -- indistinguishable from a clean sample.

Run with:  pytest tests/test_validate_watchlist.py
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "bin" / "validate_watchlist.py"


def _load_module():
    sys.path.insert(0, str(REPO_ROOT / "bin"))
    spec = importlib.util.spec_from_file_location("validate_watchlist", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


vw = _load_module()


def write_list(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def run_cli(path: Path, *extra: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), str(path), *extra],
        capture_output=True, text=True,
    )


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def test_a_minimal_watchlist_loads(tmp_path):
    path = write_list(tmp_path / "w.csv", "gene,mutation,note\nPB2,E627K,mammalian\n")
    entries = vw.load_watchlist(path)
    assert len(entries) == 1
    assert entries[0].gene == "PB2"
    assert (entries[0].ref_aa, entries[0].position, entries[0].alt_aa) == ("E", 627, "K")
    assert entries[0].note == "mammalian"


def test_the_note_column_is_optional(tmp_path):
    path = write_list(tmp_path / "w.csv", "gene,mutation\nPB2,E627K\n")
    assert vw.load_watchlist(path)[0].note == ""


def test_three_letter_spelling_is_accepted(tmp_path):
    """The same substitution either way -- snpEff's default output is 3-letter."""
    one = vw.load_watchlist(write_list(tmp_path / "a.csv", "gene,mutation\nPB2,E627K\n"))[0]
    three = vw.load_watchlist(write_list(tmp_path / "b.csv", "gene,mutation\nPB2,Glu627Lys\n"))[0]
    assert one.key == three.key


def test_gene_names_are_canonicalised(tmp_path):
    """NS-1/NS1 and NS-2/NEP name the same protein and must key the same."""
    assert vw.load_watchlist(
        write_list(tmp_path / "a.csv", "gene,mutation\nNS-1,A53D\n"))[0].gene == "NS1"
    assert vw.load_watchlist(
        write_list(tmp_path / "b.csv", "gene,mutation\nNS-2,T48A\n"))[0].gene == "NEP"
    assert vw.load_watchlist(
        write_list(tmp_path / "c.csv", "gene,mutation\nNS2,T48A\n"))[0].gene == "NEP"


def test_comment_lines_are_ignored(tmp_path):
    path = write_list(
        tmp_path / "w.csv",
        "# where this came from\n#\ngene,mutation\nPB2,E627K\n\nPB2,S715N\n",
    )
    entries = vw.load_watchlist(path)
    assert [e.mutation for e in entries] == ["E627K", "S715N"]


def test_line_numbers_survive_comment_stripping(tmp_path):
    """Messages must point at the line the user has to edit, not a post-filter index."""
    path = write_list(tmp_path / "w.csv", "# a\n# b\ngene,mutation\nPB2,E627K\nPB2,S715N\n")
    assert [e.line for e in vw.load_watchlist(path)] == [4, 5]


# ---------------------------------------------------------------------------
# Hard failures
# ---------------------------------------------------------------------------

def test_a_segment_column_is_rejected_and_names_the_rename(tmp_path):
    """The old schema is refused rather than read, because the values differ too.

    'segment' held HA where the matching key is HA1: reading one as the other is what made
    88 rows of the previous list unmatchable while looking perfectly valid.
    """
    path = write_list(tmp_path / "w.csv", "segment,mutation\nHA,T195I\n")
    with pytest.raises(vw.WatchlistError) as excinfo:
        vw.load_watchlist(path)
    message = str(excinfo.value)
    assert "gene" in message and "segment" in message
    assert "HA1" in message


def test_a_missing_required_column_is_rejected(tmp_path):
    path = write_list(tmp_path / "w.csv", "gene,note\nPB2,x\n")
    with pytest.raises(vw.WatchlistError, match="mutation"):
        vw.load_watchlist(path)


def test_an_unparseable_mutation_is_rejected(tmp_path):
    path = write_list(tmp_path / "w.csv", "gene,mutation\nNA,247-250del\n")
    with pytest.raises(vw.WatchlistError, match="247-250del"):
        vw.load_watchlist(path)


def test_a_non_ascii_mutation_is_rejected_naming_the_character(tmp_path):
    """An en dash renders almost identically to a hyphen, so the message must say which."""
    path = write_list(tmp_path / "w.csv", "gene,mutation\nNA,245–248del\n")
    with pytest.raises(vw.WatchlistError) as excinfo:
        vw.load_watchlist(path)
    assert "U+2013" in str(excinfo.value)


def test_a_synonymous_entry_is_rejected(tmp_path):
    """It could never match: aa_changes_from_info drops synonymous changes."""
    path = write_list(tmp_path / "w.csv", "gene,mutation\nPB2,E627E\n")
    with pytest.raises(vw.WatchlistError, match="synonymous"):
        vw.load_watchlist(path)


def test_an_unknown_gene_is_rejected(tmp_path):
    path = write_list(tmp_path / "w.csv", "gene,mutation\nPB3,E627K\n")
    with pytest.raises(vw.WatchlistError, match="PB3"):
        vw.load_watchlist(path)


def test_ha_is_a_known_gene(tmp_path):
    """HA is not a typo for HA1. Datasets with a bare HA CDS exist and their numbering is
    genuinely HA's; whether a list should say HA or HA1 is the watcher's coverage report to
    answer, not this one's."""
    assert vw.load_watchlist(write_list(tmp_path / "w.csv", "gene,mutation\nHA,T195I\n"))


def test_a_missing_value_is_rejected(tmp_path):
    path = write_list(tmp_path / "w.csv", "gene,mutation\nPB2,\n")
    with pytest.raises(vw.WatchlistError, match="required"):
        vw.load_watchlist(path)


def test_an_empty_file_is_rejected(tmp_path):
    path = write_list(tmp_path / "w.csv", "# only comments\n")
    with pytest.raises(vw.WatchlistError):
        vw.load_watchlist(path)


# ---------------------------------------------------------------------------
# Duplicates and conflicts
# ---------------------------------------------------------------------------

def test_exact_duplicates_are_dropped_with_a_warning(tmp_path, caplog):
    path = write_list(tmp_path / "w.csv", "gene,mutation\nPB2,E627K\nPB2,E627K\n")
    entries = vw.load_watchlist(path)
    assert len(entries) == 1
    assert "duplicate" in caplog.text


def test_duplicates_are_detected_across_spellings(tmp_path):
    """Dedup is on the parsed substitution, so E627K and Glu627Lys are one entry."""
    path = write_list(tmp_path / "w.csv", "gene,mutation\nPB2,E627K\nPB2,Glu627Lys\n")
    assert len(vw.load_watchlist(path)) == 1


def test_two_reference_residues_at_one_position_are_reported(tmp_path):
    """A position has one reference residue. Two means two numbering conventions."""
    path = write_list(tmp_path / "w.csv", "gene,mutation\nNA,H275Y\nNA,N275S\n")
    messages = vw.report_conflicts(vw.load_watchlist(path))
    assert any("275" in m and "more than one reference residue" in m for m in messages)


def test_the_same_substitution_one_position_apart_is_reported(tmp_path):
    """I222L and I223L are one marker under N2 and N1 numbering, not two markers."""
    path = write_list(tmp_path / "w.csv", "gene,mutation\nNA,I222L\nNA,I223L\n")
    messages = vw.report_conflicts(vw.load_watchlist(path))
    assert any("I222L" in m and "I223L" in m for m in messages)


def test_conflicts_are_reported_not_resolved(tmp_path):
    """Both entries survive. Which convention is wanted depends on where each came from,
    and a tool that guessed would silently delete a marker somebody added on purpose."""
    path = write_list(tmp_path / "w.csv", "gene,mutation\nNA,I222L\nNA,I223L\n")
    entries = vw.load_watchlist(path)
    assert {e.mutation for e in entries} == {"I222L", "I223L"}


def test_unrelated_neighbouring_positions_are_not_reported(tmp_path):
    """Only the *same* substitution one apart looks like a convention clash."""
    path = write_list(tmp_path / "w.csv", "gene,mutation\nNA,G147R\nNA,T148I\n")
    assert vw.report_conflicts(vw.load_watchlist(path)) == []


# ---------------------------------------------------------------------------
# Reference-residue validation
# ---------------------------------------------------------------------------

def write_aa_tsv(path: Path, rows: list[tuple[str, int, str]]) -> Path:
    lines = ["source_reference\tsegment\tgene\tsource_aa_pos\tsource_residue\t"
             "target_aa_pos\ttarget_residue\tstatus"]
    for gene, pos, residue in rows:
        lines.append(f"SRC\tNA\t{gene}\t{pos}\t{residue}\t{pos}\t{residue}\taligned")
    path.write_text("\n".join(lines) + "\n")
    return path


def test_a_matching_reference_residue_is_silent(tmp_path):
    entries = vw.load_watchlist(write_list(tmp_path / "w.csv", "gene,mutation\nNA,H275Y\n"))
    residues = vw.load_reference_residues(write_aa_tsv(tmp_path / "m.tsv", [("NA", 275, "H")]))
    assert vw.check_reference_residues(entries, residues) == []


def test_a_mismatching_reference_residue_is_reported(tmp_path):
    """This is what catches a wrong numbering convention at load time rather than after
    it has quietly mis-scanned every sample."""
    entries = vw.load_watchlist(write_list(tmp_path / "w.csv", "gene,mutation\nNA,R292K\n"))
    residues = vw.load_reference_residues(write_aa_tsv(tmp_path / "m.tsv", [("NA", 292, "C")]))
    messages = vw.check_reference_residues(entries, residues)
    assert len(messages) == 1
    assert "R292K" in messages[0] and "'C'" in messages[0]


def test_a_position_outside_the_reference_is_reported(tmp_path):
    entries = vw.load_watchlist(write_list(tmp_path / "w.csv", "gene,mutation\nNA,H999Y\n"))
    residues = vw.load_reference_residues(write_aa_tsv(tmp_path / "m.tsv", [("NA", 275, "H")]))
    assert "not present" in vw.check_reference_residues(entries, residues)[0]


def test_a_reference_mismatch_does_not_remove_the_entry(tmp_path):
    """Matching is on position and ALT. A reference that differs is a reason to warn, not
    a reason to stop looking for the residue."""
    entries = vw.load_watchlist(write_list(tmp_path / "w.csv", "gene,mutation\nNA,R292K\n"))
    residues = vw.load_reference_residues(write_aa_tsv(tmp_path / "m.tsv", [("NA", 292, "C")]))
    vw.check_reference_residues(entries, residues)
    assert entries[0].key == ("NA", 292, "K")


def test_reference_residue_lookup_canonicalises_gene_names(tmp_path):
    """The map may spell a gene NS-1 where the list says NS1; both key as NS1."""
    entries = vw.load_watchlist(write_list(tmp_path / "w.csv", "gene,mutation\nNS1,A53D\n"))
    residues = vw.load_reference_residues(write_aa_tsv(tmp_path / "m.tsv", [("NS-1", 53, "A")]))
    assert vw.check_reference_residues(entries, residues) == []


def test_a_file_that_is_not_an_aa_map_is_rejected(tmp_path):
    path = tmp_path / "not_a_map.tsv"
    path.write_text("a\tb\n1\t2\n")
    with pytest.raises(vw.WatchlistError, match="coordinate_map.aa.tsv"):
        vw.load_reference_residues(path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def test_cli_exits_zero_on_a_valid_list(tmp_path):
    path = write_list(tmp_path / "w.csv", "gene,mutation\nPB2,E627K\n")
    assert run_cli(path).returncode == 0


def test_cli_exits_nonzero_on_an_invalid_list(tmp_path):
    """The whole point: this used to be `return dict()` and a clean-looking scan."""
    path = write_list(tmp_path / "w.csv", "segment,mutation\nHA,T195I\n")
    proc = run_cli(path)
    assert proc.returncode == 1
    assert "gene" in proc.stderr


def test_cli_warnings_alone_do_not_fail(tmp_path):
    path = write_list(tmp_path / "w.csv", "gene,mutation\nNA,I222L\nNA,I223L\n")
    proc = run_cli(path)
    assert proc.returncode == 0
    assert "numbering conventions" in proc.stderr


def test_cli_strict_fails_on_warnings(tmp_path):
    path = write_list(tmp_path / "w.csv", "gene,mutation\nNA,I222L\nNA,I223L\n")
    assert run_cli(path, "--strict").returncode == 1


def test_cli_reports_a_version(tmp_path):
    """The provenance block shells out to this rather than repeating the literal."""
    proc = subprocess.run([sys.executable, str(SCRIPT_PATH), "--version"],
                          capture_output=True, text=True)
    assert proc.returncode == 0
    assert vw.SCRIPT_VERSION in proc.stdout


# ---------------------------------------------------------------------------
# The shipped watchlists
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "path", sorted((REPO_ROOT / "assets" / "watchlists").glob("*watchlist.csv")),
    ids=lambda p: p.name,
)
def test_the_shipped_watchlists_are_valid(path):
    """--watchlist_dir defaults here, so a broken file would break every default run."""
    assert vw.load_watchlist(path)
