"""Unit tests for bin/mutation_watcher.py (custom watchlist scanning of a final VCF).

Several cases here are regressions against the previous TSV-reading implementation, which
had no gene canonicalisation, no mature-peptide handling and no MNV_UNRESOLVED skip. Those
are marked in their docstrings, because each one failed silently: the tool reported zero
hits, which is also what a clean sample reports.

Run with:  pytest tests/test_mutation_watcher.py
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
    spec.loader.exec_module(module)
    return module


watcher = _load_module()


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def ann(gene: str, hgvs_p: str, feature_id: str | None = None,
        annotation: str = "missense_variant") -> str:
    """A 16-field snpEff-style ANN entry, matching what the converter emits.

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
    """Run the scan and return (hit dicts, coverage rows, mnv_unresolved count)."""
    from validate_watchlist import load_watchlist

    vcf = write_vcf(tmp_path / "in.vcf", rows)
    entries = load_watchlist(write_watchlist(tmp_path / "w.csv", watchlist))
    hits, genes_seen, mnv = watcher.scan_vcf(vcf, entries, sample, segment)
    coverage = watcher.report_coverage(entries, genes_seen, hits)
    return hits, coverage, mnv


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------

def test_a_watched_substitution_is_a_hit(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")])],
                      "gene,mutation,note\nPB2,E627K,mammalian\n")
    assert len(hits) == 1
    assert hits[0]["gene"] == "PB2"
    assert hits[0]["position"] == 627
    assert hits[0]["alt_aa"] == "K"
    assert hits[0]["note"] == "mammalian"


def test_the_sample_id_is_reported_verbatim(tmp_path):
    """It used to be args.input.split('_')[0], which truncated at the first underscore."""
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")])],
                      "gene,mutation\nPB2,E627K\n", sample="WIN_AH_2024_FAV_0415")
    assert hits[0]["sample_id"] == "WIN_AH_2024_FAV_0415"


def test_the_three_letter_spelling_gives_the_identical_hit(tmp_path):
    """snpEff's default output. Matching only the 1-letter form yields zero hits and no
    warning -- indistinguishable from a VCF with no coding changes."""
    one, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")])],
                     "gene,mutation\nPB2,E627K\n")
    three, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.Glu627Lys")])],
                       "gene,mutation\nPB2,E627K\n")
    assert len(three) == 1
    assert three[0]["position"] == one[0]["position"]
    assert three[0]["alt_aa"] == one[0]["alt_aa"]


def test_a_gene_alias_matches(tmp_path):
    """Fails on the previous implementation: it compared gene names verbatim, so a list
    saying NS-1 never matched an NS1 annotation. ISSUES.md I-2 in a second location."""
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("NS1", "p.A53D")])],
                      "gene,mutation\nNS-1,A53D\n", segment="NS")
    assert len(hits) == 1
    assert hits[0]["gene"] == "NS1"


def test_the_reverse_alias_direction_also_matches(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("NS-1", "p.A53D")])],
                      "gene,mutation\nNS1,A53D\n", segment="NS")
    assert len(hits) == 1


def test_a_wrong_alt_residue_is_not_a_hit(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E627R")])],
                      "gene,mutation\nPB2,E627K\n")
    assert hits == []


def test_a_wrong_position_is_not_a_hit(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E628K")])],
                      "gene,mutation\nPB2,E627K\n")
    assert hits == []


def test_a_synonymous_change_is_not_a_hit(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E627E")])],
                      "gene,mutation\nPB2,E627K\n")
    assert hits == []


@pytest.mark.parametrize("hgvs_p", ["p.Ala100fs", "p.Lys50del", "p.Glu627dup",
                                    "p.Met1?", "p.Glu627=", "p.Ter100ext"])
def test_non_substitutions_are_never_hits(tmp_path, hgvs_p):
    """Replaces the old `"fs" in mut or "*" in mut` string test."""
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", hgvs_p)])],
                      "gene,mutation\nPB2,E627K\n")
    assert hits == []


def test_an_intergenic_pseudo_gene_is_ignored(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("CHR_START-PB2", "p.E627K")])],
                      "gene,mutation\nPB2,E627K\n")
    assert hits == []


# ---------------------------------------------------------------------------
# The REF residue: reported, not matched on
# ---------------------------------------------------------------------------

def test_a_differing_reference_residue_is_still_a_hit(tmp_path, caplog):
    """What makes a sample resistant is ending up with K at 627, not what it started from.
    This mirrors the merge, whose MARKER_QUERY omits ref_seq for the same reason."""
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.D627K")])],
                      "gene,mutation\nPB2,E627K\n")
    assert len(hits) == 1
    assert hits[0]["ref_aa"] == "D"
    assert hits[0]["mutation"] == "D627K"
    assert "reference residue" in caplog.text


def test_the_observed_reference_residue_wins_over_the_watchlist(tmp_path):
    """The row is the fact; the watchlist's REF is the claim."""
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("NA", "p.N275Y")])],
                      "gene,mutation\nNA,H275Y\n", segment="NA")
    assert hits[0]["ref_aa"] == "N"


def test_the_hgvs_p_is_carried_through_to_the_output(tmp_path):
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")])],
                      "gene,mutation\nPB2,E627K\n")
    assert hits[0]["hgvs_p"] == "p.E627K"
    assert hits[0]["chrom"] == "CHROM_X"
    assert hits[0]["pos"] == "100"


# ---------------------------------------------------------------------------
# Mature peptides and the HA trap
# ---------------------------------------------------------------------------

def test_a_mature_peptide_resolves_through_feature_id(tmp_path):
    """NCBI-style HA: Gene_Name is the parent (HA), only Feature_ID names the subunit.
    Fails on the previous implementation, which only ever saw the gene name."""
    hits, _, _ = scan(tmp_path, [vcf_row(100, [ann("HA", "p.T195I", feature_id="HA1")])],
                      "gene,mutation\nHA1,T195I\n", segment="HA")
    assert len(hits) == 1
    assert hits[0]["gene"] == "HA1"


def test_a_parent_gene_never_matches_a_subunit_watchlist_entry(tmp_path):
    """The HA1 numbering trap. HA1 starts after HA's 16-residue signal peptide, so keying
    an HA1 marker on a bare HA CDS returns the wrong residue at 93% of positions. Reporting
    nothing is correct here; reporting a hit would be a wrong answer, silently."""
    hits, coverage, _ = scan(tmp_path, [vcf_row(100, [ann("HA", "p.T195I")])],
                             "gene,mutation\nHA1,T195I\n", segment="HA")
    assert hits == []
    assert [r for r in coverage if r["gene"] == "HA1"][0]["in_vcf"] == "no"


def test_a_subunit_annotation_never_matches_a_parent_watchlist_entry(tmp_path):
    """The same trap in the other direction, which is how 88 rows of the previous list
    silently matched nothing on the production datasets."""
    hits, coverage, _ = scan(tmp_path, [vcf_row(100, [ann("HA1", "p.T195I")])],
                             "gene,mutation\nHA,T195I\n", segment="HA")
    assert hits == []
    assert [r for r in coverage if r["gene"] == "HA"][0]["in_vcf"] == "no"


def test_a_row_carrying_both_the_parent_and_the_subunit_hits_once(tmp_path):
    """Every in-HA1 row on an NCBI-style dataset carries two AA-changing ANN entries."""
    rows = [vcf_row(100, [ann("HA", "p.T211I"), ann("HA", "p.T195I", feature_id="HA1")])]
    hits, _, _ = scan(tmp_path, rows, "gene,mutation\nHA1,T195I\n", segment="HA")
    assert len(hits) == 1


# ---------------------------------------------------------------------------
# MNV_UNRESOLVED
# ---------------------------------------------------------------------------

def test_an_mnv_unresolved_row_is_skipped(tmp_path):
    """fix_mnv_codons.py sets this on rows whose HGVS.p it has determined to be wrong. A
    hit there would be a false positive. Fails on the previous implementation."""
    hits, _, mnv = scan(tmp_path,
                        [vcf_row(100, [ann("PB2", "p.E627K")], extra_info="MNV_UNRESOLVED")],
                        "gene,mutation\nPB2,E627K\n")
    assert hits == []
    assert mnv == 1


def test_the_same_row_without_the_flag_is_a_hit(tmp_path):
    hits, _, mnv = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")])],
                        "gene,mutation\nPB2,E627K\n")
    assert len(hits) == 1
    assert mnv == 0


def test_the_flag_is_matched_as_a_token_not_a_substring(tmp_path):
    """A key must never be matched inside another key's name or value."""
    hits, _, mnv = scan(
        tmp_path,
        [vcf_row(100, [ann("PB2", "p.E627K")], extra_info="NOTE=NOT_MNV_UNRESOLVED_REALLY")],
        "gene,mutation\nPB2,E627K\n",
    )
    assert len(hits) == 1
    assert mnv == 0


# ---------------------------------------------------------------------------
# Coverage reporting -- the safeguard
# ---------------------------------------------------------------------------

def test_a_gene_absent_from_the_vcf_is_reported_as_unmatched(tmp_path):
    _, coverage, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.S715N")])],
                          "gene,mutation\nPB2,E627K\nNA,H275Y\n")
    by_gene = {r["gene"]: r for r in coverage}
    assert by_gene["PB2"]["in_vcf"] == "yes"
    assert by_gene["NA"]["in_vcf"] == "no"


def test_a_gene_present_but_unmutated_is_distinct_from_a_gene_absent(tmp_path):
    """Zero hits is what a clean sample looks like. The two must not print the same."""
    _, coverage, _ = scan(tmp_path, [vcf_row(100, [ann("PB2", "p.S715N")])],
                          "gene,mutation\nPB2,E627K\nNA,H275Y\n")
    by_gene = {r["gene"]: r for r in coverage}
    assert by_gene["PB2"]["hits"] == 0 and by_gene["PB2"]["in_vcf"] == "yes"
    assert by_gene["NA"]["hits"] == 0 and by_gene["NA"]["in_vcf"] == "no"


def test_the_unmatched_gene_warning_suggests_the_subunit(tmp_path, caplog):
    """The HA -> HA1 case, named outright: it is the one that has actually bitten."""
    watcher.warn_unmatched_genes(
        [{"gene": "HA", "in_vcf": "no", "watchlist_entries": 88, "hits": 0}],
        {"HA1", "HA2"},
    )
    assert "HA1" in caplog.text
    assert "Did you mean" in caplog.text


def test_coverage_counts_entries_and_hits_per_gene(tmp_path):
    rows = [vcf_row(100, [ann("PB2", "p.E627K")]), vcf_row(200, [ann("PB2", "p.S715N")])]
    _, coverage, _ = scan(tmp_path, rows, "gene,mutation\nPB2,E627K\nPB2,S715N\nPB2,A274T\n")
    assert coverage[0]["watchlist_entries"] == 3
    assert coverage[0]["hits"] == 2


# ---------------------------------------------------------------------------
# Output files
# ---------------------------------------------------------------------------

def run_cli(tmp_path, rows, watchlist, *extra, sample="SAMPLE-1", segment="PB2"):
    vcf = write_vcf(tmp_path / "in.vcf", rows)
    wl = write_watchlist(tmp_path / "w.csv", watchlist)
    out = tmp_path / "out_watchlist_hits.tsv"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--vcf", str(vcf), "--watchlist", str(wl),
         "--sample", sample, "--segment", segment, "-o", str(out), *extra],
        capture_output=True, text=True,
    )
    return proc, out


def test_the_output_has_the_documented_columns(tmp_path):
    proc, out = run_cli(tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")])],
                        "gene,mutation,note\nPB2,E627K,mammalian\n")
    assert proc.returncode == 0
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert list(rows[0]) == list(watcher.OUTPUT_COLUMNS)
    assert rows[0]["sample_id"] == "SAMPLE-1"
    assert rows[0]["segment"] == "PB2"
    assert rows[0]["note"] == "mammalian"


def test_no_hits_writes_a_header_only_file(tmp_path):
    """So 'scanned, found nothing' is distinguishable from 'did not run'."""
    proc, out = run_cli(tmp_path, [vcf_row(100, [ann("PB2", "p.S715N")])],
                        "gene,mutation\nPB2,E627K\n")
    assert proc.returncode == 0
    lines = out.read_text().splitlines()
    assert len(lines) == 1
    assert lines[0].split("\t") == list(watcher.OUTPUT_COLUMNS)


def test_the_coverage_sidecar_is_written_alongside(tmp_path):
    proc, out = run_cli(tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")])],
                        "gene,mutation\nPB2,E627K\nNA,H275Y\n")
    sidecar = out.with_name(out.name.replace("_hits.tsv", "_coverage.tsv"))
    assert sidecar.exists()
    rows = {r["gene"]: r for r in csv.DictReader(sidecar.open(), delimiter="\t")}
    assert rows["NA"]["in_vcf"] == "no"


def test_a_malformed_watchlist_exits_nonzero(tmp_path):
    """It used to `return dict()` and report no mutations of concern."""
    proc, _ = run_cli(tmp_path, [vcf_row(100, [ann("PB2", "p.E627K")])],
                      "segment,mutation\nHA,T195I\n")
    assert proc.returncode == 1
    assert "gene" in proc.stderr


def test_an_empty_vcf_is_not_an_error(tmp_path):
    proc, out = run_cli(tmp_path, [], "gene,mutation\nPB2,E627K\n")
    assert proc.returncode == 0
    assert out.read_text().splitlines() == ["\t".join(watcher.OUTPUT_COLUMNS)]


def test_the_cli_reports_a_version():
    """The provenance block shells out to this rather than repeating the literal."""
    proc = subprocess.run([sys.executable, str(SCRIPT_PATH), "--version"],
                          capture_output=True, text=True)
    assert proc.returncode == 0
    assert watcher.SCRIPT_VERSION in proc.stdout


# ---------------------------------------------------------------------------
# Shared-module identity
# ---------------------------------------------------------------------------

def test_the_watcher_uses_the_shared_parsers():
    """A second copy of these would be free to drift over which HGVS.p spellings parse,
    and the disagreement would be silent. Mirrors the merge's equivalent assertion."""
    from lib import aa_markers, flu_gene_names

    assert watcher.aa_changes_from_info is aa_markers.aa_changes_from_info
    assert watcher.has_info_flag is aa_markers.has_info_flag
    assert watcher.parse_hgvs_p is aa_markers.parse_hgvs_p
    assert watcher.canonical_gene is flu_gene_names.canonical_gene
