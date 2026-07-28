#!/usr/bin/env python3
"""Annotate an existing fluscan/Nextclade VCF with FluMutDB marker matches.

Version 2 of the flumut merge. Where v1 (merge_flumut_mutations.py) consumed flumut's
per-sample melted TSV and *inserted* synthetic rows for markers fluscan had not called,
this version never adds, removes or reorders a row. It walks the VCF, and for every amino
-acid change already annotated in a row's ANN= field it queries FluMutDB for a curated
marker matching that (protein, AA position, alt AA). Matching rows gain one extra INFO
field:

    FLUMUT=<Protein>|<AA_Pos>|<Alt_AA>|<Mutation_Name>|<DirectEffect>|<AssociatedEffect>

one comma-separated entry per matching ANN annotation. ANN= itself is left untouched, so
it stays strictly snpEff 16-field compatible.

A FluMutDB *marker* may require several mutations, and this script matches one mutation at
a time, so the two effect fields say how far each published effect can be attributed to the
residue actually observed:

  DirectEffect      effects of single-mutation markers -- attributable to this mutation.
  AssociatedEffect  effects of markers that need more than one mutation. The complete
                    marker is NOT asserted to be present; this residue is one of its
                    constituents. 118 of FluMutDB 6.6's 333 mutations have no
                    single-mutation marker at all, so dropping these would lose their only
                    published evidence -- they are labelled, not hidden.

Both are 'Subtype:Effect' pairs, so an effect can be traced to the experimental system it
came from; without that, effects reported for one mutation across several subtypes and
papers read as self-contradictory (PB2:E627K carries both decreased and increased avian
polymerase activity). '.' means the class is empty.

The marker set comes from the database rather than from a flumut run, so no per-sample
flumut output is needed. FluMutDB positions live in flumut's own reference coordinates and
are translated onto the Nextclade dataset via the coordinate_conversion map, the same
aa_source_to_target block v1 uses (applied in reverse here, since we start from the VCF).

FluMutDB annotates HA as HA1/HA2 subunits, so HA markers only match when the Nextclade
dataset exposes that split too -- as two CDS features, or as mature_protein_region_of_CDS
children of one HA CDS (the NCBI/RefSeq model). A mature peptide's ANN entry names the
parent gene, so gene resolution falls back to ANN Feature_ID against the aliases the
coordinate map records. Datasets with a bare HA CDS and no subunit features produce a
warning and pass through unannotated.

Rows flagged MNV_UNRESOLVED by fix_mnv_codons.py are skipped: their HGVS.p is one the
pipeline has already determined to be wrong, so a marker match there would be a false
positive.

This script does not modify any existing file; it writes a new VCF.

Pipeline position: after post_process_vcf_batch (nextclade workflow), consuming
coordinate_conversion output.
"""

import argparse
import json
import logging
import os
import re
import sqlite3
import sys
from functools import lru_cache

# bin/ holds entry points and is on PATH for every task; the modules they import live in
# bin/lib/. Python already puts a script's own directory on sys.path, but this insert uses
# realpath so the import still resolves when bin/ is reached through a symlink.
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from lib.flu_gene_names import canonical_gene  # noqa: E402
from lib.tools import read_vcf, write_vcf  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

SCRIPT_VERSION = "2.2.0"

FLUMUT_HEADER = (
    '##INFO=<ID=FLUMUT,Number=.,Type=String,'
    'Description="FluMutDB marker matches: '
    "'Protein | AA_Pos | Alt_AA | Mutation_Name | DirectEffect | AssociatedEffect'. "
    "DirectEffect lists Subtype:Effect pairs from markers this mutation alone defines; "
    "AssociatedEffect lists them for markers needing further mutations, which are not "
    "asserted to be present. '.' means the class is empty. "
    "'&' separates multiple values within a field; ',' separates one entry per "
    'matching ANN annotation. AA_Pos is in this VCF\'s coordinates.">'
)

# FluMutDB names the M segment reference 'MP'; every other segment label matches verbatim.
SEGMENT_TO_DB_REFERENCE = {"M": "MP"}

# FluMutDB protein names differ from some Nextclade annotations (NS-1/NS-2 vs NS1/NEP).
# Rather than guess, gene resolution generates candidate names and keeps whichever exists
# as a gene key in the coordinate map for the segment being processed. Carried over from
# merge_flumut_mutations.py (v1) rather than imported, so v1 stays independently removable.
GENE_ALIASES = {
    "NS2": ["NEP"],
    "NEP": ["NS2"],
    "M": ["MP"],
    "MP": ["M"],
}

# Two HGVS.p substitution spellings are accepted. bin/nextclade_to_freebayes_vcf.py emits the
# 1-letter form (p.E627K); snpEff -- which annotates the short-read and paftools branches this
# process is meant to be reusable from (modules/nextclade/merge_flumut_mutations.nf's
# task.ext.workflow) -- defaults to the 3-letter form (p.Glu627Lys). Matching only one of them
# fails silently: zero annotations, indistinguishable from a VCF with no coding changes.
HGVS_P_RE = re.compile(r"p\.([A-Za-z*])(\d+)([A-Za-z*])")
HGVS_P_RE_3LETTER = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")

# Local rather than Bio.SeqUtils.seq1: this script deliberately avoids a Biopython dependency
# (see bin/lib/flu_gene_names.py for the same reasoning), and a 23-entry table does not justify one.
AA3_TO_AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*",
}


# Forms that legitimately are not substitutions: indel/frameshift notations and snpEff's
# 'p.Glu627=' synonymous shorthand. FluMutDB 6.6 holds only SNP-type mutations, so skipping
# these loses nothing curated -- but they must not be reported as parse failures either.
HGVS_P_NON_SUB_RE = re.compile(r"del|ins|dup|fs|ext|=|\?")


def parse_hgvs_p(hgvs_p: str) -> tuple[str, int, str] | None:
    """(ref_aa, position, alt_aa) for an HGVS.p substitution, or None if unparseable.

    Non-substitution forms (del/ins/fs) return None. FluMutDB 6.6 holds only SNP-type
    mutations, so nothing curated is lost by ignoring them.
    """
    hit = HGVS_P_RE.search(hgvs_p)
    if hit:
        return hit.group(1), int(hit.group(2)), hit.group(3)
    hit = HGVS_P_RE_3LETTER.search(hgvs_p)
    if hit:
        ref_aa = AA3_TO_AA1.get(hit.group(1))
        alt_aa = AA3_TO_AA1.get(hit.group(3))
        if ref_aa and alt_aa:
            return ref_aa, int(hit.group(2)), alt_aa
    return None

# One parameterised lookup per (protein, position, alt AA). ref_seq is deliberately absent
# from the WHERE clause: a marker is defined by its alt allele, and the sample's reference
# amino acid may legitimately differ from FluMutDB's.
#
# One row per marker/effect association rather than a group_concat, so each effect keeps the
# marker's mutation count and the paper's subtype and can be classified and sanitised
# individually. Every join stays LEFT: detecting a mutation must not depend on it having any
# published effect (3 of FluMutDB 6.6's mutations have none at all).
MARKER_QUERY = """
SELECT mm.mutation_name, s.all_mutations_count, me.subtype, me.effect_name
FROM mutation_mappings mm
JOIN mutations m               ON m.name = mm.mutation_name
LEFT JOIN markers_mutations mk ON mk.mutation_name = mm.mutation_name
LEFT JOIN markers_summary s    ON s.marker_id = mk.marker_id
LEFT JOIN markers_effects me   ON me.marker_id = mk.marker_id
WHERE mm.reference_name = ? AND m.protein_name = ? AND mm.position = ? AND mm.alt_seq = ?
ORDER BY mm.mutation_name
"""

PROTEIN_QUERY = """
SELECT DISTINCT m.protein_name
FROM mutations m
JOIN mutation_mappings mm ON mm.mutation_name = m.name
WHERE mm.reference_name = ?
ORDER BY m.protein_name
"""


# ---------------------------------------------------------------------------
# Gene-name resolution
# ---------------------------------------------------------------------------

def resolve_gene(raw: str, available_genes: set) -> str | None:
    """Resolve a FluMutDB protein name to a coordinate-map gene key for this segment.

    Returns the matching key, or None if the protein does not belong to this segment
    (or cannot be reconciled with the annotation).

    Candidates are tried most-specific first, and the bare stem last. That ordering is
    load-bearing: a target dataset naming its NS1 gene 'NS' (common in GenBank-derived
    annotations) offers the stem as a match for *both* NS-1 and NS-2, and NS1 (230 aa) and
    NEP (121 aa) are unrelated in numbering -- resolving NS-2 to 'NS' would query FluMutDB
    at NEP positions through NS1's coordinate map, i.e. produce false positives rather than
    silence.
    """
    specific: list[str] = [raw]
    stem: list[str] = []
    m = re.match(r"^(.*?)-(\d+)$", raw)
    if m:
        base, suffix = m.group(1), m.group(2)
        specific.append(base + suffix)     # NS-1 -> NS1 ; NS-2 -> NS2
        stem.append(base)                  # HA1-5 -> HA1 ; NS-1 -> NS
    candidates: list[str] = []
    for group in (specific, stem):
        for c in group:
            for candidate in (c, *GENE_ALIASES.get(c, [])):
                if candidate not in candidates:
                    candidates.append(candidate)
    for c in candidates:
        if c in available_genes:
            return c
    return None


# ---------------------------------------------------------------------------
# VCF annotation helpers
# ---------------------------------------------------------------------------

def ann_entries(info: str) -> list[list[str]]:
    """Extract ANN entries (each a list of fields) from a VCF INFO string."""
    for token in info.split(";"):
        if token.startswith("ANN="):
            return [entry.split("|") for entry in token[4:].split(",") if entry]
    return []


def has_info_flag(info: str, key: str) -> bool:
    """Whether a VCF INFO string carries `key`, as a flag or a key=value pair.

    Token-wise rather than a substring test, so a key is never matched inside another key's
    name or value.
    """
    return any(token == key or token.startswith(f"{key}=") for token in info.split(";"))


def aa_changes_from_info(info: str, unparsed: list | None = None) -> list[tuple[str, str, int, str]]:
    """Non-synonymous (gene, feature_id, aa_pos, alt_aa) changes annotated on a VCF row.

    One item per ANN entry that carries a parseable HGVS.p substitution, in ANN order and
    de-duplicated, so a row hitting several features yields one candidate per feature.

    ``feature_id`` (ANN field 7) is carried alongside ``gene`` (field 4) because a gene name
    alone cannot key a mature peptide: NCBI-style HA annotations put the *parent* gene there
    for the HA1/HA2 entries, and only the feature id names the subunit.

    Non-empty HGVS.p values that parse as neither substitution spelling are appended to
    ``unparsed`` so the caller can say so rather than reporting nothing at all.
    """
    changes: list[tuple[str, str, int, str]] = []
    seen: set = set()
    for fields in ann_entries(info):
        if len(fields) < 11:
            continue
        gene, feature_id, hgvs_p = fields[3], fields[6], fields[10]
        parsed = parse_hgvs_p(hgvs_p)
        if parsed is None:
            if hgvs_p and unparsed is not None and not HGVS_P_NON_SUB_RE.search(hgvs_p):
                unparsed.append(hgvs_p)
            continue
        ref_aa, pos, alt_aa = parsed
        if ref_aa == alt_aa:  # synonymous
            continue
        key = (gene, feature_id, pos, alt_aa)
        if key not in seen:
            seen.add(key)
            changes.append(key)
    return changes


def sanitize(text: str) -> str:
    """Make a free-text database value safe for a VCF INFO value.

    VCF INFO values may not contain whitespace, ';' or ','; '|' and '=' are reserved by
    this field's own structure. FluMutDB effect names contain spaces, parentheses and an
    en dash.
    """
    text = text.replace("–", "-").replace("—", "-")
    text = re.sub(r"\s+", "_", text.strip())
    return re.sub(r"[;,=|]", "", text)


# ---------------------------------------------------------------------------
# Inputs (FluMutDB, coordinate map)
# ---------------------------------------------------------------------------

def locate_db(explicit: str | None) -> str:
    """Path to flumut_db.sqlite: --flumut-db when given, else the installed flumutdb."""
    if explicit:
        if not os.path.exists(explicit):
            log.error("FluMutDB not found at %s", explicit)
            sys.exit(1)
        return explicit
    try:
        from flumut.DbReader import get_db_file
    except ImportError:
        log.error("flumut is not importable and --flumut-db was not given; "
                  "install flumut/flumutdb (see environments/fluscan.yml) or pass a path")
        sys.exit(1)
    return str(get_db_file())


def db_version(conn: sqlite3.Connection) -> str:
    """FluMutDB version string for logging ('unknown' if the table is missing)."""
    try:
        major, minor, date = conn.execute("SELECT major, minor, date FROM db_version").fetchone()
        return f"{major}.{minor} ({date})"
    except (sqlite3.Error, TypeError):
        return "unknown"


def load_segment_map(coord_map_path: str, segment: str) -> dict:
    """The coordinate map's block for one segment ({} if the segment is absent)."""
    with open(coord_map_path) as fh:
        coord = json.load(fh)
    seg_map = coord.get("segments", {}).get(segment)
    if not seg_map:
        log.warning("Segment %s not present in coordinate map; copying VCF unannotated", segment)
        return {}
    return seg_map


def build_alias_index(genes_map: dict) -> dict[str, str]:
    """{ANN Feature_ID alias -> coordinate-map gene key} for this segment.

    Only mature-peptide gene keys carry aliases. Their ANN entries name the *parent* gene in
    Gene_Name (an NCBI-style HA1 entry says gene 'HA'), so the subunit is identifiable only
    through Feature_ID -- the peptide's protein_id or GFF3 ID, which is what these aliases
    hold. Empty for datasets whose subunits are real CDS features.
    """
    index: dict[str, str] = {}
    for gene_key, gene_doc in genes_map.items():
        for alias in gene_doc.get("aliases", []):
            index.setdefault(alias, gene_key)
    return index


def build_gene_protein_maps(conn: sqlite3.Connection, db_ref: str, genes_map: dict) -> dict:
    """Map each coordinate-map gene key to its FluMutDB protein and a target->source index.

    Returns {gene_key: (protein_name, {target_aa_pos: source_aa_pos})}. The coordinate map
    stores aa_source_to_target (flumut reference -> Nextclade); inverting it once per gene
    is what lets the main loop start from a VCF row and query the database directly.
    """
    gene_maps: dict[str, tuple[str, dict]] = {}
    available = set(genes_map)
    unresolved: list[str] = []

    for (protein,) in conn.execute(PROTEIN_QUERY, (db_ref,)):
        gene_key = resolve_gene(protein, available)
        if gene_key is None:
            unresolved.append(protein)
            continue
        if gene_key in gene_maps:
            # Two FluMutDB proteins claiming one gene key means the annotation is too coarse
            # to tell them apart (a single 'NS' gene for NS1 and NEP, say). Overwriting would
            # silently query one protein's markers through the other's coordinate map; refuse
            # instead, and keep the first -- PROTEIN_QUERY is ORDER BY'd, so that is stable.
            log.warning("FluMutDB proteins %s and %s both resolve to coordinate-map gene key "
                        "%s; keeping %s and skipping %s. The target annotation does not "
                        "distinguish them, so %s markers cannot be annotated.",
                        gene_maps[gene_key][0], protein, gene_key,
                        gene_maps[gene_key][0], protein, protein)
            continue
        # aa_source_to_target comes from a pairwise alignment, so aligned columns pair 1:1
        # and the inverse is well defined. Keep the lowest source position if that ever
        # stops holding, so the annotation cannot depend on dict ordering.
        source_to_target = genes_map[gene_key].get("aa_source_to_target", {})
        target_to_source: dict[int, int] = {}
        for src, tgt in sorted((int(s), t) for s, t in source_to_target.items()):
            if tgt is None:
                continue
            if int(tgt) in target_to_source:
                log.warning("Gene %s: target AA %d has several source positions (%d, %d); "
                            "keeping %d", gene_key, int(tgt), target_to_source[int(tgt)],
                            src, target_to_source[int(tgt)])
                continue
            target_to_source[int(tgt)] = src
        if not target_to_source:
            log.warning("Gene %s (FluMutDB protein %s) has an empty amino-acid map; skipping",
                        gene_key, protein)
            continue
        gene_maps[gene_key] = (protein, target_to_source)

    if unresolved:
        log.warning("FluMutDB protein(s) %s have no counterpart gene in the coordinate map "
                    "(available: %s); their markers cannot be annotated",
                    ", ".join(unresolved), ", ".join(sorted(available)) or "none")
    return gene_maps


# ---------------------------------------------------------------------------
# FluMutDB lookup
# ---------------------------------------------------------------------------

def make_marker_lookup(conn: sqlite3.Connection, db_ref: str):
    """Build a cached (protein, source_pos, alt_aa) -> marker matches lookup.

    Each match is (mutation_name, direct_effects, associated_effects), the two effect
    tuples holding sorted, deduplicated 'Subtype:Effect' values. An effect is *direct* when
    the marker carrying it is defined by this mutation alone (all_mutations_count == 1) and
    *associated* otherwise; every association is classified, and a mutation carried by both
    kinds of marker legitimately populates both tuples.
    """

    @lru_cache(maxsize=None)
    def lookup(protein: str, source_pos: int, alt_aa: str) -> tuple:
        # dict rather than a set of names: preserves mutation_name order (the query is
        # ORDER BY'd) while folding the several effect rows each name fans out to.
        matches: dict[str, tuple[set, set]] = {}
        for name, mutation_count, subtype, effect in conn.execute(
            MARKER_QUERY, (db_ref, protein, source_pos, alt_aa)
        ):
            direct, associated = matches.setdefault(sanitize(name), (set(), set()))
            if not effect:  # LEFT JOIN miss: mutation in no marker, or marker with no effect
                continue
            # Qualifying by subtype is what keeps the field honest: the same mutation is
            # reported with opposite effects in different experimental systems.
            qualified = f"{sanitize(subtype)}:{sanitize(effect)}"
            (direct if mutation_count == 1 else associated).add(qualified)
        return tuple(
            (name, tuple(sorted(direct)), tuple(sorted(associated)))
            for name, (direct, associated) in matches.items()
        )

    return lookup


def resolve_ann_gene(gene: str, feature_id: str, gene_maps: dict,
                     alias_index: dict) -> tuple[str, dict] | None:
    """The (protein, target->source map) an ANN entry's gene/feature names key onto.

    Three attempts, in order:

    1. The coordinate map's gene keys are canonicalised from the target GFF, and ANN's
       Gene_Name comes from the same annotation uncanonicalised -- so canonicalise here too,
       or a dataset spelling the gene NS2 never reaches its NEP map.
    2. Feature_ID against the map's aliases, which is how an NCBI-style HA1/HA2 mature
       peptide is identified: its Gene_Name is the parent gene 'HA', and only its
       protein_id (or GFF3 ID) names the subunit FluMutDB numbers by.
    3. Feature_ID as a gene key in its own right, for mature rows carrying no protein_id --
       load_gff3 then falls back to `product`, putting 'HA1' in Feature_ID directly.
    """
    for candidate in (canonical_gene(gene), alias_index.get(feature_id),
                      canonical_gene(feature_id) if feature_id else None):
        if candidate is None:
            continue
        mapped = gene_maps.get(candidate)
        if mapped is not None:
            return mapped
    return None


def row_flumut_field(info: str, gene_maps: dict, lookup, want_effects: bool,
                     alias_index: dict | None = None,
                     subunit_parent_genes: set | None = None,
                     unmapped_genes: set | None = None,
                     unparsed_hgvs_p: list | None = None) -> str:
    """The 'FLUMUT=...' INFO sub-field for one VCF row, or '' when nothing matches.

    Gene names carrying an AA change that this segment's map cannot place are collected in
    ``unmapped_genes`` so the caller can report them: a VCF annotated against a dataset
    that FluMutDB does not share a gene split with (a single ``HA`` gene rather than
    ``HA1``/``HA2``) would otherwise match nothing, silently.

    A name is only reported when *nothing* on the row resolved. On a mature-peptide dataset
    every subunit row also carries the parent CDS's own ANN entry, which cannot resolve by
    design; reporting that would put a warning on a row that was in fact annotated correctly.
    """
    alias_index = alias_index or {}
    subunit_parent_genes = subunit_parent_genes or set()
    entries: list[str] = []
    unresolved: list[str] = []
    resolved_any = False
    for gene, feature_id, target_pos, alt_aa in aa_changes_from_info(info, unparsed_hgvs_p):
        mapped = resolve_ann_gene(gene, feature_id, gene_maps, alias_index)
        if mapped is None:
            if gene not in subunit_parent_genes:
                unresolved.append(gene)
            continue
        resolved_any = True
        protein, target_to_source = mapped
        source_pos = target_to_source.get(target_pos)
        if source_pos is None:
            log.debug("%s target AA %d has no flumut-reference counterpart; skipping",
                      gene, target_pos)
            continue
        matches = lookup(protein, source_pos, alt_aa)
        if not matches:
            continue
        # A handful of (protein, position, alt) keys are shared by two mutation names
        # (e.g. HA1-5:N182K and HA1-5:V182K), which differ only in the reference residue
        # this VCF does not constrain. Name and effect fields both union over them.
        names = "&".join(name for name, _, _ in matches)
        fields = [protein, str(target_pos), alt_aa, names]
        if want_effects:
            direct = sorted({e for _, eff, _ in matches for e in eff})
            associated = sorted({e for _, _, eff in matches for e in eff})
            fields.append("&".join(direct) or ".")
            fields.append("&".join(associated) or ".")
        entries.append("|".join(fields))
    if unresolved and not resolved_any and unmapped_genes is not None:
        unmapped_genes.update(unresolved)
    return "FLUMUT=" + ",".join(entries) if entries else ""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def get_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # The process's provenance block reads this rather than repeating the literal, so the
    # recorded tool_version cannot drift from SCRIPT_VERSION.
    p.add_argument("--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}")
    p.add_argument("--vcf", required=True, help="Existing per-sample/per-segment VCF to annotate")
    p.add_argument("--coord-map", required=True, help="coordinate_map.<SEG>.json (flumut -> Nextclade)")
    p.add_argument("--segment", required=True, help="Segment label, e.g. HA, M, NS")
    p.add_argument("--sample", required=True, help="Sample id (used for logging only)")
    p.add_argument("--flumut-db", help="Path to flumut_db.sqlite (default: the installed flumutdb)")
    p.add_argument("--no-effects", action="store_true",
                   help="Omit both effect fields, emitting only Protein|AA_Pos|Alt_AA|Mutation_Name")
    p.add_argument("-o", "--output", required=True, help="Output VCF path")
    return p.parse_args()


def main() -> None:
    args = get_args()
    segment = args.segment.upper()
    want_effects = not args.no_effects

    df, header_lines = read_vcf(args.vcf)

    seg_map = load_segment_map(args.coord_map, segment)
    genes_map = seg_map.get("genes", {})
    alias_index = build_alias_index(genes_map)
    subunit_parent_genes = set(seg_map.get("subunit_parent_genes", []))
    if seg_map.get("stale_annotation"):
        log.warning("Coordinate map flags segment %s as stale (the Nextclade dataset does not "
                    "expose the HA1/HA2 split FluMutDB requires); HA markers cannot be annotated",
                    segment)

    db_path = locate_db(args.flumut_db)
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        version = db_version(conn)
        db_ref = SEGMENT_TO_DB_REFERENCE.get(segment, segment)
        gene_maps = build_gene_protein_maps(conn, db_ref, genes_map)
        lookup = make_marker_lookup(conn, db_ref)

        # Annotate in place: rows, order, POS, REF/ALT and ANN are all preserved; only INFO
        # grows, and only on rows that matched.
        info_out: list[str] = []
        annotated = 0
        markers: set = set()
        unmapped_genes: set = set()
        unparsed_hgvs_p: list[str] = []
        mnv_unresolved = 0
        for info in df["INFO"].astype(str):
            if has_info_flag(info, "MNV_UNRESOLVED"):
                # fix_mnv_codons.py sets this when it identified a multi-SNV codon but could
                # not compute the combined change, leaving the per-row HGVS.p at a value the
                # pipeline has already determined to be wrong. Matching a curated marker
                # against it would manufacture a false positive.
                mnv_unresolved += 1
                info_out.append(info)
                continue
            field = row_flumut_field(info, gene_maps, lookup, want_effects,
                                     alias_index=alias_index,
                                     subunit_parent_genes=subunit_parent_genes,
                                     unmapped_genes=unmapped_genes,
                                     unparsed_hgvs_p=unparsed_hgvs_p)
            if field:
                annotated += 1
                for entry in field[len("FLUMUT="):].split(","):
                    markers.update(entry.split("|")[3].split("&"))
                info_out.append(f"{info};{field}")
            else:
                info_out.append(info)
    finally:
        conn.close()

    if mnv_unresolved:
        log.warning("%d row(s) carry the MNV_UNRESOLVED flag; their HGVS.p is known to be "
                    "wrong, so they were skipped rather than matched against FluMutDB",
                    mnv_unresolved)

    if unparsed_hgvs_p:
        # Without this the failure mode is indistinguishable from "this VCF has no coding
        # changes": an HGVS.p spelling neither regex handles yields silence, not an error.
        log.warning("%d ANN entr(ies) carry an HGVS.p substitution this script cannot parse "
                    "(e.g. %r) and were skipped; only 1-letter (p.E627K) and 3-letter "
                    "(p.Glu627Lys) forms are recognised",
                    len(unparsed_hgvs_p), unparsed_hgvs_p[0])

    if unmapped_genes:
        # State the general condition first. Unresolvable HA is only one of its causes, and
        # leading with the HA explanation misdirects every other case (an NS2/NEP spelling
        # mismatch, say) towards a segment that has nothing to do with it.
        message = ("VCF gene(s) %s carry amino-acid changes but have no FluMutDB counterpart "
                   "in the coordinate map: an ANN Gene_Name must match one of the map's gene "
                   "keys (mapped gene(s): %s).")
        if segment == "HA" or seg_map.get("stale_annotation"):
            message += (" For HA specifically, FluMutDB numbers the HA1/HA2 subunits, so the "
                        "Nextclade dataset must expose them as separate genes.")
        log.warning(message,
                    ", ".join(sorted(unmapped_genes)),
                    ", ".join(sorted(gene_maps)) or "none")

    out_df = df.copy()
    out_df["INFO"] = info_out
    out_header = header_lines[:-1] + [FLUMUT_HEADER, header_lines[-1]]
    write_vcf(out_df, out_header, args.output)

    log.info(
        "Wrote %s: %d row(s) scanned, %d annotated, %d distinct FluMutDB marker(s) matched "
        "across gene(s) %s for %s/%s (v%s, FluMutDB %s)",
        args.output, len(df), annotated, len(markers),
        ", ".join(sorted(gene_maps)) or "none", args.sample, segment, SCRIPT_VERSION, version,
    )


if __name__ == "__main__":
    main()
