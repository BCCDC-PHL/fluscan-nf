#!/usr/bin/env python3
"""Annotate an existing fluscan/Nextclade VCF with FluMutDB marker matches.

Walks the VCF and, for every amino-acid change already annotated in a row's ANN= field,
queries FluMutDB for a curated marker matching that (protein, AA position, alt AA).
Matching rows gain one extra INFO field, one comma-separated entry per matching ANN
annotation; ANN= itself is untouched. Rows are never added, removed or reordered.

    FLUMUT=<Protein>|<AA_Pos>|<Alt_AA>|<Mutation_Name>|<DirectEffect>|<AssociatedEffect>

DirectEffect carries effects of markers this mutation alone defines; AssociatedEffect
carries them for markers needing further mutations, which are NOT asserted to be present.
Both hold 'Subtype:Effect' pairs, '.' when empty.

FluMutDB positions are in flumut's reference coordinates and are translated onto the
Nextclade dataset through the coordinate_conversion map. FluMutDB numbers HA by its
HA1/HA2 subunits, so HA markers match only where the dataset exposes that split.

Pipeline position: after nextclade_run_batch, once per segment over every sample's VCF,
consuming coordinate_conversion output. See docs/vcf_annotation_and_flumut_insertion.md.
"""

import argparse
import json
import logging
import os
import re
import sqlite3
import sys
from functools import lru_cache

# realpath so `from lib...` still resolves when bin/ is reached through a symlink.
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from lib.aa_markers import (  # noqa: E402
    AA3_TO_AA1,
    HGVS_P_NON_SUB_RE,
    HGVS_P_RE,
    HGVS_P_RE_3LETTER,
    aa_changes_from_info,
    ann_entries,
    has_info_flag,
    parse_hgvs_p,
)
from lib.flu_gene_names import canonical_gene  # noqa: E402
from lib.tools import read_vcf, write_vcf  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

SCRIPT_VERSION = "2.3.0"

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
# Gene resolution generates candidates and keeps whichever exists as a gene key in the
# coordinate map for the segment being processed.
GENE_ALIASES = {
    "NS2": ["NEP"],
    "NEP": ["NS2"],
    "M": ["MP"],
    "MP": ["M"],
}

# HGVS.p parsing, ANN extraction and the INFO-flag test are imported from lib/aa_markers.py
# rather than reimplemented: mutation_watcher.py asks the same questions of the same VCFs
# and the two must not disagree about which spellings parse.

# ref_seq is deliberately absent from the WHERE clause: a marker is defined by its alt
# allele, and the sample's reference amino acid may legitimately differ from FluMutDB's.
# One row per marker/effect association rather than a group_concat, so each effect keeps its
# marker's mutation count and subtype. Every join stays LEFT: detecting a mutation must not
# depend on it having any published effect.
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
    """Resolve a FluMutDB protein name to a coordinate-map gene key, or None.

    Candidates are tried most-specific first and the bare stem last. That ordering is
    load-bearing: a dataset naming its NS1 gene 'NS' offers the stem to both NS-1 and NS-2,
    and resolving NS-2 that way would query NEP positions through NS1's map -- false
    positives rather than silence. See docs/vcf_annotation_and_flumut_insertion.md.
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

def sanitize(text: str) -> str:
    """Make a free-text database value safe for a VCF INFO value.

    INFO may not contain whitespace, ';' or ','; '|' and '=' are reserved by this field.
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

    Only mature peptides carry aliases: their ANN entry names the *parent* gene, so the
    subunit is identifiable only through Feature_ID. Empty where subunits are real CDSs.
    """
    index: dict[str, str] = {}
    for gene_key, gene_doc in genes_map.items():
        for alias in gene_doc.get("aliases", []):
            index.setdefault(alias, gene_key)
    return index


def build_gene_protein_maps(conn: sqlite3.Connection, db_ref: str, genes_map: dict) -> dict:
    """Map each coordinate-map gene key to its FluMutDB protein and a target->source index.

    Returns {gene_key: (protein_name, {target_aa_pos: source_aa_pos})}. The map stores
    aa_source_to_target; inverting it once per gene lets the main loop start from a VCF row.
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
            # Two proteins on one gene key means the annotation cannot tell them apart.
            # Overwriting would query one protein's markers through the other's map, so
            # refuse and keep the first -- PROTEIN_QUERY is ORDER BY'd, so that is stable.
            log.warning("FluMutDB proteins %s and %s both resolve to coordinate-map gene key "
                        "%s; keeping %s and skipping %s. The target annotation does not "
                        "distinguish them, so %s markers cannot be annotated.",
                        gene_maps[gene_key][0], protein, gene_key,
                        gene_maps[gene_key][0], protein, protein)
            continue
        # aa_source_to_target comes from a pairwise alignment, so columns pair 1:1 and the
        # inverse is well defined. Keep the lowest source position so the annotation can
        # never depend on dict ordering.
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

    Each match is (mutation_name, direct_effects, associated_effects), holding sorted
    'Subtype:Effect' values. An effect is *direct* when its marker is defined by this
    mutation alone; a mutation in both kinds of marker populates both tuples.
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

    Three attempts, in order: canonicalised Gene_Name (map keys are canonicalised too, so a
    dataset spelling the gene NS2 still reaches its NEP map); Feature_ID against the map's
    aliases, which is how a mature HA1/HA2 peptide is identified; then Feature_ID as a gene
    key in its own right. See docs/vcf_annotation_and_flumut_insertion.md § Gene names.
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

    Genes this segment's map cannot place go into ``unmapped_genes`` for the caller to
    report, so a dataset lacking FluMutDB's gene split does not just match nothing silently.
    A name is reported only when *nothing* on the row resolved: on a mature-peptide dataset
    every subunit row also carries the parent CDS entry, which cannot resolve by design.
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
        # A few (protein, position, alt) keys are shared by two mutation names differing
        # only in a reference residue this VCF does not constrain; union over them.
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
    # The process's provenance block reads this, so tool_version cannot drift.
    p.add_argument("--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}")
    p.add_argument("--vcf", required=True, nargs="+",
                   help="Per-sample VCFs for this segment; each is written to "
                        "<stem>_flumut_annotated.vcf in the working directory, the stem being "
                        "its filename without .vcf and a trailing _final")
    p.add_argument("--coord-map", required=True, help="coordinate_map.<SEG>.json (flumut -> Nextclade)")
    p.add_argument("--segment", required=True, help="Segment label, e.g. HA, M, NS")
    p.add_argument("--flumut-db", help="Path to flumut_db.sqlite (default: the installed flumutdb)")
    p.add_argument("--no-effects", action="store_true",
                   help="Omit both effect fields, emitting only Protein|AA_Pos|Alt_AA|Mutation_Name")
    return p.parse_args()


def annotate_vcf(vcf_path: str, segment: str, seg_map: dict, gene_maps: dict, lookup,
                 want_effects: bool, alias_index: dict, subunit_parent_genes: set,
                 version: str) -> None:
    """Write <stem>_flumut_annotated.vcf for one VCF, logging what it matched and skipped."""
    df, header_lines = read_vcf(vcf_path)
    stem = os.path.basename(vcf_path).removesuffix(".vcf").removesuffix("_final")
    output = f"{stem}_flumut_annotated.vcf"

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
            # fix_mnv_codons.py left this row's HGVS.p at a value already determined to
            # be wrong; matching a curated marker against it would be a false positive.
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

    if mnv_unresolved:
        log.warning("%d row(s) carry the MNV_UNRESOLVED flag; their HGVS.p is known to be "
                    "wrong, so they were skipped rather than matched against FluMutDB",
                    mnv_unresolved)

    if unparsed_hgvs_p:
        # Without this an unhandled HGVS.p spelling is indistinguishable from "this VCF has
        # no coding changes": silence rather than an error.
        log.warning("%d ANN entr(ies) carry an HGVS.p substitution this script cannot parse "
                    "(e.g. %r) and were skipped; only 1-letter (p.E627K) and 3-letter "
                    "(p.Glu627Lys) forms are recognised",
                    len(unparsed_hgvs_p), unparsed_hgvs_p[0])

    if unmapped_genes:
        # General condition first: HA is only one cause, and leading with it misdirects
        # every other case towards a segment that has nothing to do with it.
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
    write_vcf(out_df, out_header, output)

    log.info(
        "Wrote %s: %d row(s) scanned, %d annotated, %d distinct FluMutDB marker(s) matched "
        "across gene(s) %s for %s/%s (v%s, FluMutDB %s)",
        output, len(df), annotated, len(markers),
        ", ".join(sorted(gene_maps)) or "none", stem.split("__")[0], segment, SCRIPT_VERSION, version,
    )


def main() -> None:
    args = get_args()
    segment = args.segment.upper()

    seg_map = load_segment_map(args.coord_map, segment)
    genes_map = seg_map.get("genes", {})
    alias_index = build_alias_index(genes_map)
    subunit_parent_genes = set(seg_map.get("subunit_parent_genes", []))
    if seg_map.get("stale_annotation"):
        log.warning("Coordinate map flags segment %s as stale (the Nextclade dataset does not "
                    "expose the HA1/HA2 split FluMutDB requires); HA markers cannot be annotated",
                    segment)

    db_path = locate_db(args.flumut_db)
    # The packaged DB never changes; immutable skips the file locking that fails on network mounts.
    conn = sqlite3.connect(f"file:{db_path}?immutable=1", uri=True)
    try:
        version = db_version(conn)
        db_ref = SEGMENT_TO_DB_REFERENCE.get(segment, segment)
        gene_maps = build_gene_protein_maps(conn, db_ref, genes_map)
        # Queries conn lazily and caches, so one lookup serves every VCF in the batch.
        lookup = make_marker_lookup(conn, db_ref)
        for vcf_path in args.vcf:
            annotate_vcf(vcf_path, segment, seg_map, gene_maps, lookup, not args.no_effects,
                         alias_index, subunit_parent_genes, version)
    finally:
        conn.close()


if __name__ == "__main__":
    main()
