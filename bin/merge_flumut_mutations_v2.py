#!/usr/bin/env python3
"""Annotate an existing fluscan/Nextclade VCF with FluMutDB marker matches.

Version 2 of the flumut merge. Where v1 (merge_flumut_mutations.py) consumed flumut's
per-sample melted TSV and *inserted* synthetic rows for markers fluscan had not called,
this version never adds, removes or reorders a row. It walks the VCF, and for every amino
-acid change already annotated in a row's ANN= field it queries FluMutDB for a curated
marker matching that (protein, AA position, alt AA). Matching rows gain one extra INFO
field:

    FLUMUT=<Protein>|<AA_Pos>|<Alt_AA>|<Mutation_Name>|<Effects>

one comma-separated entry per matching ANN annotation. ANN= itself is left untouched, so
it stays strictly snpEff 16-field compatible.

The marker set comes from the database rather than from a flumut run, so no per-sample
flumut output is needed. FluMutDB positions live in flumut's own reference coordinates and
are translated onto the Nextclade dataset via the coordinate_conversion map, the same
aa_source_to_target block v1 uses (applied in reverse here, since we start from the VCF).

FluMutDB annotates HA as HA1/HA2 subunits, so HA markers only match when the Nextclade
dataset also exposes dedicated HA1/HA2 genes. Datasets carrying a single HA gene produce a
warning and pass through unannotated.

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

# read_vcf/write_vcf round-trip a VCF byte-for-byte; tools.py lives in bin/ alongside this
# script, so adding this script's directory to sys.path makes it importable.
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from tools import read_vcf, write_vcf  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

FLUMUT_HEADER = (
    '##INFO=<ID=FLUMUT,Number=.,Type=String,'
    'Description="FluMutDB marker matches: '
    "'Protein | AA_Pos | Alt_AA | Mutation_Name | Effects'. "
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

HGVS_P_RE = re.compile(r"p\.([A-Za-z*])(\d+)([A-Za-z*])")

# One parameterised lookup per (protein, position, alt AA). ref_seq is deliberately absent
# from the WHERE clause: a marker is defined by its alt allele, and the sample's reference
# amino acid may legitimately differ from FluMutDB's. group_concat() joins on ',', which is
# safe here because no effect name in the database contains a comma.
MARKER_QUERY = """
SELECT mm.mutation_name, group_concat(DISTINCT me.effect_name)
FROM mutation_mappings mm
JOIN mutations m               ON m.name = mm.mutation_name
LEFT JOIN markers_mutations mk ON mk.mutation_name = mm.mutation_name
LEFT JOIN markers_effects me   ON me.marker_id = mk.marker_id
WHERE mm.reference_name = ? AND m.protein_name = ? AND mm.position = ? AND mm.alt_seq = ?
GROUP BY mm.mutation_name
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
    """
    candidates: list[str] = [raw]
    m = re.match(r"^(.*?)-(\d+)$", raw)
    if m:
        base, suffix = m.group(1), m.group(2)
        candidates.append(base)            # HA1-5 -> HA1 ; NS-1 -> NS
        candidates.append(base + suffix)   # NS-1 -> NS1 ; NS-2 -> NS2
    for c in list(candidates):
        candidates.extend(GENE_ALIASES.get(c, []))
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


def aa_changes_from_info(info: str) -> list[tuple[str, int, str]]:
    """Non-synonymous (gene, aa_pos, alt_aa) changes annotated on a VCF row.

    One item per ANN entry that carries a parseable HGVS.p substitution, in ANN order and
    de-duplicated, so a row hitting several features yields one candidate per feature.
    """
    changes: list[tuple[str, int, str]] = []
    seen: set = set()
    for fields in ann_entries(info):
        if len(fields) < 11:
            continue
        gene = fields[3]
        hp = HGVS_P_RE.search(fields[10])
        if not hp:
            continue
        ref_aa, pos, alt_aa = hp.group(1), int(hp.group(2)), hp.group(3)
        if ref_aa == alt_aa:  # synonymous
            continue
        key = (gene, pos, alt_aa)
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

def make_marker_lookup(conn: sqlite3.Connection, db_ref: str, want_effects: bool):
    """Build a cached (protein, source_pos, alt_aa) -> [(mutation_name, effects)] lookup."""

    @lru_cache(maxsize=None)
    def lookup(protein: str, source_pos: int, alt_aa: str) -> tuple:
        rows = conn.execute(MARKER_QUERY, (db_ref, protein, source_pos, alt_aa)).fetchall()
        results = []
        for mutation_name, effects in rows:
            names = sorted({sanitize(e) for e in (effects or "").split(",") if e.strip()})
            results.append((sanitize(mutation_name), "&".join(names) if want_effects else ""))
        return tuple(results)

    return lookup


def row_flumut_field(info: str, gene_maps: dict, lookup, want_effects: bool,
                     unmapped_genes: set | None = None) -> str:
    """The 'FLUMUT=...' INFO sub-field for one VCF row, or '' when nothing matches.

    Gene names carrying an AA change that this segment's map cannot place are collected in
    ``unmapped_genes`` so the caller can report them: a VCF annotated against a dataset
    that FluMutDB does not share a gene split with (a single ``HA`` gene rather than
    ``HA1``/``HA2``) would otherwise match nothing, silently.
    """
    entries: list[str] = []
    for gene, target_pos, alt_aa in aa_changes_from_info(info):
        mapped = gene_maps.get(gene)
        if mapped is None:
            if unmapped_genes is not None:
                unmapped_genes.add(gene)
            continue
        protein, target_to_source = mapped
        source_pos = target_to_source.get(target_pos)
        if source_pos is None:
            log.debug("%s target AA %d has no flumut-reference counterpart; skipping",
                      gene, target_pos)
            continue
        matches = lookup(protein, source_pos, alt_aa)
        if not matches:
            continue
        names = "&".join(name for name, _ in matches)
        fields = [protein, str(target_pos), alt_aa, names]
        if want_effects:
            effects = sorted({e for _, eff in matches for e in eff.split("&") if e})
            fields.append("&".join(effects))
        entries.append("|".join(fields))
    return "FLUMUT=" + ",".join(entries) if entries else ""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def get_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="Existing per-sample/per-segment VCF to annotate")
    p.add_argument("--coord-map", required=True, help="coordinate_map.<SEG>.json (flumut -> Nextclade)")
    p.add_argument("--segment", required=True, help="Segment label, e.g. HA, M, NS")
    p.add_argument("--sample", required=True, help="Sample id (used for logging only)")
    p.add_argument("--flumut-db", help="Path to flumut_db.sqlite (default: the installed flumutdb)")
    p.add_argument("--no-effects", action="store_true",
                   help="Omit the effect summary, emitting only Protein|AA_Pos|Alt_AA|Mutation_Name")
    p.add_argument("-o", "--output", required=True, help="Output VCF path")
    return p.parse_args()


def main() -> None:
    args = get_args()
    segment = args.segment.upper()
    want_effects = not args.no_effects

    df, header_lines = read_vcf(args.vcf)

    seg_map = load_segment_map(args.coord_map, segment)
    genes_map = seg_map.get("genes", {})
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
        lookup = make_marker_lookup(conn, db_ref, want_effects)

        # Annotate in place: rows, order, POS, REF/ALT and ANN are all preserved; only INFO
        # grows, and only on rows that matched.
        info_out: list[str] = []
        annotated = 0
        markers: set = set()
        unmapped_genes: set = set()
        for info in df["INFO"].astype(str):
            field = row_flumut_field(info, gene_maps, lookup, want_effects, unmapped_genes)
            if field:
                annotated += 1
                for entry in field[len("FLUMUT="):].split(","):
                    markers.update(entry.split("|")[3].split("&"))
                info_out.append(f"{info};{field}")
            else:
                info_out.append(info)
    finally:
        conn.close()

    if unmapped_genes:
        log.warning("VCF gene(s) %s carry amino-acid changes but have no FluMutDB counterpart "
                    "in the coordinate map (mapped gene(s): %s). A VCF annotated against a "
                    "dataset with a single 'HA' gene lands here — FluMutDB numbers HA as the "
                    "HA1/HA2 subunits, so the dataset must expose them too.",
                    ", ".join(sorted(unmapped_genes)),
                    ", ".join(sorted(gene_maps)) or "none")

    out_df = df.copy()
    out_df["INFO"] = info_out
    out_header = header_lines[:-1] + [FLUMUT_HEADER, header_lines[-1]]
    write_vcf(out_df, out_header, args.output)

    log.info(
        "Wrote %s: %d row(s) scanned, %d annotated, %d distinct FluMutDB marker(s) matched "
        "across gene(s) %s for %s/%s (FluMutDB %s)",
        args.output, len(df), annotated, len(markers),
        ", ".join(sorted(gene_maps)) or "none", args.sample, segment, version,
    )


if __name__ == "__main__":
    main()
