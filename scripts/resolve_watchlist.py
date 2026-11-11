#!/usr/bin/env python3
"""Resolve a source watchlist's positions onto one reference bundle's numbering.

    resolve_watchlist.py --source nai_h5n1_source.csv \\
                         --db /path/to/references/h5n1/chicken2021 \\
                         --out-dir assets/watchlists

A source list is a CSV written in whatever numbering its source used. It declares what
its columns mean in a `#` header block, so the mapping travels with the file:

    # numbering: n2
    # columns: gene=Protein, mutation=Substitution
    # offsets: 0, 20
    # source: free-text provenance

Only `numbering` is required. Columns default to `gene` and `mutation`; every other
column is carried verbatim into the generated note, so a supplier can add as many as
they like without touching this script. `offsets` lists the candidate offsets between
the source frame and the target, defaulting to 0 -- positions already in the target's
numbering, validated rather than moved.

The output is a resolved list in the target's numbering, which is the only kind the
pipeline scans, plus a report carrying every entry with its status so an exclusion is
visible rather than silent.

Positions are looked up, not aligned. Resolving a marker numbered in one lineage's
protein onto another lineage's needs alignment, because the proteins differ in length;
that is out of scope here and the resolver would have to be told to do it.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from dataclasses import dataclass
from datetime import date
from pathlib import Path

from Bio import SeqIO

# scripts/ is maintainer tooling and is not on PATH in a task, so bin/ has to be put on
# sys.path to reuse the pipeline's own CDS reader rather than writing a second one.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from build_reference_coordinate_map import (  # noqa: E402
    SEGMENTS,
    read_cds_models,
    translate_model,
)
from lib.aa_markers import parse_hgvs_p  # noqa: E402
from lib.flu_gene_names import KNOWN_GENE_KEYS, canonical_gene  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

SCRIPT_VERSION = "1.0.0"

DIRECTIVE_KEYS = ("numbering", "columns", "offsets", "source")

REPORT_COLUMNS = (
    "gene", "source_mutation", "source_position", "ref_aa", "alt_aa",
    "offset", "target_position", "target_residue", "status", "note",
)

# A status that resolves to a position the pipeline can scan. Everything else is reported
# and left out of the generated list -- never dropped without a row in the report.
INCLUDED_STATUSES = ("ok", "ref_differs")


class SourceError(Exception):
    """A source list that cannot be resolved. Always fatal."""


@dataclass(frozen=True)
class Resolution:
    gene: str
    source_mutation: str
    source_position: int
    ref_aa: str
    alt_aa: str
    offset: int | None
    target_position: int | None
    target_residue: str
    status: str
    note: str

    @property
    def key(self) -> tuple[str, int, str]:
        return (self.gene, self.target_position, self.alt_aa)

    @property
    def mutation(self) -> str:
        """Reference-oriented: the residue this bundle carries, not the one the source claims.

        Writing the source's REF letter would make the watcher's REF check fire on every
        run for every entry whose wild type differs here. Written this way it stays quiet
        unless the bundle itself changes, which turns the check into a staleness detector.
        """
        return f"{self.target_residue}{self.target_position}{self.alt_aa}"


def read_directives(lines: list[str]) -> dict[str, str]:
    """The `# key: value` header block. Comment lines that are not directives are prose."""
    directives = {}
    for line in lines:
        if not line.lstrip().startswith("#"):
            continue
        body = line.lstrip().lstrip("#").strip()
        key, sep, value = body.partition(":")
        if sep and key.strip().lower() in DIRECTIVE_KEYS:
            directives[key.strip().lower()] = value.strip()
    return directives


def parse_column_map(raw: str) -> dict[str, str]:
    """`gene=Protein, mutation=Substitution` -> {'gene': 'Protein', ...}."""
    mapping = {}
    for item in raw.split(","):
        key, sep, value = item.partition("=")
        if not sep or not value.strip():
            raise SourceError(
                f"cannot parse '# columns: {raw}'. Expected comma-separated "
                "role=column pairs, e.g. 'gene=Protein, mutation=Substitution'."
            )
        mapping[key.strip().lower()] = value.strip()
    return mapping


def parse_offsets(raw: str) -> list[int]:
    try:
        offsets = [int(item) for item in raw.split(",") if item.strip()]
    except ValueError as exc:
        raise SourceError(
            f"cannot parse '# offsets: {raw}'. Expected comma-separated integers, e.g. '0, 20'."
        ) from exc
    if not offsets:
        raise SourceError("'# offsets:' is present but empty; omit it to default to 0.")
    return offsets


def load_source(path: Path) -> tuple[dict[str, str], list[dict], list[str]]:
    """(directives, rows, extra column names) from a source list."""
    lines = path.read_text(encoding="utf-8").splitlines()
    directives = read_directives(lines)
    if "numbering" not in directives:
        raise SourceError(
            f"{path}: missing the required '# numbering: <frame>' header. It names the "
            "numbering the positions are written in -- a bundle name, or a label such as n2."
        )

    data = [line for line in lines if line.strip() and not line.lstrip().startswith("#")]
    if not data:
        raise SourceError(f"{path}: no data rows (file is empty or all comments)")

    reader = csv.DictReader(data)
    columns = [c.strip() for c in (reader.fieldnames or [])]
    mapping = parse_column_map(directives.get("columns", "gene=gene, mutation=mutation"))

    missing = [role for role in ("gene", "mutation") if mapping.get(role) not in columns]
    if missing:
        named = ", ".join(f"{role}={mapping.get(role)!r}" for role in missing)
        raise SourceError(
            f"{path}: column(s) {named} named in '# columns:' are not in the file. "
            f"Found: {', '.join(columns) or 'nothing'}."
        )

    rows = [{k.strip(): (v or "").strip() for k, v in row.items() if k} for row in reader]
    extras = [c for c in columns if c not in mapping.values()]
    return {**directives, "_columns": mapping}, rows, extras


def read_db_metadata(root: Path) -> dict[str, str]:
    """`name` and `version` from db.yml. Top-level `key: value` only, which is all it holds."""
    path = root / "db.yml"
    if not path.exists():
        raise SourceError(
            f"{root}: no db.yml. This tool resolves against bundles built from GenBank "
            "records, which carry a db.yml naming and versioning the bundle; the version "
            "is what a generated list records so staleness is visible. Rebuild the bundle "
            "from GenBank, or resolve against one that was."
        )
    meta = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line[:1].isspace() or line.lstrip().startswith("#") or ":" not in line:
            continue
        key, _, value = line.partition(":")
        meta[key.strip()] = value.strip().strip("\"'")

    missing = [k for k in ("name", "version") if not meta.get(k)]
    if missing:
        raise SourceError(f"{path}: missing {' and '.join(missing)}.")
    return meta


def build_gene_index(root: Path) -> dict[str, tuple[str, str]]:
    """gene key -> (segment, protein sequence), discovered from the bundle's own annotations.

    Discovery rather than a gene-to-segment map: the bundle's annotations already say which
    protein lives where, a missing gene is then visibly absent instead of silently wrong, and
    HA1/HA2 come back as their own proteins because read_cds_models promotes mature peptides.
    """
    index: dict[str, tuple[str, str]] = {}
    for segment in SEGMENTS:
        directory = root / "nextclade" / segment
        annotation = directory / "genome_annotation.gff3"
        if not annotation.exists():
            raise SourceError(f"{annotation}: missing. The bundle is incomplete.")
        sequence = str(SeqIO.read(directory / "reference.fasta", "fasta").seq)
        for gene, model in read_cds_models(annotation, segment).items():
            if gene in index:
                raise SourceError(
                    f"{root}: gene {gene!r} is annotated on both {index[gene][0]} and "
                    f"{segment}. One of the two annotations is wrong."
                )
            index[gene] = (segment, translate_model(sequence, model))
    return index


def resolve(row: dict, mapping: dict, extras: list[str], index: dict, offsets: list[int]) -> Resolution:
    """One source row resolved against the bundle, with its status."""
    gene = canonical_gene(row.get(mapping["gene"], ""))
    mutation = row.get(mapping["mutation"], "")
    note = ";".join(
        f"{column}={row[column].replace(';', '/')}"
        for column in extras
        if row.get(column)
    )

    def result(**kwargs) -> Resolution:
        base = dict(
            gene=gene, source_mutation=mutation, source_position=0, ref_aa="", alt_aa="",
            offset=None, target_position=None, target_residue="",
            note=f"native={mutation}" + (f";{note}" if note else ""),
        )
        return Resolution(**{**base, **kwargs})

    parsed = parse_hgvs_p(f"p.{mutation}")
    if parsed is None:
        return result(status="unparseable")
    ref_aa, position, alt_aa = parsed
    if ref_aa == alt_aa:
        return result(status="synonymous", source_position=position,
                      ref_aa=ref_aa, alt_aa=alt_aa)
    if gene not in KNOWN_GENE_KEYS:
        return result(status="unknown_gene", source_position=position,
                      ref_aa=ref_aa, alt_aa=alt_aa)

    fixed = dict(source_position=position, ref_aa=ref_aa, alt_aa=alt_aa)
    if gene not in index:
        return result(status="gene_not_in_reference", **fixed)

    protein = index[gene][1]
    candidates = [
        (offset, position + offset)
        for offset in offsets
        if 1 <= position + offset <= len(protein)
    ]
    if not candidates:
        return result(status="out_of_range", **fixed)

    fits = [(offset, target) for offset, target in candidates if protein[target - 1] == ref_aa]
    if len(fits) > 1:
        return result(status="ambiguous", **fixed)
    if fits:
        offset, target = fits[0]
        status = "ok"
    elif len(offsets) == 1:
        # One declared frame and nothing to disambiguate with, so the position is taken at
        # its word and the disagreement reported: a lineage's wild type legitimately differs.
        offset, target = candidates[0]
        status = "ref_differs"
    else:
        # Several candidate offsets and no residue agreement: nothing picks between them.
        return result(status="no_fit", **fixed)

    residue = protein[target - 1]
    if residue == alt_aa:
        # The bundle already carries the ALT, so a sample matching this entry differs from
        # the reference nowhere and produces no VCF record. Undetectable here, not a hit.
        status = "ref_is_alt"
    return result(status=status, offset=offset, target_position=target,
                  target_residue=residue, **fixed)


def deduplicate(resolutions: list[Resolution]) -> list[Resolution]:
    """One entry per (gene, position, ALT), merging the provenance of the rest.

    The watcher collapses repeated keys itself via setdefault and loses the second note,
    so a list carrying duplicates would silently discard what the duplicate said.
    """
    merged: dict[tuple[str, int, str], Resolution] = {}
    duplicates: dict[tuple[str, int, str], list[str]] = {}
    for resolution in resolutions:
        if resolution.key in merged:
            duplicates.setdefault(resolution.key, []).append(resolution.source_mutation)
        else:
            merged[resolution.key] = resolution

    for key, others in duplicates.items():
        kept = merged[key]
        merged[key] = Resolution(
            **{**vars(kept), "note": f"{kept.note};also_in={'/'.join(others)}"}
        )
    return list(merged.values())


def write_watchlist(path: Path, entries: list[Resolution], header: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        handle.writelines(f"# {line}\n" for line in header)
        writer = csv.DictWriter(handle, fieldnames=("gene", "mutation", "note"))
        writer.writeheader()
        writer.writerows(
            {"gene": e.gene, "mutation": e.mutation, "note": e.note} for e in entries
        )


def write_report(path: Path, resolutions: list[Resolution]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(
            {k: ("" if v is None else v) for k, v in vars(r).items()} for r in resolutions
        )


def get_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}")
    p.add_argument("--source", required=True, type=Path, help="Source watchlist CSV")
    p.add_argument("--db", required=True, type=Path, help="Reference bundle directory")
    p.add_argument("--out-dir", required=True, type=Path, help="Where to write both outputs")
    return p.parse_args()


def main() -> int:
    args = get_args()

    try:
        directives, rows, extras = load_source(args.source)
        meta = read_db_metadata(args.db)
        index = build_gene_index(args.db)
    except SourceError as exc:
        log.error("%s", exc)
        return 1

    offsets = parse_offsets(directives["offsets"]) if "offsets" in directives else [0]
    mapping = directives["_columns"]
    resolutions = [resolve(row, mapping, extras, index, offsets) for row in rows]

    included = deduplicate([r for r in resolutions if r.status in INCLUDED_STATUSES])
    args.out_dir.mkdir(parents=True, exist_ok=True)
    stem = args.source.stem.removesuffix("_source")
    watchlist_path = args.out_dir / f"{stem}_{meta['name']}_watchlist.csv"
    report_path = args.out_dir / f"{stem}_{meta['name']}_reference_map.tsv"

    header = [
        f"numbering: {meta['name']}",
        f"db-version: {meta['version']}",
        *([f"source: {directives['source']}"] if "source" in directives else []),
        f"resolved-from: {args.source.name} (# numbering: {directives['numbering']})",
        f"generated: {date.today().isoformat()} by scripts/resolve_watchlist.py "
        f"--source {args.source.name} --db {args.db}",
    ]
    write_watchlist(watchlist_path, included, header)
    write_report(report_path, resolutions)

    counts: dict[str, int] = {}
    for resolution in resolutions:
        counts[resolution.status] = counts.get(resolution.status, 0) + 1
    for status in sorted(counts):
        level = log.info if status in INCLUDED_STATUSES else log.warning
        level("%s: %d", status, counts[status])

    log.info(
        "%d of %d row(s) resolved onto %s %s -> %s",
        len(included), len(rows), meta["name"], meta["version"], watchlist_path,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
