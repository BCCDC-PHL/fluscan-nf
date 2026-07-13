#!/usr/bin/env python3
"""Build a nucleotide/amino-acid coordinate map from a source reference onto a
target reference coordinate system.

The map lifts flumut's fixed reference (the "source", the simpler/less
reference-entangled representation) onto the selected Nextclade dataset (the
"target", the canonical coordinate system). A downstream step consumes the map
to remap flumut mutations into Nextclade coordinates before merging them into
the pipeline's VCF.

This script's only job is to produce a reusable, auditable coordinate map
between two FluScan-style reference datasets. It is idempotent -- given the
same two reference directories it always produces the same map, including the
trivial case where the source and target references are identical (a strict
1:1 map, e.g. when the Nextclade dataset is derived from flumut's reference).

Reference directory layout (both source and target must follow this):

    reference_dir/
      ha/reference.fasta
      ha/genome_annotation.gff3
      m/reference.fasta
      m/genome_annotation.gff3
      ...

An optional top-level config.json may define segment directories and an
id_prefix used as the reference's display name; see read_reference_segments().

The eight influenza A segments are fixed:

    PB2, PB1, PA, HA, NP, NA, M, NS

Usage:
    build_reference_coordinate_map.py \\
        --source-ref-dir  /path/to/flumut/nextclade \\
        --target-ref-dir  /path/to/chicken2021/nextclade \\
        --out-dir         reference_remap \\
        --source-name     flumut \\
        --target-name     chicken2021
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import unquote

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

# Reuse the converter's GFF parser so the serialized target-CDS placement is exactly what
# the downstream insert step would otherwise compute itself (guarantees identical ANN).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nextclade_to_freebayes_vcf import load_gff3  # noqa: E402

SCRIPT_VERSION = "2.1.0"
SCHEMA_VERSION = "2.1"

SEGMENTS = ("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
SEGMENT_ALIASES = {segment.lower(): segment for segment in SEGMENTS}

NT_ALIGNMENT_PARAMS = {
    "mode": "global",
    "match_score": 2,
    "mismatch_score": -1,
    "open_gap_score": -10,
    "extend_gap_score": -0.5,
}
AA_ALIGNMENT_PARAMS = {
    "mode": "global",
    "substitution_matrix": "BLOSUM62",
    "open_gap_score": -11,
    "extend_gap_score": -1,
}

# Below these thresholds a map is still written, but flagged with a QC warning.
IDENTITY_WARNING_THRESHOLD_PCT = 80.0
COVERAGE_WARNING_THRESHOLD_PCT = 95.0


@dataclass(frozen=True)
class CdsModel:
    segment: str
    gene: str
    contig: str
    strand: str
    parts: tuple[tuple[int, int], ...]
    cds_ids: tuple[str, ...]

    def coding_positions(self) -> list[int]:
        positions: list[int] = []
        for start, end in sorted(self.parts):
            positions.extend(range(start, end + 1))
        return positions


@dataclass
class SegmentBuild:
    segment: str
    source_chrom: str
    target_chrom: str
    source_seq: str
    target_seq: str
    nt_target_aligned: str
    nt_source_aligned: str
    nt_map: dict[int, int | None]
    nt_status: dict[int, str]
    nt_identity: int
    nt_aligned: int
    nt_gaps: int
    aa_maps: dict[str, dict[int, int | None]] = field(default_factory=dict)
    aa_status: dict[str, dict[int, str]] = field(default_factory=dict)
    aa_identity: dict[str, int] = field(default_factory=dict)
    aa_aligned: dict[str, int] = field(default_factory=dict)
    aa_gaps: dict[str, int] = field(default_factory=dict)
    aa_source: dict[str, str] = field(default_factory=dict)
    aa_target: dict[str, str] = field(default_factory=dict)
    aa_target_aligned: dict[str, str] = field(default_factory=dict)
    aa_source_aligned: dict[str, str] = field(default_factory=dict)
    all_source_gene_lengths: dict[str, int] = field(default_factory=dict)
    all_target_gene_lengths: dict[str, int] = field(default_factory=dict)
    missing_target_genes: tuple[str, ...] = ()
    extra_source_genes: tuple[str, ...] = ()
    stale_annotation: bool = False


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--source-ref-dir",
        required=True,
        type=Path,
        help="Source reference directory to lift FROM (flumut's Nextclade dataset: "
        "per-segment reference.fasta + genome_annotation.gff3).",
    )
    parser.add_argument(
        "--target-ref-dir",
        required=True,
        type=Path,
        help="Target reference directory to lift ONTO (the selected Nextclade dataset: "
        "per-segment reference.fasta + genome_annotation.gff3).",
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--source-name",
        default=None,
        help="Preferred source reference name. Overridden by config.json id_prefix if present.",
    )
    parser.add_argument(
        "--target-name",
        default=None,
        help="Preferred target reference name. Overridden by config.json id_prefix if present.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv if argv is not None else sys.argv[1:])
    build_map(args)
    return 0


def build_map(args: argparse.Namespace) -> None:
    source_ref_dir = args.source_ref_dir.resolve()
    target_ref_dir = args.target_ref_dir.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    source_name = resolve_reference_name(source_ref_dir, args.source_name)
    target_name = resolve_reference_name(target_ref_dir, args.target_name)

    source_segments = read_reference_segments(source_ref_dir)
    source_models = read_reference_models(source_ref_dir)
    target_segments = read_reference_segments(target_ref_dir)
    target_models = read_reference_models(target_ref_dir)

    missing = [
        segment
        for segment in SEGMENTS
        if segment not in source_segments or segment not in target_segments
    ]
    if missing:
        raise SystemExit(f"Missing required segment(s): {', '.join(missing)}")

    builds: dict[str, SegmentBuild] = {}
    for segment in SEGMENTS:
        builds[segment] = build_segment_map(
            segment=segment,
            source_record=source_segments[segment],
            target_record=target_segments[segment],
            source_models=source_models.get(segment, {}),
            target_models=target_models.get(segment, {}),
        )

    target_placements = read_target_placements(target_ref_dir)

    map_doc = make_map_doc(
        source_name=source_name,
        target_name=target_name,
        builds=builds,
        target_placements=target_placements,
    )
    write_json(out_dir / "coordinate_map.json", map_doc)
    shard_names = write_segment_shards(out_dir, map_doc)
    write_nt_tsv(out_dir / "coordinate_map.nt.tsv", source_name, builds)
    write_aa_tsv(out_dir / "coordinate_map.aa.tsv", source_name, builds)
    write_qc_tsv(out_dir / "coordinate_map.qc.tsv", source_name, builds)
    write_nt_alignment_fasta(out_dir / f"{source_name}_alignment.fasta", builds)
    write_aa_alignment_fasta(out_dir / f"{source_name}_aa_alignment.fasta", builds)
    warn_low_quality_maps(builds)

    for produced in (
        "coordinate_map.json",
        *shard_names,
        "coordinate_map.nt.tsv",
        "coordinate_map.aa.tsv",
        "coordinate_map.qc.tsv",
        f"{source_name}_alignment.fasta",
        f"{source_name}_aa_alignment.fasta",
    ):
        print(f"Wrote {out_dir / produced}")


def resolve_reference_name(ref_dir: Path, explicit_name: str | None) -> str:
    """Resolve a stable, human-meaningful name for a reference dataset.

    Priority: config.json id_prefix/display_name (curated, embedded in the
    dataset) > explicit --source-name/--target-name argument > parent
    directory name (last resort, e.g. .../gsgd/nextclade -> "gsgd").
    """
    config_path = ref_dir / "config.json"
    if config_path.exists():
        with config_path.open() as handle:
            config = json.load(handle)
        if config.get("id_prefix"):
            return str(config["id_prefix"])
        if config.get("display_name"):
            return str(config["display_name"])
    if explicit_name:
        return explicit_name
    return ref_dir.parent.name or ref_dir.name


def read_reference_segments(reference_dir: Path) -> dict[str, Any]:
    records: dict[str, Any] = {}
    config_path = reference_dir / "config.json"
    if config_path.exists():
        with config_path.open() as handle:
            config = json.load(handle)
        for item in config.get("segments", []):
            segment = normalize_segment(item.get("name") or item.get("dir"))
            if not segment:
                continue
            fasta = reference_dir / str(item["dir"]) / "reference.fasta"
            records[segment] = read_single_fasta(fasta)
    else:
        for segment_dir in sorted(reference_dir.iterdir()):
            if not segment_dir.is_dir():
                continue
            segment = normalize_segment(segment_dir.name)
            fasta = segment_dir / "reference.fasta"
            if segment and fasta.exists():
                records[segment] = read_single_fasta(fasta)
    if not records and (reference_dir / "blastn_db.fa").exists():
        records = read_combined_fasta_segments(reference_dir / "blastn_db.fa")
    return records


def read_combined_fasta_segments(combined_fasta: Path) -> dict[str, Any]:
    records: dict[str, Any] = {}
    for record in SeqIO.parse(str(combined_fasta), "fasta"):
        segment = segment_from_record(record)
        if segment:
            records[segment] = record
    return records


def read_reference_models(reference_dir: Path) -> dict[str, dict[str, CdsModel]]:
    models: dict[str, dict[str, CdsModel]] = {}
    for segment in SEGMENTS:
        gff = reference_dir / segment.lower() / "genome_annotation.gff3"
        if not gff.exists():
            gff = reference_dir / f"{segment.lower()}.gff3"
        if gff.exists():
            models[segment] = read_cds_models(gff, segment)
    return models


def read_single_fasta(path: Path) -> Any:
    try:
        return SeqIO.read(str(path), "fasta")
    except Exception as exc:
        raise SystemExit(f"Could not read single FASTA record from {path}: {exc}") from exc


def normalize_segment(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return SEGMENT_ALIASES.get(text.lower())


def segment_from_record(record: Any) -> str | None:
    for text in (record.id, record.name, record.description):
        parts = re.split(r"[|\s_/-]+", str(text).upper())
        for part in reversed(parts):
            if part in SEGMENTS:
                return part
    return None


def parse_gff_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in raw.split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = unquote(value)
    return attrs


def read_cds_models(path: Path, segment: str) -> dict[str, CdsModel]:
    grouped: dict[str, dict[str, Any]] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = parse_gff_attributes(fields[8])
            gene = attrs.get("gene") or attrs.get("Name") or attrs.get("ID")
            if not gene:
                continue
            grouped.setdefault(
                gene,
                {
                    "segment": segment,
                    "gene": gene,
                    "contig": fields[0],
                    "strand": fields[6],
                    "parts": [],
                    "cds_ids": set(),
                },
            )
            grouped[gene]["parts"].append((int(fields[3]), int(fields[4])))
            if attrs.get("ID"):
                grouped[gene]["cds_ids"].add(attrs["ID"])

    models: dict[str, CdsModel] = {}
    for gene, values in grouped.items():
        models[gene] = CdsModel(
            segment=values["segment"],
            gene=values["gene"],
            contig=values["contig"],
            strand=values["strand"],
            parts=tuple(sorted(values["parts"])),
            cds_ids=tuple(sorted(values["cds_ids"])),
        )
    return models


def detect_stale_ha_annotation(
    segment: str, source_gene_keys: set[str], target_gene_keys: set[str]
) -> bool:
    """Flag legacy single-CDS "HA" annotations on the target (Nextclade) dataset.

    The current accepted standard splits HA into HA1/HA2 CDS features (as
    flumut's source reference does). The target (the selected Nextclade dataset)
    is the flexible, user-chosen reference; if it still annotates a single "HA"
    gene (and has not also been updated to expose HA1/HA2) it is considered out
    of date: no HA amino-acid map can be built, because the flumut source uses
    the HA1/HA2 split and the shared-gene intersection is empty.
    """
    if segment != "HA" or "HA" not in target_gene_keys:
        return False
    if {"HA1", "HA2"} <= target_gene_keys:
        return False  # target already exposes both -- not stale
    if not ({"HA1", "HA2"} <= source_gene_keys):
        return False  # source itself doesn't require the split
    return True


def build_segment_map(
    segment: str,
    source_record: Any,
    target_record: Any,
    source_models: dict[str, CdsModel],
    target_models: dict[str, CdsModel],
) -> SegmentBuild:
    source_seq = str(source_record.seq).upper()
    target_seq = str(target_record.seq).upper()
    (
        nt_map,
        nt_status,
        nt_identity,
        nt_aligned,
        nt_gaps,
        nt_target_aligned,
        nt_source_aligned,
    ) = align_position_map(target_seq=target_seq, query_seq=source_seq, kind="nt")

    build = SegmentBuild(
        segment=segment,
        source_chrom=source_record.id,
        target_chrom=target_record.id,
        source_seq=source_seq,
        target_seq=target_seq,
        nt_target_aligned=nt_target_aligned,
        nt_source_aligned=nt_source_aligned,
        nt_map=nt_map,
        nt_status=nt_status,
        nt_identity=nt_identity,
        nt_aligned=nt_aligned,
        nt_gaps=nt_gaps,
    )

    source_gene_keys = set(source_models)
    target_gene_keys = set(target_models)
    build.missing_target_genes = tuple(sorted(target_gene_keys - source_gene_keys))
    build.extra_source_genes = tuple(sorted(source_gene_keys - target_gene_keys))
    build.stale_annotation = detect_stale_ha_annotation(
        segment, source_gene_keys, target_gene_keys
    )

    for gene, model in source_models.items():
        build.all_source_gene_lengths[gene] = len(translate_model(source_seq, model))
    for gene, model in target_models.items():
        build.all_target_gene_lengths[gene] = len(translate_model(target_seq, model))

    for gene in sorted(source_gene_keys & target_gene_keys):
        source_protein = translate_model(source_seq, source_models[gene])
        target_protein = translate_model(target_seq, target_models[gene])
        (
            gene_map,
            gene_status,
            ident,
            aligned,
            gaps,
            gene_target_aligned,
            gene_source_aligned,
        ) = align_position_map(target_seq=target_protein, query_seq=source_protein, kind="aa")
        build.aa_maps[gene] = gene_map
        build.aa_status[gene] = gene_status
        build.aa_identity[gene] = ident
        build.aa_aligned[gene] = aligned
        build.aa_gaps[gene] = gaps
        build.aa_source[gene] = source_protein
        build.aa_target[gene] = target_protein
        build.aa_target_aligned[gene] = gene_target_aligned
        build.aa_source_aligned[gene] = gene_source_aligned

    return build


def make_aligner(kind: str) -> PairwiseAligner:
    if kind == "aa":
        aligner = PairwiseAligner(mode="global")
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = AA_ALIGNMENT_PARAMS["open_gap_score"]
        aligner.extend_gap_score = AA_ALIGNMENT_PARAMS["extend_gap_score"]
    else:
        aligner = PairwiseAligner(
            mode="global",
            match_score=NT_ALIGNMENT_PARAMS["match_score"],
            mismatch_score=NT_ALIGNMENT_PARAMS["mismatch_score"],
            open_gap_score=NT_ALIGNMENT_PARAMS["open_gap_score"],
            extend_gap_score=NT_ALIGNMENT_PARAMS["extend_gap_score"],
        )
    return aligner


def align_position_map(
    target_seq: str,
    query_seq: str,
    kind: str,
) -> tuple[dict[int, int | None], dict[int, str], int, int, int, str, str]:
    aligner = make_aligner(kind)
    alignment = next(iter(aligner.align(target_seq, query_seq)))
    target_aligned = str(alignment[0, :])
    query_aligned = str(alignment[1, :])

    query_to_target: dict[int, int | None] = {}
    status: dict[int, str] = {}
    target_pos = 0
    query_pos = 0
    identities = 0
    aligned = 0
    gaps = 0

    for target_char, query_char in zip(target_aligned, query_aligned, strict=True):
        target_has_base = target_char != "-"
        query_has_base = query_char != "-"
        if target_has_base:
            target_pos += 1
        if query_has_base:
            query_pos += 1

        if target_has_base and query_has_base:
            query_to_target[query_pos] = target_pos
            status[query_pos] = "mapped"
            aligned += 1
            if target_char == query_char:
                identities += 1
        elif query_has_base:
            query_to_target[query_pos] = None
            status[query_pos] = "unmapped"
            gaps += 1
        elif target_has_base:
            gaps += 1

    return query_to_target, status, identities, aligned, gaps, target_aligned, query_aligned


def translate_model(sequence: str, model: CdsModel) -> str:
    positions = model.coding_positions()
    coding = "".join(sequence[pos - 1] for pos in positions)
    if model.strand == "-":
        coding = str(Seq(coding).reverse_complement())
    usable = coding[: len(coding) - (len(coding) % 3)]
    protein = str(Seq(usable).translate(to_stop=False))
    if protein.endswith("*"):
        protein = protein[:-1]
    return protein


def safe_pct(numerator: int, denominator: int) -> float:
    if denominator == 0:
        return math.nan
    return numerator / denominator * 100.0


def compute_status(identity_pct: float, coverage_pct: float) -> str:
    if math.isnan(identity_pct) or math.isnan(coverage_pct):
        return "warning"
    if identity_pct < IDENTITY_WARNING_THRESHOLD_PCT or coverage_pct < COVERAGE_WARNING_THRESHOLD_PCT:
        return "warning"
    return "ok"


def index_transcripts(transcripts: list) -> dict[str, dict]:
    """Map gene identifiers -> transcript, preferring main (non-mature) CDS transcripts.

    Mirrors insert_flumut_mutations.index_transcripts so that lookups by coordinate-map
    gene key resolve to the same transcript on both sides.
    """
    by_key: dict[str, dict] = {}
    for tx in transcripts:
        for key in (tx.get("gene_name"), tx.get("feature_id"), tx.get("protein_id")):
            if not key:
                continue
            existing = by_key.get(key)
            if existing is None or (existing.get("is_mature") and not tx.get("is_mature")):
                by_key[key] = tx
    return by_key


def read_target_placements(target_ref_dir: Path) -> dict[str, dict[str, dict]]:
    """Per-segment {gene_key -> transcript} from the target GFFs, via the converter's
    load_gff3, so serialized placement matches the insert step's own computation."""
    placements: dict[str, dict[str, dict]] = {}
    for segment in SEGMENTS:
        gff = target_ref_dir / segment.lower() / "genome_annotation.gff3"
        if not gff.exists():
            gff = target_ref_dir / f"{segment.lower()}.gff3"
        if gff.exists():
            _, transcripts = load_gff3(str(gff))
            placements[segment] = index_transcripts(transcripts)
    return placements


def target_cds_block(target_chrom: str, tx: dict) -> dict[str, Any]:
    """Serialize a target transcript's placement + identifiers for the coordinate map.

    Carries exactly the fields the insert step needs to derive a genomic POS and build the
    ANN entry without re-parsing the GFF: 0-based half-open CDS segments (as load_gff3
    yields), strand, lengths, and the GFF-derived gene_id/feature_id/protein_id.
    """
    return {
        "chrom": target_chrom,
        "strand": tx["strand"],
        "segments": [list(seg) for seg in tx["segments"]],
        "cds_length": tx["cds_length"],
        "aa_length": tx["aa_length"],
        "gene_id": tx["gene_id"],
        "feature_id": tx["feature_id"],
        "protein_id": tx["protein_id"],
    }


def make_map_doc(
    source_name: str,
    target_name: str,
    builds: dict[str, SegmentBuild],
    target_placements: dict[str, dict[str, dict]],
) -> dict[str, Any]:
    segments: dict[str, Any] = {}
    for segment, build in builds.items():
        placements = target_placements.get(segment, {})
        genes: dict[str, Any] = {}
        for gene in sorted(build.aa_maps):
            gene_doc: dict[str, Any] = {
                "source_length": len(build.aa_source[gene]),
                "target_length": len(build.aa_target[gene]),
                "aa_identity_pct": safe_pct(build.aa_identity[gene], build.aa_aligned[gene]),
                "aa_aligned_positions": build.aa_aligned[gene],
                "aa_gap_columns": build.aa_gaps[gene],
                "aa_source_to_target": {
                    str(pos): value
                    for pos, value in sorted(build.aa_maps[gene].items())
                },
            }
            tx = placements.get(gene)
            if tx is not None:
                gene_doc["target_cds"] = target_cds_block(build.target_chrom, tx)
            genes[gene] = gene_doc

        segments[segment] = {
            "source_chrom": build.source_chrom,
            "target_chrom": build.target_chrom,
            "source_length": len(build.source_seq),
            "target_length": len(build.target_seq),
            "nt_identity_pct": safe_pct(build.nt_identity, build.nt_aligned),
            "nt_aligned_positions": build.nt_aligned,
            "nt_gap_columns": build.nt_gaps,
            "nt_source_to_target": {
                str(pos): value for pos, value in sorted(build.nt_map.items())
            },
            "missing_target_genes": list(build.missing_target_genes),
            "extra_source_genes": list(build.extra_source_genes),
            "stale_annotation": build.stale_annotation,
            "genes": genes,
        }

    return {
        "schema_version": SCHEMA_VERSION,
        "script_version": SCRIPT_VERSION,
        "source_reference": source_name,
        "target_reference": target_name,
        "alignment_parameters": {
            "nt": NT_ALIGNMENT_PARAMS,
            "aa": AA_ALIGNMENT_PARAMS,
        },
        "segments": segments,
    }


def write_json(path: Path, value: dict[str, Any]) -> None:
    with path.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_segment_shards(out_dir: Path, map_doc: dict[str, Any]) -> list[str]:
    """Write one slim per-segment shard (coordinate_map.<SEG>.json) alongside the combined
    map. Shards keep the same segments[SEG] nesting so consumers work unchanged, but drop
    the large, downstream-unused nt_source_to_target map. Returns the shard filenames."""
    names: list[str] = []
    for segment, seg_doc in map_doc["segments"].items():
        slim_seg = {k: v for k, v in seg_doc.items() if k != "nt_source_to_target"}
        shard = {
            "schema_version": map_doc["schema_version"],
            "script_version": map_doc["script_version"],
            "source_reference": map_doc["source_reference"],
            "target_reference": map_doc["target_reference"],
            "segments": {segment: slim_seg},
        }
        name = f"coordinate_map.{segment}.json"
        write_json(out_dir / name, shard)
        names.append(name)
    return names


def write_nt_tsv(path: Path, source_name: str, builds: dict[str, SegmentBuild]) -> None:
    fieldnames = (
        "source_reference",
        "segment",
        "source_chrom",
        "source_pos",
        "source_base",
        "target_chrom",
        "target_pos",
        "target_base",
        "status",
    )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for segment in SEGMENTS:
            build = builds[segment]
            for source_pos, target_pos in sorted(build.nt_map.items()):
                writer.writerow(
                    {
                        "source_reference": source_name,
                        "segment": segment,
                        "source_chrom": build.source_chrom,
                        "source_pos": source_pos,
                        "source_base": build.source_seq[source_pos - 1],
                        "target_chrom": build.target_chrom,
                        "target_pos": target_pos or "",
                        "target_base": build.target_seq[target_pos - 1] if target_pos else "",
                        "status": build.nt_status[source_pos],
                    }
                )


def write_aa_tsv(path: Path, source_name: str, builds: dict[str, SegmentBuild]) -> None:
    fieldnames = (
        "source_reference",
        "segment",
        "gene",
        "source_aa_pos",
        "source_residue",
        "target_aa_pos",
        "target_residue",
        "status",
    )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for segment in SEGMENTS:
            build = builds[segment]
            for gene in sorted(build.aa_maps):
                for source_pos, target_pos in sorted(build.aa_maps[gene].items()):
                    writer.writerow(
                        {
                            "source_reference": source_name,
                            "segment": segment,
                            "gene": gene,
                            "source_aa_pos": source_pos,
                            "source_residue": build.aa_source[gene][source_pos - 1],
                            "target_aa_pos": target_pos or "",
                            "target_residue": (
                                build.aa_target[gene][target_pos - 1] if target_pos else ""
                            ),
                            "status": build.aa_status[gene][source_pos],
                        }
                    )


def write_qc_tsv(path: Path, source_name: str, builds: dict[str, SegmentBuild]) -> None:
    fieldnames = (
        "source_reference",
        "segment",
        "layer",
        "gene",
        "source_length",
        "target_length",
        "aligned_positions",
        "gap_columns",
        "identity_pct",
        "status",
    )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for segment in SEGMENTS:
            build = builds[segment]
            nt_identity_pct = safe_pct(build.nt_identity, build.nt_aligned)
            nt_coverage_pct = safe_pct(build.nt_aligned, len(build.source_seq))
            writer.writerow(
                {
                    "source_reference": source_name,
                    "segment": segment,
                    "layer": "nt",
                    "gene": "",
                    "source_length": len(build.source_seq),
                    "target_length": len(build.target_seq),
                    "aligned_positions": build.nt_aligned,
                    "gap_columns": build.nt_gaps,
                    "identity_pct": f"{nt_identity_pct:.4f}",
                    "status": compute_status(nt_identity_pct, nt_coverage_pct),
                }
            )
            for gene in sorted(build.aa_maps):
                aa_identity_pct = safe_pct(build.aa_identity[gene], build.aa_aligned[gene])
                aa_coverage_pct = safe_pct(build.aa_aligned[gene], len(build.aa_source[gene]))
                writer.writerow(
                    {
                        "source_reference": source_name,
                        "segment": segment,
                        "layer": "aa",
                        "gene": gene,
                        "source_length": len(build.aa_source[gene]),
                        "target_length": len(build.aa_target[gene]),
                        "aligned_positions": build.aa_aligned[gene],
                        "gap_columns": build.aa_gaps[gene],
                        "identity_pct": f"{aa_identity_pct:.4f}",
                        "status": compute_status(aa_identity_pct, aa_coverage_pct),
                    }
                )
            for gene in build.missing_target_genes:
                status = "stale_annotation" if build.stale_annotation else "missing_in_source"
                writer.writerow(
                    {
                        "source_reference": source_name,
                        "segment": segment,
                        "layer": "aa",
                        "gene": gene,
                        "source_length": "",
                        "target_length": build.all_target_gene_lengths.get(gene, ""),
                        "aligned_positions": 0,
                        "gap_columns": "",
                        "identity_pct": "",
                        "status": status,
                    }
                )
            for gene in build.extra_source_genes:
                writer.writerow(
                    {
                        "source_reference": source_name,
                        "segment": segment,
                        "layer": "aa",
                        "gene": gene,
                        "source_length": build.all_source_gene_lengths.get(gene, ""),
                        "target_length": "",
                        "aligned_positions": 0,
                        "gap_columns": "",
                        "identity_pct": "",
                        "status": "extra_in_source",
                    }
                )


def write_alignment_block(handle: Any, header: str, target_aligned: str, source_aligned: str) -> None:
    handle.write(f">{header}\n")
    handle.write(f"{target_aligned}\n")
    handle.write(f"{source_aligned}\n")


def write_nt_alignment_fasta(path: Path, builds: dict[str, SegmentBuild]) -> None:
    with path.open("w") as handle:
        for segment in SEGMENTS:
            build = builds[segment]
            header = (
                f"segment={segment};alphabet=nt;target={build.target_chrom};"
                f"source={build.source_chrom}"
            )
            write_alignment_block(handle, header, build.nt_target_aligned, build.nt_source_aligned)


def write_aa_alignment_fasta(path: Path, builds: dict[str, SegmentBuild]) -> None:
    with path.open("w") as handle:
        for segment in SEGMENTS:
            build = builds[segment]
            for gene in sorted(build.aa_maps):
                header = (
                    f"segment={segment};gene={gene};alphabet=aa;target={build.target_chrom};"
                    f"source={build.source_chrom}"
                )
                write_alignment_block(
                    handle,
                    header,
                    build.aa_target_aligned[gene],
                    build.aa_source_aligned[gene],
                )


def warn_low_quality_maps(builds: dict[str, SegmentBuild]) -> None:
    for segment, build in builds.items():
        nt_identity = safe_pct(build.nt_identity, build.nt_aligned)
        nt_coverage = safe_pct(build.nt_aligned, len(build.source_seq))
        if nt_identity < IDENTITY_WARNING_THRESHOLD_PCT or nt_coverage < COVERAGE_WARNING_THRESHOLD_PCT:
            print(
                f"WARNING: {segment} nucleotide map has identity={nt_identity:.2f}% "
                f"and source coverage={nt_coverage:.2f}%",
                file=sys.stderr,
            )
        if build.stale_annotation:
            print(
                f"WARNING: {segment} target (Nextclade) annotation appears out of date "
                f"(single 'HA' CDS found; the flumut source reference requires the "
                f"HA1/HA2 split). No HA amino-acid map can be built for this target dataset.",
                file=sys.stderr,
            )
        elif build.missing_target_genes:
            print(
                f"WARNING: {segment} target gene(s) not found in source annotation: "
                f"{', '.join(build.missing_target_genes)}",
                file=sys.stderr,
            )
        for gene in sorted(build.aa_maps):
            aa_identity = safe_pct(build.aa_identity[gene], build.aa_aligned[gene])
            aa_coverage = safe_pct(build.aa_aligned[gene], len(build.aa_source[gene]))
            if aa_identity < IDENTITY_WARNING_THRESHOLD_PCT or aa_coverage < COVERAGE_WARNING_THRESHOLD_PCT:
                print(
                    f"WARNING: {segment}/{gene} amino-acid map has "
                    f"identity={aa_identity:.2f}% and source coverage={aa_coverage:.2f}%",
                    file=sys.stderr,
                )


if __name__ == "__main__":
    raise SystemExit(main())
