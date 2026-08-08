#!/usr/bin/env python3
"""Map nucleotide and protein coordinates between two influenza references.

Each reference contains uppercase segment directories with one
``reference.fasta`` and ``genome_annotation.gff3`` per segment.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import unquote

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

from lib.flu_gene_names import canonical_gene  # noqa: E402,F401  (re-exported for callers)

SCRIPT_VERSION = "3.0.0"
SCHEMA_VERSION = "2.2"

SEGMENTS = ("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

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

# Below these thresholds a map is still written, but flagged in the QC table and on stderr.
# Real reference pairs legitimately sit here -- flumut -> gsgd reports NS/NS1 at 66% identity
# and PB1/PB1-F2 at 73% -- so this is a reviewer's signal, not a rejection.
IDENTITY_WARNING_THRESHOLD_PCT = 80.0
COVERAGE_WARNING_THRESHOLD_PCT = 95.0

# Co-optimal alignment counts are only ever read as "was this choice arbitrary, and how
# badly"; the exact figure can be astronomically large, so it is capped for reporting.
CO_OPTIMAL_REPORTING_CAP = 1000


@dataclass(frozen=True)
class CdsModel:
    segment: str
    gene: str
    contig: str
    strand: str
    parts: tuple[tuple[int, int], ...]
    # Set only for models recovered from mature_protein_region_of_CDS rows. `aliases` holds
    # the raw (uncanonicalised) identifiers such a row can surface in a VCF's ANN Feature_ID
    # -- its protein_id and its GFF3 ID -- so the merge can key on them when the ANN
    # Gene_Name names the parent gene instead of the subunit. `subunit_parent` is the
    # canonical gene key of that parent CDS.
    aliases: tuple[str, ...] = ()
    is_mature: bool = False
    subunit_parent: str | None = None

    def coding_positions(self) -> list[int]:
        positions: list[int] = []
        for start, end in sorted(self.parts):
            positions.extend(range(start, end + 1))
        return positions


@dataclass(frozen=True)
class AlignmentResult:
    """One pairwise alignment, and the source-to-target position map read off it.

    Nucleotide and protein layers use this identically; `source_seq`/`target_seq` are the
    ungapped sequences that were aligned, so a caller needs nothing else to report a
    residue, a length, or a quality figure.
    """

    source_seq: str
    target_seq: str
    source_aligned: str
    target_aligned: str
    # A source position maps to a target position, or to None where the target has a gap.
    # That value is the whole status: there is no second dictionary to keep in step with it.
    source_to_target: dict[int, int | None]
    identities: int
    aligned: int
    gaps: int
    co_optimal: int

    @property
    def identity_pct(self) -> float:
        return safe_pct(self.identities, self.aligned)

    @property
    def coverage_pct(self) -> float:
        return safe_pct(self.aligned, len(self.source_seq))

    @property
    def is_questionable(self) -> bool:
        """Whether this map is poor enough to warn about.

        NaN must count as failing. `safe_pct` returns NaN when nothing aligned -- the worst
        possible outcome -- and every `nan < threshold` comparison is False, so a plain
        comparison stays silent on exactly the case that most needs reporting.
        """
        identity, coverage = self.identity_pct, self.coverage_pct
        return (
            math.isnan(identity)
            or math.isnan(coverage)
            or identity < IDENTITY_WARNING_THRESHOLD_PCT
            or coverage < COVERAGE_WARNING_THRESHOLD_PCT
        )

    @property
    def status(self) -> str:
        return "warning" if self.is_questionable else "ok"


@dataclass
class SegmentBuild:
    segment: str
    source_chrom: str
    target_chrom: str
    nt: AlignmentResult
    aa: dict[str, AlignmentResult] = field(default_factory=dict)
    missing_target_genes: tuple[str, ...] = ()
    extra_source_genes: tuple[str, ...] = ()
    # Translated lengths for exactly the genes above -- the ones with no alignment of their
    # own to read a length from. The two gene sets are disjoint, so one dictionary serves.
    unmapped_gene_lengths: dict[str, int] = field(default_factory=dict)
    stale_annotation: bool = False
    gene_aliases: dict[str, tuple[str, ...]] = field(default_factory=dict)
    subunit_parent_genes: tuple[str, ...] = ()


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
        required=True,
        help="Source reference name written to output metadata and filenames.",
    )
    parser.add_argument(
        "--target-name",
        required=True,
        help="Target reference name written to output metadata.",
    )
    # The process's provenance block reads this rather than repeating the literal, so the
    # recorded tool_version cannot drift from SCRIPT_VERSION.
    parser.add_argument("--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv if argv is not None else sys.argv[1:])
    build_map(args)
    return 0


def build_map(args: argparse.Namespace) -> None:
    source_ref_dir = args.source_ref_dir
    target_ref_dir = args.target_ref_dir
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    builds: dict[str, SegmentBuild] = {}
    for segment in SEGMENTS:
        source_dir = source_ref_dir / segment
        target_dir = target_ref_dir / segment
        builds[segment] = build_segment_map(
            segment=segment,
            source_record=SeqIO.read(source_dir / "reference.fasta", "fasta"),
            target_record=SeqIO.read(target_dir / "reference.fasta", "fasta"),
            source_models=read_cds_models(source_dir / "genome_annotation.gff3", segment),
            target_models=read_cds_models(target_dir / "genome_annotation.gff3", segment),
        )

    map_doc = make_map_doc(
        source_name=args.source_name,
        target_name=args.target_name,
        builds=builds,
    )
    write_json(out_dir / "coordinate_map.json", map_doc)
    shard_names = write_segment_shards(out_dir, map_doc)
    write_nt_tsv(out_dir / "coordinate_map.nt.tsv", args.source_name, builds)
    write_aa_tsv(out_dir / "coordinate_map.aa.tsv", args.source_name, builds)
    write_qc_tsv(out_dir / "coordinate_map.qc.tsv", args.source_name, builds)
    write_nt_alignment_fasta(out_dir / f"{args.source_name}_alignment.fasta", builds)
    write_aa_alignment_fasta(out_dir / f"{args.source_name}_aa_alignment.fasta", builds)
    warn_low_quality_maps(builds)

    for produced in (
        "coordinate_map.json",
        *shard_names,
        "coordinate_map.nt.tsv",
        "coordinate_map.aa.tsv",
        "coordinate_map.qc.tsv",
        f"{args.source_name}_alignment.fasta",
        f"{args.source_name}_aa_alignment.fasta",
    ):
        print(f"Wrote {out_dir / produced}")


def parse_gff_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in raw.split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = unquote(value)
    return attrs


def mature_peptide_frame_offset(parent: CdsModel, start: int, end: int) -> int | None:
    """Nucleotide offset of a mature peptide's first base within its parent's coding sequence.

    Returns None when that base is not among the parent's coding positions at all -- a
    malformed or mismatched Parent reference.

    Measured along the parent's *spliced* coding positions rather than as a genomic
    subtraction, so a spliced parent (M2, NEP, PA-X) gives the right answer; and taken from
    the peptide's end on the minus strand, since that is where its translation begins.
    """
    positions = parent.coding_positions()
    index = {pos: i for i, pos in enumerate(positions)}
    first = end if parent.strand == "-" else start
    if first not in index:
        return None
    offset = index[first]
    return len(positions) - 1 - offset if parent.strand == "-" else offset


def read_cds_models(path: Path, segment: str) -> dict[str, CdsModel]:
    """Gene key -> CDS model for one segment's GFF3.

    Reads `CDS` rows, and additionally `mature_protein_region_of_CDS` rows: NCBI/RefSeq flu
    annotations (the GsGD dataset among them) model HA as one CDS with HA1/HA2 as mature
    peptide children, whereas FluMutDB numbers HA by those subunits. Without the mature rows
    the gene-key intersection with the flumut source reference is empty for HA and no HA
    amino-acid map is built at all.

    A mature peptide is promoted to a gene key only when no real CDS already claims that key,
    so datasets exposing genuine HA1/HA2 CDS features (chicken2021) are unaffected.
    """
    grouped: dict[str, dict[str, Any]] = {}
    mature_rows: list[dict[str, Any]] = []
    cds_key_by_id: dict[str, str] = {}   # GFF3 ID of a CDS -> its canonical gene key
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in ("CDS", "mature_protein_region_of_CDS"):
                continue
            attrs = parse_gff_attributes(fields[8])

            if fields[2] == "mature_protein_region_of_CDS":
                # `product` names the subunit (HA1/HA2); protein_id/ID are what a VCF's ANN
                # Feature_ID can carry for it, so both are kept as lookup aliases.
                name = attrs.get("product") or attrs.get("Name") or attrs.get("protein_id")
                if not name:
                    continue
                mature_rows.append({
                    "gene": canonical_gene(name),
                    "contig": fields[0],
                    "strand": fields[6],
                    "part": (int(fields[3]), int(fields[4])),
                    "aliases": [a for a in (attrs.get("protein_id"), attrs.get("ID")) if a],
                    "parent": attrs.get("Parent", ""),
                })
                continue

            gene = attrs.get("gene") or attrs.get("Name") or attrs.get("ID")
            if not gene:
                continue
            gene = canonical_gene(gene)
            grouped.setdefault(
                gene,
                {
                    "segment": segment,
                    "gene": gene,
                    "contig": fields[0],
                    "strand": fields[6],
                    "parts": [],
                },
            )
            grouped[gene]["parts"].append((int(fields[3]), int(fields[4])))
            if attrs.get("ID"):
                cds_key_by_id[attrs["ID"]] = gene

    models: dict[str, CdsModel] = {}
    for gene, values in grouped.items():
        models[gene] = CdsModel(
            segment=values["segment"],
            gene=values["gene"],
            contig=values["contig"],
            strand=values["strand"],
            parts=tuple(sorted(values["parts"])),
        )

    mature_grouped: dict[str, dict[str, Any]] = {}
    for row in mature_rows:
        gene = row["gene"]
        if gene in models:
            continue  # a real CDS already defines this key; it wins

        # Treating a mature peptide as a CDS is only valid if it begins on one of its
        # parent's codon boundaries. translate_model reads the interval from its own first
        # base, so a peptide starting mid-codon yields a frameshifted protein -- one that
        # still translates, still aligns, still produces a position map, and reports entirely
        # different residues. Nothing errors; the markers are simply wrong. That silence is
        # the whole reason this is checked, and the reason an *unverifiable* peptide is
        # refused too: accepting the one case we cannot check would undo the point of
        # checking. RefSeq flu satisfies the condition (HA1 at 70 and HA2 at 1060 are both
        # codon-aligned to the HA CDS at 22).
        parent_key = cds_key_by_id.get(row["parent"])
        parent = models.get(parent_key) if parent_key else None
        if parent is None:
            print(f"WARNING: {segment} mature peptide {gene} "
                  f"({row['part'][0]}..{row['part'][1]}) has no resolvable parent CDS "
                  f"(Parent={row['parent']!r}), so its reading frame cannot be verified; "
                  f"skipping it rather than risking a frameshifted map",
                  file=sys.stderr)
            continue
        offset = mature_peptide_frame_offset(parent, *row["part"])
        if offset is None or offset % 3 != 0:
            reason = ("lies outside its parent's coding positions" if offset is None
                      else f"starts {offset % 3} base(s) into a codon of its parent")
            print(f"WARNING: {segment} mature peptide {gene} "
                  f"({row['part'][0]}..{row['part'][1]}) {reason} (parent CDS {parent.gene}); "
                  f"skipping it rather than mapping positions through the wrong frame",
                  file=sys.stderr)
            continue

        entry = mature_grouped.setdefault(
            gene,
            {
                "contig": row["contig"],
                "strand": row["strand"],
                "parts": [],
                "aliases": [],
                "parent": parent_key,   # non-None: an unresolvable parent was refused above
            },
        )
        entry["parts"].append(row["part"])
        for alias in row["aliases"]:
            if alias not in entry["aliases"]:
                entry["aliases"].append(alias)

    for gene, values in mature_grouped.items():
        models[gene] = CdsModel(
            segment=segment,
            gene=gene,
            contig=values["contig"],
            strand=values["strand"],
            parts=tuple(sorted(values["parts"])),
            aliases=tuple(values["aliases"]),
            is_mature=True,
            subunit_parent=values["parent"],
        )
    return models


def detect_stale_ha_annotation(
    segment: str, source_gene_keys: set[str], target_gene_keys: set[str]
) -> bool:
    """Flag target (Nextclade) HA annotations that expose no HA1/HA2 split at all.

    FluMutDB numbers HA by its HA1/HA2 subunits, as flumut's source reference does. The
    target is the flexible, user-chosen reference; it may express the split either as two
    real CDS features (chicken2021) or as mature_protein_region_of_CDS children of a single
    HA CDS (the NCBI/RefSeq model, GsGD among them) -- read_cds_models accepts both, so
    `target_gene_keys` already carries HA1/HA2 in either case.

    What remains stale is an annotation offering *only* a bare "HA" CDS with no subunit
    features of any kind: the shared-gene intersection is then empty and no HA amino-acid
    map can be built.
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
    build = SegmentBuild(
        segment=segment,
        source_chrom=source_record.id,
        target_chrom=target_record.id,
        nt=align_position_map(
            source_seq=str(source_record.seq).upper(),
            target_seq=str(target_record.seq).upper(),
            kind="nt",
        ),
    )

    source_gene_keys = set(source_models)
    target_gene_keys = set(target_models)

    # A gene whose subunits were recovered from its mature peptide children (HA, when HA1/HA2
    # came from mature_protein_region_of_CDS rows) has no source counterpart by construction.
    # Reporting it as "missing in source" would warn, every run, about the very arrangement
    # this build just handled correctly; name it as a subunit parent instead.
    build.subunit_parent_genes = tuple(sorted({
        model.subunit_parent
        for model in target_models.values()
        if model.is_mature and model.subunit_parent
    }))
    build.missing_target_genes = tuple(
        sorted(target_gene_keys - source_gene_keys - set(build.subunit_parent_genes))
    )
    build.extra_source_genes = tuple(sorted(source_gene_keys - target_gene_keys))
    build.stale_annotation = detect_stale_ha_annotation(
        segment, source_gene_keys, target_gene_keys
    )
    build.gene_aliases = {
        gene: model.aliases for gene, model in target_models.items() if model.aliases
    }

    for gene in build.missing_target_genes:
        build.unmapped_gene_lengths[gene] = len(
            translate_model(build.nt.target_seq, target_models[gene])
        )
    for gene in build.extra_source_genes:
        build.unmapped_gene_lengths[gene] = len(
            translate_model(build.nt.source_seq, source_models[gene])
        )

    for gene in sorted(source_gene_keys & target_gene_keys):
        build.aa[gene] = align_position_map(
            source_seq=translate_model(build.nt.source_seq, source_models[gene]),
            target_seq=translate_model(build.nt.target_seq, target_models[gene]),
            kind="aa",
        )

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


def count_co_optimal(alignments: Any) -> int:
    """How many co-optimal alignments the aligner found, capped and never fatal.

    The count is a dynamic-programming product and can be astronomically large for divergent
    sequences, so it is capped rather than reported exactly; all the reader needs to know is
    whether the one alignment that was used was the only choice.
    """
    try:
        count = len(alignments)
    except (OverflowError, ValueError, MemoryError, TypeError):
        return CO_OPTIMAL_REPORTING_CAP
    return min(int(count), CO_OPTIMAL_REPORTING_CAP)


def align_position_map(source_seq: str, target_seq: str, kind: str) -> AlignmentResult:
    aligner = make_aligner(kind)
    alignments = aligner.align(target_seq, source_seq)
    # PairwiseAligner can return many equally optimal alignments; this takes the first. With
    # equal-length gapless proteins there is exactly one, but wherever an indel forces a gap
    # the placement among co-optimal solutions is arbitrary and the position map inherits
    # that arbitrariness -- so record how many there were, and surface it in the QC table.
    co_optimal = count_co_optimal(alignments)
    alignment = next(iter(alignments))
    target_aligned = str(alignment[0, :])
    source_aligned = str(alignment[1, :])

    source_to_target: dict[int, int | None] = {}
    target_pos = 0
    source_pos = 0
    identities = 0
    aligned = 0
    gaps = 0

    for target_char, source_char in zip(target_aligned, source_aligned, strict=True):
        target_has_base = target_char != "-"
        source_has_base = source_char != "-"
        if target_has_base:
            target_pos += 1
        if source_has_base:
            source_pos += 1

        if target_has_base and source_has_base:
            source_to_target[source_pos] = target_pos
            aligned += 1
            if target_char == source_char:
                identities += 1
        elif source_has_base:
            source_to_target[source_pos] = None
            gaps += 1
        elif target_has_base:
            gaps += 1

    return AlignmentResult(
        source_seq=source_seq,
        target_seq=target_seq,
        source_aligned=source_aligned,
        target_aligned=target_aligned,
        source_to_target=source_to_target,
        identities=identities,
        aligned=aligned,
        gaps=gaps,
        co_optimal=co_optimal,
    )


def map_status(target_pos: int | None) -> str:
    """The per-position status published in the nt/aa TSVs, read off the coordinate itself."""
    return "mapped" if target_pos else "unmapped"


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


def format_pct(value: float) -> str:
    """A percentage for the QC TSV; empty rather than the literal 'nan'.

    A NaN means nothing aligned at all. `status` carries the verdict either way, and leaving
    the column blank matches the missing_in_source / extra_in_source rows.
    """
    return "" if math.isnan(value) else f"{value:.4f}"


def make_map_doc(
    source_name: str,
    target_name: str,
    builds: dict[str, SegmentBuild],
) -> dict[str, Any]:
    segments: dict[str, Any] = {}
    for segment, build in builds.items():
        genes: dict[str, Any] = {}
        for gene, aa in sorted(build.aa.items()):
            gene_doc: dict[str, Any] = {
                "source_length": len(aa.source_seq),
                "target_length": len(aa.target_seq),
                "aa_identity_pct": aa.identity_pct,
                "aa_aligned_positions": aa.aligned,
                "aa_gap_columns": aa.gaps,
                "aa_co_optimal_alignments": aa.co_optimal,
                "aa_source_to_target": {
                    str(pos): value for pos, value in sorted(aa.source_to_target.items())
                },
            }
            # Alternative names a VCF's ANN Feature_ID may use for this gene (a mature
            # peptide's protein_id/GFF3 ID). The merge needs them because such an entry's
            # ANN Gene_Name is the parent gene, not the subunit.
            aliases = build.gene_aliases.get(gene)
            if aliases:
                gene_doc["aliases"] = list(aliases)
            genes[gene] = gene_doc

        segments[segment] = {
            "source_chrom": build.source_chrom,
            "target_chrom": build.target_chrom,
            "source_length": len(build.nt.source_seq),
            "target_length": len(build.nt.target_seq),
            "nt_identity_pct": build.nt.identity_pct,
            "nt_aligned_positions": build.nt.aligned,
            "nt_gap_columns": build.nt.gaps,
            "nt_co_optimal_alignments": build.nt.co_optimal,
            "nt_source_to_target": {
                str(pos): value for pos, value in sorted(build.nt.source_to_target.items())
            },
            "missing_target_genes": list(build.missing_target_genes),
            "extra_source_genes": list(build.extra_source_genes),
            "stale_annotation": build.stale_annotation,
            # Target genes deliberately left unmapped because their subunits carry the map
            # instead. The merge treats an ANN Gene_Name in this list as known, so it does
            # not report the parent as an unmapped gene on every subunit row.
            "subunit_parent_genes": list(build.subunit_parent_genes),
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
            nt = build.nt
            for source_pos, target_pos in sorted(nt.source_to_target.items()):
                writer.writerow(
                    {
                        "source_reference": source_name,
                        "segment": segment,
                        "source_chrom": build.source_chrom,
                        "source_pos": source_pos,
                        "source_base": nt.source_seq[source_pos - 1],
                        "target_chrom": build.target_chrom,
                        "target_pos": target_pos or "",
                        "target_base": nt.target_seq[target_pos - 1] if target_pos else "",
                        "status": map_status(target_pos),
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
            for gene, aa in sorted(build.aa.items()):
                for source_pos, target_pos in sorted(aa.source_to_target.items()):
                    writer.writerow(
                        {
                            "source_reference": source_name,
                            "segment": segment,
                            "gene": gene,
                            "source_aa_pos": source_pos,
                            "source_residue": aa.source_seq[source_pos - 1],
                            "target_aa_pos": target_pos or "",
                            "target_residue": (
                                aa.target_seq[target_pos - 1] if target_pos else ""
                            ),
                            "status": map_status(target_pos),
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
        "co_optimal_alignments",
        "status",
    )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for segment in SEGMENTS:
            build = builds[segment]
            for layer, gene, result in [
                ("nt", "", build.nt),
                *[("aa", gene, aa) for gene, aa in sorted(build.aa.items())],
            ]:
                writer.writerow(
                    {
                        "source_reference": source_name,
                        "segment": segment,
                        "layer": layer,
                        "gene": gene,
                        "source_length": len(result.source_seq),
                        "target_length": len(result.target_seq),
                        "aligned_positions": result.aligned,
                        "gap_columns": result.gaps,
                        "identity_pct": format_pct(result.identity_pct),
                        "co_optimal_alignments": result.co_optimal,
                        "status": result.status,
                    }
                )
            for gene in build.missing_target_genes:
                writer.writerow(
                    {
                        "source_reference": source_name,
                        "segment": segment,
                        "layer": "aa",
                        "gene": gene,
                        "source_length": "",
                        "target_length": build.unmapped_gene_lengths.get(gene, ""),
                        "aligned_positions": 0,
                        "gap_columns": "",
                        "identity_pct": "",
                        "status": (
                            "stale_annotation" if build.stale_annotation else "missing_in_source"
                        ),
                    }
                )
            for gene in build.extra_source_genes:
                writer.writerow(
                    {
                        "source_reference": source_name,
                        "segment": segment,
                        "layer": "aa",
                        "gene": gene,
                        "source_length": build.unmapped_gene_lengths.get(gene, ""),
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
            write_alignment_block(handle, header, build.nt.target_aligned, build.nt.source_aligned)


def write_aa_alignment_fasta(path: Path, builds: dict[str, SegmentBuild]) -> None:
    with path.open("w") as handle:
        for segment in SEGMENTS:
            build = builds[segment]
            for gene, aa in sorted(build.aa.items()):
                header = (
                    f"segment={segment};gene={gene};alphabet=aa;target={build.target_chrom};"
                    f"source={build.source_chrom}"
                )
                write_alignment_block(handle, header, aa.target_aligned, aa.source_aligned)


def pct_text(value: float) -> str:
    return "n/a (nothing aligned)" if math.isnan(value) else f"{value:.2f}%"


def co_optimal_text(count: int) -> str:
    return f">={CO_OPTIMAL_REPORTING_CAP}" if count >= CO_OPTIMAL_REPORTING_CAP else str(count)


def warn_low_quality_maps(builds: dict[str, SegmentBuild]) -> None:
    """Report questionable maps on stderr. The maps are written either way; this and
    coordinate_map.qc.tsv are how a reviewer learns which ones to distrust."""
    for segment, build in builds.items():
        for label, result in [
            (f"{segment} nucleotide", build.nt),
            *[(f"{segment}/{gene} amino-acid", aa) for gene, aa in sorted(build.aa.items())],
        ]:
            if result.is_questionable:
                print(
                    f"WARNING: {label} map has identity={pct_text(result.identity_pct)} "
                    f"and source coverage={pct_text(result.coverage_pct)}",
                    file=sys.stderr,
                )
            if result.co_optimal > 1:
                print(
                    f"WARNING: {label} alignment had {co_optimal_text(result.co_optimal)} "
                    f"co-optimal solutions; the gap placement in this map is one arbitrary "
                    f"choice among them",
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


if __name__ == "__main__":
    raise SystemExit(main())
