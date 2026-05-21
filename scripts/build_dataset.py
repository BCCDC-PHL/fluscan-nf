#!/usr/bin/env python3
"""Build Nextclade reference dataset directories from GenBank files.

Reads a config.json from the dataset directory for all strain-specific
values (display name, ID prefix, strain, QC thresholds, segment list).

Usage:
    python build_chicken2021_dataset.py <dataset_dir>
    python build_chicken2021_dataset.py .       # run from the dataset dir

Each segment subdir must already contain a <dir>/<dir>.gb GenBank file.
Outputs per-segment: reference.fasta, genome_annotation.gff3, pathogen.json
"""

import argparse
import json
import sys
import textwrap
import warnings
from datetime import datetime, timezone
from pathlib import Path

from Bio import BiopythonParserWarning, SeqIO
from Bio.SeqFeature import CompoundLocation

warnings.simplefilter("ignore", BiopythonParserWarning)

NOW = datetime.now(timezone.utc)
TIMESTAMP = NOW.strftime("%Y-%m-%dT%H:%M:%SZ")
TAG = NOW.strftime("%Y-%m-%d--%H-%M-%SZ")

_PATHOGEN_SCHEMA = (
    "https://raw.githubusercontent.com/nextstrain/nextclade/refs/heads/release"
    "/packages/nextclade-schemas/input-pathogen-json.schema.json"
)


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def load_config(dataset_dir: Path) -> dict:
    cfg_path = dataset_dir / "config.json"
    if not cfg_path.exists():
        sys.exit(f"config.json not found in {dataset_dir}")
    with open(cfg_path) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# GFF3 helpers
# ---------------------------------------------------------------------------

def gff_row(seqid, source, ftype, start1, end1, score, strand, phase, attrs):
    attr_str = ";".join(f"{k}={v}" for k, v in attrs.items())
    return f"{seqid}\t{source}\t{ftype}\t{start1}\t{end1}\t{score}\t{strand}\t{phase}\t{attr_str}"


def compute_phases(parts):
    """Return GFF3 phase for each part of a compound CDS (0, 1, or 2)."""
    phases = []
    cumlen = 0
    for part in parts:
        phases.append((3 - (cumlen % 3)) % 3)
        cumlen += len(part)
    return phases


def build_gff3(record, seqid, seg_num, cfg):
    prefix = cfg["id_prefix"]
    strain = cfg["strain"]

    lines = [
        "##gff-version 3",
        f"##sequence-region {seqid} 1 {len(record.seq)}",
    ]
    length = len(record.seq)
    mat_pep_counter = {}
    gene_counter = 0
    current_locus_tag = ""
    current_protein_id = ""

    for feat in record.features:
        ft = feat.type
        loc = feat.location
        q = feat.qualifiers
        gene_name = q.get("gene", [""])[0]

        if ft == "source":
            attrs = {
                "ID": f"{seqid}:1..{length}",
                "Name": seg_num,
                "gbkey": "Src",
                "mol_type": "viral cRNA",
                "strain": strain,
            }
            lines.append(gff_row(seqid, ".", "region", 1, length, ".", "+", ".", attrs))

        elif ft == "gene":
            gene_counter += 1
            locus_tag = f"{prefix}_s{seg_num}gp{gene_counter}"
            protein_id = f"{prefix}_{gene_name}_p1"
            current_locus_tag = locus_tag
            current_protein_id = protein_id
            start1 = int(loc.start) + 1
            end1 = int(loc.end)
            attrs = {
                "ID": f"gene-{locus_tag}",
                "Name": gene_name,
                "gbkey": "Gene",
                "gene": gene_name,
                "gene_biotype": "protein_coding",
                "locus_tag": locus_tag,
            }
            lines.append(gff_row(seqid, ".", "gene", start1, end1, ".", "+", ".", attrs))

        elif ft == "CDS":
            product = q.get("product", [""])[0]
            protein_id = current_protein_id
            locus_tag = current_locus_tag
            cds_attrs = {
                "ID": f"cds-{protein_id}",
                "Parent": f"gene-{locus_tag}",
                "Name": protein_id,
                "gbkey": "CDS",
                "gene": gene_name,
                "locus_tag": locus_tag,
                "product": product,
                "protein_id": protein_id,
            }
            if isinstance(loc, CompoundLocation):
                parts = list(loc.parts)
                phases = compute_phases(parts)
                for part, phase in zip(parts, phases):
                    lines.append(gff_row(
                        seqid, ".", "CDS",
                        int(part.start) + 1, int(part.end),
                        ".", "+", phase, cds_attrs,
                    ))
            else:
                lines.append(gff_row(
                    seqid, ".", "CDS",
                    int(loc.start) + 1, int(loc.end),
                    ".", "+", "0", cds_attrs,
                ))

        elif ft == "sig_peptide":
            attrs = {
                "ID": f"sigpep-{current_protein_id}",
                "Parent": f"cds-{current_protein_id}",
                "gbkey": "Prot",
            }
            lines.append(gff_row(
                seqid, ".", "signal_peptide_region_of_CDS",
                int(loc.start) + 1, int(loc.end),
                ".", "+", ".", attrs,
            ))

        elif ft == "mat_peptide":
            product = q.get("product", [""])[0]
            mat_pep_counter[gene_name] = mat_pep_counter.get(gene_name, 0) + 1
            n = mat_pep_counter[gene_name]
            mp_protein_id = f"{prefix}_{gene_name}_mp{n}"
            attrs = {
                "ID": f"matpep-{current_protein_id}-{n}",
                "Parent": f"cds-{current_protein_id}",
                "gbkey": "Prot",
                "product": product,
                "protein_id": mp_protein_id,
            }
            lines.append(gff_row(
                seqid, ".", "mature_protein_region_of_CDS",
                int(loc.start) + 1, int(loc.end),
                ".", "+", ".", attrs,
            ))

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# pathogen.json
# ---------------------------------------------------------------------------

def build_pathogen_json(seg_cfg, seqid, cfg):
    compat = cfg["nextclade_compatibility"]
    data = {
        "$schema": _PATHOGEN_SCHEMA,
        "attributes": {
            "name": cfg["display_name"],
            "reference accession": seqid,
            "reference name": cfg["strain"],
            "segment": seg_cfg["dir"],
        },
        "compatibility": {"cli": compat, "web": compat},
        "defaultCds": seg_cfg["default_cds"],
        "files": {
            "genomeAnnotation": "genome_annotation.gff3",
            "pathogenJson": "pathogen.json",
            "reference": "reference.fasta",
        },
        "qc": cfg["qc"],
        "schemaVersion": "3.0.0",
        "version": {
            "updatedAt": TIMESTAMP,
            "tag": TAG,
            "compatibility": {"cli": compat, "web": compat},
        },
    }
    return json.dumps(data, indent=2) + "\n"


# ---------------------------------------------------------------------------
# Per-segment processing
# ---------------------------------------------------------------------------

def process_segment(seg_cfg, dataset_dir, cfg):
    prefix = cfg["id_prefix"]
    seqid = f"{prefix}_{seg_cfg['name']}"
    outdir = dataset_dir / seg_cfg["dir"]
    outdir.mkdir(exist_ok=True)

    gb_in = dataset_dir / f"{seg_cfg['dir']}.gb"
    if not gb_in.exists():
        sys.exit(f"GenBank file not found: {gb_in}")

    gb_in.rename(outdir / f"{seg_cfg['dir']}.gb")
    gb_in = outdir / f"{seg_cfg['dir']}.gb"

    record = next(SeqIO.parse(gb_in, "genbank"))

    seq = str(record.seq).upper()
    with open(outdir / "reference.fasta", "w") as fh:
        fh.write(f">{seqid}\n")
        for chunk in textwrap.wrap(seq, 70):
            fh.write(chunk + "\n")

    gff_text = build_gff3(record, seqid, seg_cfg["segment_num"], cfg)
    (outdir / "genome_annotation.gff3").write_text(gff_text)

    pj_text = build_pathogen_json(seg_cfg, seqid, cfg)
    (outdir / "pathogen.json").write_text(pj_text)

    print(f"  {seg_cfg['dir']}/  ({len(record.seq)} bp, {len(record.features)} features)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "dataset_dir", nargs="?", default=".",
        help="Path to the dataset directory containing config.json and segment subdirs "
             "(default: current working directory)",
    )
    return p.parse_args()


def main():
    args = get_args()
    dataset_dir = Path(args.dataset_dir).resolve()

    if not dataset_dir.is_dir():
        sys.exit(f"Dataset directory not found: {dataset_dir}")

    cfg = load_config(dataset_dir)
    print(f"Building {cfg['display_name']} Nextclade datasets — {TIMESTAMP}")

    for seg_cfg in cfg["segments"]:
        process_segment(seg_cfg, dataset_dir, cfg)

    print("Done.")


if __name__ == "__main__":
    main()
