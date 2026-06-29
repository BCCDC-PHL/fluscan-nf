#!/usr/bin/env python3
"""Report consensus-to-reference identity QC for a single influenza segment.

Globally aligns one cleaned consensus segment to its workflow reference segment
and writes a single-row TSV with the segment's identity and PASS/WARN/FAIL
status. One invocation handles one (sample, segment); cross-segment and
cross-sample aggregation is left to the caller.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import Align, SeqIO

OUTPUT_COLUMNS = (
    "sample_id",
    "workflow",
    "segment",
    "status",
    "identity_pct",
    "consensus_length",
    "reference_length",
    "reference_id",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Globally align one cleaned consensus segment to its workflow "
        "reference and write a single-row identity QC TSV."
    )
    parser.add_argument("--consensus", required=True, type=Path, help="Single-segment consensus FASTA")
    parser.add_argument("--reference", required=True, type=Path, help="Single-segment reference FASTA")
    parser.add_argument("--sample-id", required=True, help="Sample/isolate identifier")
    parser.add_argument("--segment", required=True, help="Influenza segment label (e.g. HA)")
    parser.add_argument("--workflow", required=True, help="Workflow name (e.g. paftools, nextclade)")
    parser.add_argument("--pass-pct", required=True, type=float, help="PASS identity threshold")
    parser.add_argument("--warn-pct", required=True, type=float, help="WARN identity threshold")
    parser.add_argument("--out-tsv", required=True, type=Path, help="Output single-row TSV")
    return parser.parse_args()


def load_single_sequence(fasta_path):
    """Return (record_id, sequence) of the first record, or (None, "") if none."""
    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        return None, ""
    for record in SeqIO.parse(fasta_path, "fasta"):
        return record.id, str(record.seq).upper()
    return None, ""


def build_aligner():
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    return aligner


def alignment_identity_pct(query, reference):
    """Return exact-match identity across global alignment columns."""
    if not query or not reference:
        return 0.0

    aligner = build_aligner()
    alignment = aligner.align(reference, query)[0]
    coords = alignment.coordinates

    matches = 0
    columns = 0
    for idx in range(coords.shape[1] - 1):
        ref_start, ref_end = coords[0][idx], coords[0][idx + 1]
        qry_start, qry_end = coords[1][idx], coords[1][idx + 1]
        ref_span = ref_end - ref_start
        qry_span = qry_end - qry_start
        block_columns = max(ref_span, qry_span)
        columns += block_columns

        if ref_span and qry_span:
            for offset in range(block_columns):
                if reference[ref_start + offset] == query[qry_start + offset]:
                    matches += 1

    if columns == 0:
        return 0.0
    return matches / columns * 100


def segment_status(identity_pct, pass_pct, warn_pct):
    if identity_pct >= pass_pct:
        return "PASS"
    if identity_pct >= warn_pct:
        return "WARN"
    return "FAIL"


def build_qc_row(sample_id, workflow, segment, consensus, reference, pass_pct, warn_pct):
    consensus_id, consensus_seq = consensus
    reference_id, reference_seq = reference

    if not consensus_seq or not reference_seq:
        identity = None
        status = "FAIL"
    else:
        identity = alignment_identity_pct(consensus_seq, reference_seq)
        status = segment_status(identity, pass_pct, warn_pct)

    return {
        "sample_id": sample_id,
        "workflow": workflow,
        "segment": segment,
        "status": status,
        "identity_pct": identity,               # float, or None when not evaluated
        "consensus_length": len(consensus_seq),
        "reference_length": len(reference_seq),
        "reference_id": reference_id,           # str, or None when reference missing
    }


def write_tsv(row, out_tsv):
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame([row], columns=OUTPUT_COLUMNS)
    df.to_csv(out_tsv, sep="\t", index=False, float_format="%.3f", na_rep="NA")


def main():
    args = parse_args()
    consensus = load_single_sequence(args.consensus)
    reference = load_single_sequence(args.reference)
    if consensus[1] == "":
        print(f"WARNING: no consensus sequence for {args.sample_id} {args.segment}", file=sys.stderr)
    if reference[1] == "":
        print(f"WARNING: no reference sequence for {args.sample_id} {args.segment}", file=sys.stderr)
    row = build_qc_row(
        sample_id=args.sample_id,
        workflow=args.workflow,
        segment=args.segment,
        consensus=consensus,
        reference=reference,
        pass_pct=args.pass_pct,
        warn_pct=args.warn_pct,
    )
    write_tsv(row, args.out_tsv)


if __name__ == "__main__":
    main()
