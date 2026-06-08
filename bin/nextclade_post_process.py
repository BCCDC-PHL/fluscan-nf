#!/usr/bin/env python3
"""Split Nextclade batch JSON and TSV outputs into per-sample files.

Optionally writes a per-batch status summary (one row per sample) tracking
whether each sample passed or failed Nextclade alignment.

Usage:
    nextclade_post_process.py \\
        --json   batch_nextclade.json \\
        --tsv    batch_nextclade.tsv \\
        --segment HA \\
        --status-summary HA_nextclade_status.tsv
"""

import argparse
import csv
import json
from pathlib import Path


def split_json(batch_json_path: Path, segment: str) -> None:
    """Write one {seqName}_{segment}_nextclade.json per entry in results."""
    with open(batch_json_path) as fh:
        batch = json.load(fh)

    for result in batch.get("results", []):
        sample_id = result["seqName"]
        sample_data = dict(batch)
        sample_data["results"] = [result]
        sample_data["errors"] = [
            e for e in batch.get("errors", []) if e.get("seqName") == sample_id
        ]
        out_path = Path(f"{sample_id}_{segment}_nextclade.json")
        with open(out_path, "w") as out_fh:
            json.dump(sample_data, out_fh)


def split_tsv(batch_tsv_path: Path, segment: str) -> None:
    """Write one {seqName}_{segment}_nextclade.tsv per row (header + data row)."""
    with open(batch_tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample_id = row["seqName"]
            out_path = Path(f"{sample_id}_{segment}_nextclade.tsv")
            with open(out_path, "w", newline="") as out_fh:
                writer = csv.DictWriter(out_fh, fieldnames=reader.fieldnames, delimiter="\t")
                writer.writeheader()
                writer.writerow(row)


def write_status_summary(batch_tsv_path: Path, segment: str, summary_path: Path) -> None:
    """Write a per-sample status summary using the TSV as the authoritative roster.

    Columns: sample, segment, status, error
    The TSV contains both passing and failing samples; the JSON omits failures.
    A sample is considered failed when its 'errors' column is non-empty.
    """
    with open(batch_tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)

    with open(summary_path, "w", newline="") as out_fh:
        writer = csv.DictWriter(
            out_fh,
            fieldnames=["sample", "segment", "status", "error"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            error_text = row.get("errors", "").strip()
            writer.writerow(
                {
                    "sample": row["seqName"],
                    "segment": segment,
                    "status": "fail" if error_text else "pass",
                    "error": error_text,
                }
            )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--json",    required=True,  type=Path, dest="json_path",  help="Batch Nextclade JSON")
    p.add_argument("--tsv",     required=True,  type=Path, dest="tsv_path",   help="Batch Nextclade TSV")
    p.add_argument("--segment", required=True,             dest="segment",    help="Segment name (e.g. HA)")
    p.add_argument("--status-summary", type=Path, default=None,              help="Path for per-sample status summary TSV")
    args = p.parse_args()

    split_json(args.json_path, args.segment)
    split_tsv(args.tsv_path, args.segment)

    if args.status_summary:
        write_status_summary(args.tsv_path, args.segment, args.status_summary)


if __name__ == "__main__":
    main()
