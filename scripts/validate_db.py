#!/usr/bin/env python3
"""Check a reference database bundle before a run depends on it.

A bundle is a directory holding every reference the pipeline needs for one reference
set, in a fixed layout:

    <db>/
      db.yml                                  optional metadata; every key optional
      blastn_db.fa                            headers accession|subtype|segment
      nextclade/{pb2,pb1,pa,ha,np,na,m,ns}/   reference.fasta + genome_annotation.gff3
      snpeff/snpEff.config
      snpeff/data/<genome>/snpEffectPredictor.bin

`modules/shared/database.nf` validates the same layout at the start of every run, but
only structurally and only for the components that run's workflow consumes. This script
is the deeper check, meant to be run once after building or moving a bundle: it opens
the files and confirms their *contents* line up -- that BLAST headers parse, that every
genome named in snpEff.config has a built database next to it, that the Nextclade
datasets carry the reference each segment claims.

The failure this exists to catch is a bundle that resolves fine and then annotates
nothing: a snpEff.config listing influenza_h5n1 with no influenza_h5n1/ under data/
produces an empty ANN= on every record rather than an error.

Bundles may legitimately be partial -- the Nextclade workflow needs no BLAST database.
Missing components are reported as SKIP, not failure; only a component that is present
and broken fails. Use --require to demand specific components.

Exit codes: 0 all present components valid, 1 something is wrong.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

SCRIPT_VERSION = "1.0.0"

SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
COMPONENTS = ["blastn", "nextclade", "snpeff"]

# >accession|subtype|segment, e.g. >NC_026438.1|H1N1|PB2
BLAST_HEADER = re.compile(r"^>([^|\s]+)\|([^|\s]+)\|([^|\s]+)\s*$")

# `influenza_h5n1.genome : Influenza_A_H5N1` -- the genome key is what snpEff is invoked with
SNPEFF_GENOME = re.compile(r"^\s*([A-Za-z0-9_.\-]+)\s*\.genome\s*:")

DB_YML_KEYS = {
    "name", "version", "description", "reference_strain",
    "nextclade_ref_name", "snpeff_genomes", "source",
}


@dataclass
class Report:
    """Findings for one bundle. `problems` are failures; `notes` are informational."""

    root: Path
    problems: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    components: dict[str, str] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return not self.problems


def load_db_yml(root: Path, report: Report) -> dict:
    """Parse db.yml if present. Absent is fine; unparseable or non-mapping is not."""
    manifest = root / "db.yml"
    if not manifest.is_file():
        report.warnings.append("no db.yml; name falls back to the directory name and version to 'unknown'")
        return {}

    try:
        import yaml
    except ImportError:  # pragma: no cover - yaml ships with the runtime env
        report.warnings.append("PyYAML unavailable; db.yml not checked")
        return {}

    try:
        parsed = yaml.safe_load(manifest.read_text())
    except yaml.YAMLError as exc:
        report.problems.append(f"db.yml is not valid YAML: {exc}")
        return {}

    if parsed is None:
        return {}
    if not isinstance(parsed, dict):
        report.problems.append(f"db.yml must contain a mapping, found {type(parsed).__name__}")
        return {}

    for key in sorted(set(parsed) - DB_YML_KEYS):
        report.warnings.append(f"db.yml: unrecognised key '{key}' (ignored by the pipeline)")

    # snakeyaml on the Nextflow side turns an unquoted 2025-07-07 into a date; the
    # pipeline renders it back, but quoting keeps the value the string it was written as.
    version = parsed.get("version")
    if version is not None and not isinstance(version, str):
        report.warnings.append(
            f"db.yml: version {version!r} parses as {type(version).__name__}, not a string; quote it"
        )

    genomes = parsed.get("snpeff_genomes")
    if genomes is not None and not isinstance(genomes, list):
        report.problems.append("db.yml: snpeff_genomes must be a list")

    return parsed


def check_blastn(root: Path, report: Report) -> None:
    fasta = root / "blastn_db.fa"
    if not fasta.is_file():
        report.components["blastn"] = "SKIP (no blastn_db.fa)"
        return

    total = 0
    malformed: list[str] = []
    segments: set[str] = set()
    subtypes: set[str] = set()

    with fasta.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            total += 1
            match = BLAST_HEADER.match(line.rstrip("\n"))
            if match is None:
                if len(malformed) < 5:
                    malformed.append(line.strip())
                continue
            _, subtype, segment = match.groups()
            segments.add(segment.upper())
            subtypes.add(subtype)

    if total == 0:
        report.problems.append("blastn_db.fa contains no sequences")
        report.components["blastn"] = "FAIL"
        return

    if malformed:
        report.problems.append(
            "blastn_db.fa headers must be accession|subtype|segment; "
            f"{len(malformed)} or more do not parse, e.g. {malformed[0]!r}"
        )

    missing = [s for s in SEGMENTS if s not in segments]
    if missing:
        # Not fatal: a subtype-specific database need not carry all eight, but a
        # consensus segment with no reference silently drops out of the BLAST step.
        report.warnings.append(f"blastn_db.fa has no sequence for segment(s): {', '.join(missing)}")

    report.components["blastn"] = (
        "FAIL" if malformed else f"OK ({total} sequences, {len(subtypes)} subtype(s))"
    )


def check_nextclade(root: Path, report: Report) -> None:
    datasets = root / "nextclade"
    if not datasets.is_dir():
        report.components["nextclade"] = "SKIP (no nextclade/)"
        return

    incomplete: list[str] = []
    for segment in SEGMENTS:
        directory = datasets / segment.lower()
        if not directory.is_dir():
            incomplete.append(f"{segment} (no {segment.lower()}/)")
            continue

        missing = [
            name for name in ("reference.fasta", "genome_annotation.gff3")
            if not (directory / name).is_file()
        ]
        if missing:
            incomplete.append(f"{segment} (missing {', '.join(missing)})")
            continue

        reference = directory / "reference.fasta"
        first = reference.read_text().splitlines()[:1]
        if not first or not first[0].startswith(">"):
            incomplete.append(f"{segment} (reference.fasta is not FASTA)")
            continue

        # pathogen.json is what `nextclade run --input-dataset` reads for QC config;
        # nextclade will run without it but silently skips the dataset's QC rules.
        if not (directory / "pathogen.json").is_file():
            report.warnings.append(f"nextclade/{segment.lower()}/ has no pathogen.json; dataset QC will not run")

    if incomplete:
        report.problems.append("incomplete Nextclade datasets: " + "; ".join(incomplete))
        report.components["nextclade"] = "FAIL"
    else:
        report.components["nextclade"] = f"OK ({len(SEGMENTS)} segment datasets)"


def parse_snpeff_genomes(config: Path) -> list[str]:
    """Genome keys declared in a snpEff config, in file order."""
    genomes = []
    for line in config.read_text().splitlines():
        if line.lstrip().startswith("#"):
            continue
        match = SNPEFF_GENOME.match(line)
        if match:
            genomes.append(match.group(1))
    return genomes


def check_snpeff(root: Path, report: Report, declared: list | None) -> None:
    directory = root / "snpeff"
    config = directory / "snpEff.config"
    data = directory / "data"

    if not directory.is_dir():
        report.components["snpeff"] = "SKIP (no snpeff/)"
        return

    if not config.is_file():
        report.problems.append(f"snpeff/ exists but {config.name} is missing")
        report.components["snpeff"] = "FAIL"
        return
    if not data.is_dir():
        report.problems.append("snpeff/ exists but snpeff/data/ is missing")
        report.components["snpeff"] = "FAIL"
        return

    genomes = parse_snpeff_genomes(config)
    if not genomes:
        report.problems.append("snpeff/snpEff.config declares no <genome>.genome entry")
        report.components["snpeff"] = "FAIL"
        return

    # data.dir is overridden by -dataDir at runtime, so its value cannot break a run --
    # but a stale absolute path is a reliable sign the bundle was moved without being
    # rebuilt, and it is what someone running snpEff by hand will hit.
    for line in config.read_text().splitlines():
        if line.strip().startswith("data.dir"):
            value = line.split("=", 1)[-1].strip()
            if value.startswith("/") and not Path(value).is_dir():
                report.warnings.append(
                    f"snpEff.config data.dir points outside the bundle at a path that does not exist "
                    f"({value}); the pipeline overrides it with -dataDir, but set it to ./data"
                )
            break

    unbuilt = []
    for genome in genomes:
        predictor = data / genome / "snpEffectPredictor.bin"
        if not predictor.is_file():
            unbuilt.append(genome)

    if unbuilt:
        report.problems.append(
            "snpEff.config declares genome(s) with no built database under snpeff/data/: "
            + ", ".join(unbuilt)
            + " -- annotation would produce empty ANN= rather than fail"
        )
        report.components["snpeff"] = "FAIL"
        return

    if declared:
        missing = [g for g in declared if g not in genomes]
        if missing:
            report.warnings.append(
                "db.yml lists snpeff_genomes not declared in snpEff.config: " + ", ".join(missing)
            )

    report.components["snpeff"] = f"OK ({', '.join(genomes)})"


def validate(root: Path, require: list[str]) -> Report:
    report = Report(root=root)

    if not root.is_dir():
        report.problems.append(f"{root} does not exist or is not a directory")
        return report

    report.metadata = load_db_yml(root, report)

    check_blastn(root, report)
    check_nextclade(root, report)
    check_snpeff(root, report, report.metadata.get("snpeff_genomes"))

    for component in require:
        if report.components.get(component, "").startswith("SKIP"):
            report.problems.append(f"--require {component}: component is absent from this bundle")

    return report


def render(report: Report) -> str:
    name = report.metadata.get("name") or report.root.name
    version = report.metadata.get("version", "unknown")
    lines = [f"{name} (version: {version})", f"  {report.root}", ""]
    for component in COMPONENTS:
        lines.append(f"  {component:<10} {report.components.get(component, 'SKIP')}")
    if report.warnings:
        lines.append("")
        lines.extend(f"  WARNING: {w}" for w in report.warnings)
    if report.problems:
        lines.append("")
        lines.extend(f"  PROBLEM: {p}" for p in report.problems)
    lines.append("")
    lines.append("  VALID" if report.ok else "  INVALID")
    return "\n".join(lines)


def get_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Validate a fluscan-nf reference database bundle.",
        epilog="Components absent from the bundle are reported as SKIP, not failure.",
    )
    parser.add_argument("db", type=Path, help="Path to the bundle root")
    parser.add_argument(
        "-r", "--require", action="append", choices=COMPONENTS, default=[],
        help="Fail if this component is absent (repeatable)",
    )
    parser.add_argument("--json", action="store_true", help="Emit findings as JSON")
    parser.add_argument("--version", action="version", version=f"validate_db.py {SCRIPT_VERSION}")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = get_args(argv)
    report = validate(args.db, args.require)

    if args.json:
        print(json.dumps({
            "root": str(report.root),
            "name": report.metadata.get("name") or report.root.name,
            "version": str(report.metadata.get("version", "unknown")),
            "components": report.components,
            "warnings": report.warnings,
            "problems": report.problems,
            "valid": report.ok,
        }, indent=2))
    else:
        print(render(report))

    return 0 if report.ok else 1


if __name__ == "__main__":
    sys.exit(main())
