# Reference database bundles

Every reference the pipeline needs for one reference set lives in a single directory — a
*bundle* — passed as `--db`. The FluMut coordinate reference is a second, independent
bundle passed as `--flumut_db`, because it changes far less often and only the optional
merge step reads it.

```
nextflow run BCCDC-PHL/fluscan-nf \
  -profile conda \
  --fasta_input /path/to/consensus_fastas \
  --db /path/to/refs/gsgd \
  --outdir /path/to/output
```

## Layout

```
<db>/
├── db.yml                                      optional metadata; every key optional
├── blastn_db.fa                                nt reference; headers accession|subtype|segment
├── nextclade/
│   ├── PB2/  reference.fasta, genome_annotation.gff3, pathogen.json
│   ├── PB1/  …  and PA, HA, NP, NA, M, NS  (uppercase segment names)
│   └── config.json                             optional; written by scripts/build_dataset.py
└── snpeff/
    ├── snpEff.config
    └── data/
        └── influenza_h5n1/
            ├── genes.gbk
            └── snpEffectPredictor.bin
```

The filenames are fixed and not configurable. `modules/shared/database.nf` is the only
place that knows this layout; `--db` is the only path you supply.

### Bundles may be partial

Each workflow validates only what it consumes, so a bundle need not be complete:

| Workflow | Requires | Optional |
| --- | --- | --- |
| `SHORT_READ` | `blastn_db.fa` | `snpeff/` — when present, SnpEff annotation, MNV correction and the watchlist scan run; when absent they are skipped and consensus/FluMut still complete |
| `PAFTOOLS` | `blastn_db.fa`, `snpeff/` | — |
| `NEXTCLADE` | `nextclade/` (all eight segments) | — |
| `--merge_flumut_mutations` | `--flumut_db` with `nextclade/` | — |

SnpEff annotation in the short-read workflow is now a property of the bundle rather than
of a flag. Where you previously omitted `--snpeff_config` to skip it, use a bundle with
no `snpeff/` directory.

Missing components fail at the start of the run, before any task is submitted, naming
the path that was expected.

## `db.yml`

Optional, and every key within it is optional. It exists so the reference set is
identified by what it *is* rather than by where it happens to sit on disk, and so the
database version reaches the provenance record.

```yaml
name: gsgd                       # default: the bundle directory's name
version: "2025-07-07"            # default: "unknown"
description: GsGD/1996 H5N1 reference set
reference_strain: A/goose/Guangdong/1/1996
nextclade_ref_name: gsgd         # label used in the coordinate map; default: `name`
snpeff_genomes: [influenza_h5n1] # recorded for provenance; not used for dispatch
source: NCBI RefSeq NC_007357.1–NC_007364.1
```

Quote `version`. Unquoted `2025-07-07` parses as a *date*, not a string, in both the
YAML reader the pipeline uses and the one that reads the provenance file back.

`nextclade_ref_name` is what appears as `source_reference` / `target_reference` in
`aggregate/reference_remap/coordinate_map.json`. It is separate from `name` so a bundle
can be renamed or relocated on disk without changing labels already written into
published coordinate maps.

## Why SnpEff data lives inside the bundle

A SnpEff config carries `data.dir`, an absolute path to wherever the database was built.
The pipeline passes `-dataDir <staged snpeff/data>` on every invocation, which overrides
it — so the config's `data.dir` is inert at runtime and the bundle is self-contained:
copy it, hand it to a colleague, or mount it in a container and it still works.

Set `data.dir = ./data` in the shipped config anyway, as documentation of where the data
is meant to be. `scripts/validate_db.py` warns if it points at an absolute path that does
not exist, which is a reliable sign a bundle was moved without being rebuilt.

## Building a bundle

**Nextclade datasets** — `scripts/build_dataset.py` builds per-segment datasets from
GenBank records; see its `--help`. The result goes in `<db>/nextclade/`.

**BLAST reference** — a plain nucleotide FASTA with `accession|subtype|segment` headers.
The pipeline runs `makeblastdb` itself, so no index files are needed alongside it.

```
>NC_007362.1|H5N1|HA
ACGT…
```

**SnpEff database** — one `genes.gbk` per genome, then build in place:

```bash
mkdir -p <db>/snpeff/data/influenza_h5n1
cat pb2.gb pb1.gb pa.gb ha.gb np.gb na.gb m.gb ns.gb \
    > <db>/snpeff/data/influenza_h5n1/genes.gbk
# <db>/snpeff/snpEff.config declares influenza_h5n1.genome and its chromosomes
snpEff build -genbank -v influenza_h5n1 \
    -config  <db>/snpeff/snpEff.config \
    -dataDir <db>/snpeff/data
```

The genome key must be `influenza_<subtype>` in lowercase — both annotation processes
derive it from the `subtype` field of the selected BLAST reference's header, so a
database built as `influenza_H5N1` is never found.

## Validating a bundle

```bash
scripts/validate_db.py /path/to/refs/gsgd
scripts/validate_db.py /path/to/refs/gsgd --require snpeff --json
```

It opens the files rather than just checking that they exist: BLAST headers parse, every
genome named in `snpEff.config` has a built `snpEffectPredictor.bin`, each segment
dataset carries its reference and annotation, `db.yml` is well formed. Components absent
from the bundle are reported `SKIP`, not as failures; `--require` makes a named component
mandatory. Exit code 0 is valid, 1 is not.

The case worth running it for is a `snpEff.config` that declares a genome with no built
database under `data/`: SnpEff then emits an empty `ANN=` on every record instead of
failing, and nothing downstream notices.

## Migrating an existing reference tree

If your references are still laid out as `<db>/snpEff.config` beside a separate SnpEff
data tree, per bundle:

```bash
mkdir -p <db>/snpeff/data
mv <db>/snpEff.config <db>/snpeff/snpEff.config
cp -r /path/to/snpeff/<name>/influenza_h5n1 <db>/snpeff/data/
sed -i 's|^data.dir.*|data.dir = ./data|' <db>/snpeff/snpEff.config
scripts/validate_db.py <db>
```

`nextclade/` and `blastn_db.fa` are already where they need to be. `db.yml` can be added
whenever convenient.

Point `--flumut_db` at the FluMut *bundle root*, not at its `nextclade/` subdirectory —
the reference label comes from the bundle's metadata rather than from the parent
directory's name. Passing the dataset directory reports the parent to use instead.
