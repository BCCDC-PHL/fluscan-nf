# fluscan-nf agent guide

This file is persistent, repository-scoped context for coding agents. It applies
to the entire repository. Keep it concise and update it whenever a change makes
one of these durable facts false.

## Project purpose

`fluscan-nf` is a Nextflow DSL2 pipeline for Influenza A mutation analysis. It
accepts either paired-end Illumina reads or one multi-segment consensus FASTA
per sample. The public output layout, sample identifiers, segment identifiers,
VCF fields, and QC tables are interfaces; preserve them unless a change is
explicitly intended and documented.

The canonical segment set is `PB2`, `PB1`, `PA`, `HA`, `NP`, `NA`, `M`, and
`NS`.

## Mental model

- `main.nf` dispatches by sentinel-valued parameters, and exactly one of
  `--fastq_input`, `--fasta_input` and `--samplesheet_input` may be set.
  `--fastq_input` runs `SHORT_READ` and `--fasta_input` launches both `PAFTOOLS`
  and `NEXTCLADE`. For `--samplesheet_input` the header decides: `ID,R1,R2` runs
  `SHORT_READ`, `ID,FASTA` runs the two FASTA branches.
- `SHORT_READ` performs read QC/primer trimming, downsampling, assembly,
  BLAST-based reference selection, mapping, LoFreq calling, SnpEff annotation,
  MNV correction, consensus generation, and FluMut screening.
- In FASTA mode, `FASTA_PREP` (`workflows/fasta-prep.nf`) hashes, cleans and
  QC-aggregates the consensus input and runs FluMut once, then hands the cleaned
  FASTAs and their provenance records to both branches.
  - `PAFTOOLS` selects references with BLAST, aligns each segment with
    minimap2, calls differences with `paftools.js`, sorts and annotates VCFs,
    converts them to the common VCF shape, and corrects MNV codons.
  - `NEXTCLADE` batches samples by segment, runs Nextclade once per segment,
    splits per-sample results, and converts JSON to the common VCF shape,
    annotating substitutions that share a codon with the combined change.
- `--merge_flumut_mutations` is off by default. When enabled in the Nextclade
  branch, a run-level map lifts coordinates from FluMut's own reference
  (source, the `references`/`annotations` tables of the installed `flumutdb`
  sqlite) onto the selected bundle's Nextclade datasets (`--db`, target), then
  `bin/merge_flumut_mutations.py` adds `FLUMUT=` annotations to matching
  existing VCF rows. The markers and the source coordinates come from the same
  package, so they cannot drift apart.
- `mutation_watch` scans each final VCF against the custom watchlists in
  `--watchlist_dir` (default `assets/watchlists/`) and publishes exactly one
  file, `*_watchlist_hits.tsv` (header-only when nothing matches). It runs in all
  three workflows and is independent of `--merge_flumut_mutations`.

## Repository map

| Path | Role |
| --- | --- |
| `main.nf` | Entry point and workflow dispatch |
| `nextflow.config` | Parameters, Conda profile, resources, reports |
| `workflows/` | Channel wiring for the three top-level workflows and the `FASTA_PREP` subworkflow both FASTA branches take their input from |
| `modules/shared/` | Reused BLAST, consensus, FluMut, post-processing, and provenance processes |
| `modules/{short_read,paftools,nextclade}/` | Branch-specific processes |
| `bin/` | Python CLI entry points *only* — Nextflow puts this directory on `PATH` for every task, so nothing that is merely importable belongs here |
| `bin/lib/` | Modules the entry points import: `tools.py` (VCF I/O), `flu_gene_names.py` (gene-key aliases and the known-gene set), and `aa_markers.py` (`ANN=` extraction, HGVS.p parsing, the `MNV_UNRESOLVED` test). Reached as `from lib.… import …`; Python already puts a script's own directory on `sys.path`, so no packaging or `__init__.py` is involved |
| `environments/` | Pinned Conda runtime definitions |
| `tests/` | Python unit and CLI tests |
| `docs/vcf_annotation_and_flumut_insertion.md` | Detailed Nextclade VCF and FluMutDB semantics |
| `docs/watchlists.md` | Custom watchlist schema, matching rules, and the HA1/HA numbering trap |
| `assets/watchlists/` | The shipped `*watchlist.csv` files; `--watchlist_dir` defaults here |
| `modules/shared/sample_id.nf` | The only place a consensus or sample sheet `sample_id` is derived; holds the FASTA discovery channel `FASTA_PREP` calls and `samplesheet_channel`, which all three workflows call for `--samplesheet_input`. Changing `clean_sample_id` renames every output |
| `modules/shared/database.nf` | The only place that knows the `--db` bundle layout; `resolve_db` validates it up front. Also holds the canonical segment list, `@Field DB_SEGMENTS`, which it uses to check each segment's dataset |
| `scripts/validate_db.py` | Deeper, contents-level bundle check; run after building or moving a bundle |
| `docs/database_bundles.md` | Bundle layout, `db.yml` schema, how to build and validate one |
| `README.md` | User-facing inputs, commands, outputs, and architecture |
| `assets/` | Primers and reference/database assets; do not rewrite casually |

## Data and interface contracts

- The consensus `sample_id` is derived once, from the input filename or the
  sample sheet `ID` (in either mode), by
  `clean_sample_id` in `modules/shared/sample_id.nf`, and is final from there
  on: it keys every channel and names every file. Nothing downstream may
  re-derive or edit it. `consensus_fasta_channel` and `samplesheet_channel` abort
  the run if two inputs collapse to one ID; `samplesheet_channel` also rejects an
  ID that still holds anything but letters, digits, `.` and `-` after cleaning
  (`|` would break the SHORT_READ `segment_regex`). The underscores-to-dashes rule is load-bearing:
  the `saveAs` in `modules/nextclade/nextclade.nf` and the `sample_name` column in
  `bin/filter_alignments.py` both read the ID back as everything before the first
  underscore.
- Accepted consensus headers include the legacy `SAMPLE_segmentN_SEGMENT`,
  canonical `SAMPLE_SEGMENT`, and GISAID-like pipe forms. Headers supply only
  the segment label, read as the rightmost `|`/`_` token naming a segment
  (`MP` maps to `M`). Preparation rewrites every header to
  `{sample_id}_{SEGMENT}`. One file must hold one sample, and two things are
  fatal because they mean it does not: two records claiming one segment, and
  leading header fields (first pipe field, else first underscore token) that name
  more than one sample. Both exit 2 under `errorStrategy 'ignore'`, so the sample
  leaves the run entirely — no QC row in the aggregate summary, no provenance.
- Directory-discovered FASTQ sample IDs currently stop at the first
  underscore, and colliding pairs are dropped silently by `.unique`; sample sheet
  FASTQ IDs go through `clean_sample_id` instead. Avoid introducing or silently
  accepting colliding names.
- All reference data comes from one bundle directory, `--db`, in the fixed
  layout `blastn_db.fa` + `nextclade/` + `snpeff/{snpEff.config,<genome>/}` plus an
  optional `db.yml`. Resolve it through `resolve_db()` in
  `modules/shared/database.nf` rather than composing paths at the call site;
  that function is the only thing that should know the layout.
- The nucleotide reference FASTA used for BLAST has headers
  `accession|subtype|segment`.
- Nextclade datasets live in uppercase per-segment subdirectories (`HA/`, `PB2/`,
  …) and must contain the expected dataset reference and annotation files.
  Resolve them with the segment label directly; never lowercase it. macOS is
  case-insensitive, so a case mistake here fails only on Linux.
- A dataset `config.json` is the **input manifest for
  `scripts/build_dataset.py`** and nothing else. That script is its only
  consumer; the pipeline, `bin/build_reference_coordinate_map.py`, and
  `scripts/biopython_sanity_check.py` all key off the fixed uppercase layout
  plus each segment's `pathogen.json`. Do not reintroduce a runtime read of it —
  a bundle without one must still work.
- SnpEff is always invoked with `-dataDir` pointing at the staged `snpeff/`,
  which overrides the config's `data.dir`. The `<genome>/` directories sit
  directly under `snpeff/` — a standalone SnpEff install nests them under
  `data/`, a bundle does not. Stage the whole `snpeff/` directory into a
  process, never the config alone, or the task reads a database that only
  exists on the machine the bundle was built on.
- Bundles may be partial and each workflow declares what it requires:
  SHORT_READ `blastn` (SnpEff and everything downstream of it runs only when
  `db.has_snpeff`), PAFTOOLS `blastn` + `snpeff`, NEXTCLADE `nextclade`. FASTA
  dispatch launches both branches, so it needs all three. The optional FluMut
  merge needs nothing beyond `nextclade`.
- Most unset path parameters use the literal `NO_FILE`. Do not treat that
  string as a real path. `nextflow.config` leaves `db` unset, so pass it
  explicitly unless a local config supplies it.
- `--samplesheet_input` is parsed eagerly by `samplesheet_channel`, so column,
  blank-cell, ID and file-existence errors all fail before any task is submitted.
  `FASTA_PREP` chooses between it and `consensus_fasta_channel`. `run_name` is the sheet's filename without `.csv`.
- Channel tuple keys are part of the design. Preserve `sample_id` and `segment`
  through `join`, `combine`, `transpose`, and batching operations, and comment
  non-obvious tuple shapes next to the transformation.
- Preserve `publishDir` destinations and filenames when refactoring. Outputs
  are primarily `${outdir}/${sample_id}/...`, with run-level data under
  `${outdir}/aggregate`. Several joins and parsers depend on double-underscore
  VCF stems and the `_paf`, `_nc`, and `_final.vcf` suffixes.
- Identity and coverage QC are report-only. Some processes intentionally
  ignore errors, so a successful run status does not prove every expected
  sample and segment was published; check completeness explicitly.

### VCF invariants

- Consensus-derived Nextclade VCFs are deliberately “freebayes-shaped.”
  Depth, quality, and genotype fields are placeholders; biological meaning is
  carried by coordinates, alleles, and the 16-field SnpEff-compatible `ANN=`.
- In the SnpEff-annotated PAFTOOLS and SHORT_READ VCFs, individual SNVs in one
  codon can have incorrect amino-acid annotations until `bin/fix_mnv_codons.py`
  runs; their final VCFs must remain downstream of that step. The Nextclade
  converter translates each codon with all its substitutions applied, keeping
  Nextclade's own `nucToAaMuts` call only for a codon with one substitution, and
  writes the same `MNV_CORRECTED`/`MNV_ORIGINAL_HGVS_P` fields itself. Nextclade's
  peptide alignment can gap a codon its nucleotide alignment keeps, so its call
  must not override a shared codon.
- `bin/merge_flumut_mutations.py` is the only implementation. It may add
  `FLUMUT=` to matching rows, but must not add, remove, reorder, or otherwise
  rewrite variant rows. The earlier implementation that inserted synthetic rows
  is deleted, and this script has taken back the plain name (it was briefly
  `merge_flumut_mutations_v2.py`); no `_v2` suffix should reappear. The
  coordinate map's `target_cds` block existed only for the old script and is
  gone too.
- The merge reads only Nextclade VCFs, which never carry `MNV_UNRESOLVED`, so it
  has no skip for that flag. Every row of a codon holding more than one substitution
  carries the same combined HGVS.p and so the same `FLUMUT=` entry — downstream
  marker counts must de-duplicate.
- **Watchlist hits never enter the VCF.** `FLUMUT=` comes from a versioned
  database, so the same consensus always yields the same annotation; a watchlist
  is a CSV we edit, and merging its hits in would make the same consensus
  produce different VCFs depending on that file. Hits go to
  `*_watchlist_hits.tsv`. The shared part — `ANN=` parsing, both HGVS.p
  spellings, gene canonicalisation, the watchlist's `MNV_UNRESOLVED` skip — lives once in
  `bin/lib/aa_markers.py`; do not let a consumer grow its own copy.
- A watchlist `gene` is a **protein**, not a segment, and gene resolution falls
  back from ANN `Gene_Name` to `Feature_ID` but **never** from a parent gene to
  one of its subunits. `HA` and `HA1` are not interchangeable in either
  direction. A watchlist gene that matches no annotated gene must be reported,
  not ignored — that silence is what hid 88 unmatchable `HA` rows in an earlier
  list.
- `FLUMUT=` describes curated markers among differences from the selected
  Nextclade reference. A marker already fixed in that reference has no VCF row
  to annotate; the independent per-sample FluMut reports answer a different
  question.
- `FLUMUT=` is `Protein|AA_Pos|Alt_AA|Mutation_Name|DirectEffect|AssociatedEffect`.
  `DirectEffect` carries effects of single-mutation markers; `AssociatedEffect`
  carries effects of markers needing further mutations and does **not** assert
  that the whole marker is present. Both hold `Subtype:Effect` values, `.` when
  empty. Do not merge the two classes or drop the subtype prefix — that is the
  over-claim the split exists to prevent.
- FluMutDB positions are translated source-to-target through the coordinate
  map. Never reverse the meanings of the FluMut source reference and the
  Nextclade target reference.
- Gene-key aliases (`NS-1`→`NS1`, `NS2`/`NS-2`→`NEP`) live only in
  `bin/lib/flu_gene_names.py`. The coordinate map builder and the merge must both
  canonicalise through it; a second copy lets the map be keyed under a name the
  merge never looks up, silently dropping a whole gene's markers.
- HA marker annotation needs the Nextclade annotation to expose the HA1/HA2
  split, either as two real `CDS` features (chicken2021) or as
  `mature_protein_region_of_CDS` children of one HA CDS (the NCBI/RefSeq model,
  GsGD included). A dataset offering only a bare `HA` CDS is still unusable for
  HA markers and is flagged `stale_annotation`.
- A mature-peptide ANN entry names the *parent* gene in `Gene_Name`, so the
  merge falls back to `Feature_ID` and the map records the peptide's
  `protein_id`/GFF3 `ID` as gene `aliases`. The parent is listed in the
  segment's `subunit_parent_genes` so it is not reported as an unmapped gene on
  every subunit row. Never make a parent gene resolve to a subunit's FluMutDB
  protein: HA1 numbering starts after HA's 16-residue signal peptide, so that
  returns a wrong residue at 93% of positions, silently.
- A mature peptide is only usable as a gene key when it starts on a codon
  boundary of its parent CDS — it is translated from its own first base, so a
  mid-codon start yields a frameshifted protein that still aligns and still maps
  positions, reporting wrong residues with nothing raised. A peptide that fails
  that check, or whose `Parent` cannot be resolved so it cannot be checked, must
  be skipped rather than used.

## Development conventions

- Keep orchestration in DSL2 workflows/modules and put substantial parsing or
  transformation logic in a testable Python CLI under `bin/`. `bin/` is on
  `PATH` in every task, so it holds entry points only; anything shared between
  scripts goes in `bin/lib/` and is imported as `from lib.<module> import …`.
- Follow existing process conventions: lower-snake-case names, a useful `tag`,
  explicit named `emit` outputs, the appropriate file in `environments/`, and
  narrowly scoped publishing.
- Prefer typed Python, `argparse`, `pathlib`, deterministic ordering, clear
  stderr logging, and a `main()` entry point. Preserve VCF/TSV values as
  strings when round-trip fidelity matters.
- Add or update focused tests for parsing, coordinate conversion, annotations,
  QC thresholds, empty inputs, and failure modes. Tests should not require
  network access.
- Pin runtime dependencies in the relevant Conda YAML rather than installing
  them ad hoc inside a process.
- When changing inputs, parameters, output paths, or scientific semantics,
  update `README.md` and the relevant file in `docs/` in the same change.
- Inspect `git status` before editing and preserve unrelated user work. Do not
  commit, push, or remove work/results unless explicitly requested.

## Validation

Run the smallest relevant checks first:

```bash
git diff --check
nextflow config .
python -m pytest -q
```

The FluMut merge tests use the real `flumutdb` package. They skip if the
database is unavailable; set `FLUMUT_DB=/path/to/flumut_db.sqlite` to test with
an explicit database. `pytest` is a development dependency and is not
currently declared in the runtime Conda files, so no single environment can run
the suite as written: `fluscan` ships `flumutdb` but no `pytest`, and `server`
ships `pytest` but no `flumut`. Bridge them rather than installing anything —
without `FLUMUT_DB` the 71 merge tests silently skip:

```bash
FLUMUT_DB=~/miniforge/envs/fluscan/lib/python3.12/site-packages/flumutdb/flumut_db.sqlite \
  ~/miniforge/envs/server/bin/python -m pytest -q
```

There is no self-contained, lightweight integration fixture committed for all
three branches. For workflow changes, run a targeted local smoke test with the
required reference assets and parameters, then inspect both task success and
published filenames/content. Use `-profile conda`; `conda.useMamba` is enabled.

## Known implementation cautions

- The short-read workflow currently maps the original FASTQ channel, while
  trimmed reads feed QC and assembly. Treat changing that as a scientific
  behavior change, not a cosmetic refactor.
- In short-read mode, SnpEff, MNV correction, and watchlist scanning run only
  when the bundle carries a `snpeff/` component; consensus and FluMut can still
  complete without those annotated final VCF steps.
- Sample-level provenance collection is wired in all three workflows, and each
  starts with `hash_files` (`modules/shared/hash_files.nf`, aliased to
  `short_read_hash_files` in SHORT_READ, run by `FASTA_PREP` for the FASTA
  branches). The collection chains are inner joins,
  and `prepare_consensus`, `blastn` and `paftools_align` all run
  `errorStrategy 'ignore'` — so a sample failing any of them publishes no
  provenance file, losing its input hash along with everything else. Use
  `join(..., remainder: true)` with a null guard
  (the pattern in `workflows/nextclade.nf`) if that record must survive.
- `collect_provenance` puts `task.ext.workflow` in the published filename.
  FASTA dispatch runs PAFTOOLS and NEXTCLADE over the same samples into the same
  `${outdir}/${sample_id}/`, and their timestamps collide at second resolution;
  without the label one branch silently overwrites the other. Any new workflow
  that collects provenance needs a `withName: '<NAME>:.*' { ext.workflow = … }`
  block or its records fall back to the unlabelled, collision-prone name.
- `mutation_watch` runs once per sample over all of its final VCFs and every
  watchlist, and emits one provenance YAML per sample. SHORT_READ does not join
  it in, so its watchlist scan is absent from the collected record. Both FASTA
  branches join it with `remainder: true`.
- The two FASTA branches share their input through `FASTA_PREP`, and each
  opens its provenance chain with its own `pipeline_provenance` record ahead of
  the prep records. Their watchlist blocks and provenance joins are still
  mirrored; a change to one belongs in both. NEXTCLADE also combines each
  per-segment record (`nextclade_run_batch`, `merge_flumut_mutations`) onto
  every sample.
- `tests/conftest.py` puts `bin/` on `sys.path` so the entry points can resolve
  their `from lib.…` imports under pytest. Fix an import failure there, not by
  adding path manipulation to a production script.
- `short_read_sample_reads` invokes `seqtk`, but the main environment does not
  explicitly declare it; do not rely on a transitive installation.
- The GitHub Actions test path is legacy and short-read-only. It references an
  absent output-checking script, and `--db` now needs a bundle directory that
  nothing in CI provisions; reconcile it with current code before relying on CI.
- The pipeline's process scripts assume GNU userland: `grep -oP`, `sed -i`
  without a suffix, `sed` with `\S`, and bash 4+ `${VAR,,}`. They do not run on
  a stock macOS shell; test workflow changes on Linux or in a container.
- Known documentation drift includes the PAFTOOLS aggregate filename. Executable
  code writes `aggregate/paftools_qc_summary.tsv`; fix code/docs together when
  touching these areas.
