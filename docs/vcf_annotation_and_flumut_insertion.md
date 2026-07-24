# VCF annotation & FluMutDB marker annotation

How the Nextclade-branch scripts turn consensus differences into an annotated,
freebayes-style VCF and tag the curated FluMutDB markers among them:

| Script | Role |
| --- | --- |
| [`bin/nextclade_to_freebayes_vcf.py`](../bin/nextclade_to_freebayes_vcf.py) | Convert a Nextclade JSON into a freebayes-format VCF with snpEff-style `ANN=` annotations. |
| [`bin/merge_flumut_mutations_v2.py`](../bin/merge_flumut_mutations_v2.py) | Add a `FLUMUT=` INFO field to rows whose AA change matches a FluMutDB marker. **Current version** — run by the `merge_flumut_mutations` process. |
| [`bin/merge_flumut_mutations.py`](../bin/merge_flumut_mutations.py) | v1, superseded. Merged flumut's per-sample calls and *inserted* rows for markers fluscan had not called. Retained for reference; no process invokes it. |

All three share the VCF I/O helpers in [`bin/tools.py`](../bin/tools.py).

> Both `coordinate_conversion` and `merge_flumut_mutations` are optional NEXTCLADE-workflow
> steps, gated behind `--merge_flumut_mutations` (default `false`). When disabled, the
> workflow stops at `*_final.vcf` and `--flumut_reference_dir` is not required. When enabled,
> `--flumut_reference_dir` becomes mandatory and both steps run as described below.

```
consensus.fasta ─► nextclade_run_batch ─► nextclade_to_freebayes_vcf.py ─► *_nc.vcf
                                                                             │
                                              fix_mnv_codons.py (post_process) ─► *_final.vcf
                                                                             │
      coordinate_conversion (flumut→Nextclade map) ─┐                        │
      FluMutDB (flumutdb conda package) ────────────┼─► merge_flumut_mutations_v2.py ─► *_flumut_merged.vcf
                                                     │
                              *_final.vcf ───────────┘
```

---

## Shared output format

Every data row is a **freebayes-shaped** VCF record. Because the calls come from a
consensus sequence (not read pileups), the depth/quality INFO and FORMAT fields are
**fixed placeholders** — the biology lives entirely in `REF`/`ALT`/`POS` and the `ANN=`
annotation.

- `QUAL` is a constant (default `60`).
- `INFO` carries the full freebayes field set as filler (`DP=1;AF=1.0;…`) plus `TYPE`,
  `CIGAR`, `LEN`, and the `ANN=` annotation (`bin/nextclade_to_freebayes_vcf.py:_info`).
- `FORMAT`/sample is a haploid alt genotype (`GT:DP:AD:… → 1:1:0,1:…`, `_sample`).

### The `ANN=` field

snpEff-compatible, pipe-delimited, comma-separated when a variant hits multiple features
(overlapping influenza ORFs). Sixteen fields:

```
Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type |
Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos/len |
CDS.pos/len | AA.pos/len | Distance | ERRORS/WARNINGS/INFO
```

Example (HA missense):

```
ANN=C|missense_variant|MODERATE|HA|FLUAVH5N1_s4gp1|transcript|FLUAVH5N1_s4gp1|protein_coding|1/1|c.9A>C|p.K3N|9/1707|9/1707|3/568||
```

---

## `nextclade_to_freebayes_vcf.py`

### Purpose

Translate Nextclade's per-sequence alignment result (its JSON) into a VCF whose annotations
match the snpEff output the paftools branch produces — so both branches yield comparable
VCFs. It resolves CDS structure directly from the dataset's GFF3, handling overlapping ORFs
(M1/M2, NS1/NEP, PB1/PB1-F2, PA/PA-X) and spliced CDSes (M2, NEP, PA-X) natively.

### CLI

```
nextclade_to_freebayes_vcf.py -i <nextclade.json> -d <dataset_dir> -o <out.vcf>
                              [--sample NAME] [--qual 60]
```

`-d` must contain `reference.fasta` and `genome_annotation.gff3` (the same files Nextclade
aligned against). `-o` is used as a filename prefix when the JSON holds multiple sequences.

### How it runs (`main`)

1. **Load reference** (`load_reference`) — single-record FASTA → `(chrom_id, sequence)`.
   `chrom_id` becomes the VCF `#CHROM`.
2. **Load annotation** (`load_gff3`) → `(genes, transcripts)`:
   - `genes`: name, id, half-open `[begin_0, end_0)`.
   - `transcripts`: one per CDS **and** per mature-peptide region. Each carries `segments`
     (a sorted list of half-open exon intervals — this is what makes spliced/overlapping
     ORFs work), `strand`, `cds_length`, `aa_length` (`cds_len//3 − 1`, i.e. excluding the
     stop codon), `feature_id`, `protein_id` (matches Nextclade's `cdsName`), `is_mature`,
     and any `warnings`.
3. **Iterate JSON `results`** — one sequence at a time. Sequences whose result carries
   `errors`, or that appear under the top-level `errors` array (e.g. *"Unable to align"*),
   produce an **empty VCF** (header only) with a warning rather than failing the run.
4. **Emit variants** (`variants_from_json_result`) and **write** (`write_vcf`). One output
   file per sequence; single-sequence input writes directly to `-o`.

### Turning JSON deltas into rows (`variants_from_json_result`)

Consumes three delta lists from the Nextclade result, plus `nucToAaMuts`, `aaDeletions`,
`aaInsertions`, and `frameShifts` for annotation:

| JSON source | VCF handling | Anchor / allele helper |
| --- | --- | --- |
| `substitutions` (`pos`, `refNuc`, `qryNuc`) | one SNV row, `POS = pos+1` | direct `REF/ALT` |
| `deletions` (`range.begin/end`) | left-anchored to the preceding base | `_deletion_alleles` |
| `insertions` (`pos`, `ins`) | anchored to `pos`; terminal overhangs (`pos<0` or `pos ≥ len−1`) are **dropped** to avoid edge false-positives; `FILTER=AmbiguousInsert` if the inserted bases contain `N` | `_insertion_alleles` |

Rows are sorted by `POS` before return.

### Substitution annotation decision tree (`_ann_field_for_sub`)

For each substitution the script walks every transcript and appends an `ANN` entry per
feature it touches (so a base in an overlap region yields several entries):

1. **Inside a main CDS** — take the AA change from `nucToAaMuts` when present, otherwise
   recompute it by translating the reference codon (`_synonymous_aa_change`). Classify the
   consequence (`_classify_aa_consequence`: `synonymous_variant` / `missense_variant` /
   `stop_gained` / `stop_lost` / `start_lost`, with LOW/MODERATE/HIGH impact) and format
   `HGVS.c` (`c.<n><ref>><alt>`) and `HGVS.p`. `Rank` = exon/exon-count.
2. **Intronic** (a gap between the exons of a spliced transcript) — `intron_variant`,
   MODIFIER, `HGVS.c` in `c.<cds_end>+<offset><ref>><alt>` form.
3. **Mature-peptide region** — same AA logic, tagged with
   `WARNING_TRANSCRIPT_NO_START_CODON`; positions 5′ of the region become
   `upstream_gene_variant` with a `Distance`.
4. **Upstream of a main CDS** (`pos < CDS start`) — `upstream_gene_variant`, MODIFIER,
   `c.-<distance>…`.
5. **Intergenic** (outside every gene) — `intergenic_region`, MODIFIER, `n.<pos>…`.

Deletions and insertions have analogous builders (`_ann_field_for_del` / `_ann_field_for_ins`)
that distinguish in-frame indels (`conservative_inframe_deletion/insertion`) from
`frameshift_variant` (HIGH) using `aaDeletions`/`aaInsertions`/`frameShifts`.

### Header (`build_header`)

`##fileformat=VCFv4.2`, `PASS` + `AmbiguousInsert` filters, `fileDate`, `source`,
`reference`, a single `##contig`, the complete freebayes `##INFO`/`##FORMAT` definitions,
and the `##INFO=<ID=ANN…>` definition.

### Known caveat — multi-nucleotide codons

Each SNV is annotated independently, assuming the rest of its codon matches the reference.
When two consensus SNVs fall in the same codon, their individual `HGVS.p` values can be
wrong. This is corrected downstream by `fix_mnv_codons.py` (the `post_process_vcf*`
processes), which recomputes the real codon before the `*_final.vcf` is produced.

---

## `merge_flumut_mutations_v2.py`

### Purpose

FluMutDB is a curated catalogue of influenza amino-acid markers and their published
phenotypic effects. This script tags the fluscan calls that *are* such markers: it walks the
VCF and, for every AA change already annotated in a row's `ANN=`, asks the database whether a
marker exists for that (protein, AA position, alt AA). It **writes a new file**
(`*_flumut_merged.vcf`) and never touches the input `*_final.vcf`.

**No row is added, removed or reordered**, and no existing field is rewritten — matching rows
simply gain one extra INFO key. `ANN=` stays strictly snpEff 16-field compatible.

> **Consequence worth knowing:** a marker is only reported where fluscan called a *difference
> from the Nextclade reference*. Markers already fixed in that reference produce no VCF row
> and so cannot be annotated. Against `chicken2021`, 10 of FluMutDB's 71 HA1 markers fall in
> this category — `flumut` itself still reports them, because it compares the consensus
> against FluMutDB's own (older) reference. The two tools answer different questions.

### The `FLUMUT=` field

```
##INFO=<ID=FLUMUT,Number=.,Type=String,Description="FluMutDB marker matches: 'Protein |
AA_Pos | Alt_AA | Mutation_Name | Effects'. '&' separates multiple values within a field;
',' separates one entry per matching ANN annotation. AA_Pos is in this VCF's coordinates.">
```

```
FLUMUT=PB2|715|N|PB2:S715N|Decreased_virulence_in_mice
FLUMUT=HA1|222|L|HA1-5:Q222L|Decreased_antiviral_response_in_host&Enhanced_replication_in_ferrets&…
```

One comma-separated entry per matching `ANN` annotation, so a row hitting two features can
carry two entries. `AA_Pos` is the position in *this* VCF's coordinates; the FluMutDB
reference position is carried inside `Mutation_Name`. Effect names are sanitised (whitespace
→ `_`, `;,=|` stripped, en dash → `-`) because VCF INFO values may not contain whitespace,
`;` or `,`. `--no-effects` emits the first four sub-fields only. Rows with no match get **no**
`FLUMUT` key.

### CLI

```
merge_flumut_mutations_v2.py --vcf <final.vcf> --coord-map <coordinate_map.SEG.json>
                             --segment <HA|M|…> --sample <id> -o <merged.vcf>
                             [--flumut-db <path>] [--no-effects]
```

### Inputs and where they come from

| Input | Producer | Used for |
| --- | --- | --- |
| `--vcf` (`*_final.vcf`) | `post_process_vcf_batch` | the AA changes to look up, and the rows to annotate |
| `--coord-map` (`coordinate_map.<SEG>.json`) | `coordinate_conversion` | flumut→Nextclade **AA position** mapping, per gene |
| FluMutDB sqlite | `flumutdb` conda package | the marker catalogue; located via `flumut.DbReader.get_db_file()`, overridable with `--flumut-db` |

Wired in `workflows/nextclade.nf` per `(sample, segment)`: the per-segment VCF batches are
flattened and combined with the matching coordinate-map shard. The per-sample `flumut` TSV is
**not** an input — the marker set comes from the database, so the `flumut` process and this
one are independent.

### How it runs (`main`)

1. **Read the VCF** with `tools.read_vcf` → DataFrame + header lines.
2. **Load this segment's block** from the coordinate map (`load_segment_map`), warning if it
   carries `stale_annotation`.
3. **Pair genes with FluMutDB proteins** (`build_gene_protein_maps`): query the distinct
   `protein_name`s the database maps onto this segment's reference, resolve each to a
   coordinate-map gene key (`resolve_gene`, see below), and **invert** that gene's
   `aa_source_to_target` once into `{target_pos: source_pos}`. Inverting is what lets the
   main loop start from a VCF row and query the database directly. Alignment-derived maps
   pair columns 1:1, so the inverse is well defined; a collision would be warned about and
   resolved to the lowest source position.
4. **Per row, per `ANN` entry** (`row_flumut_field` → `aa_changes_from_info`): take
   `Gene_Name` and the `HGVS.p` substitution, skipping synonymous changes and entries with no
   `HGVS.p`. Translate the position target→source, then run one parameterised query
   (`MARKER_QUERY`, `lru_cache`d) for `(reference, protein, position, alt)`.
5. **Write** with `tools.write_vcf`, inserting the `##INFO=<ID=FLUMUT…>` definition
   immediately before the `#CHROM` line.

`ref_seq` is deliberately **not** part of the match: a marker is defined by its alt allele,
and the sample's reference amino acid may legitimately differ from FluMutDB's.

### Gene-name resolution (`resolve_gene`)

FluMutDB protein names differ from some Nextclade annotations. Rather than hard-code
assumptions, the resolver generates candidates and keeps whichever exists as a real gene key
in the coordinate map for that segment:

| FluMutDB `protein_name` | candidates tried | resolves to (if present) |
| --- | --- | --- |
| `HA1` / `HA2` | `HA1` / `HA2` | `HA1` / `HA2` |
| `NS-1` | `NS-1`, `NS`, `NS1` | `NS1` |
| `NS-2` | `NS-2`, `NS`, `NS2`, `NEP` | `NEP` |
| `M1` / `M2` | `M1` / `M2` | `M1` / `M2` |
| `PB1-F2` | `PB1-F2` | `PB1-F2` |

The segment label maps to FluMutDB's `reference_name` verbatim except **M → `MP`**.

### Edge cases (all non-fatal)

- **HA without the subunit split** — FluMutDB numbers HA by its HA1/HA2 subunits, so the
  Nextclade dataset must expose dedicated `HA1`/`HA2` genes. A dataset with a single `HA`
  gene yields VCF rows whose `Gene_Name` is `HA`, which resolves to nothing; the script warns
  (`VCF gene(s) HA … have no FluMutDB counterpart`) and passes the VCF through unannotated.
- **Unmapped AA position** — a VCF position with no counterpart in flumut's reference
  (alignment gap) is debug-logged and skipped.
- **Missing gene maps** — a segment absent from the map, or one whose `genes` block is empty,
  warns and passes through.
- **Empty-variant VCF** — handled; the `FLUMUT` header is still added.
- FluMutDB holds **substitutions only** (`mutations.type` is `SNP` throughout), so indel and
  non-coding rows never match.

### v1 (`merge_flumut_mutations.py`, superseded)

v1 read flumut's melted per-sample TSV, labelled every row `DETECTED_BY=fluscan|flumut|both`,
and **inserted** synthetic rows (filler `REF=N`/`ALT=N`, `POS` derived from the target CDS)
for markers fluscan had not called. Those rows carried no nucleotide evidence. v2 drops both
the insertion behaviour and the `DETECTED_BY` field: with the per-sample TSV gone, the script
cannot say what flumut detected — only that a fluscan call coincides with a curated marker,
which is what `FLUMUT=` states. The file is retained but no process invokes it.

---

## Shared VCF I/O — `tools.py`

`read_vcf(path)` → `(DataFrame, header_lines)` and `write_vcf(df, header_lines, path)` are
the inverse of one another and round-trip a VCF **byte-for-byte**. Key properties:

- Transparent `.gz` read **and** write.
- Header scan raises `ValueError` if there is no `#CHROM` line, and tolerates very large
  headers (many `##contig`/`##INFO` lines).
- Data columns are read as **strings with NA-detection off**, so `.`, `NA`-like tokens, and
  empty fields survive unchanged and there are no `NaN` surprises. `POS` is therefore a
  string — cast with `.astype(int)` when you need integers. This verbatim round-trip is what
  lets `merge_flumut_mutations_v2.py` guarantee that only `INFO` changes.
