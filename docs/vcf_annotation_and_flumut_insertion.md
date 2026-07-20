# VCF annotation & flumut insertion

How two of the Nextclade-branch scripts turn consensus differences and flumut markers into
an annotated, freebayes-style VCF:

| Script | Role |
| --- | --- |
| [`bin/nextclade_to_freebayes_vcf.py`](../bin/nextclade_to_freebayes_vcf.py) | Convert a Nextclade JSON into a freebayes-format VCF with snpEff-style `ANN=` annotations. |
| [`bin/merge_flumut_mutations.py`](../bin/merge_flumut_mutations.py) | Merge flumut amino-acid markers into that VCF and label every AA change by the tool(s) that found it. |

Both share the VCF I/O helpers in [`bin/tools.py`](../bin/tools.py), and the second script
reuses the annotation helpers of the first, so the two produce identically-shaped rows.

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
      flumut (per-sample AA markers) ───────────────┼─► merge_flumut_mutations.py ─► *_flumut_merged.vcf
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

## `merge_flumut_mutations.py`

### Purpose

flumut screens the consensus for curated influenza amino-acid markers, but its calls never
reached the VCF. This script folds them in: it labels every AA change in the VCF with the
tool(s) that detected it and appends rows for flumut markers Nextclade/fluscan did not call.
It **writes a new file** (`*_flumut_merged.vcf`) and never touches the input `*_final.vcf`.

Every row gains one INFO field:

| `DETECTED_BY` | Meaning |
| --- | --- |
| `fluscan` | called only by the Nextclade branch (includes all non-coding/synonymous rows) |
| `flumut` | a flumut marker not called by fluscan — a **newly inserted** row |
| `both` | the *same exact* AA substitution (gene + position + alt AA) called by both tools |

### CLI

```
merge_flumut_mutations.py --vcf <final.vcf> --coord-map <coordinate_map.json>
                           --flumut-mutations <melted.tsv> --dataset-dir <dataset>
                           --segment <HA|M|…> --sample <id> -o <merged.vcf>
```

### Inputs and where they come from

| Input | Producer | Used for |
| --- | --- | --- |
| `--vcf` (`*_final.vcf`) | `post_process_vcf_batch` | existing fluscan AA changes + rows to relabel |
| `--coord-map` (`coordinate_map.json`) | `coordinate_conversion` | flumut→Nextclade **AA position** mapping, per gene |
| `--flumut-mutations` (melted TSV) | `flumut` + `melt_flumut.py` | the marker calls (`sample`,`gene`,`mutation`,`allele`) |
| `--dataset-dir` | Nextclade dataset | target `genome_annotation.gff3` → CDS start for placing new rows |

Wired in `workflows/nextclade.nf` per `(sample, segment)`: the per-segment VCF batches are
flattened, joined to flumut mutations and the consensus by sample, and combined with the
single coordinate map.

### How it runs (`main`)

1. **Read the VCF** with `tools.read_vcf` → DataFrame + header lines. `#CHROM` comes from
   the data (or the `##contig` header for an empty VCF).
2. **Load this segment's gene maps** from the coordinate map (`load_segment_genes`) and the
   **target transcripts** from the dataset GFF (`load_target_transcripts` → `load_gff3`).
3. **Collect fluscan AA changes** already in the VCF — parse each row's `ANN=` for
   `(gene, aa_pos, alt_aa)` (`row_aa_changes`, skipping synonymous).
4. **Read flumut markers** for this sample (`read_flumut_mutations`, a pandas read of the
   melted TSV filtered to `--sample`).
5. **Map each marker into Nextclade coordinates** (`build_flumut_inserts`):
   - **Resolve the gene name** to a coordinate-map key (`resolve_gene`, see table below).
     A marker that doesn't resolve to this segment is skipped (it belongs to another VCF).
   - **Map the AA position** flumut→Nextclade via `genes[gene].aa_source_to_target`. An
     unmapped (gapped) position is logged and skipped.
   - Form the key `(gene, target_pos, alt_aa)`. If it matches a fluscan change → the marker
     is **both** (the existing row is relabelled; no duplicate is inserted). Otherwise it's
     a **flumut-only** insert.
6. **Build flumut-only rows** (`flumut_only_row`): derive a real `POS` from the target CDS
   (`POS = cds_start + (target_pos−1)·3`, via the reused `_cds_to_genomic`), set filler
   `REF=N`/`ALT=N`, and carry the AA change in `ANN=` (reusing `_build_ann_entry` and
   `variant_row` so the row is shaped exactly like a converter row). A marker that can't be
   placed on the CDS is logged and skipped.
7. **Assemble output** (`assemble_output`): relabel every existing row `both`/`fluscan`,
   append the flumut-only rows, cast `POS` to int, and stable-sort by `POS`.
8. **Write** with `tools.write_vcf`, inserting the `##INFO=<ID=DETECTED_BY…>` definition
   immediately before the `#CHROM` line.

### Gene-name resolution (`resolve_gene`)

flumut names carry a trailing `-<digits>` that is sometimes a FluMutDB version and
sometimes a protein number, and some names differ from the Nextclade annotation. Rather
than hard-code assumptions, the resolver generates candidates and keeps whichever exists as
a real gene key in the coordinate map for that segment:

| flumut name | candidates tried | resolves to (if present) |
| --- | --- | --- |
| `HA1-5` | `HA1-5`, `HA1`, `HA15` | `HA1` |
| `NS-1` | `NS-1`, `NS`, `NS1` | `NS1` |
| `NS-2` | `NS-2`, `NS`, `NS2`, `NEP` | `NEP` |
| `M1` / `M2` | `M1` / `M2` | `M1` / `M2` |
| `PB1-F2` | `PB1-F2` | `PB1-F2` |

### Edge cases (all non-fatal)

- **Stale/mismatched target annotation** — if the Nextclade dataset lacks a split flumut
  needs (e.g. GsGD has a single `HA` gene, not HA1/HA2; or lacks NS1/NEP), the coordinate
  map has no genes for that segment, so those markers are skipped and the VCF passes through
  with `fluscan` labels only. Use a dataset that splits them (e.g. chicken2021) to insert.
- **Unmapped AA position** — a flumut position with no Nextclade counterpart (alignment gap)
  is logged and skipped.
- **Empty-variant VCF** — handled; `#CHROM` is read from `##contig` and the `DETECTED_BY`
  header is still added.
- flumut reports **substitutions only**, so only missense rows are ever inserted.

---

## Shared VCF I/O — `tools.py`

`read_vcf(path)` → `(DataFrame, header_lines)` and `write_vcf(df, header_lines, path)` are
the inverse of one another and round-trip a VCF **byte-for-byte**. Key properties:

- Transparent `.gz` read **and** write.
- Header scan raises `ValueError` if there is no `#CHROM` line, and tolerates very large
  headers (many `##contig`/`##INFO` lines).
- Data columns are read as **strings with NA-detection off**, so `.`, `NA`-like tokens, and
  empty fields survive unchanged and there are no `NaN` surprises. `POS` is therefore a
  string — cast with `.astype(int)` when you need integers (as `assemble_output` does before
  sorting).
