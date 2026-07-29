# Custom mutation watchlists

A watchlist is a CSV of amino-acid substitutions the pipeline scans every sample's final
VCF for. It exists for mutations *we* care about on *our* terms — the curated FluMutDB
markers are handled separately and automatically, and do not belong in a watchlist.

Lists live in `assets/watchlists/`. `--watchlist_dir` points at that directory by default;
set it to `NO_FILE` to skip the scan entirely.

## Schema

```csv
gene,mutation,note
PB2,E627K,mammalian adaptation
HA1,T195I,
NA,H275Y,oseltamivir resistance
```

| Column | Required | Meaning |
|---|---|---|
| `gene` | yes | The protein the position is numbered in. Canonicalised, so `NS-1`/`NS1` and `NS-2`/`NEP` are interchangeable. |
| `mutation` | yes | `[REF][POS][ALT]`, 1-letter (`E627K`) or 3-letter (`Glu627Lys`). |
| `note` | no | Free text, passed through to the output. |

Lines beginning with `#` are comments, so a list can record where it came from without a
second file.

Only single substitutions are supported. Insertions, deletions and frameshifts are
rejected at load time rather than silently ignored.

The file is named `*watchlist.csv`; the part before `_watchlist` becomes part of the output
filename, so several lists can be scanned in one run without colliding.

## It is `gene`, not `segment`

The column names a **protein**, because it is matched against a VCF annotation's gene name.
`HA1` is not a segment. This is not a cosmetic distinction — see *The HA trap* below, which
is what a `segment` column cost the previous list. A file still carrying a `segment` header
is rejected outright rather than read, because the values need checking too.

## What matching actually keys on

**`(gene, position, ALT)`. The `REF` letter is not part of the key.**

What makes a sample resistant is ending up with `Y` at 275, not which residue it started
from. If a dataset's reference carries `N` there and a sample reads `p.N275Y`, that sample
is a hit. Requiring `REF` to agree would miss it. This mirrors the FluMut merge, whose
`MARKER_QUERY` deliberately omits `ref_seq` for the same reason.

The residue actually observed is written to the output's `ref_aa` column, and a
disagreement with the watchlist is warned about, so the difference is visible rather than
absorbed.

### So what is `REF` for?

Sanity-checking the list. `validate_watchlist.py --reference-aa coordinate_map.aa.tsv`
compares every entry's `REF` against the reference's actual residue at that gene and
position. That is what catches a wrong numbering convention **before** it mis-scans every
sample:

```
$ validate_watchlist.py assets/watchlists/nai_watchlist.csv \
      --reference-aa results/aggregate/reference_remap/coordinate_map.aa.tsv
WARNING: NA:R292K (line 54): reference residue at 292 is 'C', not 'R'. ...
```

A mismatch warns rather than fails, because two causes produce it and nothing in the data
separates them:

1. the entry is wrong — wrong numbering convention, or a typo;
2. the entry is fine and this reference genome simply differs at that position.

Hard-failing would block a legitimate run for cause 2; silence would hide cause 1. In
practice the two look nothing alike: a systematic numbering error lights up dozens of
entries at once, a genuine reference difference shows up as one or two.

## Mixed numbering conventions

N1 and N2 number NA differently, and the offset is not constant — `H275Y` in N1 is `H274Y`
in N2, but `R292K` in N2 is `R293K` in N1. A list mixing them will silently miss some
markers and, worse, match a *different* residue that happens to sit at the quoted position.

`validate_watchlist.py` reports the two shapes this makes:

- **one position given two reference residues** — a position has one reference residue, so
  two means two coordinate systems;
- **the same substitution one position apart** — `I222L` and `I223L` is one marker written
  twice, not two markers.

It reports them and resolves none of them. Which convention is wanted depends on what this
lab reports in and where each entry came from, and a tool that guessed would silently
delete a marker somebody added on purpose. Curating to one convention is a one-time edit.

`assets/watchlists/nai_watchlist.csv` currently has seven such conflicts, flagged in a
comment at the top of the file.

## The HA trap

HA1 numbering starts *after* HA's 16-residue signal peptide. Keying an HA1 marker on a bare
`HA` CDS returns the wrong residue at **93% of positions**, silently. So:

| Dataset | What the VCF annotates | Correct watchlist key |
|---|---|---|
| `chicken2021` | `HA1`, `HA2` as two CDS features | `HA1` / `HA2` |
| `gsgd` | `HA` CDS with `HA1`/`HA2` mature peptides | `HA1` / `HA2` |
| a bare `HA` CDS | `HA` | `HA` |

The scanner resolves an annotation's gene name first and its `Feature_ID` second, which is
what makes the mature-peptide case work: an NCBI-style HA annotation emits `Gene_Name=HA`
and `Feature_ID=HA1`, and only the second names the protein the numbering belongs to.

**The fallback never runs in reverse.** A watchlist saying `HA` will not match an `HA1`
annotation, and vice versa. Reporting nothing is correct; reporting a hit would be a wrong
answer with nothing to show it was wrong.

An earlier list keyed 88 HA rows — 24.5% of the whole file — on `HA` while the production
datasets annotate `HA1`/`HA2`. None of them could ever match, and nothing said so, because
zero hits is also what a clean sample looks like. Which is why:

## The coverage report

Every scan writes a `*_watchlist_coverage.tsv` alongside the hits:

```
sample_id     segment  gene  in_vcf  watchlist_entries  hits
SAMPLE-1      HA       HA1   yes     12                 1
SAMPLE-1      HA       HA    no      88                 0
```

and warns on stderr:

```
WARNING: watchlist gene 'HA' matched no annotated gene in this VCF (88 entries unused).
         Genes present: HA1, HA2. Did you mean HA1?
```

"This gene has no mutations" and "this gene does not exist here" are different statements,
and the point of the sidecar is that they no longer print the same.

## Output

`${outdir}/${sample_id}/${workflow}/${sample_id}_[segment_]${list}_watchlist_hits.tsv`:

```
sample_id  segment  gene  position  ref_aa  alt_aa  mutation  chrom  pos  hgvs_p  note
```

A scan that finds nothing writes a header-only file, so "scanned, found nothing" is
distinguishable from "did not run".

## Watchlists are never written into the VCF

`FLUMUT=` goes into the VCF; watchlist hits do not. This is deliberate and should stay that
way.

`FLUMUT=` is a curated, always-on annotation from a versioned database — it is a property
of the variant, and the same consensus always produces the same annotation. A watchlist is
a file we edit. Merging watchlist hits into the VCF would make the same consensus produce
different VCFs depending on the contents of a CSV, which breaks every downstream assumption
about a VCF being reproducible from its inputs.

The two share the mechanics of reading `ANN=` — HGVS.p parsing, gene canonicalisation, the
`MNV_UNRESOLVED` skip — and those live in exactly one place, `bin/lib/aa_markers.py`.

## Running it outside the pipeline

The scanner reads published VCFs, so a new list can be tried against existing results
without a pipeline run:

```bash
mutation_watcher.py \
    --vcf results/SAMPLE-1/nextclade/SAMPLE-1__REF__H5N1__NA_nc_final.vcf \
    --watchlist assets/watchlists/nai_watchlist.csv \
    --sample SAMPLE-1 --segment NA \
    -o SAMPLE-1_nai_watchlist_hits.tsv
```

Validate a list before committing it:

```bash
validate_watchlist.py assets/watchlists/nai_watchlist.csv --strict
```

`--strict` turns warnings into a non-zero exit, which is what you want while curating.

## Rows that are skipped

Rows flagged `MNV_UNRESOLVED` by `bin/fix_mnv_codons.py` are not scanned. The pipeline has
already determined their HGVS.p to be wrong, so a hit there would be a false positive. The
count is reported.
