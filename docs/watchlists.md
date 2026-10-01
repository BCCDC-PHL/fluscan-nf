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
| `gene` | yes | The protein the position is numbered in. Canonicalised, so `NS-1`/`NS1` and `NS-2`/`NS2`/`NEP` are interchangeable. |
| `mutation` | yes | `[REF][POS][ALT]`, 1-letter (`E627K`) or 3-letter (`Glu627Lys`). |
| `note` | no | Free text, passed through to the output. |

Lines beginning with `#` are comments, so a list can record where it came from without a
second file. Only single substitutions are supported; insertions, deletions and
frameshifts are rejected at load time rather than silently ignored.

Anything malformed — a missing column, an unparseable or synonymous mutation, an unknown
gene — is a fatal error naming the line to edit. A watchlist that failed to load quietly
would report no mutations of concern, which is exactly what a clean sample reports.

The file is named `*watchlist.csv`; the part before `_watchlist` becomes part of the output
filename, so several lists can be scanned in one run without colliding.

## It is `gene`, not `segment`

The column names a **protein**, because it is matched against a VCF annotation's gene name.
`HA1` is not a segment. See *The HA trap* below for what the difference costs.

## What matching keys on

**`(gene, position, ALT)`. The `REF` letter is not part of the key.**

What makes a sample resistant is ending up with `Y` at 275, not which residue it started
from. If a dataset's reference carries `N` there and a sample reads `p.N275Y`, that sample
is a hit. Requiring `REF` to agree would miss it. This mirrors the FluMut merge, whose
`MARKER_QUERY` deliberately omits `ref_seq` for the same reason.

The residue actually observed is what reaches the output's `ref_aa` and `mutation` columns.

## The REF residue is checked against the VCF

`REF` is there to sanity-check the list. Every watchlist entry at a position the VCF
annotates has its `REF` compared against the residue the VCF reports, and a disagreement
warns:

```
WARNING: NA H275Y (nai_watchlist.csv:19): VCF reports reference N at NA:275, not H
```

This is what catches a wrong numbering convention. It fires whether or not the entry
produced a hit, so a list keyed in N2 numbering scanned against N1 data lights up even
when nothing matches.

Two things to know about it:

- **It warns rather than fails**, because two causes produce it and nothing in the data
  separates them: the entry is wrong, or this reference genuinely differs there. A
  systematic numbering error lights up many entries at once; a real difference shows up as
  one or two. The entry is still scanned either way — matching is on position and ALT.
- **It only covers positions the VCF annotates.** A VCF carries variants, not every
  residue, so an entry at a position where the sample has no variant cannot be checked.

### Mixed numbering conventions

N1 and N2 number NA differently, and the offset is not constant — `H275Y` in N1 is `H274Y`
in N2, but `R292K` in N2 is `R293K` in N1. A list mixing them will silently miss some
markers and, worse, match a *different* residue that happens to sit at the quoted position.
The REF warning above is how a mix surfaces. Curating a list to one convention is a
one-time edit; `assets/watchlists/nai_watchlist.csv` has not had it yet, and says so in a
comment at the top of the file.

## The HA trap

HA1 numbering starts *after* HA's 16-residue signal peptide. Keying an HA1 marker on a bare
`HA` CDS returns the wrong residue at **93% of positions**, silently. So:

| Dataset | What the VCF annotates | Correct watchlist key |
|---|---|---|
| `chicken2021` | `HA1`, `HA2` as two CDS features | `HA1` / `HA2` |
| `gsgd` | `HA` CDS with `HA1`/`HA2` mature peptides | `HA1` / `HA2` |
| a bare `HA` CDS | `HA` | `HA` |

The scanner resolves an annotation's `Feature_ID` first and its gene name second, because
the feature is what the numbering belongs to: an NCBI-style HA annotation emits
`Gene_Name=HA` with `Feature_ID=HA1`, and only the second names the protein the position is
numbered in.

**There is no path from a parent gene to one of its subunits.** A watchlist saying `HA` will
not match an `HA1` annotation, and vice versa. Reporting nothing is correct; reporting a hit
would be a wrong answer with nothing to show it was wrong.

## Output

`${outdir}/${sample_id}/${workflow}/${sample_id}_[segment_]${list}_watchlist_hits.tsv`:

```
sample_id  segment  gene  position  ref_aa  alt_aa  mutation  chrom  pos  hgvs_p  note
```

One row per watchlist entry the VCF annotates, in the order the entries appear in the list.
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

The two share the mechanics of reading `ANN=` — HGVS.p parsing and gene canonicalisation —
and those live in exactly one place, `bin/lib/aa_markers.py`.

## Rows that are skipped

Rows flagged `MNV_UNRESOLVED` by `bin/fix_mnv_codons.py` are not scanned. The pipeline has
already determined their HGVS.p to be wrong, so a hit there would be a false positive. The
count is reported.

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

`--segment` is reported verbatim in the output and is optional; it labels the scan, and
nothing is inferred from it.
