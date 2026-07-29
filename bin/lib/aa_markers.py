#!/usr/bin/env python3
"""Amino-acid marker mechanics shared by every consumer of a fluscan VCF's ANN= field.

Two tools ask the same question of a VCF — "which amino-acid changes does this row
annotate, and where?" — from different curated sources:

  `merge_flumut_mutations.py`  matches them against FluMutDB and writes FLUMUT= back.
  `mutation_watcher.py`        matches them against a custom watchlist CSV and writes a
                               report; it never touches the VCF.

The sources differ; the mechanics must not. This module holds those mechanics once: ANN
field extraction, both HGVS.p substitution spellings, and the INFO-flag test used to skip
rows the pipeline has already flagged as wrong. `bin/lib/flu_gene_names.py` records what a
second copy of shared logic costs — a target dataset naming its gene NS2 built the map
under one key and looked it up under another, and every NEP marker vanished silently. The
same failure is available here for any tool that reimplements HGVS.p parsing and quietly
disagrees about which spellings, or which non-substitution forms, count.

Deliberately dependency-free, for the same reason as `flu_gene_names.py`: importers must
not be forced to pay for pandas/numpy or Biopython to parse a string.

Gene-name canonicalisation is *not* here. It lives in `flu_gene_names.py` and callers
import it from there, so there stays exactly one alias table.
"""

from __future__ import annotations

import re

# Two HGVS.p substitution spellings are accepted. bin/nextclade_to_freebayes_vcf.py emits the
# 1-letter form (p.E627K); snpEff -- which annotates the short-read and paftools branches --
# defaults to the 3-letter form (p.Glu627Lys). Matching only one of them fails silently: zero
# annotations, indistinguishable from a VCF with no coding changes.
HGVS_P_RE = re.compile(r"p\.([A-Za-z*])(\d+)([A-Za-z*])")
HGVS_P_RE_3LETTER = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")

# Local rather than Bio.SeqUtils.seq1: this module is deliberately dependency-free (see the
# module docstring), and a 23-entry table does not justify a Biopython dependency.
AA3_TO_AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*",
}


# Forms that legitimately are not substitutions: indel/frameshift notations and snpEff's
# 'p.Glu627=' synonymous shorthand. FluMutDB 6.6 holds only SNP-type mutations, so skipping
# these loses nothing curated -- but they must not be reported as parse failures either.
HGVS_P_NON_SUB_RE = re.compile(r"del|ins|dup|fs|ext|=|\?")


def parse_hgvs_p(hgvs_p: str) -> tuple[str, int, str] | None:
    """(ref_aa, position, alt_aa) for an HGVS.p substitution, or None if unparseable.

    Non-substitution forms (del/ins/fs) return None. FluMutDB 6.6 holds only SNP-type
    mutations, so nothing curated is lost by ignoring them.
    """
    hit = HGVS_P_RE.search(hgvs_p)
    if hit:
        return hit.group(1), int(hit.group(2)), hit.group(3)
    hit = HGVS_P_RE_3LETTER.search(hgvs_p)
    if hit:
        ref_aa = AA3_TO_AA1.get(hit.group(1))
        alt_aa = AA3_TO_AA1.get(hit.group(3))
        if ref_aa and alt_aa:
            return ref_aa, int(hit.group(2)), alt_aa
    return None


def ann_entries(info: str) -> list[list[str]]:
    """Extract ANN entries (each a list of fields) from a VCF INFO string."""
    for token in info.split(";"):
        if token.startswith("ANN="):
            return [entry.split("|") for entry in token[4:].split(",") if entry]
    return []


def has_info_flag(info: str, key: str) -> bool:
    """Whether a VCF INFO string carries `key`, as a flag or a key=value pair.

    Token-wise rather than a substring test, so a key is never matched inside another key's
    name or value.
    """
    return any(token == key or token.startswith(f"{key}=") for token in info.split(";"))


def aa_changes_from_info(info: str, unparsed: list | None = None) -> list[tuple[str, str, int, str]]:
    """Non-synonymous (gene, feature_id, aa_pos, alt_aa) changes annotated on a VCF row.

    One item per ANN entry that carries a parseable HGVS.p substitution, in ANN order and
    de-duplicated, so a row hitting several features yields one candidate per feature.

    ``feature_id`` (ANN field 7) is carried alongside ``gene`` (field 4) because a gene name
    alone cannot key a mature peptide: NCBI-style HA annotations put the *parent* gene there
    for the HA1/HA2 entries, and only the feature id names the subunit.

    Non-empty HGVS.p values that parse as neither substitution spelling are appended to
    ``unparsed`` so the caller can say so rather than reporting nothing at all.
    """
    changes: list[tuple[str, str, int, str]] = []
    seen: set = set()
    for fields in ann_entries(info):
        if len(fields) < 11:
            continue
        gene, feature_id, hgvs_p = fields[3], fields[6], fields[10]
        parsed = parse_hgvs_p(hgvs_p)
        if parsed is None:
            if hgvs_p and unparsed is not None and not HGVS_P_NON_SUB_RE.search(hgvs_p):
                unparsed.append(hgvs_p)
            continue
        ref_aa, pos, alt_aa = parsed
        if ref_aa == alt_aa:  # synonymous
            continue
        key = (gene, feature_id, pos, alt_aa)
        if key not in seen:
            seen.add(key)
            changes.append(key)
    return changes
