#!/usr/bin/env python3
"""Reading amino-acid changes out of a fluscan VCF's ANN= field.

Shared by `merge_flumut_mutations.py` and `mutation_watcher.py` so the two cannot
disagree about which HGVS.p spellings parse, or which forms count as substitutions.

Gene-name canonicalisation is not here; it lives in `flu_gene_names.py`.

Deliberately dependency-free -- importers must not pay for pandas or Biopython to parse
a string.
"""

from __future__ import annotations

import re
from typing import NamedTuple

# Both spellings are accepted because both are produced: bin/nextclade_to_freebayes_vcf.py
# emits p.E627K, and snpEff -- which annotates the short-read and paftools branches --
# defaults to p.Glu627Lys. Matching only one yields zero annotations, which is
# indistinguishable from a VCF with no coding changes.
HGVS_P_RE = re.compile(r"p\.([A-Za-z*])(\d+)([A-Za-z*])")
HGVS_P_RE_3LETTER = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")

AA3_TO_AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*",
}

# Forms that legitimately are not substitutions: indel/frameshift notations and snpEff's
# 'p.Glu627=' synonymous shorthand. They must not be reported as parse failures.
HGVS_P_NON_SUB_RE = re.compile(r"del|ins|dup|fs|ext|=|\?")


class Change(NamedTuple):
    """One non-synonymous amino-acid change annotated on a VCF row.

    ``feature_id`` (ANN field 7) is carried alongside ``gene`` (field 4) because a gene name
    alone cannot key a mature peptide: NCBI-style HA annotations put the parent gene there
    for HA1/HA2 entries, and only the feature id names the subunit the numbering belongs to.
    """

    gene: str
    feature_id: str
    ref_aa: str
    position: int
    alt_aa: str
    hgvs_p: str


def parse_hgvs_p(hgvs_p: str) -> tuple[str, int, str] | None:
    """(ref_aa, position, alt_aa) for an HGVS.p substitution, or None if unparseable."""
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
    """ANN entries, each a list of fields, from a VCF INFO string."""
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


def annotated_changes(info: str, unparsed: list | None = None) -> list[Change]:
    """Non-synonymous changes on a VCF row, in ANN order and de-duplicated.

    A row hitting several features yields one Change per feature. Non-empty HGVS.p values
    that parse as neither spelling are appended to ``unparsed`` so a caller can say so
    rather than reporting nothing at all.
    """
    changes: list[Change] = []
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
        ref_aa, position, alt_aa = parsed
        if ref_aa == alt_aa:  # synonymous
            continue
        key = (gene, feature_id, position, alt_aa)
        if key not in seen:
            seen.add(key)
            changes.append(Change(gene, feature_id, ref_aa, position, alt_aa, hgvs_p))
    return changes


def aa_changes_from_info(info: str, unparsed: list | None = None) -> list[tuple[str, str, int, str]]:
    """(gene, feature_id, aa_pos, alt_aa) per non-synonymous change on a VCF row."""
    return [
        (c.gene, c.feature_id, c.position, c.alt_aa)
        for c in annotated_changes(info, unparsed)
    ]
