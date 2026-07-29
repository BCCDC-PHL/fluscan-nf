#!/usr/bin/env python3
"""Canonical influenza gene keys, shared by the coordinate map builder and the VCF merge.

The same protein is spelled differently across annotations: flumut's reference uses the
FluMutDB protein names NS-1/NS-2, Nextclade datasets use NS1/NEP, and NCBI-derived GFF3s
often use NS2. `build_reference_coordinate_map.py` canonicalises the target GFF's gene
attribute when it builds a segment's gene keys; `merge_flumut_mutations.py` must apply
the *same* function to the ANN Gene_Name it looks those keys up with.

That symmetry is the point of this module. Both sides ultimately read the same target GFF
gene/Name attribute — the map builder from the CDS feature, the ANN producer
(`nextclade_to_freebayes_vcf.py`) from its parent gene feature — so putting one
canonicalisation on both sides makes them agree by construction, for every alias here and
every alias added later. When the table lived in the map builder alone, a target dataset
naming the gene NS2 produced a map key of NEP and a VCF gene of NS2, and every NEP marker
was silently dropped.

Deliberately dependency-free: `build_reference_coordinate_map.py` does not import
`tools.py` (and so does not pay for pandas/numpy), and a three-entry lookup table should
not change that.
"""

from __future__ import annotations

# Only aliases that name the *same* protein belong here. Anything that merely looks similar
# (M vs MP, HA vs HA1) is a different question and is resolved elsewhere, against the gene
# keys a given segment actually has.
GENE_KEY_ALIASES = {
    "NS-1": "NS1",
    "NS-2": "NEP",
    "NS2": "NEP",
}


# The influenza A proteins a gene key may name, after canonicalisation. Used to reject
# typos in hand-written inputs (bin/validate_watchlist.py); it is deliberately not consulted
# when reading an annotation, where an unrecognised gene is the dataset's business and not
# an error.
#
# 'HA' belongs here and is not a synonym for HA1: datasets exist with a bare HA CDS, and
# their coordinates really are HA's. Whether a given watchlist should say HA or HA1 depends
# on the dataset, and is answered by the watcher's gene-coverage report rather than here.
KNOWN_GENE_KEYS = frozenset({
    "PB2",
    "PB1", "PB1-F2",
    "PA", "PA-X",
    "HA", "HA1", "HA2",
    "NP",
    "NA",
    "M1", "M2",
    "NS1", "NEP",
})


def canonical_gene(name: str) -> str:
    """Canonical gene key used to pair source and target CDSes.

    Without canonicalisation the gene-key intersection in build_segment_map is empty for NS
    and no amino-acid map is built for that segment; applying it on only one side is worse
    still, because the map is then built under a key the merge never looks up.
    """
    return GENE_KEY_ALIASES.get(name, name)
