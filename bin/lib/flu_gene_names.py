#!/usr/bin/env python3
"""Canonical influenza gene keys.

The same protein is spelled differently across annotations: FluMutDB uses NS-1/NS-2,
Nextclade datasets use NS1/NEP, NCBI-derived GFF3s often use NS2. Every side of a lookup
must apply `canonical_gene` to the name it keys on -- the coordinate map builder to the
target GFF's gene attribute, the VCF merge and the watchlist scanner to the ANN Gene_Name --
so that the keys agree by construction. Applied to one side only, the map is built under
one key and read under another, and matches vanish silently.

Deliberately dependency-free, so importers do not pay for pandas to read a lookup table.
"""

from __future__ import annotations

# Only aliases naming the *same* protein belong here. Anything that merely looks similar
# (M vs MP, HA vs HA1) is a different question, resolved against the gene keys a given
# segment actually has.
GENE_KEY_ALIASES = {
    "NS-1": "NS1",
    "NS-2": "NEP",
    "NS2": "NEP",
}


# The proteins a gene key may name, after canonicalisation. Used to reject typos in
# hand-written inputs; deliberately not consulted when reading an annotation, where an
# unrecognised gene is the dataset's business and not an error.
#
# 'HA' belongs here and is not a synonym for HA1: datasets exist with a bare HA CDS, and
# their coordinates really are HA's.
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
    """The canonical key for a gene name."""
    return GENE_KEY_ALIASES.get(name, name)
