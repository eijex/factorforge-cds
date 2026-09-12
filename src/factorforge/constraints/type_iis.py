"""Shared Type IIS restriction enzyme motifs and constraint registries."""

from __future__ import annotations
from typing import Dict, List, Optional, Set

TYPE_IIS_MOTIFS: Dict[str, List[str]] = {
    "BsaI": ["GGTCTC", "GAGACC"],
    "BsmBI": ["CGTCTC", "GAGACG"],
    "BbsI": ["GAAGAC", "GTCTTC"],
    "SapI": ["GCTCTTC", "GAAGAGC"],
    "BpiI": ["GAAGAC", "GTCTTC"],
}


def get_canonical_forbidden_motifs(forbidden_enzymes: Optional[Set[str]] = None) -> List[str]:
    """Extract canonical forbidden motif sequences from the registry."""
    enzymes = forbidden_enzymes if forbidden_enzymes is not None else {"BsaI", "BsmBI", "BbsI", "BpiI"}
    motifs: Set[str] = set()
    for enz in enzymes:
        if enz in TYPE_IIS_MOTIFS:
            motifs.update(TYPE_IIS_MOTIFS[enz])
    return sorted(motifs)
