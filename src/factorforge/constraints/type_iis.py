"""Shared Type IIS restriction enzyme motifs and constraint registries."""

from __future__ import annotations
from collections.abc import Iterable
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
    enzymes = (
        forbidden_enzymes
        if forbidden_enzymes is not None
        else {"BsaI", "BsmBI", "BbsI", "BpiI", "SapI"}
    )
    motifs: Set[str] = set()
    for enz in enzymes:
        if enz in TYPE_IIS_MOTIFS:
            motifs.update(TYPE_IIS_MOTIFS[enz])
    return sorted(motifs)


_RESTRICTION_SITE_MOTIFS = {
    **{name.casefold(): motifs for name, motifs in TYPE_IIS_MOTIFS.items()},
    "noti": ("GCGGCCGC",),
    "noti-hf": ("GCGGCCGC",),
}


def resolve_restriction_site_motifs(enzymes: Iterable[str] | None) -> list[str]:
    """Resolve names for the full-construct policy, rejecting unknown names."""
    motifs: set[str] = set()
    for enzyme in TYPE_IIS_MOTIFS if enzymes is None else enzymes:
        sites = _RESTRICTION_SITE_MOTIFS.get(enzyme.strip().casefold())
        if sites is None:
            raise ValueError(f"Unknown restriction enzyme: {enzyme!r}")
        motifs.update(sites)
    return sorted(motifs)
