"""Machine-learning support utilities for FactorForge."""

from .metrics import (
    amino_acid_identity,
    calculate_cai,
    calculate_first_region_gc,
    calculate_gc,
    calculate_gc_windows,
    codon_usage_profile,
    count_internal_stops,
    detect_forbidden_motifs,
    detect_homopolymers,
    detect_invalid_codons,
    detect_repeats,
    translate_dna,
)
from .feasibility import analyze_feasibility

__all__ = [
    "analyze_feasibility",
    "amino_acid_identity",
    "calculate_cai",
    "calculate_first_region_gc",
    "calculate_gc",
    "calculate_gc_windows",
    "codon_usage_profile",
    "count_internal_stops",
    "detect_forbidden_motifs",
    "detect_homopolymers",
    "detect_invalid_codons",
    "detect_repeats",
    "translate_dna",
]
