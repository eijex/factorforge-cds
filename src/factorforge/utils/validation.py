"""Structured candidate DNA validation for v3-alpha."""

from __future__ import annotations

from typing import Any

from factorforge.ml.metrics import (
    amino_acid_identity,
    calculate_first_region_gc,
    calculate_gc,
    calculate_gc_windows,
    count_internal_stops,
    detect_forbidden_motifs,
    detect_homopolymers,
    detect_invalid_codons,
    detect_repeats,
)


DEFAULT_CONFIG: dict[str, Any] = {
    "gc_window_low": 30.0,
    "gc_window_high": 70.0,
    "gc_window_size": 60,
    "gc_window_step": 30,
    "forbidden_motifs": [],
    "fail_forbidden_motifs": False,
    "homopolymer_max_run": 6,
}


def validate_candidate_sequence(
    input_protein: str,
    dna_sequence: str,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Validate a candidate CDS and return a machine-readable result."""
    cfg = {**DEFAULT_CONFIG, **(config or {})}
    seq = "".join(dna_sequence.upper().replace("U", "T").split())
    errors: list[str] = []
    warnings: list[str] = []

    if not seq:
        errors.append("empty sequence")
    if len(seq) % 3 != 0:
        errors.append("sequence length is not divisible by 3")

    invalid_codons = detect_invalid_codons(seq)
    internal_stop_count = count_internal_stops(seq)
    identity = amino_acid_identity(input_protein, seq) if seq else 0.0
    gc = calculate_gc(seq)
    windows = calculate_gc_windows(
        seq,
        window_size=int(cfg["gc_window_size"]),
        step=int(cfg["gc_window_step"]),
    )
    window_values = [float(window["gc"]) for window in windows]
    gc_window_min = min(window_values) if window_values else 0.0
    gc_window_max = max(window_values) if window_values else 0.0
    low = float(cfg["gc_window_low"])
    high = float(cfg["gc_window_high"])
    gc_window_outlier_count = sum(1 for value in window_values if value < low or value > high)

    first_region = calculate_first_region_gc(seq)
    forbidden = detect_forbidden_motifs(seq, list(cfg.get("forbidden_motifs", [])))
    homopolymers = detect_homopolymers(seq, max_run=int(cfg["homopolymer_max_run"]))
    repeats = detect_repeats(seq)

    if identity < 1.0:
        errors.append("amino acid identity is below 1.0")
    if invalid_codons:
        errors.append("invalid codons detected")
    if internal_stop_count:
        errors.append("internal stop codons detected")
    if forbidden and bool(cfg.get("fail_forbidden_motifs")):
        errors.append("forbidden motifs detected")

    if gc_window_outlier_count:
        warnings.append("local GC window outliers detected")
    if forbidden and not bool(cfg.get("fail_forbidden_motifs")):
        warnings.append("forbidden motifs detected")
    if homopolymers:
        warnings.append("homopolymers detected")
    if repeats:
        warnings.append("simple repeats detected")

    return {
        "passed": not errors,
        "amino_acid_identity": identity,
        "gc": gc,
        "gc_window_min": gc_window_min,
        "gc_window_max": gc_window_max,
        "gc_window_outlier_count": gc_window_outlier_count,
        "first_30nt_gc": first_region["first_30nt_gc"],
        "first_60nt_gc": first_region["first_60nt_gc"],
        "first_90nt_gc": first_region["first_90nt_gc"],
        "internal_stop_count": internal_stop_count,
        "invalid_codon_count": len(invalid_codons),
        "forbidden_motif_count": len(forbidden),
        "homopolymer_count": len(homopolymers),
        "repeat_count": len(repeats),
        "warnings": warnings,
        "errors": errors,
    }

