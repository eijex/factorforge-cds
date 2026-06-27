"""Pareto feasibility analysis for synonymous codon constraints."""

from __future__ import annotations

import math
from typing import Any

from factorforge.analysis.metrics import (
    STANDARD_GENETIC_CODE,
    calculate_cai,
    calculate_first_region_gc,
    calculate_gc,
    calculate_gc_windows,
)


# DEFAULT_CAI_TARGET=0.82 aligns with industry practice (>0.8) and is achievable
# (internal benchmark, n=49, avg CAI=0.76).
# DEFAULT_GC_LOW/HIGH = native genome-composition anchor for N. benthamiana CDS
# (_analysis/025 STEP 2: 004 endogenous n=10 measured range 40-47%, cross-checked
# against nbev11_cds_hc/all and qld183_v103 derived-asset GC ~42.8-43.1% and
# external ground truth ~44%). NOT an empirically validated expression optimum —
# this is a composition anchor, not a target to maximize toward.
# Exported as named constants so tests/test_registry_production_sync.py can
# strictly compare them against the registry (single source of truth).
DEFAULT_CAI_TARGET: float = 0.82
DEFAULT_GC_LOW: float = 40.0
DEFAULT_GC_HIGH: float = 47.0


AA_TO_CODONS: dict[str, list[str]] = {}
for _codon, _aa in STANDARD_GENETIC_CODE.items():
    if _aa == "*":
        continue
    AA_TO_CODONS.setdefault(_aa, []).append(_codon)


def _normalize_gc_bound(value: float) -> float:
    """Accept either fraction 0-1 or percentage 0-100 and return percentage."""
    return value * 100.0 if 0.0 <= value <= 1.0 else value


def _gc_count(codon: str) -> int:
    return codon.count("G") + codon.count("C")


def _cai_from_log(log_sum: float, codon_count: int) -> float:
    if codon_count <= 0 or not math.isfinite(log_sum):
        return 0.0
    return math.exp(log_sum / codon_count)


def _candidate_summary(
    dna_sequence: str | None,
    codon_weights: dict[str, float],
) -> dict[str, Any] | None:
    if dna_sequence is None:
        return None
    windows = calculate_gc_windows(dna_sequence)
    window_values = [float(window["gc"]) for window in windows]
    return {
        "dna_sequence": dna_sequence,
        "cai": calculate_cai(dna_sequence, codon_weights),
        "gc": calculate_gc(dna_sequence),
        "first_region_gc": calculate_first_region_gc(dna_sequence),
        "gc_window_min": min(window_values) if window_values else 0.0,
        "gc_window_max": max(window_values) if window_values else 0.0,
    }


def _best_gc_under_gc_range(
    states: dict[int, float],
    protein_length: int,
    gc_low: float,
    gc_high: float,
) -> int | None:
    total_bases = protein_length * 3
    best_gc: int | None = None
    best_log_sum = float("-inf")
    for gc_count, log_sum in states.items():
        gc_percent = (gc_count / total_bases) * 100.0 if total_bases else 0.0
        if gc_low <= gc_percent <= gc_high and log_sum > best_log_sum:
            best_gc = gc_count
            best_log_sum = log_sum
    return best_gc


def _reconstruct_sequence(
    backrefs: list[dict[int, tuple[int, str]]],
    final_gc_count: int,
) -> str:
    codons: list[str] = []
    gc_count = final_gc_count
    for position_refs in reversed(backrefs):
        previous_gc, codon = position_refs[gc_count]
        codons.append(codon)
        gc_count = previous_gc
    codons.reverse()
    return "".join(codons)


def analyze_feasibility(
    protein_sequence: str,
    codon_weights: dict[str, float],
    target_cai: float = DEFAULT_CAI_TARGET,
    target_gc_low: float = DEFAULT_GC_LOW,
    target_gc_high: float = DEFAULT_GC_HIGH,
    gc_ranges: list[tuple[float, float]] | None = None,
) -> dict[str, Any]:
    """Compute exact CAI/GC feasibility over synonymous codon choices.

    The dynamic program keeps the best log-CAI sequence for each reachable
    global GC count. This is exact for global GC and CAI under the supplied
    codon weights.

    See module-level DEFAULT_CAI_TARGET / DEFAULT_GC_LOW / DEFAULT_GC_HIGH for
    the calibration rationale (internal benchmark, n=49).
    """
    protein = "".join(protein_sequence.upper().split()).rstrip("*")
    if not protein:
        raise ValueError("protein_sequence must not be empty")

    # Default exploration ranges: genome-grounded native anchor first, then
    # progressively wider windows. (55.0, 65.0) is retained ONLY as an explicit
    # non-native/high-GC option (matches legacy engine-output-calibrated band) —
    # callers must opt in explicitly to it, it is not the production default.
    ranges = gc_ranges or [(40.0, 47.0), (35.0, 50.0), (55.0, 65.0)]
    normalized_ranges = [
        (_normalize_gc_bound(low), _normalize_gc_bound(high)) for low, high in ranges
    ]
    target_low = _normalize_gc_bound(target_gc_low)
    target_high = _normalize_gc_bound(target_gc_high)
    if target_low > target_high:
        raise ValueError("target_gc_low must be <= target_gc_high")

    states: dict[int, float] = {0: 0.0}
    backrefs: list[dict[int, tuple[int, str]]] = []
    min_gc_count = 0
    max_gc_count = 0

    for aa in protein:
        codons = AA_TO_CODONS.get(aa)
        if not codons:
            raise ValueError(f"No synonymous codons for amino acid: {aa}")

        min_gc_count += min(_gc_count(codon) for codon in codons)
        max_gc_count += max(_gc_count(codon) for codon in codons)

        next_states: dict[int, float] = {}
        next_backrefs: dict[int, tuple[int, str]] = {}
        for current_gc, log_sum in states.items():
            for codon in codons:
                weight = codon_weights.get(codon, 0.0)
                if weight <= 0.0:
                    continue
                new_gc = current_gc + _gc_count(codon)
                new_log_sum = log_sum + math.log(weight)
                old_log_sum = next_states.get(new_gc)
                if old_log_sum is None or new_log_sum > old_log_sum:
                    next_states[new_gc] = new_log_sum
                    next_backrefs[new_gc] = (current_gc, codon)
        if not next_states:
            raise ValueError(f"No codons with positive CAI weight for amino acid: {aa}")
        states = next_states
        backrefs.append(next_backrefs)

    protein_length = len(protein)
    total_bases = protein_length * 3
    best_any_gc = max(states, key=lambda gc_count: states[gc_count])
    best_any_log_sum = states[best_any_gc]
    best_any_sequence = _reconstruct_sequence(backrefs, best_any_gc)
    target_gc_count = _best_gc_under_gc_range(states, protein_length, target_low, target_high)

    range_results: dict[str, dict[str, Any]] = {}
    for low, high in normalized_ranges:
        gc_count = _best_gc_under_gc_range(states, protein_length, low, high)
        key = f"{low:g}-{high:g}"
        if gc_count is None:
            range_results[key] = {
                "feasible": False,
                "max_cai": None,
                "best_candidate": None,
            }
        else:
            sequence = _reconstruct_sequence(backrefs, gc_count)
            range_results[key] = {
                "feasible": True,
                "max_cai": _cai_from_log(states[gc_count], protein_length),
                "best_candidate": _candidate_summary(sequence, codon_weights),
            }

    target_cai_possible = (
        target_gc_count is not None
        and _cai_from_log(states[target_gc_count], protein_length) >= target_cai
    )

    return {
        "protein_length": protein_length,
        "minimum_possible_gc": (min_gc_count / total_bases) * 100.0,
        "maximum_possible_gc": (max_gc_count / total_bases) * 100.0,
        "maximum_achievable_cai_without_gc": _cai_from_log(best_any_log_sum, protein_length),
        "best_candidate_without_gc": _candidate_summary(best_any_sequence, codon_weights),
        "ranges": range_results,
        "target": {
            "cai": target_cai,
            "gc_low": target_low,
            "gc_high": target_high,
            "feasible": target_cai_possible,
            "max_cai_under_gc": (
                _cai_from_log(states[target_gc_count], protein_length)
                if target_gc_count is not None
                else None
            ),
            "best_candidate": (
                _candidate_summary(_reconstruct_sequence(backrefs, target_gc_count), codon_weights)
                if target_gc_count is not None
                else None
            ),
        },
    }
