"""Formal v2 adapter for v3-alpha baseline, teacher, and fallback use."""

from __future__ import annotations

from typing import Any

from factorforge.engines.v2.rules.reverse_translator import OptimizationProfile, ReverseTranslator
from factorforge.engines.v2.rules.rule_engine import RuleEngine
from factorforge.engines.v2.scoring import calculate_composite_score
from factorforge.engines.v3.metrics import load_codon_usage_table
from factorforge.ml.metrics import calculate_cai, calculate_gc
from factorforge.utils.validation import validate_candidate_sequence


def optimize_with_v2(
    protein_sequence: str,
    options: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Optimize a protein sequence with v2 semantics and return v3-alpha metadata."""
    opts = options or {}
    profile_name = str(opts.get("profile", "high_cai")).lower()
    scan_mode = str(opts.get("scan_mode", "fast"))
    try:
        profile = OptimizationProfile(profile_name)
    except ValueError as exc:
        supported = ", ".join(item.value for item in OptimizationProfile)
        raise ValueError(f"Unknown v2 profile: {profile_name}. Supported profiles: {supported}") from exc

    protein = "".join(protein_sequence.upper().split()).rstrip("*")
    if not protein:
        raise ValueError("protein_sequence must not be empty")

    translator = ReverseTranslator()
    candidates = translator.generate_candidates(protein, profile=profile, n=1)
    dna_sequence = candidates[0]["sequence"]
    table = load_codon_usage_table()
    metrics = {
        "cai": calculate_cai(dna_sequence, table.codon_weights),
        "gc": calculate_gc(dna_sequence),
        "gc_content": calculate_gc(dna_sequence),
        "score": calculate_composite_score(
            cai=candidates[0]["cai"],
            gc=candidates[0]["gc"],
            sequence=dna_sequence,
            profile=profile.value,
        ),
    }

    rule_engine = RuleEngine()
    scan_results = rule_engine.scan_all(dna_sequence, mode=scan_mode)
    validator = validate_candidate_sequence(protein, dna_sequence)
    warnings = list(validator["warnings"])
    errors = list(validator["errors"])
    violation_count = sum(len(value) for value in scan_results.values())
    if violation_count:
        warnings.append(f"v2 rule scan reported {violation_count} violation(s)")

    return {
        "engine": "v2",
        "protein_sequence": protein,
        "dna_sequence": dna_sequence,
        "metrics": metrics,
        "validator": validator,
        "warnings": warnings,
        "errors": errors,
        "metadata": {
            "profile": profile.value,
            "scan_mode": scan_mode,
            "scan_results": scan_results,
        },
    }

