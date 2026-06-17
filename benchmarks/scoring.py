"""Shared per-sequence scorer. Reuses production functions only; reads thresholds
from the registry-resolved BenchmarkConfig (no hardcoded values)."""
from __future__ import annotations
import pandas as pd
from factorforge.analysis.metrics import calculate_cai, calculate_gc, translate_dna, load_codon_usage_table
from factorforge.utils.sequence_validator import validate_cds_output, detect_invalid_codons
from factorforge.engines.profile.rules.domesticator import Domesticator
from benchmarks.config import BenchmarkConfig

_WEIGHTS = load_codon_usage_table().codon_weights
_DOM = Domesticator()


def canonical_multi_constraint_pass(df: pd.DataFrame, gc_min: float, gc_max: float) -> pd.Series:
    """Recompute multi_constraint_pass from primitive columns.
    Used when reading historical CSVs that may have stale multi_constraint_pass values.
    Definition: biological_pass AND assembly_pass AND gc_in_target_range (scoring_contract v1.1)
    """
    return (
        df["biological_pass"].fillna(False).astype(bool)
        & df["assembly_pass"].fillna(False).astype(bool)
        & df["gc_in_target_range"].fillna(False).astype(bool)
    )


def score_cds(method: str, method_type: str, sequence_id: str, protein: str, cds: str,
              config: BenchmarkConfig, runtime_seconds: float,
              codon_weights: dict[str, float] | None = None) -> dict:
    weights = codon_weights if codon_weights is not None else _WEIGHTS
    val = validate_cds_output(protein, cds)
    translated = translate_dna(cds)
    internal_stop = translated[:-1].count("*")
    invalid = len(detect_invalid_codons(cds))
    length_ok = len(cds) % 3 == 0
    biological_pass = bool(val["passed"]) and internal_stop == 0 and invalid == 0 and length_ok
    type_iis = _DOM.scan_restriction_sites(cds, "golden_gate")
    assembly_pass = len(type_iis) == 0
    gc = round(calculate_gc(cds), 2)
    gc_in_range = config.gc_min <= gc <= config.gc_max
    return {
        "method": method, "method_type": method_type, "sequence_id": sequence_id,
        "aa_identity": round(float(val["aa_identity"]), 4),
        "internal_stop_count": internal_stop,
        "invalid_codon_count": invalid,
        "length_multiple_of_three": length_ok,
        "cai": round(calculate_cai(cds, weights), 4),
        "gc_percent": gc,
        "gc_in_target_range": gc_in_range,
        "forbidden_type_iis_site_count": len(type_iis),
        "biological_pass": biological_pass,
        "assembly_pass": assembly_pass,
        "multi_constraint_pass": biological_pass and assembly_pass and gc_in_range,
        "runtime_seconds": round(runtime_seconds, 6),
    }
