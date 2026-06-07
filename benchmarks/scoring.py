"""Shared per-sequence scorer. Reuses production functions only; reads thresholds
from the registry-resolved BenchmarkConfig (no hardcoded values)."""
from __future__ import annotations
from factorforge.analysis.metrics import calculate_cai, calculate_gc, translate_dna, load_codon_usage_table
from factorforge.utils.sequence_validator import validate_cds_output, detect_invalid_codons
from factorforge.engines.profile.rules.domesticator import Domesticator
from benchmarks.config import BenchmarkConfig

_WEIGHTS = load_codon_usage_table().codon_weights
_DOM = Domesticator()


def score_cds(method: str, method_type: str, sequence_id: str, protein: str, cds: str,
              config: BenchmarkConfig, runtime_seconds: float) -> dict:
    val = validate_cds_output(protein, cds)
    translated = translate_dna(cds)
    internal_stop = translated[:-1].count("*")
    invalid = len(detect_invalid_codons(cds))
    length_ok = len(cds) % 3 == 0
    biological_pass = bool(val["passed"]) and internal_stop == 0 and invalid == 0 and length_ok
    type_iis = _DOM.scan_restriction_sites(cds, "golden_gate")
    assembly_pass = len(type_iis) == 0
    gc = round(calculate_gc(cds), 2)
    return {
        "method": method, "method_type": method_type, "sequence_id": sequence_id,
        "aa_identity": round(float(val["aa_identity"]), 4),
        "internal_stop_count": internal_stop,
        "invalid_codon_count": invalid,
        "length_multiple_of_three": length_ok,
        "cai": round(calculate_cai(cds, _WEIGHTS), 4),
        "gc_percent": gc,
        "gc_in_target_range": config.gc_min <= gc <= config.gc_max,
        "forbidden_type_iis_site_count": len(type_iis),
        "biological_pass": biological_pass,
        "assembly_pass": assembly_pass,
        "multi_constraint_pass": biological_pass and assembly_pass,
        "runtime_seconds": round(runtime_seconds, 6),
    }
