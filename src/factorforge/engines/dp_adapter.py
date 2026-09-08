from typing import Any, Optional
from factorforge.core.interfaces import OptimizationResult, OptimizerEngine
from factorforge.evaluation.evaluator import SharedEvaluator
from factorforge.analysis.metrics import load_codon_usage_table
from factorforge.analysis.feasibility import analyze_feasibility

class DPEngineAdapter(OptimizerEngine):
    """Deterministic Constrained Optimizer (DP) wrapped for Benchmark."""

    def __init__(self) -> None:
        pass # Evaluator is now passed in or handled by runner

    @property
    def name(self) -> str:
        return "FactorForge-DP"

    @property
    def version(self) -> str:
        return "1.2.0"

    def optimize(
        self,
        sequence: str,
        profile: str | None = None,
        host: str = "nbenthamiana",
        **kwargs: Any,
    ) -> OptimizationResult:
        
        protein = sequence.upper().strip().rstrip("*")
        target_gc_min = kwargs.get("target_gc_min", 0.40)
        target_gc_max = kwargs.get("target_gc_max", 0.47)
        
        # DP engine internally uses percentages (0-100), not fractions (0-1)
        gc_low = target_gc_min * 100 if target_gc_min <= 1.0 else target_gc_min
        gc_high = target_gc_max * 100 if target_gc_max <= 1.0 else target_gc_max

        # Load host table (e.g. from built-in standard table)
        from factorforge.engines.profile.utils import load_golden_set
        from factorforge.engines.profile.rules.reverse_translator import ReverseTranslator
        
        golden_table = load_golden_set()
        codon_weights = ReverseTranslator._build_ref_weights(golden_table)
        
        res = analyze_feasibility(
            protein_sequence=protein,
            codon_weights=codon_weights,
            target_gc_low=gc_low,
            target_gc_high=gc_high,
            codon_reference_id=f"host_{host}"
        )
        
        target_info = res["target"]
        best_cand = target_info.get("best_candidate")
        target_intersection_exists = best_cand is not None
        
        # Fallback if unfeasible under target GC: use best without GC constraints
        if best_cand is None:
            best_cand = res.get("best_candidate_without_gc")
            
        if best_cand is None:
            cds = "ATG" * len(protein) # Ultimate fallback, shouldn't happen for valid proteins
        else:
            cds = best_cand["dna_sequence"]
        
        terminal_stop_policy = kwargs.get("terminal_stop_policy", "preserve")
        if terminal_stop_policy == "append" or (terminal_stop_policy == "preserve" and sequence.endswith("*")):
            cds += "TAA"

        metrics = {
            "score": 0.0,
        }

        return OptimizationResult(
            sequence=cds,
            metrics=metrics,
            metadata={
                "engine": "dp",
                "version": self.version,
                "host": host,
                "inference_mode": "deterministic_solver",
                "min_achievable_gc": res.get("minimum_possible_gc"),
                "max_achievable_gc": res.get("maximum_possible_gc"),
                "target_intersection_exists": target_intersection_exists,
                "feasible": target_info.get("feasible", False),
            },
        )

    def validate(self, sequence: str) -> bool:
        return True
