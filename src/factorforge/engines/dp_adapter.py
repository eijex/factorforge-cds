from typing import Any
from factorforge.core.interfaces import OptimizationResult, OptimizerEngine
from factorforge.registry.versioning import engine_version

class DPEngineAdapter(OptimizerEngine):
    """Deterministic Constrained Optimizer (DP) wrapped for Benchmark."""

    def __init__(self) -> None:
        pass # Evaluator is now passed in or handled by runner

    @property
    def name(self) -> str:
        return "FactorForge-DP"

    @property
    def version(self) -> str:
        return engine_version("dp")

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
        
        # Load host table (e.g. from built-in standard table)
        from factorforge.engines.profile.utils import load_golden_set
        from factorforge.engines.profile.rules.reverse_translator import ReverseTranslator
        from factorforge.engines.dp_v2 import DPV2Optimizer
        
        golden_table = load_golden_set()
        codon_weights = ReverseTranslator._build_ref_weights(golden_table)
        
        forbidden_motifs = kwargs.get("forbidden_motifs", ["GGTCTC", "CGTCTC", "GAAGAC"])
        left_flank = kwargs.get("left_flank", "")
        right_flank = kwargs.get("right_flank", "")

        optimizer_v2 = DPV2Optimizer(forbidden_motifs=forbidden_motifs)
        dp_result = optimizer_v2.optimize(
            protein_sequence=protein,
            codon_weights=codon_weights,
            target_gc_min=target_gc_min,
            target_gc_max=target_gc_max,
            left_flank=left_flank,
            right_flank=right_flank,
        )
        
        cds = dp_result["sequence"]
        
        terminal_stop_policy = kwargs.get("terminal_stop_policy", "preserve")
        if terminal_stop_policy == "append" or (terminal_stop_policy == "preserve" and sequence.endswith("*")):
            cds += "TAA"

        metrics = {
            "score": dp_result["score"],
            "cai": dp_result["cai"],
            "gc_percent": dp_result["gc_percent"],
        }

        return OptimizationResult(
            sequence=cds,
            metrics=metrics,
            metadata={
                "engine": "dp_v2",
                "version": self.version,
                "host": host,
                "inference_mode": "exact_3d_automaton_dp",
                "gc_feasible": dp_result["gc_feasible"],
                "target_intersection_exists": dp_result["gc_feasible"],
                "feasible": True,
            },
        )

    def validate(self, sequence: str) -> bool:
        return True
