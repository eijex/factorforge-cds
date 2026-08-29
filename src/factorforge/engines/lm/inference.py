"""Experimental FactorForge-LM constrained decoding scaffold.

This module deliberately exposes a small deterministic research scaffold rather
than a trained production LM. It preserves amino-acid identity through the shared
FactorForge genetic-code table and reports validation-derived status fields.
"""

from __future__ import annotations

from typing import Any, Optional

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE
from factorforge.core.interfaces import OptimizationResult, OptimizerEngine
from factorforge.engines.lm.adapter import FactorForgeLogitMasker
from factorforge.engines.lm.models.mbart_codon import TORCH_AVAILABLE
from factorforge.engines.lm.tokenizer.control_tokenizer import FactorForgeControlTokenizer
from factorforge.evaluation.evaluator import SharedEvaluator

if TORCH_AVAILABLE:
    import torch  # noqa: F401

AA_TO_CODONS: dict[str, tuple[str, ...]] = {}
for _codon, _aa in STANDARD_GENETIC_CODE.items():
    if _aa != "*":
        AA_TO_CODONS.setdefault(_aa, tuple())
        AA_TO_CODONS[_aa] = (*AA_TO_CODONS[_aa], _codon)


class ConstrainedBeamSearchEngine:
    """Experimental constrained codon selector for LM integration tests.

    The class is not a trained expression optimizer. If ``model`` is omitted, it
    uses deterministic synonymous codon selection plus the same mask helpers that
    a future LM decoder would call. Results are computational review artifacts.
    """

    def __init__(
        self,
        tokenizer: Optional[FactorForgeControlTokenizer] = None,
        model: Optional[Any] = None,
        masker: Optional[FactorForgeLogitMasker] = None,
        beam_width: int = 5,
        target_gc_min: float = 0.40,
        target_gc_max: float = 0.47,
    ) -> None:
        self.tokenizer = tokenizer or FactorForgeControlTokenizer()
        self.target_gc_min = target_gc_min
        self.target_gc_max = target_gc_max
        self.masker = masker or FactorForgeLogitMasker(
            tokenizer=self.tokenizer,
            target_gc_min=self.target_gc_min,
            target_gc_max=self.target_gc_max,
        )
        self.model = model
        self.beam_width = beam_width
        self.evaluator = SharedEvaluator(version="1.0.0")

    def optimize_cds(
        self,
        amino_acids: str,
        host: str = "nbenthamiana",
        gc_band: str = "40-47",
        type2is_clean: bool = True,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Generate a computational CDS candidate for a protein sequence.

        Unknown/unsupported residues fail closed with ``ValueError`` instead of
        silently substituting an unrelated codon.
        """
        protein = "".join(amino_acids.upper().split()).rstrip("*")
        if not protein:
            raise ValueError("amino_acids must contain at least one residue")

        unsupported = sorted({aa for aa in protein if aa not in AA_TO_CODONS})
        if unsupported:
            raise ValueError(f"Unsupported amino-acid residues: {', '.join(unsupported)}")

        # Exercise control-token encoding so invalid control-token regressions are visible.
        self.tokenizer.encode_encoder_input(
            amino_acids=protein,
            host=host,
            gc_band=gc_band,
            type2is_clean=type2is_clean,
        )

        current_cds = ""
        for aa in protein:
            synonym_codons = AA_TO_CODONS[aa]
            logits = [-float("inf")] * self.tokenizer.vocab_size
            for codon in synonym_codons:
                token_id = self.tokenizer.token_to_id[codon]
                # Simple deterministic baseline preference. A future trained model
                # may replace these scores before the same masks are applied.
                logits[token_id] = 0.0

            masked_logits = self.masker.apply_logit_masks(
                current_cds,
                logits,
                expected_next_aa=aa,
            )
            ranked = sorted(
                synonym_codons,
                key=lambda codon: masked_logits[self.tokenizer.token_to_id[codon]],
                reverse=True,
            )
            best = ranked[0]
            if masked_logits[self.tokenizer.token_to_id[best]] == -float("inf"):
                raise ValueError(f"No valid synonymous codon remained for residue {aa!r}")
            current_cds += best

        # Handle terminal stop based on policy
        terminal_stop_policy = kwargs.get("terminal_stop_policy", "preserve")
        
        has_stop = amino_acids.strip().endswith("*")
        needs_stop = False
        if terminal_stop_policy == "append":
            needs_stop = True
        elif terminal_stop_policy == "preserve":
            needs_stop = has_stop
        
        if needs_stop:
            current_cds += "TAA"
        
        # Set up forbidden Type IIS based on type2is_clean flag
        forbidden_sites = {"BsaI", "BsmBI", "BpiI"} if type2is_clean else set()
        
        # Delegate to SharedEvaluator instead of ad-hoc checking
        eval_result = self.evaluator.evaluate_candidate(
            candidate_dna=current_cds,
            expected_protein=protein,
            candidate_id="lm-candidate-01",
            target_gc_min=self.target_gc_min,
            target_gc_max=self.target_gc_max,
            forbidden_type_iis=forbidden_sites,
        )

        return {
            "generation_engine": "lm",
            "inference_mode": "deterministic_scaffold",
            "trained_model_used": False,
            "optimized_sequence": current_cds,
            "sequence_length": len(current_cds),
            "gc_percent": eval_result.metrics.gc_percent,
            "validator_passed": eval_result.passed,
            "evaluation_report": eval_result.model_dump(),
        }


class LMEngineAdapter(OptimizerEngine):
    """Adapter wrapping ConstrainedBeamSearchEngine for standard OptimizerEngine interface."""

    def __init__(self, beam_engine: Optional[ConstrainedBeamSearchEngine] = None) -> None:
        self.beam_engine = beam_engine or ConstrainedBeamSearchEngine()

    @property
    def name(self) -> str:
        return "FactorForge-LM"

    @property
    def version(self) -> str:
        return "3.6.0-scaffold"

    def optimize(
        self,
        sequence: str,
        profile: str | None = None,
        host: str = "nbenthamiana",
        **kwargs: Any,
    ) -> OptimizationResult:
        res = self.beam_engine.optimize_cds(sequence, host=host, **kwargs)
        metrics = {
            "cai": res.get("evaluation_report", {}).get("metrics", {}).get("cai", 0.0),
            "gc_percent": res.get("gc_percent", 0.0),
            "score": 0.0,
        }
        return OptimizationResult(
            sequence=res["optimized_sequence"],
            metrics=metrics,
            metadata={
                "engine": "lm", 
                "version": self.version, 
                "host": host,
                "inference_mode": res["inference_mode"],
                "validator_passed": res["validator_passed"],
                "evaluation_report": res.get("evaluation_report", {}),
            },
        )

    def validate(self, sequence: str) -> bool:
        return bool(sequence and isinstance(sequence, str))
