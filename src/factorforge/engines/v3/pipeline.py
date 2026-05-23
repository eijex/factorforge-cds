"""FactorForge v3 pipeline scaffolding."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from factorforge.core.interfaces import OptimizationResult, OptimizerEngine
from factorforge.engines.registry import EngineRegistry
from factorforge.engines.v2.rules.domesticator import Domesticator
from factorforge.engines.v2.rules.rule_engine import RuleEngine

from .explain import ExplainabilityInputs, build_fda_report
from .metrics import CodonUsageTable, compute_cai, compute_gc, load_codon_usage_table
from .tokenizer import AA_TOKENS, AATokenizer, CodonTokenizer

if TYPE_CHECKING:
    from .modeling_bart_decoder import BartDecoderSkeleton


VALID_AA = set(AA_TOKENS)


@dataclass
class V3Result:
    sequence: str
    metrics: dict[str, float]
    metadata: dict[str, Any] = field(default_factory=dict)
    report: dict[str, Any] | None = None


class V3Pipeline:
    """Minimal v3 pipeline scaffold (decoder + post-guard hook)."""

    def __init__(
        self,
        aa_tokenizer: AATokenizer | None = None,
        codon_tokenizer: CodonTokenizer | None = None,
        codon_usage: CodonUsageTable | None = None,
        decoder: "BartDecoderSkeleton | None" = None,
        model_id: str = "v3-bart-decoder-skeleton",
    ) -> None:
        self.aa_tokenizer = aa_tokenizer or AATokenizer.default()
        self.codon_tokenizer = codon_tokenizer or CodonTokenizer.default()
        self.codon_usage = codon_usage or load_codon_usage_table()
        self.decoder = decoder
        self.model_id = model_id

    def run(
        self,
        sequence: str,
        encoder_embeddings: Any | None = None,
        max_new_tokens: int = 24,
        apply_post_guard: bool = True,
        seed: int = 0,
        config: dict[str, Any] | None = None,
    ) -> V3Result:
        aa_sequence = _normalize_aa(sequence)

        if encoder_embeddings is not None:
            if self.decoder is None:
                raise ValueError("Decoder not configured. Provide a BartDecoderSkeleton instance.")
            from factorforge.engines.v3.inference.constrained_decoder import (
                constrained_greedy_decode,
                validate_candidate_or_fallback,
            )

            token_ids = constrained_greedy_decode(
                self.decoder,
                encoder_embeddings,
                aa_sequence,
                self.codon_tokenizer,
            )
            dna_sequence = self.codon_tokenizer.decode(token_ids[0].tolist())
            fallback = validate_candidate_or_fallback(aa_sequence, dna_sequence)
            dna_sequence = fallback["dna_sequence"]
        else:
            dna_sequence = self._fallback_gc_target(aa_sequence)

        post_guard: dict[str, Any] | None = None
        if apply_post_guard:
            dna_sequence, post_guard = self._apply_post_guard(dna_sequence)

        metrics = {
            "cai": compute_cai(dna_sequence, self.codon_usage),
            "gc_content": compute_gc(dna_sequence),
        }

        report_inputs = ExplainabilityInputs(
            aa_sequence=aa_sequence,
            dna_sequence=dna_sequence,
            metrics=metrics,
            model_id=self.model_id,
            tokenizer_hash=self.codon_tokenizer.mapping_hash(),
            seed=seed,
            config=config or {},
            post_guard=post_guard,
        )
        report = build_fda_report(report_inputs)

        metadata = {
            "model_id": self.model_id,
            "post_guard": post_guard or {},
        }

        return V3Result(sequence=dna_sequence, metrics=metrics, metadata=metadata, report=report)

    def _reverse_translate_best(self, aa_sequence: str) -> str:
        codons: list[str] = []
        for aa in aa_sequence:
            codon = self.codon_usage.best_codon_for_aa.get(aa)
            if codon is None:
                raise ValueError(f"No codon mapping for amino acid: {aa}")
            codons.append(codon)
        return "".join(codons)

    def _fallback_gc_target(self, aa_sequence: str) -> str:
        from factorforge.engines.v3.inference.v2_adapter import optimize_with_v2

        return optimize_with_v2(aa_sequence, options={"profile": "gc_target"})["dna_sequence"]

    def _apply_post_guard(self, dna_sequence: str) -> tuple[str, dict[str, Any]]:
        rule_engine = RuleEngine()
        domesticator = Domesticator()

        # Step 1: PolyA iterative fix (v2 pipeline과 동일한 방식)
        polya_fix: dict[str, Any] | None = None
        current_seq = dna_sequence
        has_polya = any(p in dna_sequence for p in rule_engine.POLYA_PATTERNS)
        if has_polya:
            polya_result = rule_engine.fix_polya_iterative(current_seq)
            if polya_result["success"]:
                current_seq = polya_result["modified_seq"]
            polya_fix = polya_result

        # Step 2: Full scan on fixed sequence
        scan_results = rule_engine.scan_all(current_seq)

        # Step 3: Domestication
        domestication = domesticator.domesticate(current_seq, standard="golden_gate")
        domesticated = domestication.get("domesticated_seq", current_seq)

        post_guard = {
            "scan_results": scan_results,
            "domestication": domestication,
            "polya_fix": polya_fix,
            "edited": domesticated != dna_sequence,
        }
        return domesticated, post_guard


class V3Optimizer(OptimizerEngine):
    """Optimizer wrapper for v3 pipeline."""

    name = "v3 BART Decoder"
    version = "3.0.0"

    def __init__(self, pipeline: V3Pipeline | None = None) -> None:
        self.pipeline = pipeline or V3Pipeline()

    def optimize(
        self,
        sequence: str,
        profile: str | None = None,
        **kwargs: Any,
    ) -> OptimizationResult:
        result = self.pipeline.run(sequence, **kwargs)
        return OptimizationResult(
            sequence=result.sequence,
            metrics=result.metrics,
            metadata=result.metadata,
        )

    def validate(self, sequence: str) -> bool:
        try:
            _normalize_aa(sequence)
            return True
        except ValueError:
            return False

    def get_supported_profiles(self) -> list[str]:
        return ["default"]


def _normalize_aa(sequence: str) -> str:
    seq = sequence.strip().replace(" ", "").replace("\n", "").upper()
    if not seq:
        raise ValueError("Protein sequence is empty.")
    invalid = {aa for aa in seq if aa not in VALID_AA}
    if invalid:
        raise ValueError(f"Invalid amino acids found: {''.join(sorted(invalid))}")
    return seq
