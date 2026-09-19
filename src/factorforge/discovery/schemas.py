"""Data contracts and schemas for FactorForge Discovery Slate (Job 283A / 284B)."""

from __future__ import annotations
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional


@dataclass
class TargetMetadata:
    """Metadata describing the target protein construct and novelty status."""

    target_name: str
    mature_protein_aa_length: int
    construct_aa_length: int
    signal_peptide_included: bool = False
    novelty_class: str = "Class_A_InDistribution"
    uncertainty_status: str = "not_calibrated"
    evidence_tier: str = "COMPUTATIONAL_HEURISTIC"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class TraitVector:
    """Multi-dimensional evaluated biophysical and sequence traits."""

    cai_golden_set: float
    global_gc_percent: float
    initiation_mfe_kcal_mol: Optional[float]
    local_50bp_gc_min: float
    local_50bp_gc_max: float
    homopolymer_max_run: int
    synthesis_penalty: float = 0.0
    codon_context_prior: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class SlateCandidate:
    """A single candidate design within the Top-K slate."""

    rank: int
    candidate_id: str
    generation_contract: str
    strategy_cluster: str
    utility_score: float
    trait_vector: TraitVector
    sequence_dna: str
    sequence_digest: str
    rationale: str
    is_pareto_optimal: bool = True

    def to_dict(self) -> Dict[str, Any]:
        res = asdict(self)
        res["trait_vector"] = self.trait_vector.to_dict()
        return res


@dataclass
class SlateSummary:
    """Statistical summary of the discovery pipeline execution."""

    generated_pool_size: int = 0
    unique_pool_size: int = 0
    hard_feasible_pool_size: int = 0
    selected_slate_size: int = 0
    requested_top_k: int = 0
    slate_complete: bool = True
    incomplete_reason: Optional[str] = None
    pareto_front_size: int = 0
    diversity_index: float = 0.0
    feasible_pool_size: int = 0  # Backward-compatibility alias for hard_feasible_pool_size
    top_k_count: int = 0  # Backward-compatibility alias for selected_slate_size

    def __post_init__(self) -> None:
        if self.unique_pool_size == 0 and self.generated_pool_size > 0:
            self.unique_pool_size = self.generated_pool_size
        if self.hard_feasible_pool_size == 0 and self.feasible_pool_size > 0:
            self.hard_feasible_pool_size = self.feasible_pool_size
        if self.feasible_pool_size == 0 and self.hard_feasible_pool_size > 0:
            self.feasible_pool_size = self.hard_feasible_pool_size
        if self.selected_slate_size == 0 and self.top_k_count > 0:
            self.selected_slate_size = self.top_k_count
        if self.top_k_count == 0 and self.selected_slate_size > 0:
            self.top_k_count = self.selected_slate_size
        if self.requested_top_k == 0 and self.selected_slate_size > 0:
            self.requested_top_k = self.selected_slate_size

    def to_dict(self) -> Dict[str, Any]:
        res = asdict(self)
        if not res.get("feasible_pool_size"):
            res["feasible_pool_size"] = self.hard_feasible_pool_size
        if not res.get("top_k_count"):
            res["top_k_count"] = self.selected_slate_size
        return res


@dataclass
class CandidateSlate:
    """Full Candidate Slate payload for API response and report generation."""

    run_id: str
    target_metadata: TargetMetadata
    slate_summary: SlateSummary
    candidates: List[SlateCandidate]
    provenance: Dict[str, Any] = field(
        default_factory=lambda: {
            "schema_version": "0.1.0",
            "created_at": datetime.now(timezone.utc).isoformat(),
        }
    )

    @property
    def summary(self) -> SlateSummary:
        return self.slate_summary

    def to_dict(self) -> Dict[str, Any]:
        return {
            "$schema": "https://eijex.com/schemas/factorforge/candidate-slate-v0.1.json",
            "run_id": self.run_id,
            "target_metadata": self.target_metadata.to_dict(),
            "slate_summary": self.slate_summary.to_dict(),
            "candidates": [c.to_dict() for c in self.candidates],
            "provenance": self.provenance,
        }
