"""FactorForge Slate v2 Engine Orchestrator (Job 293A).

Coordinates the 5-Stage Scientific Pipeline:
Stage 0: Diverse Candidate Generation (Heuristic Beam v1 / 512+ candidates)
Stage 1: Hard Invariant Gate (AA Invariance, 0 Forbidden Sites, Extreme GC Guard)
Stage 2: Multi-Factor Scoring & NSGA-II Pareto Ranking (Real 7-D Trait Vector)
Stage 3: Phenotype Diversity Preselection (L ≈ 75 candidates)
Stage 4: Multi-Resolution ValidationHub (Integrity + True RNAfold or NOT_AVAILABLE + Cross-Host Diagnostic)
Stage 5: Exact Top-K Slate Assembly (Successive Pareto Fronts + Boundary Furthest-First Diversity)
"""

from __future__ import annotations

from dataclasses import dataclass, field
import datetime
import logging
import time
from typing import Any, Dict, List, Optional

from factorforge.core.host_model import HostModel
from factorforge.core.reranker import (
    MultiFactorSlateReranker,
    PhenotypeDiversitySelector,
    ScoredCandidate,
)
from factorforge.core.slate_generator import (
    DEFAULT_FORBIDDEN_SITES,
    DiverseCandidateGenerator,
    HardInvariantGate,
    SUPPORTED_GENERATION_MODES,
    normalize_protein_sequence,
)
from factorforge.validation.hub import MultiResolutionValidationHub

logger = logging.getLogger(__name__)


class SlateV2Engine:
    """Orchestrates the hardened 5-stage diverse slate generation and validation pipeline."""

    def __init__(
        self,
        host: str = "nbenthamiana",
        target_gc: Optional[float] = None,
    ) -> None:
        self.host_id = host
        self.host_model = HostModel.load(host)
        self.target_gc = target_gc if target_gc is not None else self.host_model.default_target_gc

    def generate_slate(
        self,
        protein_sequence: str,
        source_cds: Optional[str] = None,
        slate_size: int = 25,
        stop_policy: str = "append_preferred",
        generation_mode: str = "heuristic_beam_v1",
        weights: Optional[Dict[str, float]] = None,
        forbidden_sites: Optional[Dict[str, Any]] = None,
        seed: int = 42,
    ) -> Dict[str, Any]:
        """Execute the end-to-end 5-stage slate generation pipeline."""
        t0 = time.perf_counter()

        # 1. Canonical Protein Normalization (Single entry point)
        canonical_protein = normalize_protein_sequence(protein_sequence)

        # 2. Strict Input Validation (Fail Closed)
        if not isinstance(slate_size, int) or slate_size < 1 or slate_size > 50:
            raise ValueError(f"slate_size must be an integer between 1 and 50. Got {slate_size}")

        if not (0.20 <= self.target_gc <= 0.80):
            raise ValueError(f"target_gc must be between 0.20 and 0.80. Got {self.target_gc}")

        if weights is not None:
            if not isinstance(weights, dict):
                raise ValueError("weights must be a dictionary")
            for k, v in weights.items():
                if not isinstance(v, (int, float)) or v < 0:
                    raise ValueError(f"Weight '{k}' must be a non-negative number. Got {v}")
            if sum(weights.values()) <= 0:
                raise ValueError("Sum of weights must be strictly positive.")

        # Canonical generation mode resolution
        canonical_mode = "heuristic_beam_v1"
        if generation_mode in SUPPORTED_GENERATION_MODES:
            canonical_mode = "heuristic_beam_v1"
        else:
            raise ValueError(
                f"Unsupported generation_mode '{generation_mode}'. "
                f"Supported modes: {sorted(SUPPORTED_GENERATION_MODES)}"
            )

        # Parse forbidden sites
        active_forbidden_sites = list(DEFAULT_FORBIDDEN_SITES)
        if forbidden_sites:
            custom_list = []
            if isinstance(forbidden_sites, dict):
                names = []
                for v in forbidden_sites.values():
                    if isinstance(v, list):
                        names.extend(v)
                    elif isinstance(v, str):
                        names.append(v)
                for name in names:
                    found = False
                    for s in DEFAULT_FORBIDDEN_SITES:
                        if s["name"].lower() == name.lower():
                            custom_list.append(s)
                            found = True
                            break
                    if not found:
                        custom_list.append({"name": name, "sequence": name, "scan_rc": True})
            elif isinstance(forbidden_sites, list):
                for item in forbidden_sites:
                    if isinstance(item, dict):
                        custom_list.append(item)
                    elif isinstance(item, str):
                        custom_list.append({"name": item, "sequence": item, "scan_rc": True})
            if custom_list:
                active_forbidden_sites = custom_list

        # --- Stage 0: Diverse Generation ---
        t_gen_start = time.perf_counter()
        generator = DiverseCandidateGenerator(
            host_model=self.host_model,
            target_gc=self.target_gc,
        )
        target_pool_size = max(512, slate_size * 20)
        raw_pool = generator.generate_pool(
            protein_sequence=canonical_protein,
            target_pool_size=target_pool_size,
            source_cds=source_cds,
            stop_policy=stop_policy,
            seed=seed,
            generation_mode=canonical_mode,
        )
        t_gen_end = time.perf_counter()

        # --- Stage 1: Hard Invariant Gate ---
        gate = HardInvariantGate(
            forbidden_sites=active_forbidden_sites,
            extreme_gc_guard=self.host_model.extreme_gc_guard,
        )
        feasible_pool = gate.filter_candidates(raw_pool)
        if not feasible_pool:
            raise RuntimeError("No candidate survived Hard Invariant Gating.")

        # --- Stage 2: Scoring & NSGA-II Pareto Sorting ---
        reranker = MultiFactorSlateReranker(
            host_model=self.host_model,
            target_gc=self.target_gc,
            weights=weights,
        )
        scored_pool = reranker.score_candidates(feasible_pool)

        # --- Stage 3: Phenotype Diversity Preselection (L ≈ 75) ---
        selector = PhenotypeDiversitySelector(slate_size=slate_size)
        target_l = min(len(scored_pool), max(slate_size * 3, 75))
        preselected_pool = selector.preselect(scored_pool, target_l=target_l)

        # --- Stage 4: MultiResolutionValidationHub on Preselected Pool ---
        t_val_start = time.perf_counter()
        val_hub = MultiResolutionValidationHub(
            primary_host_model=self.host_model,
            forbidden_sites=active_forbidden_sites,
        )
        annotated_preselected, val_summary = val_hub.validate_and_annotate_pool(
            scored_candidates=preselected_pool,
            expected_canonical_protein=canonical_protein,
        )
        t_val_end = time.perf_counter()

        # Map candidate SHA-256 to annotated dict
        annotated_map = {item["sha256"]: item for item in annotated_preselected}

        # --- Stage 5: Final Top-K Slate Assembly ---
        final_selected_scored = selector.select_final_slate(preselected_pool, slate_size=slate_size)
        final_slate: List[Dict[str, Any]] = []

        for idx, sc in enumerate(final_selected_scored):
            ann = annotated_map.get(sc.candidate.sha256)
            if ann:
                ann_copy = dict(ann)
                ann_copy["rank"] = idx + 1
                final_slate.append(ann_copy)
            else:
                # Fallback if somehow not in map
                final_slate.append({
                    "rank": idx + 1,
                    "pareto_front": sc.pareto_front,
                    "profile_name": sc.profile_name,
                    "dna_sequence": sc.candidate.dna_sequence,
                    "coding_sequence": sc.candidate.coding_sequence,
                    "terminal_stop": sc.candidate.terminal_stop,
                    "scores": {
                        "composite_utility": round(sc.composite_utility, 3),
                        "cai": round(sc.candidate.cai, 3),
                        "gc_global": round(sc.candidate.gc_global, 3),
                        "five_prime_structure_proxy_score": sc.five_prime_structure_proxy_score,
                        "rare_codon_count": sc.candidate.rare_codon_count,
                        "rare_codon_guard_score": round(sc.rare_codon_guard_score, 3),
                    },
                    "sha256": sc.candidate.sha256,
                    "trait_vector_normalized_7d": sc.trait_vector_7d,
                })

        t_total = time.perf_counter() - t0

        now_str = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        job_id = f"slate_{now_str}"

        response = {
            "status": "success",
            "job_id": job_id,
            "metadata": {
                "protein_length_aa": len(canonical_protein),
                "coding_codons": len(canonical_protein),
                "stop_codon_included": True if stop_policy != "exclude" else False,
                "host": self.host_model.host_id,
                "host_display_name": self.host_model.display_name,
                "codon_model_version": self.host_model.model_version,
                "target_gc": self.target_gc,
                "extreme_gc_guard": list(self.host_model.extreme_gc_guard),
                "total_candidates_generated": len(raw_pool),
                "unique_feasible_candidates": len(feasible_pool),
                "preselected_candidates_evaluated": len(preselected_pool),
                "slate_size": len(final_slate),
                "generation_mode": canonical_mode,
                "mfe_5p_proxy_method": "first_48_nt_gc_density_surrogate_v1",
                "mfe_5p_evidence_class": "SURROGATE",
                "timing_benchmark": {
                    "candidate_generation_ms": round((t_gen_end - t_gen_start) * 1000.0, 2),
                    "validation_hub_diagnostic_ms": round((t_val_end - t_val_start) * 1000.0, 2),
                    "total_pipeline_ms": round(t_total * 1000.0, 2),
                },
            },
            "slate": final_slate,
            "validation_summary": val_summary.to_dict(),
        }

        return response
