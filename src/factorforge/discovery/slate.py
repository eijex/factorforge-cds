"""Discovery Slate Engine for FactorForge (Job 283A / 284B).

Orchestrates multi-contract candidate generation, hard constraint filtering,
trait vector extraction, Pareto non-dominated filtering, diversity-aware
clustering, and Top-K candidate slate assembly with research-preview sLLM support.
"""

from __future__ import annotations
from datetime import datetime, timezone
import hashlib
from typing import Dict, List, Optional, Tuple

from factorforge.discovery.filter import HardConstraintConfig, HardConstraintFilter
from factorforge.discovery.generator import MultiContractCandidateGenerator, RawCandidate
from factorforge.discovery.schemas import (
    CandidateSlate,
    SlateCandidate,
    SlateSummary,
    TargetMetadata,
    TraitVector,
)
from factorforge.discovery.traits import TraitVectorExtractor
from factorforge.engines.sllm.adapters import SLMModelAdapter
from factorforge.registry.versioning import engine_version, product_version


class DiscoverySlateEngine:
    """Engine for generating diverse, Pareto-evaluated Top-K candidate slates."""

    def __init__(
        self,
        host: str = "nbenthamiana",
        hard_filter_config: Optional[HardConstraintConfig] = None,
        enable_sllm_preview: bool = False,
        sllm_model_adapter: Optional[SLMModelAdapter] = None,
    ) -> None:
        self.host = host
        self.hard_filter_config = hard_filter_config or HardConstraintConfig()
        self.enable_sllm_preview = enable_sllm_preview
        self.sllm_model_adapter = sllm_model_adapter
        self.generator = MultiContractCandidateGenerator(
            host=self.host,
            enable_sllm_preview=self.enable_sllm_preview,
            sllm_model_adapter=self.sllm_model_adapter,
            hard_filter_config=self.hard_filter_config,
        )
        self.hard_filter = HardConstraintFilter(config=self.hard_filter_config)
        self.trait_extractor = TraitVectorExtractor(host=self.host)

    def compute_utility(
        self,
        trait: TraitVector,
        weights: Optional[Dict[str, float]] = None,
    ) -> float:
        """Computes transparent scalar utility score from trait vector."""
        w = weights or {
            "w_cai": 0.40,
            "w_gc": 0.20,
            "w_mfe": 0.15,
            "w_prior": 0.15,
            "w_syn": 0.10,
        }

        # CAI component (0.0 to 1.0)
        s_cai = trait.cai_golden_set

        # GC component: peak at 43.5%, decays outside 40-47%
        gc = trait.global_gc_percent
        if 40.0 <= gc <= 47.0:
            s_gc = 1.0
        else:
            dev = min(abs(gc - 40.0), abs(gc - 47.0))
            s_gc = max(0.0, 1.0 - (dev / 10.0))

        # 5' MFE component (0.0 to 1.0, higher is less secondary structure)
        if trait.initiation_mfe_kcal_mol is not None:
            # Typical 5' 60-nt MFE ranges -25.0 to 0.0 kcal/mol
            mfe_norm = max(0.0, min(1.0, 1.0 + (trait.initiation_mfe_kcal_mol / 25.0)))
        else:
            mfe_norm = 0.5  # Neutral fallback when MFE not calculated

        # Prior component
        s_prior = trait.codon_context_prior

        # Synthesis penalty component
        s_syn = 1.0 - trait.synthesis_penalty

        utility = (
            w.get("w_cai", 0.40) * s_cai
            + w.get("w_gc", 0.20) * s_gc
            + w.get("w_mfe", 0.15) * mfe_norm
            + w.get("w_prior", 0.15) * s_prior
            + w.get("w_syn", 0.10) * s_syn
        )
        return round(max(0.0, min(1.0, utility)), 4)

    def filter_pareto_front(
        self,
        items: List[Tuple[RawCandidate, TraitVector, float]],
    ) -> List[bool]:
        """Identifies non-dominated candidates across the multi-objective trait space."""
        n = len(items)
        is_pareto = [True] * n

        # Extract objective vectors: (CAI, -abs(GC-43.5), MFE, Prior, -SynPenalty)
        vectors = []
        for raw, trait, util in items:
            mfe_val = (
                trait.initiation_mfe_kcal_mol
                if trait.initiation_mfe_kcal_mol is not None
                else -15.0
            )
            gc_dev = -abs(trait.global_gc_percent - 43.5)
            syn_risk = -trait.synthesis_penalty
            vectors.append(
                [
                    trait.cai_golden_set,
                    gc_dev,
                    mfe_val,
                    trait.codon_context_prior,
                    syn_risk,
                ]
            )

        for i in range(n):
            for j in range(n):
                if i != j:
                    # Check if candidate j strictly dominates candidate i
                    j_dominates_i = all(vectors[j][k] >= vectors[i][k] for k in range(5)) and any(
                        vectors[j][k] > vectors[i][k] for k in range(5)
                    )
                    if j_dominates_i:
                        is_pareto[i] = False
                        break

        return is_pareto

    def calculate_pairwise_diversity(self, sequences: List[str]) -> float:
        """Calculates normalized average Hamming distance across candidate sequences."""
        if len(sequences) < 2:
            return 0.0

        total_dist = 0
        pairs_count = 0

        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                s1 = sequences[i]
                s2 = sequences[j]
                min_len = min(len(s1), len(s2))
                dist = sum(1 for k in range(min_len) if s1[k] != s2[k])
                dist += abs(len(s1) - len(s2))
                total_dist += dist / max(1, min_len)
                pairs_count += 1

        return round(total_dist / pairs_count, 4) if pairs_count > 0 else 0.0

    def generate_slate(
        self,
        target_aa: str,
        target_name: str = "Target-Protein",
        mature_protein_aa_length: Optional[int] = None,
        construct_aa_length: Optional[int] = None,
        signal_peptide_included: bool = False,
        novelty_class: str = "Class_A_InDistribution",
        top_k: int = 3,
        stop_codon: str = "TAA",
        utility_weights: Optional[Dict[str, float]] = None,
        enable_sllm_preview: Optional[bool] = None,
    ) -> CandidateSlate:
        """Executes full discovery pipeline to produce a diverse, Pareto-filtered Top-K slate."""
        clean_aa = "".join(target_aa.upper().split()).rstrip("*")
        run_id = f"FF-SLATE-{datetime.now(timezone.utc).strftime('%Y%m%d%H%M%S')}"

        # 1. Target Metadata
        mat_len = mature_protein_aa_length or len(clean_aa)
        cons_len = construct_aa_length or len(clean_aa)
        target_meta = TargetMetadata(
            target_name=target_name,
            mature_protein_aa_length=mat_len,
            construct_aa_length=cons_len,
            signal_peptide_included=signal_peptide_included,
            novelty_class=novelty_class,
            uncertainty_status="not_calibrated",
            evidence_tier="COMPUTATIONAL_HEURISTIC",
        )

        # 2. Candidate Pool Generation
        raw_pool = self.generator.generate_pool(
            protein_sequence=clean_aa,
            stop_codon=stop_codon,
            enable_sllm_override=enable_sllm_preview,
        )
        generated_pool_size = len(raw_pool)

        # 3. Hard Constraint Filtering
        feasible_items: List[Tuple[RawCandidate, TraitVector, float]] = []
        for raw in raw_pool:
            filter_res = self.hard_filter.evaluate(
                sequence_dna=raw.sequence_dna,
                expected_protein_aa=clean_aa,
                expected_stop=stop_codon,
            )
            if filter_res.is_feasible:
                trait = self.trait_extractor.extract(raw.sequence_dna)
                utility = self.compute_utility(trait, utility_weights)
                feasible_items.append((raw, trait, utility))

        feasible_pool_size = len(feasible_items)

        # Fallback if hard filter eliminated everything: take raw candidates with trait extraction
        if not feasible_items:
            for raw in raw_pool:
                trait = self.trait_extractor.extract(raw.sequence_dna)
                utility = self.compute_utility(trait, utility_weights)
                feasible_items.append((raw, trait, utility))

        # 4. Pareto Frontier Identification
        pareto_flags = self.filter_pareto_front(feasible_items)
        pareto_front_size = sum(1 for f in pareto_flags if f)
        evaluated_items: List[Tuple[RawCandidate, TraitVector, float, bool]] = [
            (item[0], item[1], item[2], is_p) for item, is_p in zip(feasible_items, pareto_flags)
        ]

        # 5. Diversity-Aware Cluster Selection
        # Group by strategy_cluster
        clusters: Dict[str, List[Tuple[RawCandidate, TraitVector, float, bool]]] = {}
        for raw, trait, util, is_par in evaluated_items:
            clusters.setdefault(raw.strategy_cluster, []).append((raw, trait, util, is_par))

        # From each cluster, pick the best representative (prefer Pareto optimal, then highest utility)
        cluster_representatives: List[Tuple[RawCandidate, TraitVector, float, bool]] = []
        for cluster_name, cand_list in clusters.items():
            sorted_cands = sorted(
                cand_list,
                key=lambda x: (1 if x[3] else 0, x[2]),  # (is_pareto, utility)
                reverse=True,
            )
            cluster_representatives.append(sorted_cands[0])

        # Sort all selected cluster representatives by utility
        cluster_representatives.sort(key=lambda x: (1 if x[3] else 0, x[2]), reverse=True)

        # Select Top-K (preferring distinct cluster representatives, then remaining high-utility candidates)
        selected_items = list(cluster_representatives[:top_k])
        if len(selected_items) < top_k and len(evaluated_items) > len(selected_items):
            chosen_seqs = {item[0].sequence_dna for item in selected_items}
            remaining = [
                item for item in evaluated_items if item[0].sequence_dna not in chosen_seqs
            ]
            remaining.sort(key=lambda x: (1 if x[3] else 0, x[2]), reverse=True)
            selected_items.extend(remaining[: top_k - len(selected_items)])
        top_k_items = selected_items

        # Assemble SlateCandidate objects
        candidates: List[SlateCandidate] = []
        for rank, (raw, trait, util, is_par) in enumerate(top_k_items, start=1):
            seq_digest = f"sha256:{hashlib.sha256(raw.sequence_dna.encode('utf-8')).hexdigest()}"
            candidates.append(
                SlateCandidate(
                    rank=rank,
                    candidate_id=raw.candidate_id,
                    generation_contract=raw.generation_contract,
                    strategy_cluster=raw.strategy_cluster,
                    utility_score=util,
                    trait_vector=trait,
                    sequence_dna=raw.sequence_dna,
                    sequence_digest=seq_digest,
                    rationale=raw.rationale,
                    is_pareto_optimal=is_par,
                )
            )

        # Compute Slate Summary & Diversity Index
        selected_sequences = [c.sequence_dna for c in candidates]
        diversity_index = self.calculate_pairwise_diversity(selected_sequences)
        unique_pool_size = len({r.sequence_dna for r in raw_pool})
        slate_complete = len(candidates) == top_k
        incomplete_reason = None if slate_complete else "insufficient_distinct_feasible_candidates"

        summary = SlateSummary(
            generated_pool_size=generated_pool_size,
            unique_pool_size=unique_pool_size,
            hard_feasible_pool_size=feasible_pool_size,
            selected_slate_size=len(candidates),
            requested_top_k=top_k,
            slate_complete=slate_complete,
            incomplete_reason=incomplete_reason,
            pareto_front_size=pareto_front_size,
            diversity_index=diversity_index,
        )

        provenance = {
            "engine": "FactorForge Discovery Slate",
            "version": product_version(),
            "engine_version": engine_version("dp_v2_1_1"),
            "host": self.host,
            "sllm_preview_enabled": (
                self.enable_sllm_preview if enable_sllm_preview is None else enable_sllm_preview
            ),
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
        }

        return CandidateSlate(
            run_id=run_id,
            target_metadata=target_meta,
            slate_summary=summary,
            candidates=candidates,
            provenance=provenance,
        )
