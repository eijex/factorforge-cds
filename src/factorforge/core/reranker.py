"""Multi-Factor Slate Reranker & Phenotype Diversity Selector for FactorForge Slate v2 (Job 293A).

Implements Stage 2 (Multi-Factor Scoring & NSGA-II Non-Dominated Sorting) and
Stage 3 (Real 7-D Normalized Trait Space Furthest-First Diversity Selection).
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

from factorforge.analysis.metrics import calculate_gc_windows, detect_homopolymers
from factorforge.core.host_model import HostModel
from factorforge.core.slate_generator import CandidateDesign


@dataclass
class ScoredCandidate:
    """A candidate with normalized multi-objective features, real 7-D trait vector, and Pareto ranking."""

    candidate: CandidateDesign
    cai_norm: float
    gc_fidelity: float
    mfe_proxy_norm: float
    rare_codon_guard_score: float
    composite_utility: float
    five_prime_structure_proxy_score: float
    trait_vector_7d: List[float]
    pareto_front: int = 1
    rank: int = 0
    profile_name: str = "Balanced Optimal"


class MultiFactorSlateReranker:
    """Reranks feasible candidates using HostModel calibration curves and non-dominated sorting."""

    def __init__(
        self,
        host_model: HostModel,
        target_gc: Optional[float] = None,
        weights: Optional[Dict[str, float]] = None,
        cai_baseline: float = 0.70,
        cai_ceiling: float = 0.98,
        mfe_optimal_target_proxy: float = -2.0,
    ) -> None:
        self.host_model = host_model
        self.target_gc = target_gc if target_gc is not None else host_model.default_target_gc
        self.weights = weights or {
            "cai": 0.40,
            "gc_fidelity": 0.25,
            "mfe_initiation": 0.20,
            "rare_codon_guard": 0.15,
        }
        self.cai_baseline = cai_baseline
        self.cai_ceiling = cai_ceiling
        self.mfe_optimal_target_proxy = mfe_optimal_target_proxy

    def compute_5p_structure_proxy(self, coding_sequence: str) -> float:
        """Fast surrogate proxy score of 5' mRNA initiation secondary structure (first 48 nt).
        
        NOTE: This is a heuristic surrogate (evidence_class: SURROGATE), NOT a thermodynamic RNAfold calculation.
        """
        ramp_len = min(48, len(coding_sequence))
        ramp_seq = coding_sequence[:ramp_len]
        gc_count = ramp_seq.count("G") + ramp_seq.count("C")
        gc_pct = gc_count / float(max(1, ramp_len))

        # Proxy score: -0.5 - (gc_pct * 28.0)
        proxy_score = -0.5 - (gc_pct * 28.0)
        return round(float(proxy_score), 2)

    def score_candidates(
        self, candidates: List[CandidateDesign]
    ) -> List[ScoredCandidate]:
        """Score each candidate using fixed host calibration curves and construct real 7-D trait vectors."""
        if not candidates:
            return []

        scored: List[ScoredCandidate] = []
        top_cand_seq = candidates[0].coding_sequence

        for cand in candidates:
            # 1. Fixed Calibration: CAI norm
            cai_range = max(0.01, self.cai_ceiling - self.cai_baseline)
            cai_n = max(0.0, min(1.0, (cand.cai - self.cai_baseline) / cai_range))

            # 2. GC fidelity
            gc_diff = abs(cand.gc_global - self.target_gc)
            gc_fid = max(0.0, 1.0 - min(1.0, gc_diff / 0.15))

            # 3. 5' Initiation proxy accessibility norm
            mfe_proxy = self.compute_5p_structure_proxy(cand.coding_sequence)
            mfe_diff = abs(mfe_proxy - self.mfe_optimal_target_proxy)
            mfe_n = max(0.0, 1.0 - min(1.0, mfe_diff / 15.0))

            # 4. Rare Codon Guard Score (1.0 = best, 0 rare codons)
            rare_guard = max(0.0, 1.0 - min(1.0, cand.rare_codon_count / 5.0))

            # Composite utility
            w = self.weights
            composite = (
                w.get("cai", 0.40) * cai_n
                + w.get("gc_fidelity", 0.25) * gc_fid
                + w.get("mfe_initiation", 0.20) * mfe_n
                + w.get("rare_codon_guard", 0.15) * rare_guard
            )

            # --- REAL 7-DIMENSIONAL NORMALIZED TRAIT VECTOR T(x) in [0, 1] ---
            t_cai = cai_n
            t_gc = max(0.0, min(1.0, cand.gc_global))

            # Local GC variance
            gc_wins = calculate_gc_windows(cand.coding_sequence, window_size=50, step=10) if len(cand.coding_sequence) >= 50 else []
            if gc_wins:
                gc_vals = [float(win.get("gc", win.get("gc_percent", 50.0))) / 100.0 for win in gc_wins]
                mean_gc = sum(gc_vals) / len(gc_vals)
                var_gc = sum((v - mean_gc) ** 2 for v in gc_vals) / len(gc_vals)
            else:
                var_gc = 0.0
            t_var_gc = max(0.0, min(1.0, var_gc / 0.05))

            t_mfe_proxy = max(0.0, min(1.0, (mfe_proxy + 20.0) / 20.0))
            t_rare = rare_guard

            # Homopolymers
            hp = detect_homopolymers(cand.coding_sequence, max_run=4)
            max_hp = max([len(h.get("sequence", "")) for h in hp], default=3)
            t_hp = max(0.0, min(1.0, max_hp / 6.0))

            # Normalized sequence edit distance to top candidate
            diff_nt = sum(c1 != c2 for c1, c2 in zip(cand.coding_sequence, top_cand_seq))
            t_edit = diff_nt / float(max(1, len(cand.coding_sequence)))

            trait_vector_7d = [
                round(t_cai, 4),
                round(t_gc, 4),
                round(t_var_gc, 4),
                round(t_mfe_proxy, 4),
                round(t_rare, 4),
                round(t_hp, 4),
                round(t_edit, 4),
            ]

            # Label heuristic archetype
            if mfe_proxy >= -4.0:
                profile = "Open 5' Ramp"
            elif cand.cai >= 0.95:
                profile = "High-CAI Champion"
            elif abs(cand.gc_global - self.target_gc) <= 0.01:
                profile = "GC-Envelope Robust"
            elif cand.rare_codon_count == 0:
                profile = "Zero-Rare Defended"
            else:
                profile = "Balanced Optimal"

            sc = ScoredCandidate(
                candidate=cand,
                cai_norm=round(float(cai_n), 4),
                gc_fidelity=round(float(gc_fid), 4),
                mfe_proxy_norm=round(float(mfe_n), 4),
                rare_codon_guard_score=round(float(rare_guard), 4),
                composite_utility=round(float(composite), 4),
                five_prime_structure_proxy_score=mfe_proxy,
                trait_vector_7d=trait_vector_7d,
                profile_name=profile,
            )
            scored.append(sc)

        # Run Non-Dominated Sorting
        self._non_dominated_sort(scored)
        return scored

    def non_dominated_sort(self, candidates: List[ScoredCandidate]) -> List[ScoredCandidate]:
        """Public alias to partition candidates into Pareto Fronts (NSGA-II)."""
        self._non_dominated_sort(candidates)
        return candidates

    def _non_dominated_sort(self, candidates: List[ScoredCandidate]) -> None:
        """Partition candidates into Pareto Fronts (NSGA-II)."""
        n_candidates = len(candidates)
        domination_counts = [0] * n_candidates
        dominated_sets: List[List[int]] = [[] for _ in range(n_candidates)]
        fronts: List[List[int]] = [[]]

        for p in range(n_candidates):
            p_obj = (
                candidates[p].cai_norm,
                candidates[p].gc_fidelity,
                candidates[p].mfe_proxy_norm,
                candidates[p].rare_codon_guard_score,
            )
            for q in range(n_candidates):
                if p == q:
                    continue
                q_obj = (
                    candidates[q].cai_norm,
                    candidates[q].gc_fidelity,
                    candidates[q].mfe_proxy_norm,
                    candidates[q].rare_codon_guard_score,
                )
                p_dominates = all(pv >= qv for pv, qv in zip(p_obj, q_obj)) and any(
                    pv > qv for pv, qv in zip(p_obj, q_obj)
                )
                q_dominates = all(qv >= pv for pv, qv in zip(q_obj, p_obj)) and any(
                    qv > pv for qv, pv in zip(q_obj, p_obj)
                )

                if p_dominates:
                    dominated_sets[p].append(q)
                elif q_dominates:
                    domination_counts[p] += 1

            if domination_counts[p] == 0:
                candidates[p].pareto_front = 1
                fronts[0].append(p)

        curr_front = 0
        while len(fronts[curr_front]) > 0:
            next_front: List[int] = []
            for p in fronts[curr_front]:
                for q in dominated_sets[p]:
                    domination_counts[q] -= 1
                    if domination_counts[q] == 0:
                        candidates[q].pareto_front = curr_front + 2
                        next_front.append(q)
            curr_front += 1
            fronts.append(next_front)

        candidates.sort(key=lambda sc: (sc.pareto_front, -sc.composite_utility, sc.candidate.candidate_id))
        for idx, sc in enumerate(candidates):
            sc.rank = idx + 1


class PhenotypeDiversitySelector:
    """Selects diverse slates using weighted Furthest-First clustering in 7-D trait space."""

    def __init__(
        self,
        slate_size: int = 25,
        diversity_weights: Optional[List[float]] = None,
    ) -> None:
        self.slate_size = slate_size
        # 7-Dimensional Trait Weights:
        # [CAI, Global GC, Local GC Var, 5' MFE Proxy, Rare Guard, Homopolymer, Edit Dist]
        self.diversity_weights = diversity_weights or [
            1.5, # CAI
            1.2, # Global GC
            1.0, # Local GC variance
            1.4, # 5' MFE proxy
            0.8, # Rare codon guard
            0.8, # Homopolymer
            1.2, # Edit distance
        ]

    def weighted_trait_distance(self, v1: List[float], v2: List[float]) -> float:
        """Calculate weighted Euclidean distance between two normalized 7-D trait vectors."""
        dist_sq = sum(
            w * ((a - b) ** 2)
            for w, a, b in zip(self.diversity_weights, v1, v2)
        )
        return math.sqrt(dist_sq)

    def preselect(
        self, scored_candidates: List[ScoredCandidate], target_l: int = 75
    ) -> List[ScoredCandidate]:
        """Preselect Top-L candidates for deep ValidationHub evaluation."""
        if len(scored_candidates) <= target_l:
            return list(scored_candidates)

        selected: List[ScoredCandidate] = [scored_candidates[0]]
        remaining = list(scored_candidates[1:])

        while len(selected) < target_l and remaining:
            best_cand_idx = -1
            max_min_dist = -1.0

            for idx, cand in enumerate(remaining):
                min_dist_to_sel = min(
                    self.weighted_trait_distance(cand.trait_vector_7d, sel.trait_vector_7d)
                    for sel in selected
                )
                front_penalty = (cand.pareto_front - 1) * 0.15
                combined_score = min_dist_to_sel - front_penalty

                if combined_score > max_min_dist:
                    max_min_dist = combined_score
                    best_cand_idx = idx

            if best_cand_idx >= 0:
                selected.append(remaining.pop(best_cand_idx))
            else:
                break

        return selected

    def select_final_slate(
        self, candidates: List[ScoredCandidate], slate_size: int = 25
    ) -> List[ScoredCandidate]:
        """Assemble final Top-K slate guaranteeing exact K candidates from Pareto fronts."""
        if len(candidates) <= slate_size:
            for idx, sc in enumerate(candidates):
                sc.rank = idx + 1
            return candidates

        fronts_map: Dict[int, List[ScoredCandidate]] = {}
        for sc in candidates:
            fronts_map.setdefault(sc.pareto_front, []).append(sc)

        final_slate: List[ScoredCandidate] = []

        for front_num in sorted(fronts_map.keys()):
            front_members = fronts_map[front_num]
            needed = slate_size - len(final_slate)

            if needed <= 0:
                break

            if len(front_members) <= needed:
                final_slate.extend(front_members)
            else:
                boundary_selected: List[ScoredCandidate] = []
                remaining_members = list(front_members)

                while len(boundary_selected) < needed and remaining_members:
                    best_idx = -1
                    max_dist = -1.0

                    for idx, cand in enumerate(remaining_members):
                        if not final_slate and not boundary_selected:
                            dist = cand.composite_utility
                        else:
                            all_curr = final_slate + boundary_selected
                            dist = min(
                                self.weighted_trait_distance(cand.trait_vector_7d, s.trait_vector_7d)
                                for s in all_curr
                            )
                        if dist > max_dist:
                            max_dist = dist
                            best_idx = idx

                    if best_idx >= 0:
                        boundary_selected.append(remaining_members.pop(best_idx))
                    else:
                        boundary_selected.extend(remaining_members[: (needed - len(boundary_selected))])
                        break

                final_slate.extend(boundary_selected)

        for idx, sc in enumerate(final_slate):
            sc.rank = idx + 1

        return final_slate
