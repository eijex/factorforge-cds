"""Experimental position-conditioned DP. Available through bounded direct web testing.

Adaptiveness is not measured elongation velocity. Manual structure annotations
and uncalibrated penalties do not establish folding or expression benefits.
The existing three DP coordinates are computational, not protein geometry.
"""

from __future__ import annotations
import math
from typing import Any, Dict, List, Optional, Tuple

from factorforge.engines.dp_v2 import DPV2Optimizer
from factorforge.analysis.metrics import STANDARD_GENETIC_CODE

AA_TO_CODONS: Dict[str, List[str]] = {}
for codon, aa in STANDARD_GENETIC_CODE.items():
    if aa != "*":
        AA_TO_CODONS.setdefault(aa, []).append(codon)


class StructuralPauseDPOptimizer:
    """Experimental DP over position, GC count and automaton state.

    These are computational state coordinates, not protein geometry.
    """

    def __init__(
        self, *, positions, target_adaptiveness, penalty_weight, structure_reference, **kwargs
    ):
        if not structure_reference or not isinstance(structure_reference, str):
            raise ValueError("An explicit structure reference is required")
        self.positions = frozenset(positions)
        if any(type(p) is not int or p < 1 for p in self.positions):
            raise ValueError("Positions must be one-based positive integers")
        if not math.isfinite(target_adaptiveness) or not 0 < target_adaptiveness <= 1:
            raise ValueError("Target adaptiveness must be finite and in (0, 1]")
        if not math.isfinite(penalty_weight) or penalty_weight < 0:
            raise ValueError("Penalty weight must be finite and nonnegative")
        self.target_adaptiveness = target_adaptiveness
        self.penalty_weight = penalty_weight
        self.structure_reference = structure_reference
        self.automaton = DPV2Optimizer(**kwargs).automaton

    def optimize(
        self,
        protein_sequence: str,
        codon_weights: Dict[str, float],
        target_gc_min: float = 0.40,
        target_gc_max: float = 0.47,
        left_flank: str = "",
        right_flank: str = "",
    ) -> Dict[str, Any]:
        """Optimize the documented position-conditioned objective over feasible CDS."""
        protein = "".join(protein_sequence.upper().split()).rstrip("*")
        if any(p > len(protein) for p in self.positions):
            raise ValueError("Structure position is outside the protein")
        if any(not math.isfinite(w) or not 0 < w <= 1 for w in codon_weights.values()):
            raise ValueError("Codon adaptiveness must be finite and in (0, 1]")
        if not protein:
            raise ValueError("protein_sequence must not be empty")

        for aa in protein:
            if aa not in AA_TO_CODONS:
                raise ValueError(f"Unsupported amino acid residue: {aa}")

        # Compute initial automaton state from left_flank
        s0 = 0
        for char in left_flank.upper():
            s0 = self.automaton.step_nucleotide(s0, char)
            if s0 in self.automaton.terminal_states:
                raise ValueError(f"Left flank '{left_flank}' contains a forbidden motif.")

        n_codons = len(protein)
        total_nt = n_codons * 3

        # Convert target GC percentages/fractions to integer base counts
        gc_min_frac = target_gc_min if target_gc_min <= 1.0 else target_gc_min / 100.0
        gc_max_frac = target_gc_max if target_gc_max <= 1.0 else target_gc_max / 100.0

        min_gc_count = int(math.ceil(gc_min_frac * total_nt))
        max_gc_count = int(math.floor(gc_max_frac * total_nt))

        # DP layer structure: dict (gc_count, automaton_state) -> (log_score, codon, prev_gc, prev_state)
        # dp_layers[i] contains reachable states after emitting i codons.
        dp_layers: List[
            Dict[Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int]]]
        ] = []
        dp_layers.append({(0, s0): (0.0, None, None, None)})
        step_codon = self.automaton.step_codon
        step_cache: Dict[int, Dict[str, Tuple[int, bool]]] = {}

        for i, aa in enumerate(protein):
            current_layer = dp_layers[-1]
            next_layer: Dict[
                Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int]]
            ] = {}
            in_pause_pos = i + 1 in self.positions
            codon_info = []
            for codon in AA_TO_CODONS[aa]:
                w = codon_weights.get(codon, 1e-4)
                log_w = math.log(max(w, 1e-6))
                penalty = self.penalty_weight * abs(w - self.target_adaptiveness) if in_pause_pos else 0.0
                codon_info.append((codon, codon.count("G") + codon.count("C"), log_w - penalty))

            for (g, s), (score, _, _, _) in current_layer.items():
                s_trans = step_cache.get(s)
                if s_trans is None:
                    s_trans = {}
                    step_cache[s] = s_trans
                for codon, codon_gc, step_score in codon_info:
                    # 1. Check Automaton Veto (cached per automaton state & codon)
                    transition = s_trans.get(codon)
                    if transition is None:
                        transition = step_codon(s, codon)
                        s_trans[codon] = transition
                    next_s, is_forbidden = transition
                    if is_forbidden:
                        continue  # Hard Constraint: Skip invalid transition

                    # 2. Update GC Count & Position-Conditioned Score
                    next_g = g + codon_gc
                    next_score = score + step_score

                    # 3. DP State Relaxation (Maximize Score with deterministic tie-breaking)
                    state_key = (next_g, next_s)
                    existing = next_layer.get(state_key)
                    if existing is None:
                        next_layer[state_key] = (next_score, codon, g, s)
                    else:
                        existing_score, existing_codon, _, _ = existing
                        if next_score > existing_score + 1e-12:
                            next_layer[state_key] = (next_score, codon, g, s)
                        elif abs(next_score - existing_score) <= 1e-12:
                            # Deterministic Tie-Break: Lexicographically smaller codon
                            if existing_codon is not None and codon < existing_codon:
                                next_layer[state_key] = (next_score, codon, g, s)

            if not next_layer:
                raise RuntimeError(
                    f"DP v2: State space exhausted at position {i + 1} ({aa}). "
                    f"No valid synonymous codon path satisfies the automaton constraints."
                )

            dp_layers.append(next_layer)

        final_layer = dp_layers[-1]

        # Filter states compatible with right_flank
        valid_final_states: List[
            Tuple[Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int]]]
        ] = []
        for (g, s), val in final_layer.items():
            curr = s
            right_valid = True
            for char in right_flank.upper():
                curr = self.automaton.step_nucleotide(curr, char)
                if curr in self.automaton.terminal_states:
                    right_valid = False
                    break
            if right_valid:
                valid_final_states.append(((g, s), val))

        if not valid_final_states:
            raise RuntimeError("DP v2: No candidate satisfies right_flank construct constraints.")

        # Candidate Selection with Deterministic Multi-Tier Tie-Breaking:
        # Tier 1: Primary Score (log-CAI)
        # Tier 2: Proximity to GC midpoint
        # Tier 3: Deterministic state order
        target_mid_gc_count = (min_gc_count + max_gc_count) / 2.0
        in_band_candidates = [
            item for item in valid_final_states if min_gc_count <= item[0][0] <= max_gc_count
        ]

        if in_band_candidates:
            best_item = max(
                in_band_candidates,
                key=lambda item: (
                    item[1][0],  # Primary: log-CAI score
                    -abs(item[0][0] - target_mid_gc_count),  # Tie 1: Distance to target GC midpoint
                    -item[0][0],  # Tie 2: GC count index
                    -item[0][1],  # Tie 3: Automaton state index
                ),
            )
            gc_feasible = True
        else:
            raise RuntimeError("No admissible sequence satisfies the requested GC band")

        # Backtracking
        best_state, (best_score, _, _, _) = best_item
        chosen_codons: List[str] = []
        curr_g, curr_s = best_state

        for layer_idx in range(n_codons, 0, -1):
            entry = dp_layers[layer_idx][(curr_g, curr_s)]
            _, codon, prev_g, prev_s = entry
            assert codon is not None
            assert prev_g is not None
            assert prev_s is not None
            chosen_codons.append(codon)
            curr_g, curr_s = prev_g, prev_s

        chosen_codons.reverse()
        optimized_dna = "".join(chosen_codons)

        final_gc_percent = (optimized_dna.count("G") + optimized_dna.count("C")) / total_nt * 100.0
        final_cai = math.exp(
            sum(math.log(max(codon_weights.get(c, 1e-4), 1e-6)) for c in chosen_codons) / n_codons
        )

        return {
            "experimental": True,
            "structure_reference": self.structure_reference,
            "structure_positions": sorted(self.positions),
            "target_adaptiveness": self.target_adaptiveness,
            "penalty_weight": self.penalty_weight,
            "objective": "log_adaptiveness_minus_position_penalty",
            "sequence": optimized_dna,
            "cai": final_cai,
            "gc_percent": final_gc_percent,
            "gc_feasible": gc_feasible,
            "constraint_scope": "FULL_CONSTRUCT" if (left_flank or right_flank) else "CDS_ONLY",
            "score": best_score,
            "length_nt": total_nt,
            "length_aa": n_codons,
        }
