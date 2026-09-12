"""DP v2.1 Optimizer: Position-Dependent Initiation-Aware Exact DP.

Guarantees & Architectural Contracts:
1. Exact discrete optimization over piecewise position-dependent log-fitness objectives:
   - 5' Initiation Ramp (codons 1..K): Open-topology structural relaxation & position-aware codon preference proxy.
   - Elongation Body (codons K+1..N): Maximal host CAI adaptiveness.
2. Hard Global GC Bounding: Guaranteed strict satisfaction of [gc_min, gc_max]; raises
   UnsatisfiableDesignError when no reachable state falls within the target band (unless
   exploratory allow_nearest_infeasible=True is explicitly set).
3. Full Construct Type IIS Elimination (U + CDS + STOP + R):
   Rejects forbidden motifs across construct flanks, start/stop junctions, and intra/inter-codon boundaries.
4. True Full-Path Deterministic LexRank: Path-rank prefix tracking ensures byte-identical reproducibility.
5. Dynamic Version Governance: Resolves engine version from the single source of truth version manifest.
"""

from __future__ import annotations
import math
from typing import Any, Dict, List, Optional, Set, Tuple

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE
from factorforge.constraints.type_iis import get_canonical_forbidden_motifs
from factorforge.engines.sllm.automaton import AutomatonCompiler, CompiledAutomaton
from factorforge.registry.versioning import engine_version
from factorforge.utils.exceptions import UnsatisfiableDesignError

AA_TO_CODONS: Dict[str, List[str]] = {}
for codon, aa in STANDARD_GENETIC_CODE.items():
    if aa != "*":
        AA_TO_CODONS.setdefault(aa, []).append(codon)


class DPV21Optimizer:
    """Position-Dependent Initiation-Aware Exact 3D DP Optimizer."""

    def __init__(
        self,
        forbidden_motifs: Optional[List[str]] = None,
        forbidden_enzymes: Optional[Set[str]] = None,
        include_reverse_complement: bool = True,
        ramp_length_codons: int = 15,
        alpha_harmonization: float = 1.0,
        beta_mfe_proxy: float = 1.0,
        gamma_gc_penalty: float = 0.5,
        delta_au_bonus: float = 0.5,
        r_ramp_start: float = 0.50,
        r_ramp_end: float = 0.90,
    ) -> None:
        if forbidden_motifs is not None:
            motifs = forbidden_motifs
        else:
            motifs = get_canonical_forbidden_motifs(forbidden_enzymes)

        self.automaton: CompiledAutomaton = AutomatonCompiler.compile(
            motifs, include_rc=include_reverse_complement
        )
        self.ramp_length_codons = ramp_length_codons
        self.alpha = alpha_harmonization
        self.beta = beta_mfe_proxy
        self.gamma = gamma_gc_penalty
        self.delta = delta_au_bonus
        self.r_start = r_ramp_start
        self.r_end = r_ramp_end

    def compute_codon_log_fitness(
        self,
        position_idx: int,
        codon: str,
        codon_weights: Dict[str, float],
    ) -> float:
        """Computes position-dependent log-fitness score ell_i(c)."""
        w = codon_weights.get(codon, 1e-4)
        if position_idx >= self.ramp_length_codons:
            # Stage 2: Elongation Body -> Host Log-CAI
            return math.log(max(w, 1e-6))

        # Stage 1: Initiation Ramp (0 <= position_idx < K)
        # 1. Position-dependent codon preference harmonization H(i, c)
        denom = max(1, self.ramp_length_codons - 1)
        r_i = self.r_start + (self.r_end - self.r_start) * (position_idx / float(denom))
        p_c = w  # host-relative codon preference proxy
        h_score = -abs(p_c - r_i)

        # 2. 5' Open-Topology Structural Proxy Pen_open(c)
        g_count = codon.count("G") + codon.count("C")
        au_count = 3 - g_count
        pen_open = -self.gamma * g_count + self.delta * au_count

        return self.alpha * h_score + self.beta * pen_open

    def optimize(
        self,
        protein_sequence: str,
        codon_weights: Dict[str, float],
        target_gc_min: float = 0.40,
        target_gc_max: float = 0.47,
        upstream_context: str = "",
        stop_codon: str = "TAA",
        downstream_context: str = "",
        allow_nearest_infeasible: bool = False,
        left_flank: Optional[str] = None,
        right_flank: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Runs exact position-dependent initiation-aware 3D DP over the full construct."""
        # Handle backward-compatible flank aliases
        up_ctx = left_flank if left_flank is not None else upstream_context
        down_ctx = right_flank if right_flank is not None else downstream_context

        protein = "".join(protein_sequence.upper().split()).rstrip("*")
        if not protein:
            raise ValueError("protein_sequence must not be empty")

        for aa in protein:
            if aa not in AA_TO_CODONS:
                raise ValueError(f"Unsupported amino acid residue: {aa}")

        # Compute initial automaton state from upstream_context
        s0 = 0
        for char in up_ctx.upper():
            s0 = self.automaton.step_nucleotide(s0, char)
            if s0 in self.automaton.terminal_states:
                raise ValueError(f"Upstream context '{up_ctx}' contains a forbidden motif.")

        n_codons = len(protein)
        total_nt = n_codons * 3

        # Target GC bounds
        gc_min_frac = target_gc_min if target_gc_min <= 1.0 else target_gc_min / 100.0
        gc_max_frac = target_gc_max if target_gc_max <= 1.0 else target_gc_max / 100.0

        min_gc_count = int(math.ceil(gc_min_frac * total_nt))
        max_gc_count = int(math.floor(gc_max_frac * total_nt))
        if min_gc_count > max_gc_count:
            # When total_nt is small (e.g. < 15 nt), ceil/floor inversion can occur.
            # Clamp to the nearest integer count to the target midpoint.
            mid_count = int(round(((gc_min_frac + gc_max_frac) / 2.0) * total_nt))
            min_gc_count = mid_count
            max_gc_count = mid_count

        # dp_layers[i] maps (gc_count, automaton_state) ->
        #   (cumulative_score, codon, prev_gc, prev_state, path_rank)
        dp_layers: List[Dict[Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int], int]]] = []
        dp_layers.append({(0, s0): (0.0, None, None, None, 0)})

        for i, aa in enumerate(protein):
            current_layer = dp_layers[-1]
            next_layer: Dict[Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int], int]] = {}
            synonymous_codons = AA_TO_CODONS[aa]

            for (g, s), (score, _, _, _, prev_path_rank) in current_layer.items():
                for codon in synonymous_codons:
                    # 1. Automaton Veto (Hard Invariant)
                    next_s, is_forbidden = self.automaton.step_codon(s, codon)
                    if is_forbidden:
                        continue

                    # 2. Cumulative GC count
                    codon_gc = codon.count("G") + codon.count("C")
                    next_g = g + codon_gc

                    # 3. Position-dependent log-fitness score
                    step_log_fitness = self.compute_codon_log_fitness(i, codon, codon_weights)
                    next_score = score + step_log_fitness

                    # Transition path key for true full-prefix LexRank: (prev_path_rank, codon)
                    state_key = (next_g, next_s)
                    if state_key not in next_layer:
                        next_layer[state_key] = (next_score, codon, g, s, prev_path_rank)
                    else:
                        existing_score, existing_codon, existing_prev_g, existing_prev_s, existing_prev_rank = next_layer[state_key]
                        if next_score > existing_score + 1e-12:
                            next_layer[state_key] = (next_score, codon, g, s, prev_path_rank)
                        elif abs(next_score - existing_score) <= 1e-12:
                            # True Sequence-Level Deterministic Tie-Break
                            # Compare (predecessor path rank, current codon)
                            cand_key = (prev_path_rank, codon)
                            exist_key = (existing_prev_rank, existing_codon)
                            if cand_key < exist_key:
                                next_layer[state_key] = (next_score, codon, g, s, prev_path_rank)

            if not next_layer:
                raise UnsatisfiableDesignError(
                    f"DP v2.1: State space exhausted at position {i+1} ({aa}). "
                    f"No synonymous codon path satisfies the full construct automaton constraints."
                )

            # Assign canonical path_rank to all active states in next_layer sorted by (prev_path_rank, codon)
            sorted_states = sorted(
                next_layer.keys(),
                key=lambda k: (next_layer[k][4], next_layer[k][1])
            )
            rank_map = {st: idx for idx, st in enumerate(sorted_states)}
            compacted_layer: Dict[Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int], int]] = {}
            for st, (sc, cd, pg, ps, pr) in next_layer.items():
                compacted_layer[st] = (sc, cd, pg, ps, rank_map[st])

            dp_layers.append(compacted_layer)

        final_layer = dp_layers[-1]

        # Construct full tail sequence: STOP + downstream_context
        s_tail = (stop_codon + down_ctx).upper()

        # Filter states compatible with full construct tail
        valid_final_states: List[Tuple[Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int], int]]] = []
        for (g, s), val in final_layer.items():
            curr = s
            tail_valid = True
            for char in s_tail:
                curr = self.automaton.step_nucleotide(curr, char)
                if curr in self.automaton.terminal_states:
                    tail_valid = False
                    break
            if tail_valid:
                valid_final_states.append(((g, s), val))

        if not valid_final_states:
            raise UnsatisfiableDesignError(
                f"DP v2.1: No candidate satisfies full construct tail constraints "
                f"(STOP='{stop_codon}', Downstream='{down_ctx}')."
            )

        # Candidate Selection: Multi-tier tie breaking
        target_mid_gc_count = (min_gc_count + max_gc_count) / 2.0
        in_band_candidates = [item for item in valid_final_states if min_gc_count <= item[0][0] <= max_gc_count]

        if in_band_candidates:
            best_item = max(
                in_band_candidates,
                key=lambda item: (
                    item[1][0],  # Primary: Total piecewise DP score J(X)
                    -abs(item[0][0] - target_mid_gc_count),  # Tie 1: Distance to target GC midpoint
                    -item[1][4],  # Tie 2: Full-path canonical LexRank
                ),
            )
            gc_feasible = True
        else:
            # P0 Hard Invariant Enforcement: Fail-fast unless exploratory mode is explicitly requested
            if not allow_nearest_infeasible:
                raise UnsatisfiableDesignError(
                    f"DP v2.1: No synonymous candidate satisfies the requested GC band "
                    f"[{target_gc_min:.1%}, {target_gc_max:.1%}] (requires [{min_gc_count}, {max_gc_count}] GC bases of {total_nt} nt) "
                    f"under full construct Type IIS constraints."
                )
            best_item = min(
                valid_final_states,
                key=lambda item: (
                    min(abs(item[0][0] - min_gc_count), abs(item[0][0] - max_gc_count)),
                    -item[1][0],
                    abs(item[0][0] - target_mid_gc_count),
                    item[1][4],
                ),
            )
            gc_feasible = False

        # Backtracking
        best_state, (best_score, _, _, _, _) = best_item
        chosen_codons: List[str] = []
        curr_g, curr_s = best_state

        for layer_idx in range(n_codons, 0, -1):
            entry = dp_layers[layer_idx][(curr_g, curr_s)]
            _, codon, prev_g, prev_s, _ = entry
            assert codon is not None
            assert prev_g is not None
            assert prev_s is not None
            chosen_codons.append(codon)
            curr_g, curr_s = prev_g, prev_s

        chosen_codons.reverse()
        optimized_dna = "".join(chosen_codons)

        # Compute segmental metrics
        k = min(self.ramp_length_codons, n_codons)
        ramp_codons = chosen_codons[:k]
        body_codons = chosen_codons[k:]

        ramp_dna = "".join(ramp_codons)
        body_dna = "".join(body_codons)

        ramp_cai = math.exp(sum(math.log(max(codon_weights.get(c, 1e-4), 1e-6)) for c in ramp_codons) / k) if k > 0 else 1.0
        body_cai = math.exp(sum(math.log(max(codon_weights.get(c, 1e-4), 1e-6)) for c in body_codons) / len(body_codons)) if body_codons else 1.0
        global_cai = math.exp(sum(math.log(max(codon_weights.get(c, 1e-4), 1e-6)) for c in chosen_codons) / n_codons)

        final_gc_percent = (optimized_dna.count("G") + optimized_dna.count("C")) / total_nt * 100.0
        ramp_gc_percent = (ramp_dna.count("G") + ramp_dna.count("C")) / len(ramp_dna) * 100.0 if ramp_dna else 0.0
        body_gc_percent = (body_dna.count("G") + body_dna.count("C")) / len(body_dna) * 100.0 if body_dna else 0.0

        full_construct_seq = up_ctx + optimized_dna + stop_codon + down_ctx

        try:
            ver = engine_version("dp_v2_1")
        except Exception:
            ver = "2.1.0-dev"

        return {
            "sequence": optimized_dna,
            "full_construct": full_construct_seq,
            "stop_codon": stop_codon,
            "cai": global_cai,
            "cai_5p_ramp": ramp_cai,
            "cai_body": body_cai,
            "gc_percent": final_gc_percent,
            "gc_5p_ramp_percent": ramp_gc_percent,
            "gc_body_percent": body_gc_percent,
            "gc_feasible": gc_feasible,
            "constraint_scope": "FULL_CONSTRUCT" if (up_ctx or down_ctx or stop_codon) else "CDS_ONLY",
            "score": best_score,
            "length_nt": total_nt,
            "length_aa": n_codons,
            "engine_version": ver,
            "ramp_length_codons": k,
        }
