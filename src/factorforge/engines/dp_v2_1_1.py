"""DP v2.1.1 Optimizer: Local Composition Guard & Initiation-Aware Exact DP.

Guarantees & Architectural Contracts:
1. Exact Hard Initiation GC Band Pruning (i = K = 15 codons / 45 nt):
   Strictly enforces ceil(45 * gc_min) <= h_15 <= floor(45 * gc_max) in the DP state space
   (for default 20%-30% band, exactly 9 <= h_15 <= 13 GC bases) clamped to the protein's
   achievable synonymous GC envelope [GC_min_syn, GC_max_syn] with ZERO additional state dimensions.
2. Homopolymer Aho-Corasick Automaton Rejection:
   Rejects >= 6-mer homopolymers (AAAAAA, TTTTTT, GGGGGG, CCCCCC) via exact automaton compilation.
3. Configured Multi-Enzyme Type IIS Invariants:
   Full construct elimination across upstream flank, CDS, stop codon, and downstream flank.
4. Deterministic LexRank Sequence Tie-Breaking:
   Guarantees byte-identical, reproducible path selection.
5. Local Composition & Biophysical Metrics:
   Emits 50-bp sliding GC window statistics, context-complete 60-nt 5' MFE
   (-15 nt upstream + 45 nt CDS), and GC constraint active-bound status.
"""

from __future__ import annotations
import math
from typing import Any, Dict, List, Optional, Set, Tuple

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE, detect_homopolymers
from factorforge.constraints.type_iis import get_canonical_forbidden_motifs
from factorforge.engines.profile.scoring import compute_5p_mfe_evidence
from factorforge.engines.sllm.automaton import AutomatonCompiler, CompiledAutomaton
from factorforge.registry.versioning import engine_version
from factorforge.utils.exceptions import UnsatisfiableDesignError

AA_TO_CODONS: Dict[str, List[str]] = {}
for codon, aa in STANDARD_GENETIC_CODE.items():
    if aa != "*":
        AA_TO_CODONS.setdefault(aa, []).append(codon)


class DPV211Optimizer:
    """Generation 2.1.1 Initiation-Aware Exact DP Optimizer with Local Composition Guard."""

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
        homopolymer_max_run: int = 5,
        initiation_gc_min: float = 0.20,
        initiation_gc_max: float = 0.30,
        nominal_ramp_gc: float = 0.25,
    ) -> None:
        if forbidden_motifs is not None:
            motifs = list(forbidden_motifs)
        else:
            motifs = get_canonical_forbidden_motifs(forbidden_enzymes)

        # Homopolymer Policy: Add forbidden runs >= (homopolymer_max_run + 1)
        self.homopolymer_max_run = homopolymer_max_run
        if self.homopolymer_max_run > 0:
            for base in ("A", "T", "G", "C"):
                poly = base * (self.homopolymer_max_run + 1)
                if poly not in motifs:
                    motifs.append(poly)

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
        self.initiation_gc_min = initiation_gc_min
        self.initiation_gc_max = initiation_gc_max
        self.nominal_ramp_gc = nominal_ramp_gc

        if self.ramp_length_codons < 1:
            raise ValueError("ramp_length_codons must be at least 1")
        if self.homopolymer_max_run < 1:
            raise ValueError("homopolymer_max_run must be at least 1")
        if not 0.0 <= self.initiation_gc_min <= self.initiation_gc_max <= 1.0:
            raise ValueError("initiation GC bounds must satisfy 0 <= min <= max <= 1")

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
        """Runs exact position-dependent initiation-aware 3D DP with local composition guard."""
        up_ctx = left_flank if left_flank is not None else upstream_context
        down_ctx = right_flank if right_flank is not None else downstream_context

        protein = "".join(protein_sequence.upper().split()).rstrip("*")
        if not protein:
            raise ValueError("protein_sequence must not be empty")

        for aa in protein:
            if aa not in AA_TO_CODONS:
                raise ValueError(f"Unsupported amino acid residue: {aa}")

        # Initial automaton state from upstream_context
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
            mid_count = int(round(((gc_min_frac + gc_max_frac) / 2.0) * total_nt))
            min_gc_count = mid_count
            max_gc_count = mid_count

        # 5' Initiation Ramp GC Constraints at codon i = K = ramp_length_codons
        ramp_k = min(self.ramp_length_codons, n_codons)
        ramp_nt = ramp_k * 3
        apply_initiation_gc_guard = n_codons >= self.ramp_length_codons
        nom_init_min = int(math.ceil(self.initiation_gc_min * ramp_nt))
        nom_init_max = int(math.floor(self.initiation_gc_max * ramp_nt))
        if nom_init_min > nom_init_max:
            mid_init = int(round(self.nominal_ramp_gc * ramp_nt))
            nom_init_min = mid_init
            nom_init_max = mid_init

        # Reachable synonymous GC range for the initiation ramp
        ramp_min_possible_gc = sum(
            min(c.count("G") + c.count("C") for c in AA_TO_CODONS[aa]) for aa in protein[:ramp_k]
        )
        ramp_max_possible_gc = sum(
            max(c.count("G") + c.count("C") for c in AA_TO_CODONS[aa]) for aa in protein[:ramp_k]
        )

        # Intersect requested initiation band with the synonymous GC envelope.
        # This clamping handles intrinsically low/high-GC prefixes without adding
        # another DP state dimension.
        init_min_gc_count = max(nom_init_min, ramp_min_possible_gc)
        init_max_gc_count = min(nom_init_max, ramp_max_possible_gc)
        initiation_gc_status = (
            "REQUESTED_BAND" if apply_initiation_gc_guard else "NOT_APPLIED_SHORT_SEQUENCE"
        )
        if apply_initiation_gc_guard and (
            init_min_gc_count != nom_init_min or init_max_gc_count != nom_init_max
        ):
            initiation_gc_status = "SYNONYMOUS_ENVELOPE_CLAMPED"
        if init_min_gc_count > init_max_gc_count:
            if apply_initiation_gc_guard:
                initiation_gc_status = "SYNONYMOUS_ENVELOPE_CLAMPED"
            if ramp_min_possible_gc > nom_init_max:
                init_min_gc_count = ramp_min_possible_gc
                init_max_gc_count = min(ramp_min_possible_gc + 2, ramp_max_possible_gc)
            else:
                init_max_gc_count = ramp_max_possible_gc
                init_min_gc_count = max(ramp_max_possible_gc - 2, ramp_min_possible_gc)

        # dp_layers[i] maps (gc_count, automaton_state) ->
        #   (cumulative_score, codon, prev_gc, prev_state, path_rank)
        dp_layers: List[
            Dict[Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int], int]]
        ] = []
        dp_layers.append({(0, s0): (0.0, None, None, None, 0)})

        for i, aa in enumerate(protein):
            current_layer = dp_layers[-1]
            next_layer: Dict[
                Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int], int]
            ] = {}
            synonymous_codons = AA_TO_CODONS[aa]

            is_ramp_boundary = i == ramp_k - 1

            for (g, s), (score, _, _, _, prev_path_rank) in current_layer.items():
                for codon in synonymous_codons:
                    # 1. Automaton Veto (Hard Invariant: Type IIS + Homopolymers)
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
                        (
                            existing_score,
                            existing_codon,
                            existing_prev_g,
                            existing_prev_s,
                            existing_prev_rank,
                        ) = next_layer[state_key]
                        if next_score > existing_score + 1e-12:
                            next_layer[state_key] = (next_score, codon, g, s, prev_path_rank)
                        elif abs(next_score - existing_score) <= 1e-12:
                            # Deterministic Tie-Break: (predecessor path rank, current codon)
                            cand_key = (prev_path_rank, codon)
                            exist_key = (existing_prev_rank, existing_codon)
                            if cand_key < exist_key:
                                next_layer[state_key] = (next_score, codon, g, s, prev_path_rank)

            if not next_layer:
                raise UnsatisfiableDesignError(
                    f"DP v2.1.1: State space exhausted at position {i + 1} ({aa}). "
                    f"No synonymous codon path satisfies the full construct automaton constraints."
                )

            # Exact 5' Initiation GC Band Hard Pruning at i = ramp_k - 1
            if is_ramp_boundary and apply_initiation_gc_guard:
                filtered_layer = {
                    k_st: v
                    for k_st, v in next_layer.items()
                    if init_min_gc_count <= k_st[0] <= init_max_gc_count
                }
                if filtered_layer:
                    next_layer = filtered_layer
                else:
                    # If automaton vetoes make the active band unreachable, retain
                    # only the nearest reachable GC layer and report the relaxation.
                    # The selected active layer is still enforced exactly.
                    initiation_gc_status = "AUTOMATON_REACHABILITY_RELAXED"
                    min_dist_to_band = min(
                        0
                        if init_min_gc_count <= k_st[0] <= init_max_gc_count
                        else min(abs(k_st[0] - init_min_gc_count), abs(k_st[0] - init_max_gc_count))
                        for k_st in next_layer.keys()
                    )
                    next_layer = {
                        k_st: v
                        for k_st, v in next_layer.items()
                        if (
                            0
                            if init_min_gc_count <= k_st[0] <= init_max_gc_count
                            else min(
                                abs(k_st[0] - init_min_gc_count), abs(k_st[0] - init_max_gc_count)
                            )
                        )
                        == min_dist_to_band
                    }
                    reachable_counts = [state[0] for state in next_layer]
                    init_min_gc_count = min(reachable_counts)
                    init_max_gc_count = max(reachable_counts)

            # Assign canonical path_rank to all active states in next_layer sorted by (prev_path_rank, codon)
            sorted_states = sorted(
                next_layer.keys(), key=lambda k_st: (next_layer[k_st][4], next_layer[k_st][1])
            )
            rank_map = {st: idx for idx, st in enumerate(sorted_states)}
            compacted_layer: Dict[
                Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int], int]
            ] = {}
            for st, (sc, cd, pg, ps, pr) in next_layer.items():
                compacted_layer[st] = (sc, cd, pg, ps, rank_map[st])

            dp_layers.append(compacted_layer)

        final_layer = dp_layers[-1]

        # Construct full tail sequence: STOP + downstream_context
        s_tail = (stop_codon + down_ctx).upper()

        # Filter states compatible with full construct tail
        valid_final_states: List[
            Tuple[Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int], int]]
        ] = []
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
                f"DP v2.1.1: No candidate satisfies full construct tail constraints "
                f"(STOP='{stop_codon}', Downstream='{down_ctx}')."
            )

        # Candidate Selection: Multi-tier tie breaking
        target_mid_gc_count = (min_gc_count + max_gc_count) / 2.0
        in_band_candidates = [
            item for item in valid_final_states if min_gc_count <= item[0][0] <= max_gc_count
        ]

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
            if not allow_nearest_infeasible:
                raise UnsatisfiableDesignError(
                    f"DP v2.1.1: No synonymous candidate satisfies the requested GC band "
                    f"[{gc_min_frac:.1%}, {gc_max_frac:.1%}] (requires [{min_gc_count}, {max_gc_count}] GC bases of {total_nt} nt) "
                    f"under full construct Type IIS and initiation constraints."
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

        # Compute segmental and local metrics
        ramp_codons = chosen_codons[:ramp_k]
        body_codons = chosen_codons[ramp_k:]

        ramp_dna = "".join(ramp_codons)
        body_dna = "".join(body_codons)

        ramp_cai = (
            math.exp(
                sum(math.log(max(codon_weights.get(c, 1e-4), 1e-6)) for c in ramp_codons) / ramp_k
            )
            if ramp_k > 0
            else 1.0
        )
        body_cai = (
            math.exp(
                sum(math.log(max(codon_weights.get(c, 1e-4), 1e-6)) for c in body_codons)
                / len(body_codons)
            )
            if body_codons
            else 1.0
        )
        global_cai = math.exp(
            sum(math.log(max(codon_weights.get(c, 1e-4), 1e-6)) for c in chosen_codons) / n_codons
        )

        final_gc_percent = (optimized_dna.count("G") + optimized_dna.count("C")) / total_nt * 100.0
        ramp_gc_percent = (
            (ramp_dna.count("G") + ramp_dna.count("C")) / len(ramp_dna) * 100.0 if ramp_dna else 0.0
        )
        body_gc_percent = (
            (body_dna.count("G") + body_dna.count("C")) / len(body_dna) * 100.0 if body_dna else 0.0
        )

        gc_5p_30nt_percent = (
            (optimized_dna[:30].count("G") + optimized_dna[:30].count("C"))
            / min(30, len(optimized_dna))
            * 100.0
            if optimized_dna
            else 0.0
        )
        gc_5p_45nt_percent = (
            (optimized_dna[:45].count("G") + optimized_dna[:45].count("C"))
            / min(45, len(optimized_dna))
            * 100.0
            if optimized_dna
            else 0.0
        )

        # GC Constraint Status
        target_min_p = target_gc_min if target_gc_min > 1.0 else target_gc_min * 100.0
        target_max_p = target_gc_max if target_gc_max > 1.0 else target_gc_max * 100.0
        dist_to_min = round(final_gc_percent - target_min_p, 2)
        dist_to_max = round(target_max_p - final_gc_percent, 2)

        if abs(final_gc_percent - target_min_p) < 0.2:
            gc_constraint_status = "LOWER_BOUND_ACTIVE"
        elif abs(final_gc_percent - target_max_p) < 0.2:
            gc_constraint_status = "UPPER_BOUND_ACTIVE"
        else:
            gc_constraint_status = "INTERIOR"

        # 50-bp Sliding Window Evaluation
        window_size = 50
        min_50bp_gc = 100.0
        outlier_50bp_count = 0
        if len(optimized_dna) >= window_size:
            for st_idx in range(len(optimized_dna) - window_size + 1):
                win = optimized_dna[st_idx : st_idx + window_size]
                win_gc = (win.count("G") + win.count("C")) / float(window_size) * 100.0
                if win_gc < min_50bp_gc:
                    min_50bp_gc = win_gc
                if win_gc < 25.0 or win_gc > 75.0:
                    outlier_50bp_count += 1
        else:
            min_50bp_gc = final_gc_percent

        # Homopolymer scanning (expression stability threshold >= 6)
        homopolymer_findings = detect_homopolymers(optimized_dna, max_run=6)
        max_homopolymer_run = 0
        if optimized_dna:
            curr_run = 1
            curr_base = optimized_dna[0]
            for ch in optimized_dna[1:]:
                if ch == curr_base:
                    curr_run += 1
                else:
                    if curr_run > max_homopolymer_run:
                        max_homopolymer_run = curr_run
                    curr_base = ch
                    curr_run = 1
            if curr_run > max_homopolymer_run:
                max_homopolymer_run = curr_run

        # A "true 60-nt" initiation-window MFE requires the complete 15-nt
        # upstream context.  Do not silently substitute a CDS-only window.
        if len(up_ctx) < 15:
            mfe_5p_result = {
                "mfe_5p_window_kcal_mol": None,
                "mfe_5p_status": "not_computed",
                "mfe_5p_reason": "insufficient_upstream_context",
                "mfe_5p_warning": "15 nt of upstream context is required for the 60-nt initiation window.",
                "fold_window_nt": 0,
            }
        else:
            mfe_5p_result = compute_5p_mfe_evidence(
                optimized_dna,
                upstream_context=up_ctx,
                window_nt=60,
            )
        mfe_5p_val = mfe_5p_result.get("mfe_5p_window_kcal_mol")

        full_construct_seq = up_ctx + optimized_dna + stop_codon + down_ctx

        try:
            ver = engine_version("dp_v2_1_1")
        except Exception:
            ver = "2.1.1-dev"

        return {
            "sequence": optimized_dna,
            "full_construct": full_construct_seq,
            "stop_codon": stop_codon,
            "cai": global_cai,
            "cai_5p_ramp": ramp_cai,
            "cai_body": body_cai,
            "gc_percent": final_gc_percent,
            "gc_5p_30nt_percent": gc_5p_30nt_percent,
            "gc_5p_45nt_percent": gc_5p_45nt_percent,
            "gc_5p_ramp_percent": ramp_gc_percent,
            "initiation_gc_status": initiation_gc_status,
            "initiation_gc_requested_min_count": nom_init_min,
            "initiation_gc_requested_max_count": nom_init_max,
            "initiation_gc_active_min_count": init_min_gc_count,
            "initiation_gc_active_max_count": init_max_gc_count,
            "initiation_gc_synonymous_min_count": ramp_min_possible_gc,
            "initiation_gc_synonymous_max_count": ramp_max_possible_gc,
            "gc_body_percent": body_gc_percent,
            "gc_feasible": gc_feasible,
            "gc_constraint_status": gc_constraint_status,
            "dist_to_gc_min": dist_to_min,
            "dist_to_gc_max": dist_to_max,
            "min_50bp_gc_percent": min_50bp_gc,
            "outlier_50bp_gc_count": outlier_50bp_count,
            "mfe_5p_window_kcal_mol": mfe_5p_val,
            "mfe_5p_status": mfe_5p_result.get("mfe_5p_status"),
            "mfe_5p_reason": mfe_5p_result.get("mfe_5p_reason"),
            "mfe_5p_warning": mfe_5p_result.get("mfe_5p_warning"),
            "mfe_5p_fold_window_nt": mfe_5p_result.get("fold_window_nt"),
            "mfe_5p_upstream_context_nt": min(15, len(up_ctx)),
            "homopolymer_findings": homopolymer_findings,
            "max_homopolymer_run": max_homopolymer_run,
            "constraint_scope": "FULL_CONSTRUCT"
            if (up_ctx or down_ctx or stop_codon)
            else "CDS_ONLY",
            "score": best_score,
            "length_nt": total_nt,
            "length_aa": n_codons,
            "engine_version": ver,
            "ramp_length_codons": ramp_k,
        }
