"""DP v2.2.0 Optimizer: Unified 3D Dynamic Programming with Initiation Ramp and Structural Pause.

Guarantees & Architectural Contracts:
1. Exact Hard Initiation GC Band Pruning (i = K = 15 codons / 45 nt):
   Enforces ceil(45 * gc_min) <= h_15 <= floor(45 * gc_max) clamped to the achievable
   synonymous GC envelope [GC_min_syn, GC_max_syn] with ZERO additional state dimensions.
2. Position-Conditioned Structural Pause Optimization:
   Evaluates w_struct(i, c) = ell_i(c) - lambda * I(i+1 in L) * |w(c) - tau|.
   When lambda = 0 or L = empty, guarantees exact bit-identical parity with DPV211Optimizer.
3. Homopolymer Aho-Corasick Automaton Rejection:
   Rejects >= 6-mer homopolymers via compiled DFA transitions.
4. Configured Multi-Enzyme Type IIS and Assembly Invariants:
   Full construct elimination across upstream flank, CDS, stop codon, and downstream flank.
5. Deterministic LexRank Sequence Tie-Breaking:
   Guarantees canonical byte-identical, reproducible path selection.
"""

from __future__ import annotations
import math
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE, detect_homopolymers
from factorforge.constraints.type_iis import get_canonical_forbidden_motifs
from factorforge.engines.profile.scoring import compute_5p_mfe_evidence
from factorforge.engines.sllm.automaton import AutomatonCompiler, CompiledAutomaton
from factorforge.utils.exceptions import UnsatisfiableDesignError

AA_TO_CODONS: Dict[str, List[str]] = {}
for codon, aa in STANDARD_GENETIC_CODE.items():
    if aa != "*":
        AA_TO_CODONS.setdefault(aa, []).append(codon)


class DPV220Optimizer:
    """Exact Multi-Constraint Dynamic Programming CDS Compiler (Engine v2.2.0).

    Unifies 5' initiation ramp harmonization, synonymous envelope clamping,
    automaton reachability relaxation, construct tail validation, LexRank tie-breaking,
    and position-conditioned structural pause optimization.
    """

    def __init__(
        self,
        forbidden_enzymes: Optional[Set[str]] = None,
        forbidden_motifs: Optional[Iterable[str]] = None,
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
        pause_positions: Optional[Iterable[int]] = None,
        pause_target_adaptiveness: float = 0.45,
        pause_penalty_weight: float = 0.0,
        structure_reference: Optional[str] = None,
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

        # Structural Pause Conditioning (3DP)
        self.pause_positions = (
            frozenset(pause_positions) if pause_positions is not None else frozenset()
        )
        if any(type(p) is not int or p < 1 for p in self.pause_positions):
            raise ValueError("pause_positions must be positive 1-based integers")
        if (
            not math.isfinite(pause_target_adaptiveness)
            or not 0.0 < pause_target_adaptiveness <= 1.0
        ):
            raise ValueError("pause_target_adaptiveness must be in (0, 1]")
        if not math.isfinite(pause_penalty_weight) or pause_penalty_weight < 0.0:
            raise ValueError("pause_penalty_weight must be a finite non-negative number")
        self.pause_target_adaptiveness = pause_target_adaptiveness
        self.pause_penalty_weight = pause_penalty_weight
        self.structure_reference = structure_reference or "none"

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
        """Computes position-dependent structural log-fitness score w_struct(i, c)."""
        w = codon_weights.get(codon, 1e-4)

        # 1. Base piecewise fitness ell_i(c)
        if position_idx >= self.ramp_length_codons:
            # Stage 2: Elongation Body -> Host Log-CAI
            base_fitness = math.log(max(w, 1e-6))
        else:
            # Stage 1: Initiation Ramp (0 <= position_idx < K)
            denom = max(1, self.ramp_length_codons - 1)
            r_i = self.r_start + (self.r_end - self.r_start) * (position_idx / float(denom))
            h_score = -abs(w - r_i)
            g_count = codon.count("G") + codon.count("C")
            au_count = 3 - g_count
            pen_open = -self.gamma * g_count + self.delta * au_count
            base_fitness = self.alpha * h_score + self.beta * pen_open

        # 2. Position-conditioned structural pause penalty
        one_based_pos = position_idx + 1
        if self.pause_penalty_weight > 0.0 and one_based_pos in self.pause_positions:
            pause_penalty = self.pause_penalty_weight * abs(w - self.pause_target_adaptiveness)
            return base_fitness - pause_penalty

        return base_fitness

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
        """Runs exact position-dependent initiation- and pause-aware 3D DP."""
        up_ctx = left_flank if left_flank is not None else upstream_context
        down_ctx = right_flank if right_flank is not None else downstream_context

        protein = "".join(protein_sequence.upper().split()).rstrip("*")
        if not protein:
            raise ValueError("protein_sequence must not be empty")

        for aa in protein:
            if aa not in AA_TO_CODONS:
                raise ValueError(f"Unsupported amino acid residue: {aa}")

        if any(p > len(protein) for p in self.pause_positions):
            raise ValueError("Configured pause position exceeds protein sequence length")

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

                    # 3. Position-dependent log-fitness score (with structural pause)
                    step_log_fitness = self.compute_codon_log_fitness(i, codon, codon_weights)
                    next_score = score + step_log_fitness

                    # Transition path key for true full-prefix LexRank
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
                            cand_key = (prev_path_rank, codon)
                            exist_key = (existing_prev_rank, existing_codon)
                            if cand_key < exist_key:
                                next_layer[state_key] = (next_score, codon, g, s, prev_path_rank)

            if not next_layer:
                raise UnsatisfiableDesignError(
                    f"DP v2.2.0: State space exhausted at position {i + 1} ({aa}). "
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
                f"DP v2.2.0: No candidate satisfies full construct tail constraints "
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
                    f"DP v2.2.0: No synonymous candidate satisfies the requested GC band "
                    f"[{gc_min_frac:.1%}, {gc_max_frac:.1%}] (requires [{min_gc_count}, {max_gc_count}] GC bases of {total_nt} nt) "
                    f"under full construct Type IIS, pause, and initiation constraints."
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

        global_gc = (optimized_dna.count("G") + optimized_dna.count("C")) / float(total_nt)
        ramp_gc = (
            (ramp_dna.count("G") + ramp_dna.count("C")) / float(ramp_nt) if ramp_nt > 0 else 0.0
        )
        body_gc = (
            (body_dna.count("G") + body_dna.count("C")) / float(len(body_dna))
            if len(body_dna) > 0
            else 0.0
        )

        # 50-bp sliding GC window audit
        window_size = 50
        step = 10
        gc_windows = []
        if total_nt >= window_size:
            for w_start in range(0, total_nt - window_size + 1, step):
                w_seq = optimized_dna[w_start : w_start + window_size]
                w_gc = (w_seq.count("G") + w_seq.count("C")) / float(window_size)
                gc_windows.append(w_gc)
        else:
            gc_windows.append(global_gc)

        # Compute 5' MFE proxy evidence
        mfe_evidence = compute_5p_mfe_evidence(
            cds_sequence=optimized_dna,
            upstream_context=up_ctx,
        )

        # Homopolymer runs
        homopolymers = detect_homopolymers(optimized_dna, max_run=self.homopolymer_max_run + 1)

        # Pause-specific reporting
        pause_report = {}
        for p in sorted(self.pause_positions):
            idx = p - 1
            if idx < len(chosen_codons):
                cd = chosen_codons[idx]
                w_cd = codon_weights.get(cd, 1e-4)
                pause_report[p] = {
                    "amino_acid": protein[idx],
                    "chosen_codon": cd,
                    "codon_weight": w_cd,
                    "target_adaptiveness": self.pause_target_adaptiveness,
                    "delta": abs(w_cd - self.pause_target_adaptiveness),
                }

        return {
            "sequence": optimized_dna,
            "score": best_score,
            "protein_length": n_codons,
            "total_nt": total_nt,
            "gc_percent": round(global_gc * 100.0, 2),
            "gc_feasible": gc_feasible,
            "gc_min_target": gc_min_frac,
            "gc_max_target": gc_max_frac,
            "initiation_gc_status": initiation_gc_status,
            "initiation_gc_requested_min_count": nom_init_min,
            "initiation_gc_requested_max_count": nom_init_max,
            "initiation_gc_active_min_count": init_min_gc_count,
            "initiation_gc_active_max_count": init_max_gc_count,
            "initiation_gc_synonymous_min_count": ramp_min_possible_gc,
            "initiation_gc_synonymous_max_count": ramp_max_possible_gc,
            "gc_5p_45nt_percent": round(ramp_gc * 100.0, 2),
            "gc_body_percent": round(body_gc * 100.0, 2),
            "sliding_window_gc": {
                "window_size_bp": window_size,
                "step_bp": step,
                "min_window_gc_percent": round(min(gc_windows) * 100.0, 2),
                "max_window_gc_percent": round(max(gc_windows) * 100.0, 2),
                "mean_window_gc_percent": round((sum(gc_windows) / len(gc_windows)) * 100.0, 2),
            },
            "cai_global": round(global_cai, 4),
            "cai_ramp_45nt": round(ramp_cai, 4),
            "cai_body": round(body_cai, 4),
            "mfe_5p_evidence": mfe_evidence,
            "max_homopolymer_run": (
                max([hp["length"] for hp in homopolymers], default=0)
                if homopolymers
                else (self.homopolymer_max_run)
            ),
            "homopolymers_detected": homopolymers,
            "pause_conditioning": {
                "positions": sorted(list(self.pause_positions)),
                "penalty_weight": self.pause_penalty_weight,
                "target_adaptiveness": self.pause_target_adaptiveness,
                "structure_reference": self.structure_reference,
                "details": pause_report,
            },
            "engine_version": "2.2.0",
        }


from factorforge.core.interfaces import OptimizationResult, OptimizerEngine


class DPV220EngineAdapter(OptimizerEngine):
    """Deterministic Constrained Optimizer (DP v2.2.0) wrapped for Benchmark / Registry."""

    def __init__(self) -> None:
        pass

    @property
    def name(self) -> str:
        return "FactorForge-DP-v2.2.0"

    @property
    def version(self) -> str:
        return "2.2.0"

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

        from factorforge.engines.profile.utils import load_golden_set
        from factorforge.engines.profile.rules.reverse_translator import ReverseTranslator

        golden_table = load_golden_set()
        codon_weights = ReverseTranslator._build_ref_weights(golden_table)

        forbidden_motifs = kwargs.get("forbidden_motifs", ["GGTCTC", "CGTCTC", "GAAGAC"])
        left_flank = kwargs.get("left_flank", "")
        right_flank = kwargs.get("right_flank", "")
        pause_positions = kwargs.get("pause_positions", None)
        pause_target_adaptiveness = kwargs.get("pause_target_adaptiveness", 0.45)
        pause_penalty_weight = kwargs.get("pause_penalty_weight", 0.0)

        optimizer = DPV220Optimizer(
            forbidden_motifs=forbidden_motifs,
            pause_positions=pause_positions,
            pause_target_adaptiveness=pause_target_adaptiveness,
            pause_penalty_weight=pause_penalty_weight,
        )
        dp_result = optimizer.optimize(
            protein_sequence=protein,
            codon_weights=codon_weights,
            target_gc_min=target_gc_min,
            target_gc_max=target_gc_max,
            left_flank=left_flank,
            right_flank=right_flank,
        )

        cds = dp_result["sequence"]
        terminal_stop_policy = kwargs.get("terminal_stop_policy", "preserve")
        if terminal_stop_policy == "append" or (
            terminal_stop_policy == "preserve" and sequence.endswith("*")
        ):
            cds += "TAA"

        metrics = {
            "score": dp_result["score"],
            "cai": dp_result["cai_global"],
            "gc_percent": dp_result["gc_percent"],
            "cai_ramp_45nt": dp_result["cai_ramp_45nt"],
            "cai_body": dp_result["cai_body"],
        }

        return OptimizationResult(
            sequence=cds,
            metrics=metrics,
            metadata={
                "engine": "dp_v2_2",
                "version": self.version,
                "host": host,
                "inference_mode": "exact_3d_automaton_dp_v2_2",
                "gc_feasible": dp_result["gc_feasible"],
                "target_intersection_exists": dp_result["gc_feasible"],
                "feasible": True,
                "pause_conditioning": dp_result["pause_conditioning"],
            },
        )

    def validate(self, sequence: str) -> bool:
        return True
