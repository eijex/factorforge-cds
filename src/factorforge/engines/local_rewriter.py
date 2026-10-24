"""Equivalence-Preserving Constrained Local Codon Rewriter (Job 440 / Paper 5 Comparative Track).

Implements synonymous peephole and coupled two-site local rewriting under the
exact same objective function, configured sequence-level design constraints, and
tie-breaking rules as DPV230Optimizer.

Scope & Governance Boundaries:
1. Amino Acid Translation Equivalence:
   Every candidate rewrite preserves the exact translated amino acid sequence
   T(c') == T(c) == a and full-construct flanking context.
2. Configured Sequence-Level Constraint Invariance:
   Maintains compiled Aho-Corasick automaton non-violation (configured Type IIS
   recognition motifs, reverse complements, and homopolymer runs >= max_run + 1),
   5' initiation ramp GC bounds, and global coding GC bounds.
3. Identical Objective & Budget Accounting:
   Uses DPV230Optimizer's position-conditioned objective (supporting 1DP, 3DP
   static boundaries, and 4DP shifted exit-tunnel windows W_k(L, delta)) and
   tracks local codon transition evaluations; DP is an uncapped reference
   with a separately replayed transition-count estimate, not a matched end-to-end budget.
4. Scientific Disclaimer:
   Satisfaction of configured sequence-level constraints is an in-silico
   computational property only and does NOT guarantee biological safety,
   in-vivo expression, folding, or clinical suitability.
   Biological validation: NOT_AVAILABLE.
   MCTS / A* / sLLM neural guidance and GNN / GRN extensions are deferred to
   separate follow-up comparative experiments.
"""

from __future__ import annotations

import itertools
import math
import time
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from typing import Any

from factorforge.analysis.metrics import (
    STANDARD_GENETIC_CODE,
    detect_homopolymers,
    translate_dna,
)
from factorforge.engines.dp_v2_3 import AA_TO_CODONS, DPV230Optimizer
from factorforge.engines.profile.scoring import compute_5p_mfe_evidence
from factorforge.utils.exceptions import UnsatisfiableDesignError


@dataclass(frozen=True)
class RewriteStepRecord:
    """Audit record for a single accepted equivalence-preserving local rewrite."""

    step_index: int
    move_type: str
    positions: tuple[int, ...]
    old_codons: tuple[str, ...]
    new_codons: tuple[str, ...]
    score_before: float
    score_after: float
    score_delta: float
    gc_count_after: int
    init_gc_count_after: int
    cumulative_step_evaluations: int


@dataclass
class _ProblemSpec:
    """Resolved problem specification shared between Exact DP and Local Rewriter."""

    protein: str
    n_codons: int
    total_nt: int
    up_ctx: str
    stop_codon: str
    down_ctx: str
    s0: int
    s_tail: str
    active_pause_map: dict[int, list[tuple[str, float, float]]]
    gc_min_frac: float
    gc_max_frac: float
    min_gc_count: int
    max_gc_count: int
    target_mid_gc_count: float
    ramp_k: int
    ramp_nt: int
    apply_initiation_gc_guard: bool
    nom_init_min: int
    nom_init_max: int
    init_min_gc_count: int
    init_max_gc_count: int
    ramp_min_possible_gc: int
    ramp_max_possible_gc: int
    initiation_gc_status: str
    pos_codons: list[list[str]]
    pos_fitness: list[dict[str, float]]
    codon_gc: dict[str, int]


def build_4dp_pause_classes(
    boundaries: Iterable[int],
    tunnel_length: int = 35,
    half_width: int = 2,
    weight: float = 2.0,
    target_tau: float = 0.45,
    protein_length: int | None = None,
    class_name: str = "CUSTOM",
) -> dict[str, dict[str, Any]]:
    """Constructs 4DP shifted boundary windows W_k(L, delta) = [b_k + L - delta, b_k + L + delta].

    Note: L=35 is a configurable model assumption per Roadmap v1.7.0 Section 2,
    not a universally validated biological law.
    """
    if type(tunnel_length) is not int or tunnel_length < 0:
        raise ValueError("tunnel_length L must be a non-negative integer")
    if type(half_width) is not int or half_width < 0:
        raise ValueError("half_width delta must be a non-negative integer")
    if type(weight) not in (int, float) or not math.isfinite(weight) or weight < 0.0:
        raise ValueError("weight lambda_k must be a finite non-negative number")
    if (
        type(target_tau) not in (int, float)
        or not math.isfinite(target_tau)
        or not 0.0 <= target_tau <= 1.0
    ):
        raise ValueError("target_tau must be in [0.0, 1.0]")

    if protein_length is not None and (type(protein_length) is not int or protein_length < 1):
        raise ValueError("protein_length must be a positive integer or None")
    if not isinstance(class_name, str) or not class_name.strip():
        raise ValueError("class_name must be nonempty")
    shifted_positions: set[int] = set()
    for b_k in boundaries:
        if type(b_k) is not int or b_k < 1:
            raise ValueError(
                f"Boundary position b_k must be a positive 1-based integer, got {b_k!r}"
            )
        win_start = b_k + tunnel_length - half_width
        win_end = b_k + tunnel_length + half_width
        for pos in range(win_start, win_end + 1):
            if pos >= 1 and (protein_length is None or pos <= protein_length):
                shifted_positions.add(pos)

    return {
        class_name: {
            "positions": sorted(shifted_positions),
            "weight": float(weight),
            "target_tau": float(target_tau),
            "tunnel_length_assumption": tunnel_length,
            "half_width": half_width,
        }
    }


class ConstrainedLocalRewriter:
    """Equivalence-preserving constrained local codon rewriter sharing DPV230Optimizer's contract."""

    def __init__(
        self,
        dp_optimizer: DPV230Optimizer | None = None,
        forbidden_enzymes: set[str] | None = None,
        forbidden_motifs: Iterable[str] | None = None,
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
        pause_positions: Iterable[int] | None = None,
        pause_target_adaptiveness: float = 0.45,
        pause_penalty_weight: float = 0.0,
        auto_annotate_pauses: bool = False,
        pause_classes: dict[str, dict[str, Any]] | None = None,
        structure_reference: str | None = None,
    ) -> None:
        if dp_optimizer is not None:
            self.dp_optimizer = dp_optimizer
        else:
            self.dp_optimizer = DPV230Optimizer(
                forbidden_enzymes=forbidden_enzymes,
                forbidden_motifs=forbidden_motifs,
                include_reverse_complement=include_reverse_complement,
                ramp_length_codons=ramp_length_codons,
                alpha_harmonization=alpha_harmonization,
                beta_mfe_proxy=beta_mfe_proxy,
                gamma_gc_penalty=gamma_gc_penalty,
                delta_au_bonus=delta_au_bonus,
                r_ramp_start=r_ramp_start,
                r_ramp_end=r_ramp_end,
                homopolymer_max_run=homopolymer_max_run,
                initiation_gc_min=initiation_gc_min,
                initiation_gc_max=initiation_gc_max,
                nominal_ramp_gc=nominal_ramp_gc,
                pause_positions=pause_positions,
                pause_target_adaptiveness=pause_target_adaptiveness,
                pause_penalty_weight=pause_penalty_weight,
                auto_annotate_pauses=auto_annotate_pauses,
                pause_classes=pause_classes,
                structure_reference=structure_reference,
            )

    def _resolve_problem_spec(
        self,
        protein_sequence: str,
        codon_weights: dict[str, float],
        target_gc_min: float = 0.40,
        target_gc_max: float = 0.47,
        upstream_context: str = "",
        stop_codon: str = "TAA",
        downstream_context: str = "",
    ) -> _ProblemSpec:
        """Resolves objective tables and constraint bounds identically to DPV230Optimizer."""
        protein = "".join(protein_sequence.upper().split()).rstrip("*")
        if not protein:
            raise ValueError("protein_sequence must not be empty")
        for aa in protein:
            if aa not in AA_TO_CODONS:
                raise ValueError(f"Unsupported amino acid residue: {aa}")

        up_ctx = upstream_context.upper()
        down_ctx = downstream_context.upper()
        stop_up = stop_codon.upper()
        if stop_up and (len(stop_up) != 3 or STANDARD_GENETIC_CODE.get(stop_up) != "*"):
            raise ValueError(f"Invalid stop_codon: {stop_codon!r}")

        automaton = self.dp_optimizer.automaton
        s0 = 0
        for ch in up_ctx:
            s0 = automaton.step_nucleotide(s0, ch)
            if s0 in automaton.terminal_states:
                raise ValueError(
                    f"Upstream context '{upstream_context}' contains a forbidden motif."
                )

        n_codons = len(protein)
        total_nt = n_codons * 3
        active_pause_map = self.dp_optimizer._resolve_active_pause_map(protein)

        gc_min_frac = target_gc_min if target_gc_min <= 1.0 else target_gc_min / 100.0
        gc_max_frac = target_gc_max if target_gc_max <= 1.0 else target_gc_max / 100.0
        if not (0.0 <= gc_min_frac <= gc_max_frac <= 1.0):
            raise ValueError("target_gc_min and target_gc_max must satisfy 0 <= min <= max <= 1")

        min_gc_count = math.ceil(gc_min_frac * total_nt)
        max_gc_count = math.floor(gc_max_frac * total_nt)
        if min_gc_count > max_gc_count:
            mid_count = round(((gc_min_frac + gc_max_frac) / 2.0) * total_nt)
            min_gc_count = mid_count
            max_gc_count = mid_count
        target_mid_gc_count = (min_gc_count + max_gc_count) / 2.0

        ramp_k = min(self.dp_optimizer.ramp_length_codons, n_codons)
        ramp_nt = ramp_k * 3
        apply_initiation_gc_guard = n_codons >= self.dp_optimizer.ramp_length_codons
        nom_init_min = math.ceil(self.dp_optimizer.initiation_gc_min * ramp_nt)
        nom_init_max = math.floor(self.dp_optimizer.initiation_gc_max * ramp_nt)
        if nom_init_min > nom_init_max:
            mid_init = round(self.dp_optimizer.nominal_ramp_gc * ramp_nt)
            nom_init_min = mid_init
            nom_init_max = mid_init

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

        codon_gc: dict[str, int] = {
            c: c.count("G") + c.count("C")
            for c in STANDARD_GENETIC_CODE
            if STANDARD_GENETIC_CODE[c] != "*"
        }
        pos_codons: list[list[str]] = []
        pos_fitness: list[dict[str, float]] = []
        for i, aa in enumerate(protein):
            syns = sorted(AA_TO_CODONS[aa])
            pos_codons.append(syns)
            pos_fitness.append(
                {
                    c: self.dp_optimizer.compute_codon_log_fitness(
                        i, c, codon_weights, active_pause_map
                    )
                    for c in syns
                }
            )

        return _ProblemSpec(
            protein=protein,
            n_codons=n_codons,
            total_nt=total_nt,
            up_ctx=up_ctx,
            stop_codon=stop_up,
            down_ctx=down_ctx,
            s0=s0,
            s_tail=stop_up + down_ctx,
            active_pause_map=active_pause_map,
            gc_min_frac=gc_min_frac,
            gc_max_frac=gc_max_frac,
            min_gc_count=min_gc_count,
            max_gc_count=max_gc_count,
            target_mid_gc_count=target_mid_gc_count,
            ramp_k=ramp_k,
            ramp_nt=ramp_nt,
            apply_initiation_gc_guard=apply_initiation_gc_guard,
            nom_init_min=nom_init_min,
            nom_init_max=nom_init_max,
            init_min_gc_count=init_min_gc_count,
            init_max_gc_count=init_max_gc_count,
            ramp_min_possible_gc=ramp_min_possible_gc,
            ramp_max_possible_gc=ramp_max_possible_gc,
            initiation_gc_status=initiation_gc_status,
            pos_codons=pos_codons,
            pos_fitness=pos_fitness,
            codon_gc=codon_gc,
        )

    def evaluate_candidate(
        self,
        protein_sequence: str,
        coding_dna: str,
        codon_weights: dict[str, float],
        target_gc_min: float = 0.40,
        target_gc_max: float = 0.47,
        upstream_context: str = "",
        stop_codon: str = "TAA",
        downstream_context: str = "",
    ) -> dict[str, Any]:
        """Evaluates a candidate coding DNA sequence under the shared objective and constraints."""
        spec = self._resolve_problem_spec(
            protein_sequence=protein_sequence,
            codon_weights=codon_weights,
            target_gc_min=target_gc_min,
            target_gc_max=target_gc_max,
            upstream_context=upstream_context,
            stop_codon=stop_codon,
            downstream_context=downstream_context,
        )
        dna = "".join(coding_dna.upper().split())
        violations: list[str] = []

        if len(dna) != spec.total_nt or len(dna) % 3 != 0:
            violations.append(f"length_mismatch:{len(dna)}!={spec.total_nt}")
            return {
                "is_feasible": False,
                "translation_equivalent": False,
                "automaton_valid": False,
                "initiation_gc_feasible": False,
                "global_gc_feasible": False,
                "score": -math.inf,
                "gc_count": 0,
                "init_gc_count": 0,
                "violations": violations,
            }

        translated = translate_dna(dna)
        translation_equivalent = translated == spec.protein
        if not translation_equivalent:
            violations.append("amino_acid_translation_mismatch")

        codons = [dna[i * 3 : (i + 1) * 3] for i in range(spec.n_codons)]
        automaton = self.dp_optimizer.automaton
        curr_s = spec.s0
        automaton_valid = True
        for i, c in enumerate(codons):
            next_s, is_forbidden = automaton.step_codon(curr_s, c)
            if is_forbidden:
                automaton_valid = False
                violations.append(f"automaton_violation_at_codon_{i + 1}")
                break
            curr_s = next_s

        if automaton_valid:
            for ch in spec.s_tail:
                curr_s = automaton.step_nucleotide(curr_s, ch)
                if curr_s in automaton.terminal_states:
                    automaton_valid = False
                    violations.append("automaton_violation_in_construct_tail")
                    break

        init_gc_count = sum(spec.codon_gc.get(c, 0) for c in codons[: spec.ramp_k])
        total_gc_count = sum(spec.codon_gc.get(c, 0) for c in codons)

        initiation_gc_feasible = True
        if spec.apply_initiation_gc_guard and not (
            spec.init_min_gc_count <= init_gc_count <= spec.init_max_gc_count
        ):
            initiation_gc_feasible = False
            violations.append(
                f"initiation_gc_out_of_bounds:{init_gc_count}_not_in_[{spec.init_min_gc_count},{spec.init_max_gc_count}]"
            )

        global_gc_feasible = spec.min_gc_count <= total_gc_count <= spec.max_gc_count
        if not global_gc_feasible:
            violations.append(
                f"global_gc_out_of_bounds:{total_gc_count}_not_in_[{spec.min_gc_count},{spec.max_gc_count}]"
            )

        score = (
            sum(spec.pos_fitness[i][codons[i]] for i in range(spec.n_codons))
            if translation_equivalent
            else -math.inf
        )

        is_feasible = (
            translation_equivalent
            and automaton_valid
            and initiation_gc_feasible
            and global_gc_feasible
        )
        return {
            "is_feasible": is_feasible,
            "translation_equivalent": translation_equivalent,
            "automaton_valid": automaton_valid,
            "initiation_gc_feasible": initiation_gc_feasible,
            "global_gc_feasible": global_gc_feasible,
            "score": score,
            "gc_count": total_gc_count,
            "init_gc_count": init_gc_count,
            "violations": violations,
        }

    def _compute_prefix_states(self, spec: _ProblemSpec, codons: Sequence[str]) -> list[int] | None:
        """Computes automaton state trajectory [s0, s1, ..., sN] and checks tail validity."""
        automaton = self.dp_optimizer.automaton
        states = [spec.s0]
        curr_s = spec.s0
        for c in codons:
            curr_s, is_forbidden = automaton.step_codon(curr_s, c)
            if is_forbidden:
                return None
            states.append(curr_s)
        tail_s = curr_s
        for ch in spec.s_tail:
            tail_s = automaton.step_nucleotide(tail_s, ch)
            if tail_s in automaton.terminal_states:
                return None
        return states

    def _verify_suffix_from_window(
        self,
        spec: _ProblemSpec,
        codons: Sequence[str],
        prefix_states: Sequence[int],
        end_idx: int,
        state_at_end: int,
        step_counter: list[int],
        max_step_evaluations: int | None,
    ) -> tuple[bool, list[int]]:
        """Verifies automaton validity from end_idx until state convergence with prefix_states."""
        if state_at_end == prefix_states[end_idx]:
            return True, []

        automaton = self.dp_optimizer.automaton
        curr_s = state_at_end
        updated_suffix_states: list[int] = []
        k = end_idx
        while k < spec.n_codons:
            if max_step_evaluations is not None and step_counter[0] >= max_step_evaluations:
                return False, []
            step_counter[0] += 1
            curr_s, is_forbidden = automaton.step_codon(curr_s, codons[k])
            if is_forbidden:
                return False, []
            k += 1
            if curr_s == prefix_states[k]:
                return True, updated_suffix_states
            updated_suffix_states.append(curr_s)

        # Reached end of coding region with a different final state; verify s_tail
        tail_s = curr_s
        for ch in spec.s_tail:
            tail_s = automaton.step_nucleotide(tail_s, ch)
            if tail_s in automaton.terminal_states:
                return False, []
        return True, updated_suffix_states

    def _construct_feasible_seed(
        self,
        spec: _ProblemSpec,
        seed_mode: str,
        step_counter: list[int],
        max_step_evaluations: int | None,
    ) -> list[str]:
        """Constructs an initial synonymous feasible seed using a narrow GC-guided beam."""
        automaton = self.dp_optimizer.automaton
        # Narrow beam width to find a valid starting point quickly within budget
        beam_width = 6
        # State: (g, s) -> (priority_score, codons_list)
        beam: dict[tuple[int, int], tuple[float, list[str]]] = {(0, spec.s0): (0.0, [])}

        target_init_mid = (spec.init_min_gc_count + spec.init_max_gc_count) / 2.0
        target_total_mid = spec.target_mid_gc_count

        for i in range(spec.n_codons):
            next_beam: dict[tuple[int, int], tuple[float, list[str]]] = {}
            syns = spec.pos_codons[i]
            is_ramp_boundary = i == spec.ramp_k - 1

            # Target cumulative GC at step i + 1
            if i < spec.ramp_k:
                desired_g = target_init_mid * ((i + 1) / float(max(1, spec.ramp_k)))
            else:
                rem_frac = (i + 1 - spec.ramp_k) / float(max(1, spec.n_codons - spec.ramp_k))
                desired_g = target_init_mid + (target_total_mid - target_init_mid) * rem_frac

            for (g, s), (acc_score, prefix) in beam.items():
                for c in syns:
                    if max_step_evaluations is not None and step_counter[0] >= max_step_evaluations:
                        break
                    step_counter[0] += 1
                    next_s, is_forbidden = automaton.step_codon(s, c)
                    if is_forbidden:
                        continue
                    next_g = g + spec.codon_gc[c]
                    fit = spec.pos_fitness[i][c]
                    if seed_mode == "low_fitness_feasible":
                        # Deliberately prefer lower-fitness feasible codons to test rewriter ascent
                        step_pref = -fit - 2.5 * abs(next_g - desired_g)
                    else:
                        # Default: prioritize GC trajectory feasibility with mild fitness guidance
                        step_pref = 0.25 * fit - 2.5 * abs(next_g - desired_g)
                    cand_score = acc_score + step_pref
                    key = (next_g, next_s)
                    cand_seq = prefix + [c]
                    if key not in next_beam or (cand_score, [-ord(x[0]) for x in cand_seq]) > (
                        next_beam[key][0],
                        [-ord(x[0]) for x in next_beam[key][1]],
                    ):
                        next_beam[key] = (cand_score, cand_seq)

            if not next_beam:
                raise UnsatisfiableDesignError(
                    f"Local rewriter seed construction exhausted budget or state space at codon {i + 1}."
                )

            if is_ramp_boundary and spec.apply_initiation_gc_guard:
                in_init = {
                    k: v
                    for k, v in next_beam.items()
                    if spec.init_min_gc_count <= k[0] <= spec.init_max_gc_count
                }
                if in_init:
                    next_beam = in_init

            # Keep top `beam_width` states closest to desired GC trajectory
            ranked = sorted(
                next_beam.items(),
                key=lambda item: (
                    item[1][0],
                    -abs(item[0][0] - desired_g),
                ),
                reverse=True,
            )
            beam = dict(ranked[:beam_width])

        # Filter by tail validity and global GC band
        valid_candidates: list[tuple[int, float, list[str]]] = []
        for (g, s), (sc, seq) in beam.items():
            if not (spec.min_gc_count <= g <= spec.max_gc_count):
                continue
            tail_s = s
            tail_ok = True
            for ch in spec.s_tail:
                tail_s = automaton.step_nucleotide(tail_s, ch)
                if tail_s in automaton.terminal_states:
                    tail_ok = False
                    break
            if tail_ok:
                valid_candidates.append((g, sc, seq))

        if not valid_candidates:
            raise UnsatisfiableDesignError(
                "Local rewriter seed construction could not reach a feasible full-construct state within budget."
            )

        valid_candidates.sort(
            key=lambda item: (item[1], -abs(item[0] - target_total_mid)), reverse=True
        )
        return valid_candidates[0][2]

    def _repair_infeasible_seed(
        self,
        spec: _ProblemSpec,
        codons: list[str],
        step_counter: list[int],
        max_step_evaluations: int | None,
    ) -> list[str]:
        """Repairs an initially infeasible synonymous sequence via local window edits."""
        automaton = self.dp_optimizer.automaton
        current = list(codons)

        # 1. Repair automaton violations if any
        for _ in range(3):
            states = self._compute_prefix_states(spec, current)
            if states is not None:
                break
            # Locate first offending codon
            curr_s = spec.s0
            bad_idx = spec.n_codons - 1
            for i, c in enumerate(current):
                next_s, is_forbidden = automaton.step_codon(curr_s, c)
                if is_forbidden:
                    bad_idx = i
                    break
                curr_s = next_s
            win_start = max(0, bad_idx - 1)
            win_end = min(spec.n_codons, bad_idx + 2)
            choices = [spec.pos_codons[k] for k in range(win_start, win_end)]
            repaired = False
            for combo in itertools.product(*choices):
                if (
                    max_step_evaluations is not None
                    and step_counter[0] + len(combo) > max_step_evaluations
                ):
                    break
                step_counter[0] += len(combo)
                trial = current[:win_start] + list(combo) + current[win_end:]
                if self._compute_prefix_states(spec, trial) is not None:
                    current = trial
                    repaired = True
                    break
            if not repaired:
                break

        prefix_states = self._compute_prefix_states(spec, current)
        if prefix_states is None:
            # Fall back to feasible seed constructor if local motif repair could not resolve
            return self._construct_feasible_seed(
                spec, "feasible_greedy", step_counter, max_step_evaluations
            )

        # 2. Repair initiation GC and global GC bounds via single-codon synonymous steps
        for _ in range(spec.n_codons * 2):
            init_gc = sum(spec.codon_gc[c] for c in current[: spec.ramp_k])
            total_gc = sum(spec.codon_gc[c] for c in current)

            init_deficit = 0
            if spec.apply_initiation_gc_guard:
                if init_gc < spec.init_min_gc_count:
                    init_deficit = spec.init_min_gc_count - init_gc
                elif init_gc > spec.init_max_gc_count:
                    init_deficit = spec.init_max_gc_count - init_gc  # negative

            total_deficit = 0
            if total_gc < spec.min_gc_count:
                total_deficit = spec.min_gc_count - total_gc
            elif total_gc > spec.max_gc_count:
                total_deficit = spec.max_gc_count - total_gc  # negative

            if init_deficit == 0 and total_deficit == 0:
                return current

            best_move = None
            best_gain = (-1, -math.inf)
            search_range = range(spec.ramp_k) if init_deficit != 0 else range(spec.n_codons)
            for i in search_range:
                old_c = current[i]
                old_gc = spec.codon_gc[old_c]
                for c in spec.pos_codons[i]:
                    if c == old_c:
                        continue
                    if max_step_evaluations is not None and step_counter[0] >= max_step_evaluations:
                        break
                    step_counter[0] += 1
                    next_s, is_forbidden = automaton.step_codon(prefix_states[i], c)
                    if is_forbidden:
                        continue
                    dg = spec.codon_gc[c] - old_gc
                    if init_deficit > 0 and dg <= 0:
                        continue
                    if init_deficit < 0 and dg >= 0:
                        continue
                    if init_deficit == 0:
                        if (
                            i < spec.ramp_k
                            and spec.apply_initiation_gc_guard
                            and not (
                                spec.init_min_gc_count <= init_gc + dg <= spec.init_max_gc_count
                            )
                        ):
                            continue
                        if total_deficit > 0 and dg <= 0:
                            continue
                        if total_deficit < 0 and dg >= 0:
                            continue
                    ok_suf, _ = self._verify_suffix_from_window(
                        spec,
                        current,
                        prefix_states,
                        i + 1,
                        next_s,
                        step_counter,
                        max_step_evaluations,
                    )
                    if not ok_suf:
                        continue
                    dfit = spec.pos_fitness[i][c] - spec.pos_fitness[i][old_c]
                    cand_gain = (abs(dg), dfit)
                    if cand_gain > best_gain:
                        best_gain = cand_gain
                        best_move = (i, c)

            if best_move is None:
                break
            idx, new_c = best_move
            current[idx] = new_c
            recomputed = self._compute_prefix_states(spec, current)
            if recomputed is None:
                break
            prefix_states = recomputed

        return current

    def rewrite(
        self,
        protein_sequence: str,
        codon_weights: dict[str, float],
        initial_coding_dna: str | None = None,
        seed_mode: str = "feasible_greedy",
        window_sizes: Sequence[int] = (1, 2, 3),
        include_paired_gc_swaps: bool = True,
        max_step_evaluations: int | None = None,
        max_passes: int = 8,
        require_feasible_seed: bool = False,
        target_gc_min: float = 0.40,
        target_gc_max: float = 0.47,
        upstream_context: str = "",
        stop_codon: str = "TAA",
        downstream_context: str = "",
    ) -> dict[str, Any]:
        """Executes equivalence-preserving constrained local codon rewriting.

        All accepted rewrites preserve exact amino acid translation equivalence and
        configured sequence-level constraints while monotonically increasing the
        shared DPV230Optimizer objective score.
        """
        if max_step_evaluations is not None and (
            type(max_step_evaluations) is not int or max_step_evaluations < 1
        ):
            raise ValueError("max_step_evaluations must be a positive integer or None")
        if type(max_passes) is not int or max_passes < 1:
            raise ValueError("max_passes must be a positive integer")
        if not window_sizes or any(type(w) is not int or w < 1 or w > 5 for w in window_sizes):
            raise ValueError("window_sizes must be a non-empty sequence of integers in [1, 5]")

        t_start = time.perf_counter()
        spec = self._resolve_problem_spec(
            protein_sequence=protein_sequence,
            codon_weights=codon_weights,
            target_gc_min=target_gc_min,
            target_gc_max=target_gc_max,
            upstream_context=upstream_context,
            stop_codon=stop_codon,
            downstream_context=downstream_context,
        )

        step_counter = [0]
        seed_repaired = False

        if initial_coding_dna is not None:
            clean_seed = "".join(initial_coding_dna.upper().split())
            if len(clean_seed) != spec.total_nt or translate_dna(clean_seed) != spec.protein:
                raise ValueError(
                    "initial_coding_dna violates amino acid translation equivalence with protein_sequence."
                )
            codons = [clean_seed[i * 3 : (i + 1) * 3] for i in range(spec.n_codons)]
            init_eval = self.evaluate_candidate(
                protein_sequence=spec.protein,
                coding_dna=clean_seed,
                codon_weights=codon_weights,
                target_gc_min=target_gc_min,
                target_gc_max=target_gc_max,
                upstream_context=upstream_context,
                stop_codon=stop_codon,
                downstream_context=downstream_context,
            )
            if not init_eval["is_feasible"]:
                if require_feasible_seed:
                    raise UnsatisfiableDesignError(
                        f"Initial seed is infeasible ({init_eval['violations']}) and require_feasible_seed=True."
                    )
                codons = self._repair_infeasible_seed(
                    spec, codons, step_counter, max_step_evaluations
                )
                seed_repaired = True
        else:
            codons = self._construct_feasible_seed(
                spec, seed_mode, step_counter, max_step_evaluations
            )

        prefix_states = self._compute_prefix_states(spec, codons)
        init_gc = sum(spec.codon_gc[c] for c in codons[: spec.ramp_k])
        total_gc = sum(spec.codon_gc[c] for c in codons)

        if (
            prefix_states is None
            or (
                spec.apply_initiation_gc_guard
                and not (spec.init_min_gc_count <= init_gc <= spec.init_max_gc_count)
            )
            or not (spec.min_gc_count <= total_gc <= spec.max_gc_count)
        ):
            raise UnsatisfiableDesignError(
                "Could not establish a feasible initial sequence within the computational budget."
            )

        current_score = sum(spec.pos_fitness[i][codons[i]] for i in range(spec.n_codons))
        initial_feasible_dna = "".join(codons)
        initial_feasible_score = current_score

        trace: list[RewriteStepRecord] = []
        automaton = self.dp_optimizer.automaton
        passes_executed = 0
        budget_exhausted = False

        for _ in range(max_passes):
            if max_step_evaluations is not None and step_counter[0] >= max_step_evaluations:
                budget_exhausted = True
                break
            passes_executed += 1
            accepted_in_pass = 0

            # 1. Contiguous Sliding-Window Peephole Rewrites
            for w in window_sizes:
                if w > spec.n_codons:
                    continue
                for start in range(spec.n_codons - w + 1):
                    if max_step_evaluations is not None and step_counter[0] >= max_step_evaluations:
                        budget_exhausted = True
                        break
                    end = start + w
                    old_win = tuple(codons[start:end])
                    old_win_score = sum(spec.pos_fitness[k][codons[k]] for k in range(start, end))
                    old_win_total_gc = sum(spec.codon_gc[codons[k]] for k in range(start, end))
                    old_win_init_gc = sum(
                        spec.codon_gc[codons[k]] for k in range(start, min(end, spec.ramp_k))
                    )

                    best_win: tuple[str, ...] | None = None
                    best_win_delta = 0.0
                    best_win_total_gc = total_gc
                    best_win_init_gc = init_gc
                    best_suffix_states: list[int] = []
                    best_inter_states: list[int] = []

                    # Recursive depth-first window enumeration with immediate automaton pruning
                    def _search_window(
                        offset: int,
                        curr_s: int,
                        acc_fit: float,
                        acc_tot_gc: int,
                        acc_init_gc: int,
                        chosen_codons: list[str],
                        chosen_states: list[int],
                        *,
                        w=w,
                        old_win=old_win,
                        init_gc=init_gc,
                        old_win_init_gc=old_win_init_gc,
                        total_gc=total_gc,
                        old_win_total_gc=old_win_total_gc,
                        old_win_score=old_win_score,
                        codons=codons,
                        prefix_states=prefix_states,
                        end=end,
                        start=start,
                    ) -> None:
                        nonlocal best_win, best_win_delta, best_win_total_gc, best_win_init_gc
                        nonlocal best_suffix_states, best_inter_states, budget_exhausted

                        if budget_exhausted:
                            return

                        if offset == w:
                            cand_win = tuple(chosen_codons)
                            if cand_win == old_win:
                                return
                            new_init_gc = init_gc - old_win_init_gc + acc_init_gc
                            if spec.apply_initiation_gc_guard and not (
                                spec.init_min_gc_count <= new_init_gc <= spec.init_max_gc_count
                            ):
                                return
                            new_total_gc = total_gc - old_win_total_gc + acc_tot_gc
                            if not (spec.min_gc_count <= new_total_gc <= spec.max_gc_count):
                                return

                            delta_score = acc_fit - old_win_score
                            # Apply identical priority hierarchy as DPV230Optimizer:
                            # Primary: score improvement > 1e-12
                            # Tie 1: closer to target_mid_gc_count
                            # Tie 2: lexicographically smaller codon tuple
                            ref_win = best_win if best_win is not None else old_win
                            ref_delta = best_win_delta if best_win is not None else 0.0
                            ref_tot_gc = best_win_total_gc if best_win is not None else total_gc

                            is_better = False
                            if delta_score > ref_delta + 1e-12:
                                is_better = True
                            elif abs(delta_score - ref_delta) <= 1e-12:
                                cand_gc_dist = abs(new_total_gc - spec.target_mid_gc_count)
                                ref_gc_dist = abs(ref_tot_gc - spec.target_mid_gc_count)
                                if (
                                    cand_gc_dist < ref_gc_dist - 1e-12
                                    or abs(cand_gc_dist - ref_gc_dist) <= 1e-12
                                    and cand_win < ref_win
                                ):
                                    is_better = True

                            if not is_better:
                                return

                            ok_suf, suf_states = self._verify_suffix_from_window(
                                spec,
                                codons,
                                prefix_states,
                                end,
                                curr_s,
                                step_counter,
                                max_step_evaluations,
                            )
                            if not ok_suf:
                                if (
                                    max_step_evaluations is not None
                                    and step_counter[0] >= max_step_evaluations
                                ):
                                    budget_exhausted = True
                                return

                            best_win = cand_win
                            best_win_delta = delta_score
                            best_win_total_gc = new_total_gc
                            best_win_init_gc = new_init_gc
                            best_inter_states = list(chosen_states)
                            best_suffix_states = suf_states
                            return

                        pos = start + offset
                        for c in spec.pos_codons[pos]:
                            if (
                                max_step_evaluations is not None
                                and step_counter[0] >= max_step_evaluations
                            ):
                                budget_exhausted = True
                                return
                            step_counter[0] += 1
                            next_s, is_forbidden = automaton.step_codon(curr_s, c)
                            if is_forbidden:
                                continue
                            cg = spec.codon_gc[c]
                            _search_window(
                                offset + 1,
                                next_s,
                                acc_fit + spec.pos_fitness[pos][c],
                                acc_tot_gc + cg,
                                acc_init_gc + (cg if pos < spec.ramp_k else 0),
                                chosen_codons + [c],
                                chosen_states + [next_s],
                            )

                    _search_window(0, prefix_states[start], 0.0, 0, 0, [], [])

                    if best_win is not None:
                        score_before = current_score
                        for idx_off, new_c in enumerate(best_win):
                            codons[start + idx_off] = new_c
                            prefix_states[start + idx_off + 1] = best_inter_states[idx_off]
                        for suf_off, suf_s in enumerate(best_suffix_states):
                            prefix_states[end + 1 + suf_off] = suf_s
                        total_gc = best_win_total_gc
                        init_gc = best_win_init_gc
                        current_score = sum(
                            spec.pos_fitness[i][codons[i]] for i in range(spec.n_codons)
                        )
                        accepted_in_pass += 1
                        trace.append(
                            RewriteStepRecord(
                                step_index=len(trace) + 1,
                                move_type=f"peephole_w{w}",
                                positions=tuple(range(start, end)),
                                old_codons=old_win,
                                new_codons=best_win,
                                score_before=score_before,
                                score_after=current_score,
                                score_delta=current_score - score_before,
                                gc_count_after=total_gc,
                                init_gc_count_after=init_gc,
                                cumulative_step_evaluations=step_counter[0],
                            )
                        )

                if budget_exhausted:
                    break

            if budget_exhausted:
                break

            # 2. Coupled Two-Site GC-Compensating Rewrites
            if include_paired_gc_swaps and spec.n_codons >= 2:
                # Collect best single-codon alternatives per position and GC delta
                # Entry: (pos, new_codon, delta_fit, delta_total_gc, delta_init_gc)
                candidates_by_dg: dict[int, list[tuple[int, str, float, int, int]]] = {}
                for i in range(spec.n_codons):
                    if max_step_evaluations is not None and step_counter[0] >= max_step_evaluations:
                        budget_exhausted = True
                        break
                    old_c = codons[i]
                    old_fit = spec.pos_fitness[i][old_c]
                    old_g = spec.codon_gc[old_c]
                    best_per_dg: dict[int, tuple[str, float]] = {}
                    for c in spec.pos_codons[i]:
                        if c == old_c:
                            continue
                        if (
                            max_step_evaluations is not None
                            and step_counter[0] >= max_step_evaluations
                        ):
                            budget_exhausted = True
                            break
                        step_counter[0] += 1
                        next_s, is_forbidden = automaton.step_codon(prefix_states[i], c)
                        if is_forbidden:
                            continue
                        ok_suf, _ = self._verify_suffix_from_window(
                            spec,
                            codons,
                            prefix_states,
                            i + 1,
                            next_s,
                            step_counter,
                            max_step_evaluations,
                        )
                        if not ok_suf:
                            continue
                        dg = spec.codon_gc[c] - old_g
                        if dg == 0:
                            continue
                        dfit = spec.pos_fitness[i][c] - old_fit
                        if dg not in best_per_dg or dfit > best_per_dg[dg][1] + 1e-12:
                            best_per_dg[dg] = (c, dfit)
                    for dg, (c_best, dfit_best) in best_per_dg.items():
                        d_init = dg if i < spec.ramp_k else 0
                        candidates_by_dg.setdefault(dg, []).append(
                            (i, c_best, dfit_best, dg, d_init)
                        )

                if not budget_exhausted:
                    for dg_list in candidates_by_dg.values():
                        dg_list.sort(key=lambda x: x[2], reverse=True)

                    best_pair: tuple[int, str, int, str, float] | None = None
                    for dg1, list1 in candidates_by_dg.items():
                        for dg2, list2 in candidates_by_dg.items():
                            if dg1 > dg2:
                                continue
                            new_tot_gc = total_gc + dg1 + dg2
                            if not (spec.min_gc_count <= new_tot_gc <= spec.max_gc_count):
                                continue
                            for i1, c1, df1, _, di1 in list1[:6]:
                                for i2, c2, df2, _, di2 in list2[:6]:
                                    if i1 == i2:
                                        continue
                                    pair_delta = df1 + df2
                                    if pair_delta <= 1e-12:
                                        continue
                                    new_in_gc = init_gc + di1 + di2
                                    if spec.apply_initiation_gc_guard and not (
                                        spec.init_min_gc_count
                                        <= new_in_gc
                                        <= spec.init_max_gc_count
                                    ):
                                        continue
                                    if best_pair is None or pair_delta > best_pair[4] + 1e-12:
                                        pos_a, cod_a, pos_b, cod_b = (
                                            (i1, c1, i2, c2) if i1 < i2 else (i2, c2, i1, c1)
                                        )
                                        best_pair = (pos_a, cod_a, pos_b, cod_b, pair_delta)

                    if best_pair is not None:
                        pos_a, cod_a, pos_b, cod_b, _ = best_pair
                        trial = list(codons)
                        old_a, old_b = trial[pos_a], trial[pos_b]
                        trial[pos_a] = cod_a
                        trial[pos_b] = cod_b
                        new_states = self._compute_prefix_states(spec, trial)
                        if new_states is not None:
                            score_before = current_score
                            codons = trial
                            prefix_states = new_states
                            init_gc = sum(spec.codon_gc[c] for c in codons[: spec.ramp_k])
                            total_gc = sum(spec.codon_gc[c] for c in codons)
                            current_score = sum(
                                spec.pos_fitness[i][codons[i]] for i in range(spec.n_codons)
                            )
                            accepted_in_pass += 1
                            trace.append(
                                RewriteStepRecord(
                                    step_index=len(trace) + 1,
                                    move_type="coupled_gc_swap",
                                    positions=(pos_a, pos_b),
                                    old_codons=(old_a, old_b),
                                    new_codons=(cod_a, cod_b),
                                    score_before=score_before,
                                    score_after=current_score,
                                    score_delta=current_score - score_before,
                                    gc_count_after=total_gc,
                                    init_gc_count_after=init_gc,
                                    cumulative_step_evaluations=step_counter[0],
                                )
                            )

            if accepted_in_pass == 0:
                break

        wall_time_ms = (time.perf_counter() - t_start) * 1000.0
        optimized_dna = "".join(codons)
        final_audit = self.evaluate_candidate(
            protein_sequence=spec.protein,
            coding_dna=optimized_dna,
            codon_weights=codon_weights,
            target_gc_min=target_gc_min,
            target_gc_max=target_gc_max,
            upstream_context=upstream_context,
            stop_codon=stop_codon,
            downstream_context=downstream_context,
        )
        if not final_audit["is_feasible"]:
            raise UnsatisfiableDesignError(
                f"Internal invariant violation: final rewritten sequence failed verification: {final_audit['violations']}"
            )

        ramp_codons = codons[: spec.ramp_k]
        body_codons = codons[spec.ramp_k :]
        ramp_dna = "".join(ramp_codons)
        body_dna = "".join(body_codons)

        global_cai = math.exp(
            sum(math.log(max(codon_weights.get(c, 1e-4), 1e-6)) for c in codons) / spec.n_codons
        )
        ramp_cai = (
            math.exp(
                sum(math.log(max(codon_weights.get(c, 1e-4), 1e-6)) for c in ramp_codons)
                / spec.ramp_k
            )
            if spec.ramp_k > 0
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

        global_gc = total_gc / float(spec.total_nt)
        ramp_gc = (
            (ramp_dna.count("G") + ramp_dna.count("C")) / float(spec.ramp_nt)
            if spec.ramp_nt > 0
            else 0.0
        )
        body_gc = (
            (body_dna.count("G") + body_dna.count("C")) / float(len(body_dna))
            if len(body_dna) > 0
            else 0.0
        )

        mfe_evidence = compute_5p_mfe_evidence(
            cds_sequence=optimized_dna,
            upstream_context=spec.up_ctx,
        )
        homopolymers = detect_homopolymers(
            optimized_dna, max_run=self.dp_optimizer.homopolymer_max_run + 1
        )

        return {
            "sequence": optimized_dna,
            "score": current_score,
            "initial_feasible_sequence": initial_feasible_dna,
            "initial_feasible_score": initial_feasible_score,
            "score_improvement": current_score - initial_feasible_score,
            "protein_length": spec.n_codons,
            "total_nt": spec.total_nt,
            "gc_percent": round(global_gc * 100.0, 2),
            "gc_feasible": final_audit["global_gc_feasible"],
            "initiation_gc_feasible": final_audit["initiation_gc_feasible"],
            "automaton_valid": final_audit["automaton_valid"],
            "translation_equivalent": final_audit["translation_equivalent"],
            "is_feasible": final_audit["is_feasible"],
            "initiation_gc_status": spec.initiation_gc_status,
            "gc_5p_45nt_percent": round(ramp_gc * 100.0, 2),
            "gc_body_percent": round(body_gc * 100.0, 2),
            "cai_global": round(global_cai, 4),
            "cai_ramp_45nt": round(ramp_cai, 4),
            "cai_body": round(body_cai, 4),
            "mfe_5p_evidence": mfe_evidence,
            "homopolymers_detected": homopolymers,
            "rewrite_metrics": {
                "step_evaluations": step_counter[0],
                "max_step_evaluations": max_step_evaluations,
                "budget_exhausted": budget_exhausted,
                "passes_executed": passes_executed,
                "accepted_rewrites": len(trace),
                "seed_repaired": seed_repaired,
                "window_sizes": list(window_sizes),
                "include_paired_gc_swaps": include_paired_gc_swaps,
                "wall_time_ms": round(wall_time_ms, 3),
                "trace": [asdict(rec) for rec in trace],
            },
            "biological_validation": "NOT_AVAILABLE",
            "engine_id": "constrained_local_rewriter",
            "engine_version": "0.1.0",
        }


def run_exact_dp_with_budget_trace(
    dp_optimizer: DPV230Optimizer,
    protein_sequence: str,
    codon_weights: dict[str, float],
    target_gc_min: float = 0.40,
    target_gc_max: float = 0.47,
    upstream_context: str = "",
    stop_codon: str = "TAA",
    downstream_context: str = "",
) -> dict[str, Any]:
    """Run exact DP, then estimate transition counts by a separate state-graph replay."""
    rewriter = ConstrainedLocalRewriter(dp_optimizer=dp_optimizer)
    spec = rewriter._resolve_problem_spec(
        protein_sequence=protein_sequence,
        codon_weights=codon_weights,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        upstream_context=upstream_context,
        stop_codon=stop_codon,
        downstream_context=downstream_context,
    )

    t0 = time.perf_counter()
    dp_result = dp_optimizer.optimize(
        protein_sequence=protein_sequence,
        codon_weights=codon_weights,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        upstream_context=upstream_context,
        stop_codon=stop_codon,
        downstream_context=downstream_context,
        allow_nearest_infeasible=False,
    )
    wall_time_ms = (time.perf_counter() - t0) * 1000.0

    # Count exact transition evaluations performed across DP layers
    automaton = dp_optimizer.automaton
    active_states: set[tuple[int, int]] = {(0, spec.s0)}
    step_evaluations = 0
    peak_active_states = 1

    for i, aa in enumerate(spec.protein):
        syns = AA_TO_CODONS[aa]
        next_states: set[tuple[int, int]] = set()
        for g, s in active_states:
            for c in syns:
                step_evaluations += 1
                next_s, is_forbidden = automaton.step_codon(s, c)
                if not is_forbidden:
                    next_states.add((g + spec.codon_gc[c], next_s))
        if i == spec.ramp_k - 1 and spec.apply_initiation_gc_guard:
            filtered = {
                st
                for st in next_states
                if spec.init_min_gc_count <= st[0] <= spec.init_max_gc_count
            }
            if filtered:
                next_states = filtered
        active_states = next_states
        peak_active_states = max(peak_active_states, len(active_states))

    audit = rewriter.evaluate_candidate(
        protein_sequence=spec.protein,
        coding_dna=dp_result["sequence"],
        codon_weights=codon_weights,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        upstream_context=upstream_context,
        stop_codon=stop_codon,
        downstream_context=downstream_context,
    )

    return {
        "sequence": dp_result["sequence"],
        "score": dp_result["score"],
        "gc_percent": dp_result["gc_percent"],
        "cai_global": dp_result["cai_global"],
        "step_evaluations": step_evaluations,
        "step_evaluations_kind": "SEPARATE_STATE_GRAPH_REPLAY_ESTIMATE",
        "constraint_contract_matches_local": audit["is_feasible"],
        "peak_active_states": peak_active_states,
        "wall_time_ms": round(wall_time_ms, 3),
        "is_feasible": audit["is_feasible"],
        "translation_equivalent": audit["translation_equivalent"],
    }


def compare_rewriter_vs_exact_dp(
    protein_sequence: str,
    codon_weights: dict[str, float],
    dp_optimizer: DPV230Optimizer | None = None,
    initial_coding_dna: str | None = None,
    seed_mode: str = "low_fitness_feasible",
    budgets: Sequence[int | None] = (200, 500, 1500, None),
    window_sizes: Sequence[int] = (1, 2, 3),
    include_paired_gc_swaps: bool = True,
    target_gc_min: float = 0.40,
    target_gc_max: float = 0.47,
    upstream_context: str = "",
    stop_codon: str = "TAA",
    downstream_context: str = "",
) -> dict[str, Any]:
    """Compare capped local search against an uncapped DP reference, with explicit accounting limits."""
    optimizer = dp_optimizer or DPV230Optimizer()
    rewriter = ConstrainedLocalRewriter(dp_optimizer=optimizer)

    dp_summary = run_exact_dp_with_budget_trace(
        dp_optimizer=optimizer,
        protein_sequence=protein_sequence,
        codon_weights=codon_weights,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        upstream_context=upstream_context,
        stop_codon=stop_codon,
        downstream_context=downstream_context,
    )

    if not dp_summary["is_feasible"]:
        raise ValueError(
            "Exact DP reference is outside the strict local constraint contract; comparison rejected"
        )
    seed_step_evaluations = 0
    # Establish a shared initial feasible seed once so all budget tiers start from the exact same sequence
    if initial_coding_dna is not None:
        shared_seed_dna = initial_coding_dna
    else:
        spec = rewriter._resolve_problem_spec(
            protein_sequence=protein_sequence,
            codon_weights=codon_weights,
            target_gc_min=target_gc_min,
            target_gc_max=target_gc_max,
            upstream_context=upstream_context,
            stop_codon=stop_codon,
            downstream_context=downstream_context,
        )
        seed_counter = [0]
        shared_seed_codons = rewriter._construct_feasible_seed(spec, seed_mode, seed_counter, None)
        shared_seed_dna = "".join(shared_seed_codons)
        seed_step_evaluations = seed_counter[0]

    seed_eval = rewriter.evaluate_candidate(
        protein_sequence=protein_sequence,
        coding_dna=shared_seed_dna,
        codon_weights=codon_weights,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        upstream_context=upstream_context,
        stop_codon=stop_codon,
        downstream_context=downstream_context,
    )

    rewriter_runs: list[dict[str, Any]] = []
    for budget in budgets:
        res = rewriter.rewrite(
            protein_sequence=protein_sequence,
            codon_weights=codon_weights,
            initial_coding_dna=shared_seed_dna,
            window_sizes=window_sizes,
            include_paired_gc_swaps=include_paired_gc_swaps,
            max_step_evaluations=budget,
            require_feasible_seed=True,
            target_gc_min=target_gc_min,
            target_gc_max=target_gc_max,
            upstream_context=upstream_context,
            stop_codon=stop_codon,
            downstream_context=downstream_context,
        )
        trace_records = res["rewrite_metrics"]["trace"]
        monotonic = all(rec["score_after"] >= rec["score_before"] - 1e-12 for rec in trace_records)
        gap = dp_summary["score"] - res["score"]
        rewriter_runs.append(
            {
                "budget_limit": budget,
                "seed_step_evaluations": seed_step_evaluations,
                "total_step_evaluations_including_seed": seed_step_evaluations
                + res["rewrite_metrics"]["step_evaluations"],
                "step_evaluations": res["rewrite_metrics"]["step_evaluations"],
                "budget_fraction_of_dp": round(
                    res["rewrite_metrics"]["step_evaluations"]
                    / float(max(1, dp_summary["step_evaluations"])),
                    4,
                ),
                "sequence": res["sequence"],
                "score": res["score"],
                "score_improvement_over_seed": res["score_improvement"],
                "optimality_gap_vs_dp": gap,
                "matches_dp_optimum": abs(gap) <= 1e-9,
                "gc_percent": res["gc_percent"],
                "cai_global": res["cai_global"],
                "accepted_rewrites": res["rewrite_metrics"]["accepted_rewrites"],
                "passes_executed": res["rewrite_metrics"]["passes_executed"],
                "monotonic_trace_verified": monotonic,
                "is_feasible": res["is_feasible"],
                "translation_equivalent": res["translation_equivalent"],
                "wall_time_ms": res["rewrite_metrics"]["wall_time_ms"],
            }
        )

    return {
        "evidence_status": "COMPUTATIONAL_COMPARISON_ONLY",
        "biological_validation": "NOT_AVAILABLE",
        "scope_limitations": (
            "In-silico comparative computation under identical objective, configured "
            "sequence-level constraints, and local transition caps versus uncapped DP. "
            "Configured constraint satisfaction does not guarantee biological safety, "
            "expression, folding, or clinical suitability. MCTS/A*/sLLM guidance and "
            "GNN/GRN extensions remain unimplemented future comparative tasks."
        ),
        "matched_budget_comparison": False,
        "budget_accounting": {
            "exact_dp": "Uncapped reference; step count is a separate replay estimate, not measured execution instrumentation",
            "local_budget": "Codon-transition counter during rewriting; shared seed preprocessing excluded from tier limit and separately reported",
            "excluded_costs": "Setup, full-sequence audits, prefix recomputations and nucleotide-tail checks are not a common end-to-end budget",
        },
        "protein_length": len("".join(protein_sequence.upper().split()).rstrip("*")),
        "initial_seed": {
            "sequence": shared_seed_dna,
            "score": seed_eval["score"],
            "is_feasible": seed_eval["is_feasible"],
            "gc_count": seed_eval["gc_count"],
        },
        "exact_dp": dp_summary,
        "rewriter_runs": rewriter_runs,
    }


def verify_small_instance_ground_truth(
    protein_sequence: str,
    codon_weights: dict[str, float],
    dp_optimizer: DPV230Optimizer | None = None,
    target_gc_min: float = 0.40,
    target_gc_max: float = 0.47,
    upstream_context: str = "",
    stop_codon: str = "TAA",
    downstream_context: str = "",
) -> dict[str, Any]:
    """Exhaustively enumerates all synonymous sequences on a small instance (N <= 8) to verify exact ground truth."""
    optimizer = dp_optimizer or DPV230Optimizer(ramp_length_codons=3)
    rewriter = ConstrainedLocalRewriter(dp_optimizer=optimizer)
    spec = rewriter._resolve_problem_spec(
        protein_sequence=protein_sequence,
        codon_weights=codon_weights,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        upstream_context=upstream_context,
        stop_codon=stop_codon,
        downstream_context=downstream_context,
    )
    if spec.n_codons > 8:
        raise ValueError("Exhaustive ground-truth enumeration is restricted to N <= 8 codons.")

    all_choices = [spec.pos_codons[i] for i in range(spec.n_codons)]
    total_synonymous = 1
    for ch in all_choices:
        total_synonymous *= len(ch)

    feasible_count = 0
    best_key: tuple[float, float, tuple[str, ...]] | None = None
    best_dna: str | None = None
    best_score = -math.inf

    for combo in itertools.product(*all_choices):
        dna = "".join(combo)
        ev = rewriter.evaluate_candidate(
            protein_sequence=spec.protein,
            coding_dna=dna,
            codon_weights=codon_weights,
            target_gc_min=target_gc_min,
            target_gc_max=target_gc_max,
            upstream_context=upstream_context,
            stop_codon=stop_codon,
            downstream_context=downstream_context,
        )
        if not ev["is_feasible"]:
            continue
        feasible_count += 1
        sc = ev["score"]
        gc_dist = abs(ev["gc_count"] - spec.target_mid_gc_count)
        if best_key is None:
            best_key = (sc, -gc_dist, combo)
            best_dna = dna
            best_score = sc
        else:
            prev_sc, prev_neg_dist, prev_combo = best_key
            if sc > prev_sc + 1e-12:
                best_key = (sc, -gc_dist, combo)
                best_dna = dna
                best_score = sc
            elif abs(sc - prev_sc) <= 1e-12:
                if (
                    -gc_dist > prev_neg_dist + 1e-12
                    or abs(-gc_dist - prev_neg_dist) <= 1e-12
                    and combo < prev_combo
                ):
                    best_key = (sc, -gc_dist, combo)
                    best_dna = dna
                    best_score = sc

    if best_dna is None:
        raise UnsatisfiableDesignError("No feasible sequence exists in exhaustive enumeration.")

    dp_res = optimizer.optimize(
        protein_sequence=spec.protein,
        codon_weights=codon_weights,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        upstream_context=upstream_context,
        stop_codon=stop_codon,
        downstream_context=downstream_context,
    )
    rw_res = rewriter.rewrite(
        protein_sequence=spec.protein,
        codon_weights=codon_weights,
        seed_mode="low_fitness_feasible",
        window_sizes=(1, 2, 3),
        include_paired_gc_swaps=True,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        upstream_context=upstream_context,
        stop_codon=stop_codon,
        downstream_context=downstream_context,
    )

    return {
        "total_synonymous_candidates": total_synonymous,
        "feasible_candidates": feasible_count,
        "exhaustive_best_sequence": best_dna,
        "exhaustive_best_score": best_score,
        "dp_sequence": dp_res["sequence"],
        "dp_score": dp_res["score"],
        "dp_matches_exhaustive_score": abs(dp_res["score"] - best_score) <= 1e-12,
        "dp_matches_exhaustive_sequence": dp_res["sequence"] == best_dna,
        "rewriter_sequence": rw_res["sequence"],
        "rewriter_score": rw_res["score"],
        "rewriter_matches_exhaustive_score": abs(rw_res["score"] - best_score) <= 1e-12,
        "rewriter_feasible": rw_res["is_feasible"],
    }
