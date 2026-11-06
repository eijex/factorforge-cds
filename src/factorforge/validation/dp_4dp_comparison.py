"""Matched-Condition Comparison Prototype between Deterministic 1DP, 3DP, and 4DP (Job 508 / Job 509).

Evaluates 1DP (baseline ramp + host adaptation), 3DP (static boundary window D_k = [b_k - delta, b_k + delta], L=0),
and 4DP (shifted ribosome exit-tunnel latency window W_k(L, delta) = [b_k + L - delta, b_k + L + delta])
under strictly identical problem constraints, host codon tables, Aho-Corasick forbidden motif DFAs,
and lexicographic tie-breaking rules.

Scope & Methodological Boundaries (Job 509 Audit Corrections):
1. Identical Constraint Contract:
   All three modes optimize over the exact same feasible state space: amino acid translation equivalence,
   compiled Aho-Corasick forbidden motifs (Type IIS sites, reverse complements, homopolymer runs),
   5' initiation ramp GC bounds, and global coding GC bounds.
2. Cross-Objective Matrix vs. Exhaustive Ground-Truth Optimality:
   - The 3x3 cross-objective matrix evaluates the three generated candidates (CDS_1DP, CDS_3DP, CDS_4DP)
     against J_1DP, J_3DP, and J_4DP. Diagonal dominance (J_m(CDS_m) >= J_m(CDS_m') for m, m' in {1DP, 3DP, 4DP})
     confirms self-preference among the three generated candidates; by itself, a 3-candidate comparison
     does NOT prove global optimality over the entire exponential feasible space.
   - Exact global optimality of the DP recurrence and tie-breaking is separately verified on short
     instances (N <= 8 codons) via exhaustive brute-force enumeration (`verify_small_instance_exhaustive_dp_4dp`).
3. Honest Invariant Reporting (No Unexecuted True Flags):
   - When collapse invariants (lambda_k = 0, L = 0, and out-of-bounds exclusion) are skipped on longer
     fixtures to save runtime, their status fields are explicitly recorded as `None` with
     `execution_status = "NOT_RUN"`. They are NEVER marked `True` when unexecuted.
4. State-Graph Replay & Memory Measurement Scope:
   - State-space transition counts (`step_evaluations`) and `peak_active_states` can be either
     independently replayed per mode (`independent_state_replay=True`) or shared from a single
     common state-graph replay (`independent_state_replay=False`), with `state_metrics_measurement_mode`
     explicitly recorded in the output.
   - Memory instrumentation (`python_tracemalloc_peak_kb`) measures Python-traced heap allocations
     during `optimizer.optimize(...)` via `tracemalloc`, which does not capture external C/native
     process RSS (e.g., ViennaRNA C buffers) or guarantee identical OS-level resident memory.
5. Benchmark Fixture Provenance:
   - Each benchmark panel entry explicitly declares `fixture_kind` (`"SYNTHETIC_ALGORITHMIC_FIXTURE"` or
     `"PUBLIC_SEQUENCE_WITH_SYNTHETIC_BOUNDARIES"`), `sequence_provenance`, and `boundary_provenance`
     so synthetic algorithmic test cases are never conflated with experimentally validated domain boundaries.
   - Biological validation: NOT_AVAILABLE (`UNVERIFIED_PUBLIC_CORPUS`).
"""

from __future__ import annotations

import time
import tracemalloc
from collections.abc import Iterable, Sequence
from typing import Any

from factorforge.analysis.metrics import (
    calculate_cai,
    calculate_gc,
    detect_homopolymers,
    load_codon_usage_table,
)
from factorforge.engines.dp_v2_3 import AA_TO_CODONS, DPV230Optimizer
from factorforge.engines.local_rewriter import (
    ConstrainedLocalRewriter,
    build_4dp_pause_classes,
    verify_small_instance_ground_truth,
)
from factorforge.engines.profile.scoring import compute_5p_mfe_evidence
from factorforge.utils.exceptions import UnsatisfiableDesignError


def _split_codons(dna: str) -> list[str]:
    seq = "".join(dna.upper().split())
    return [seq[i : i + 3] for i in range(0, len(seq), 3)]


def _window_mean_weight(
    codons: Sequence[str],
    positions_1based: Sequence[int],
    codon_weights: dict[str, float],
) -> float | None:
    if not positions_1based:
        return None
    vals = [
        float(codon_weights.get(codons[pos - 1], 0.01))
        for pos in positions_1based
        if 1 <= pos <= len(codons)
    ]
    if not vals:
        return None
    return round(sum(vals) / len(vals), 4)


def _window_mean_abs_tau_dev(
    codons: Sequence[str],
    positions_1based: Sequence[int],
    codon_weights: dict[str, float],
    target_tau: float,
) -> float | None:
    if not positions_1based:
        return None
    vals = [
        abs(float(codon_weights.get(codons[pos - 1], 0.01)) - target_tau)
        for pos in positions_1based
        if 1 <= pos <= len(codons)
    ]
    if not vals:
        return None
    return round(sum(vals) / len(vals), 4)


def _replay_state_graph_metrics(
    optimizer: DPV230Optimizer,
    spec: Any,
) -> tuple[int, int]:
    """Counts exact DP transition evaluations and peak active states across the constraint state graph."""
    automaton = optimizer.automaton
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

    return step_evaluations, peak_active_states


def _execute_instrumented_dp(
    optimizer: DPV230Optimizer,
    protein_sequence: str,
    codon_weights: dict[str, float],
    target_gc_min: float,
    target_gc_max: float,
    upstream_context: str,
    stop_codon: str,
    downstream_context: str,
    homopolymer_max_run: int,
    precomputed_state_metrics: tuple[int, int] | None = None,
) -> tuple[dict[str, Any], ConstrainedLocalRewriter]:
    """Runs DPV230Optimizer with wall-clock, Python tracemalloc heap memory, and state-space instrumentation."""
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

    was_tracing = tracemalloc.is_tracing()
    if not was_tracing:
        tracemalloc.start()
    else:
        tracemalloc.reset_peak()

    t0 = time.perf_counter()
    dp_result = optimizer.optimize(
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
    _, peak_bytes = tracemalloc.get_traced_memory()
    if not was_tracing:
        tracemalloc.stop()
    python_tracemalloc_peak_kb = round(peak_bytes / 1024.0, 2)

    if precomputed_state_metrics is not None:
        step_evaluations, peak_active_states = precomputed_state_metrics
        state_metrics_mode = "SHARED_COMMON_GRAPH_REUSE"
    else:
        step_evaluations, peak_active_states = _replay_state_graph_metrics(optimizer, spec)
        state_metrics_mode = "INDEPENDENT_STATE_GRAPH_REPLAY"

    dna = dp_result["sequence"]
    audit = rewriter.evaluate_candidate(
        protein_sequence=spec.protein,
        coding_dna=dna,
        codon_weights=codon_weights,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        upstream_context=upstream_context,
        stop_codon=stop_codon,
        downstream_context=downstream_context,
    )
    full_construct = spec.up_ctx + dna + spec.s_tail
    homopolymers = detect_homopolymers(full_construct, max_run=homopolymer_max_run + 1)
    mfe_5p = compute_5p_mfe_evidence(dna, upstream_context=upstream_context)

    summary = {
        "sequence": dna,
        "score": round(float(dp_result["score"]), 6),
        "raw_score": float(dp_result["score"]),
        "cai_global": round(calculate_cai(dna, codon_weights), 4),
        "cai_ramp": dp_result["cai_ramp_45nt"],
        "cai_body": dp_result["cai_body"],
        "gc_percent": round(calculate_gc(dna), 2),
        "gc_count": audit["gc_count"],
        "init_gc_count": audit["init_gc_count"],
        "is_feasible": audit["is_feasible"] and len(homopolymers) == 0,
        "translation_equivalent": audit["translation_equivalent"],
        "automaton_valid": audit["automaton_valid"],
        "initiation_gc_feasible": audit["initiation_gc_feasible"],
        "global_gc_feasible": audit["global_gc_feasible"],
        "homopolymers_detected": len(homopolymers),
        "violations": audit["violations"],
        "wall_time_ms": round(wall_time_ms, 3),
        "python_tracemalloc_peak_kb": python_tracemalloc_peak_kb,
        "peak_memory_kb": python_tracemalloc_peak_kb,
        "memory_measurement_scope": "PYTHON_TRACEMALLOC_HEAP_ONLY_EXCLUDES_NATIVE_RSS",
        "step_evaluations": step_evaluations,
        "peak_active_states": peak_active_states,
        "state_metrics_measurement_mode": state_metrics_mode,
        "mfe_5p_evidence": mfe_5p,
    }
    return summary, rewriter


def _compare_candidate_pair(
    dna_a: str,
    dna_b: str,
    positions_3dp: Sequence[int],
    positions_4dp: Sequence[int],
) -> dict[str, Any]:
    """Computes codon/nucleotide divergence and spatial localization between two synonymous CDS candidates."""
    codons_a = _split_codons(dna_a)
    codons_b = _split_codons(dna_b)
    if len(codons_a) != len(codons_b):
        raise ValueError("Candidate sequences must have identical codon length")

    set_3dp = set(positions_3dp)
    set_4dp = set(positions_4dp)
    union_windows = set_3dp | set_4dp

    edited_positions_1based = [
        i + 1 for i, (ca, cb) in enumerate(zip(codons_a, codons_b)) if ca != cb
    ]
    nt_hamming = sum(1 for na, nb in zip(dna_a, dna_b) if na != nb)

    edits_in_3dp = [pos for pos in edited_positions_1based if pos in set_3dp]
    edits_in_4dp = [pos for pos in edited_positions_1based if pos in set_4dp]
    edits_outside_both = [pos for pos in edited_positions_1based if pos not in union_windows]

    total_codon_edits = len(edited_positions_1based)
    return {
        "codon_hamming_distance": total_codon_edits,
        "codon_divergence_fraction": round(total_codon_edits / max(1, len(codons_a)), 4),
        "nucleotide_hamming_distance": nt_hamming,
        "edited_positions_1based": edited_positions_1based,
        "edits_in_3dp_window": len(edits_in_3dp),
        "edits_in_4dp_window": len(edits_in_4dp),
        "compensatory_edits_outside_windows": len(edits_outside_both),
        "edited_3dp_positions": edits_in_3dp,
        "edited_4dp_positions": edits_in_4dp,
        "edited_outside_positions": edits_outside_both,
    }


def run_matched_dp_comparison(
    protein_sequence: str,
    boundaries: Iterable[int],
    codon_weights: dict[str, float] | None = None,
    tunnel_length: int = 35,
    half_width: int = 2,
    weight: float = 2.0,
    target_tau: float = 0.45,
    target_gc_min: float = 0.40,
    target_gc_max: float = 0.47,
    upstream_context: str = "",
    stop_codon: str = "TAA",
    downstream_context: str = "",
    ramp_length_codons: int = 15,
    homopolymer_max_run: int = 5,
    forbidden_enzymes: set[str] | None = None,
    forbidden_motifs: Iterable[str] | None = None,
    independent_state_replay: bool = True,
) -> dict[str, Any]:
    """Executes a matched-condition comparison across 1DP, 3DP (L=0), and 4DP (L=tunnel_length)."""
    protein = "".join(protein_sequence.upper().split()).rstrip("*")
    if not protein:
        raise ValueError("protein_sequence must not be empty")
    n_codons = len(protein)
    boundary_list = sorted({int(b) for b in boundaries})

    weights = (
        codon_weights if codon_weights is not None else load_codon_usage_table().codon_weights
    )

    pause_3dp = build_4dp_pause_classes(
        boundaries=boundary_list,
        tunnel_length=0,
        half_width=half_width,
        weight=weight,
        target_tau=target_tau,
        protein_length=n_codons,
        class_name="STATIC_3DP",
    )
    pause_4dp = build_4dp_pause_classes(
        boundaries=boundary_list,
        tunnel_length=tunnel_length,
        half_width=half_width,
        weight=weight,
        target_tau=target_tau,
        protein_length=n_codons,
        class_name="SHIFTED_4DP",
    )

    positions_3dp: list[int] = pause_3dp["STATIC_3DP"]["positions"]
    positions_4dp: list[int] = pause_4dp["SHIFTED_4DP"]["positions"]

    shared_kwargs: dict[str, Any] = {
        "forbidden_enzymes": forbidden_enzymes,
        "forbidden_motifs": forbidden_motifs,
        "ramp_length_codons": ramp_length_codons,
        "homopolymer_max_run": homopolymer_max_run,
        "auto_annotate_pauses": False,
    }

    opt_1dp = DPV230Optimizer(pause_classes=None, **shared_kwargs)
    opt_3dp = DPV230Optimizer(pause_classes=pause_3dp, **shared_kwargs)
    opt_4dp = DPV230Optimizer(pause_classes=pause_4dp, **shared_kwargs)

    res_1dp, rw_1dp = _execute_instrumented_dp(
        opt_1dp,
        protein,
        weights,
        target_gc_min,
        target_gc_max,
        upstream_context,
        stop_codon,
        downstream_context,
        homopolymer_max_run,
        precomputed_state_metrics=None,
    )
    precomputed_metrics = (
        None
        if independent_state_replay
        else (res_1dp["step_evaluations"], res_1dp["peak_active_states"])
    )
    res_3dp, rw_3dp = _execute_instrumented_dp(
        opt_3dp,
        protein,
        weights,
        target_gc_min,
        target_gc_max,
        upstream_context,
        stop_codon,
        downstream_context,
        homopolymer_max_run,
        precomputed_state_metrics=precomputed_metrics,
    )
    res_4dp, rw_4dp = _execute_instrumented_dp(
        opt_4dp,
        protein,
        weights,
        target_gc_min,
        target_gc_max,
        upstream_context,
        stop_codon,
        downstream_context,
        homopolymer_max_run,
        precomputed_state_metrics=precomputed_metrics,
    )

    mode_results = {"1dp": res_1dp, "3dp": res_3dp, "4dp": res_4dp}
    rewriters = {"1dp": rw_1dp, "3dp": rw_3dp, "4dp": rw_4dp}

    for res in mode_results.values():
        codons = _split_codons(res["sequence"])
        res["mean_w_in_3dp_window"] = _window_mean_weight(codons, positions_3dp, weights)
        res["mean_w_in_4dp_window"] = _window_mean_weight(codons, positions_4dp, weights)
        res["mean_abs_tau_dev_in_3dp_window"] = _window_mean_abs_tau_dev(
            codons, positions_3dp, weights, target_tau
        )
        res["mean_abs_tau_dev_in_4dp_window"] = _window_mean_abs_tau_dev(
            codons, positions_4dp, weights, target_tau
        )

    # 3x3 Cross-objective score matrix: evaluates the 3 generated candidates against each objective
    cross_objective_matrix: dict[str, dict[str, float]] = {}
    candidate_self_preference_checks: dict[str, bool] = {}

    for obj_mode, rw in rewriters.items():
        row: dict[str, float] = {}
        for cand_mode, cand_res in mode_results.items():
            ev = rw.evaluate_candidate(
                protein_sequence=protein,
                coding_dna=cand_res["sequence"],
                codon_weights=weights,
                target_gc_min=target_gc_min,
                target_gc_max=target_gc_max,
                upstream_context=upstream_context,
                stop_codon=stop_codon,
                downstream_context=downstream_context,
            )
            row[cand_mode] = round(float(ev["score"]), 6)
        cross_objective_matrix[obj_mode] = row
        own_score = row[obj_mode]
        candidate_self_preference_checks[obj_mode] = all(
            own_score >= other_score - 1e-6 for other_score in row.values()
        )

    pairwise_diffs = {
        "1dp_vs_3dp": _compare_candidate_pair(
            res_1dp["sequence"], res_3dp["sequence"], positions_3dp, positions_4dp
        ),
        "1dp_vs_4dp": _compare_candidate_pair(
            res_1dp["sequence"], res_4dp["sequence"], positions_3dp, positions_4dp
        ),
        "3dp_vs_4dp": _compare_candidate_pair(
            res_3dp["sequence"], res_4dp["sequence"], positions_3dp, positions_4dp
        ),
    }

    all_constraints_satisfied = all(r["is_feasible"] for r in mode_results.values())
    all_self_preferred = all(candidate_self_preference_checks.values())

    state_counts_equal = (
        res_1dp["step_evaluations"]
        == res_3dp["step_evaluations"]
        == res_4dp["step_evaluations"]
        and res_1dp["peak_active_states"]
        == res_3dp["peak_active_states"]
        == res_4dp["peak_active_states"]
    )

    return {
        "evidence_classification": "UNVERIFIED_PUBLIC_CORPUS",
        "biological_validation": "NOT_AVAILABLE",
        "epistemic_note": (
            f"Matched-condition algorithmic comparison prototype. Tunnel length L={tunnel_length} is a "
            "configured model assumption per Roadmap v1.7.0 Section 2, not a universally validated "
            "biological constant. Cross-objective matrix diagonal dominance confirms relative self-preference "
            "among the three generated candidates (CDS_1DP, CDS_3DP, CDS_4DP); exact global optimality of the "
            "DP recurrence is separately verified on short exhaustive instances (N<=8). Python tracemalloc "
            "measures Python heap allocations only and does not capture process-wide native RSS."
        ),
        "problem_setup": {
            "protein_length_aa": n_codons,
            "coding_length_nt": n_codons * 3,
            "boundaries_1based": boundary_list,
            "tunnel_length_assumption": tunnel_length,
            "half_width": half_width,
            "pause_weight_lambda": weight,
            "target_tau": target_tau,
            "target_gc_min": target_gc_min,
            "target_gc_max": target_gc_max,
            "ramp_length_codons": ramp_length_codons,
            "static_3dp_window_positions": positions_3dp,
            "shifted_4dp_window_positions": positions_4dp,
            "independent_state_replay": independent_state_replay,
        },
        "all_constraints_satisfied": all_constraints_satisfied,
        "candidate_self_preference_verified": all_self_preferred,
        "candidate_self_preference_checks": candidate_self_preference_checks,
        "cross_objective_matrix_scope": (
            "THREE_CANDIDATE_SELF_PREFERENCE_ONLY_NOT_EXHAUSTIVE_GLOBAL_OPTIMALITY_PROOF"
        ),
        "state_graph_parity_verified": state_counts_equal,
        "state_metrics_measurement_mode": (
            "INDEPENDENT_STATE_GRAPH_REPLAY"
            if independent_state_replay
            else "SHARED_COMMON_GRAPH_REUSE"
        ),
        "cross_objective_matrix": cross_objective_matrix,
        "modes": mode_results,
        "pairwise_candidate_divergence": pairwise_diffs,
    }


def verify_small_instance_exhaustive_dp_4dp(
    protein_sequence: str = "MVSKGEEL",
    boundaries: Iterable[int] = (2,),
    codon_weights: dict[str, float] | None = None,
    tunnel_length: int = 4,
    half_width: int = 1,
    weight: float = 2.5,
    target_tau: float = 0.45,
    target_gc_min: float = 0.35,
    target_gc_max: float = 0.55,
    ramp_length_codons: int = 3,
    homopolymer_max_run: int = 5,
    forbidden_enzymes: set[str] | None = None,
    forbidden_motifs: Iterable[str] | None = None,
    upstream_context: str = "",
    stop_codon: str = "TAA",
    downstream_context: str = "",
) -> dict[str, Any]:
    """Exhaustively enumerates all synonymous sequences on a short instance (N <= 8) for 1DP, 3DP, and 4DP.

    Provides an exact ground-truth check on the tested instance that DPV230Optimizer's state transitions,
    score accumulation, tie-breaking, and unsatisfiable constraint rejection match brute-force enumeration.
    """
    protein = "".join(protein_sequence.upper().split()).rstrip("*")
    n_codons = len(protein)
    if n_codons > 8:
        raise ValueError("Exhaustive ground-truth enumeration is restricted to N <= 8 codons.")

    weights = (
        codon_weights if codon_weights is not None else load_codon_usage_table().codon_weights
    )
    boundary_list = sorted({int(b) for b in boundaries})

    pause_3dp = build_4dp_pause_classes(
        boundaries=boundary_list,
        tunnel_length=0,
        half_width=half_width,
        weight=weight,
        target_tau=target_tau,
        protein_length=n_codons,
        class_name="STATIC_3DP",
    )
    pause_4dp = build_4dp_pause_classes(
        boundaries=boundary_list,
        tunnel_length=tunnel_length,
        half_width=half_width,
        weight=weight,
        target_tau=target_tau,
        protein_length=n_codons,
        class_name="SHIFTED_4DP",
    )

    shared_opt_kwargs: dict[str, Any] = {
        "ramp_length_codons": ramp_length_codons,
        "homopolymer_max_run": homopolymer_max_run,
        "forbidden_enzymes": forbidden_enzymes,
        "forbidden_motifs": forbidden_motifs,
        "auto_annotate_pauses": False,
    }

    opts = {
        "1dp": DPV230Optimizer(pause_classes=None, **shared_opt_kwargs),
        "3dp": DPV230Optimizer(pause_classes=pause_3dp, **shared_opt_kwargs),
        "4dp": DPV230Optimizer(pause_classes=pause_4dp, **shared_opt_kwargs),
    }

    mode_ground_truth: dict[str, Any] = {}
    for mode_key, opt in opts.items():
        gt = verify_small_instance_ground_truth(
            protein_sequence=protein,
            codon_weights=weights,
            dp_optimizer=opt,
            target_gc_min=target_gc_min,
            target_gc_max=target_gc_max,
            upstream_context=upstream_context,
            stop_codon=stop_codon,
            downstream_context=downstream_context,
        )
        mode_ground_truth[mode_key] = gt

    all_match_score = all(
        gt["dp_matches_exhaustive_score"] for gt in mode_ground_truth.values()
    )
    all_match_seq = all(
        gt["dp_matches_exhaustive_sequence"] for gt in mode_ground_truth.values()
    )

    return {
        "execution_status": "EXECUTED",
        "protein_sequence": protein,
        "protein_length_aa": n_codons,
        "boundaries_1based": boundary_list,
        "tunnel_length": tunnel_length,
        "half_width": half_width,
        "static_3dp_window_positions": pause_3dp["STATIC_3DP"]["positions"],
        "shifted_4dp_window_positions": pause_4dp["SHIFTED_4DP"]["positions"],
        "all_modes_match_exhaustive_global_optimum_score": all_match_score,
        "all_modes_match_exhaustive_global_optimum_sequence": all_match_seq,
        "modes": mode_ground_truth,
    }


def verify_small_instance_unsatisfiable_dp_4dp(
    protein_sequence: str,
    boundaries: Iterable[int],
    codon_weights: dict[str, float] | None = None,
    tunnel_length: int = 3,
    half_width: int = 1,
    weight: float = 2.5,
    target_tau: float = 0.45,
    target_gc_min: float = 0.40,
    target_gc_max: float = 0.47,
    ramp_length_codons: int = 3,
    homopolymer_max_run: int = 5,
    forbidden_enzymes: set[str] | None = None,
    forbidden_motifs: Iterable[str] | None = None,
    upstream_context: str = "",
    stop_codon: str = "TAA",
    downstream_context: str = "",
) -> dict[str, Any]:
    """Verifies on an infeasible short instance (N <= 8) that brute-force finds 0 feasible sequences and all DP modes raise UnsatisfiableDesignError."""
    protein = "".join(protein_sequence.upper().split()).rstrip("*")
    n_codons = len(protein)
    if n_codons > 8:
        raise ValueError("Exhaustive unsatisfiability check is restricted to N <= 8 codons.")

    weights = (
        codon_weights if codon_weights is not None else load_codon_usage_table().codon_weights
    )
    boundary_list = sorted({int(b) for b in boundaries})

    pause_3dp = build_4dp_pause_classes(
        boundaries=boundary_list,
        tunnel_length=0,
        half_width=half_width,
        weight=weight,
        target_tau=target_tau,
        protein_length=n_codons,
        class_name="STATIC_3DP",
    )
    pause_4dp = build_4dp_pause_classes(
        boundaries=boundary_list,
        tunnel_length=tunnel_length,
        half_width=half_width,
        weight=weight,
        target_tau=target_tau,
        protein_length=n_codons,
        class_name="SHIFTED_4DP",
    )

    shared_opt_kwargs: dict[str, Any] = {
        "ramp_length_codons": ramp_length_codons,
        "homopolymer_max_run": homopolymer_max_run,
        "forbidden_enzymes": forbidden_enzymes,
        "forbidden_motifs": forbidden_motifs,
        "auto_annotate_pauses": False,
    }

    opts = {
        "1dp": DPV230Optimizer(pause_classes=None, **shared_opt_kwargs),
        "3dp": DPV230Optimizer(pause_classes=pause_3dp, **shared_opt_kwargs),
        "4dp": DPV230Optimizer(pause_classes=pause_4dp, **shared_opt_kwargs),
    }

    # Check brute-force enumeration raises UnsatisfiableDesignError
    exhaustive_raised = False
    try:
        verify_small_instance_ground_truth(
            protein_sequence=protein,
            codon_weights=weights,
            dp_optimizer=opts["1dp"],
            target_gc_min=target_gc_min,
            target_gc_max=target_gc_max,
            upstream_context=upstream_context,
            stop_codon=stop_codon,
            downstream_context=downstream_context,
        )
    except UnsatisfiableDesignError:
        exhaustive_raised = True

    dp_modes_raised: dict[str, bool] = {}
    for mode_key, opt in opts.items():
        raised = False
        try:
            opt.optimize(
                protein_sequence=protein,
                codon_weights=weights,
                target_gc_min=target_gc_min,
                target_gc_max=target_gc_max,
                upstream_context=upstream_context,
                stop_codon=stop_codon,
                downstream_context=downstream_context,
                allow_nearest_infeasible=False,
            )
        except UnsatisfiableDesignError:
            raised = True
        dp_modes_raised[mode_key] = raised

    all_consistent = exhaustive_raised and all(dp_modes_raised.values())
    return {
        "execution_status": "EXECUTED",
        "protein_sequence": protein,
        "protein_length_aa": n_codons,
        "exhaustive_zero_feasible_confirmed": exhaustive_raised,
        "dp_modes_raised_unsatisfiable": dp_modes_raised,
        "unsatisfiable_rejection_consistent": all_consistent,
    }


def run_step1_exhaustive_verification_suite(
    codon_weights: dict[str, float] | None = None,
) -> dict[str, Any]:
    """STEP 1 Suite: Exhaustively verifies 1DP, 3DP, and 4DP across diverse short inputs, overlapping windows, ties, and infeasible constraints."""
    weights = (
        codon_weights if codon_weights is not None else load_codon_usage_table().codon_weights
    )
    uniform_weights = {codon: 0.75 for codon in weights}

    feasible_cases: list[dict[str, Any]] = [
        {
            "case_id": "disjoint_windows_8aa",
            "description": "Disjoint 3DP [1..3] and 4DP [5..7] windows on 8-aa mixed sequence (4,608 synonymous candidates)",
            "kwargs": {
                "protein_sequence": "MVSKGEEL",
                "boundaries": [2],
                "codon_weights": weights,
                "tunnel_length": 4,
                "half_width": 1,
                "weight": 2.5,
                "target_tau": 0.45,
                "target_gc_min": 0.35,
                "target_gc_max": 0.55,
            },
        },
        {
            "case_id": "overlapping_dual_boundaries_8aa",
            "description": "Overlapping dual boundaries b=[2, 3] with delta=1 => merged 3DP [1..4] and merged 4DP [4..7] on 8-aa sequence",
            "kwargs": {
                "protein_sequence": "MKLSTVAG",
                "boundaries": [2, 3],
                "codon_weights": weights,
                "tunnel_length": 3,
                "half_width": 1,
                "weight": 3.0,
                "target_tau": 0.40,
                "target_gc_min": 0.35,
                "target_gc_max": 0.60,
            },
        },
        {
            "case_id": "high_degeneracy_6box_7aa",
            "description": "High-degeneracy 6-box (L, R, S) + 4-box (P, V) sequence (candidate count is enumerated in the result)",
            "kwargs": {
                "protein_sequence": "MLRSPVG",
                "boundaries": [2],
                "codon_weights": weights,
                "tunnel_length": 3,
                "half_width": 1,
                "weight": 2.5,
                "target_tau": 0.45,
                "target_gc_min": 0.40,
                "target_gc_max": 0.65,
            },
        },
        {
            "case_id": "exact_score_tie_breaking_7aa",
            "description": "Uniform codon weights (all w=0.75) inducing massive exact score ties resolved by GC-midpoint distance and LexRank",
            "kwargs": {
                "protein_sequence": "MVSKGEE",
                "boundaries": [2],
                "codon_weights": uniform_weights,
                "tunnel_length": 3,
                "half_width": 1,
                "weight": 2.0,
                "target_tau": 0.75,
                "target_gc_min": 0.35,
                "target_gc_max": 0.55,
            },
        },
        {
            "case_id": "flanking_context_and_motif_avoidance_7aa",
            "description": "Upstream/downstream flanking contexts with forbidden motif avoidance and C-terminal partial window clipping",
            "kwargs": {
                "protein_sequence": "MGSKELV",
                "boundaries": [4],
                "codon_weights": weights,
                "tunnel_length": 3,
                "half_width": 1,
                "weight": 2.5,
                "target_tau": 0.45,
                "target_gc_min": 0.35,
                "target_gc_max": 0.60,
                "upstream_context": "CACC",
                "stop_codon": "TAA",
                "downstream_context": "GCTT",
                "forbidden_motifs": ["GGATCC", "GAATTC"],
            },
        },
    ]

    infeasible_cases: list[dict[str, Any]] = [
        {
            "case_id": "impossible_high_gc_band_on_at_rich_aa",
            "description": "AT-locked amino acids (M-I-F-Y-N-K) cannot reach [80%, 95%] GC band",
            "kwargs": {
                "protein_sequence": "MIFYNK",
                "boundaries": [2],
                "codon_weights": weights,
                "tunnel_length": 2,
                "half_width": 1,
                "target_gc_min": 0.80,
                "target_gc_max": 0.95,
            },
        },
        {
            "case_id": "impossible_homopolymer_on_poly_phe",
            "description": "MFFFF with homopolymer_max_run=2 and forbidden TTC motif leaves 0 valid paths",
            "kwargs": {
                "protein_sequence": "MFFFF",
                "boundaries": [2],
                "codon_weights": weights,
                "tunnel_length": 2,
                "half_width": 1,
                "target_gc_min": 0.10,
                "target_gc_max": 0.90,
                "homopolymer_max_run": 2,
                "forbidden_motifs": ["TTC"],
            },
        },
        {
            "case_id": "impossible_stop_downstream_junction_motif",
            "description": "Stop codon TAA + downstream AAAA triggers homopolymer >= 6 nt violation on every candidate",
            "kwargs": {
                "protein_sequence": "MVSKGE",
                "boundaries": [2],
                "codon_weights": weights,
                "tunnel_length": 2,
                "half_width": 1,
                "target_gc_min": 0.30,
                "target_gc_max": 0.60,
                "stop_codon": "TAA",
                "downstream_context": "AAAA",
                "homopolymer_max_run": 5,
            },
        },
    ]

    feasible_results: list[dict[str, Any]] = []
    for item in feasible_cases:
        res = verify_small_instance_exhaustive_dp_4dp(**item["kwargs"])
        feasible_results.append(
            {
                "case_id": item["case_id"],
                "description": item["description"],
                "result": res,
            }
        )

    infeasible_results: list[dict[str, Any]] = []
    for item in infeasible_cases:
        res = verify_small_instance_unsatisfiable_dp_4dp(**item["kwargs"])
        infeasible_results.append(
            {
                "case_id": item["case_id"],
                "description": item["description"],
                "result": res,
            }
        )

    all_feasible_exact = all(
        r["result"]["all_modes_match_exhaustive_global_optimum_score"]
        and r["result"]["all_modes_match_exhaustive_global_optimum_sequence"]
        for r in feasible_results
    )
    all_infeasible_exact = all(
        r["result"]["unsatisfiable_rejection_consistent"] for r in infeasible_results
    )

    return {
        "step": "STEP_1_EXHAUSTIVE_VERIFICATION",
        "scope_note": (
            "Confirms that 1DP, 3DP, and 4DP match brute-force exhaustive enumeration (score, sequence, "
            "tie-breaking, and unsatisfiable rejection) on the tested short instances (N <= 8). "
            "Passing these short-instance tests does not by itself prove correctness on arbitrary unexamined inputs "
            "or establish wet-lab biological efficacy."
        ),
        "feasible_case_count": len(feasible_results),
        "infeasible_case_count": len(infeasible_results),
        "all_feasible_cases_match_exhaustive_optimum": all_feasible_exact,
        "all_infeasible_cases_consistently_rejected": all_infeasible_exact,
        "step1_gate_passed": all_feasible_exact and all_infeasible_exact,
        "feasible_cases": feasible_results,
        "infeasible_cases": infeasible_results,
    }


def verify_dp_4dp_collapse_invariants(
    protein_sequence: str,
    boundaries: Iterable[int],
    codon_weights: dict[str, float] | None = None,
    tunnel_length: int = 35,
    half_width: int = 2,
    target_tau: float = 0.45,
    target_gc_min: float = 0.40,
    target_gc_max: float = 0.47,
) -> dict[str, Any]:
    """Verifies mathematical collapse and boundary exclusion invariants across 1DP, 3DP, and 4DP."""
    protein = "".join(protein_sequence.upper().split()).rstrip("*")
    n_codons = len(protein)
    weights = (
        codon_weights if codon_weights is not None else load_codon_usage_table().codon_weights
    )
    boundary_list = sorted({int(b) for b in boundaries})

    # 1. Invariant A: lambda = 0 => 1DP == 3DP == 4DP
    comp_lambda_zero = run_matched_dp_comparison(
        protein_sequence=protein,
        boundaries=boundary_list,
        codon_weights=weights,
        tunnel_length=tunnel_length,
        half_width=half_width,
        weight=0.0,
        target_tau=target_tau,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        independent_state_replay=False,
    )
    seq_1dp_l0 = comp_lambda_zero["modes"]["1dp"]["sequence"]
    seq_3dp_l0 = comp_lambda_zero["modes"]["3dp"]["sequence"]
    seq_4dp_l0 = comp_lambda_zero["modes"]["4dp"]["sequence"]
    sc_1dp_l0 = comp_lambda_zero["modes"]["1dp"]["raw_score"]
    sc_3dp_l0 = comp_lambda_zero["modes"]["3dp"]["raw_score"]
    sc_4dp_l0 = comp_lambda_zero["modes"]["4dp"]["raw_score"]

    lambda_zero_collapse_ok = (
        seq_1dp_l0 == seq_3dp_l0 == seq_4dp_l0
        and abs(sc_1dp_l0 - sc_3dp_l0) <= 1e-12
        and abs(sc_1dp_l0 - sc_4dp_l0) <= 1e-12
    )

    # 2. Invariant B: L = 0 => 4DP == 3DP (for positive lambda > 0)
    comp_tunnel_zero = run_matched_dp_comparison(
        protein_sequence=protein,
        boundaries=boundary_list,
        codon_weights=weights,
        tunnel_length=0,
        half_width=half_width,
        weight=2.0,
        target_tau=target_tau,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        independent_state_replay=False,
    )
    seq_3dp_t0 = comp_tunnel_zero["modes"]["3dp"]["sequence"]
    seq_4dp_t0 = comp_tunnel_zero["modes"]["4dp"]["sequence"]
    sc_3dp_t0 = comp_tunnel_zero["modes"]["3dp"]["raw_score"]
    sc_4dp_t0 = comp_tunnel_zero["modes"]["4dp"]["raw_score"]

    tunnel_zero_collapse_ok = (
        seq_3dp_t0 == seq_4dp_t0 and abs(sc_3dp_t0 - sc_4dp_t0) <= 1e-12
    )

    # 3. Invariant C: Out-of-bounds shifted window (b_k + L - delta > N) => empty 4DP window & collapse to 1DP
    oob_boundary = max(1, n_codons - max(1, tunnel_length // 2))
    if oob_boundary + tunnel_length - half_width <= n_codons:
        oob_boundary = n_codons
    comp_oob = run_matched_dp_comparison(
        protein_sequence=protein,
        boundaries=[oob_boundary],
        codon_weights=weights,
        tunnel_length=tunnel_length,
        half_width=half_width,
        weight=2.0,
        target_tau=target_tau,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        independent_state_replay=False,
    )
    oob_4dp_positions = comp_oob["problem_setup"]["shifted_4dp_window_positions"]
    oob_exclusion_ok = (
        len(oob_4dp_positions) == 0
        and comp_oob["modes"]["4dp"]["sequence"] == comp_oob["modes"]["1dp"]["sequence"]
        and abs(comp_oob["modes"]["4dp"]["raw_score"] - comp_oob["modes"]["1dp"]["raw_score"])
        <= 1e-12
    )

    return {
        "execution_status": "EXECUTED",
        "lambda_zero_collapse_verified": lambda_zero_collapse_ok,
        "tunnel_zero_collapse_verified": tunnel_zero_collapse_ok,
        "out_of_bounds_exclusion_verified": oob_exclusion_ok,
        "all_invariants_verified": (
            lambda_zero_collapse_ok and tunnel_zero_collapse_ok and oob_exclusion_ok
        ),
        "details": {
            "lambda_zero_scores": {
                "1dp": round(sc_1dp_l0, 6),
                "3dp": round(sc_3dp_l0, 6),
                "4dp": round(sc_4dp_l0, 6),
            },
            "tunnel_zero_scores": {
                "3dp": round(sc_3dp_t0, 6),
                "4dp": round(sc_4dp_t0, 6),
            },
            "oob_boundary_tested": oob_boundary,
            "oob_shifted_4dp_positions": oob_4dp_positions,
        },
    }


def make_not_run_invariants_record(reason: str) -> dict[str, Any]:
    """Returns an explicit NOT_RUN invariant record with None fields rather than unexecuted True flags."""
    return {
        "execution_status": "NOT_RUN",
        "lambda_zero_collapse_verified": None,
        "tunnel_zero_collapse_verified": None,
        "out_of_bounds_exclusion_verified": None,
        "all_invariants_verified": None,
        "reason": reason,
    }


BENCHMARK_PROTEIN_PANEL: tuple[dict[str, Any], ...] = (
    {
        "id": "synthetic_vhh_65aa",
        "name": "Synthetic VHH-like Algorithmic Test Subsequence (65 aa)",
        "fixture_kind": "SYNTHETIC_ALGORITHMIC_FIXTURE",
        "sequence_provenance": (
            "Truncated 65 aa N-terminal VHH framework subsequence with synthetic ATG prefix "
            "for fast algorithmic window-shift testing."
        ),
        "boundary_provenance": (
            "Arbitrary test index b_1=22 placed near CDR1/FR2 for disjoint D_1=[20..24] vs W_1=[55..59] verification; "
            "not an experimentally validated folding pause site."
        ),
        "protein": "ATGGVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYAD",
        "boundaries": [22],
        "tunnel_length": 35,
        "half_width": 2,
    },
    {
        "id": "two_domain_linker_95aa",
        "name": "Truncated GFP N-Terminal Subsequence with Synthetic Dual Boundaries (95 aa)",
        "fixture_kind": "SYNTHETIC_ALGORITHMIC_FIXTURE",
        "sequence_provenance": (
            "95 aa N-terminal subsequence of GFP used as a medium-length algorithmic fixture."
        ),
        "boundary_provenance": (
            "Synthetic dual test positions b=[30, 50] chosen to test two non-overlapping windows "
            "under lower-bound GC constraint coupling; not true domain boundaries."
        ),
        "protein": (
            "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPT"
            "LVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQ"
        ),
        "boundaries": [30, 50],
        "tunnel_length": 35,
        "half_width": 2,
    },
    {
        "id": "sfgfp_238aa",
        "name": "Superfolder GFP Sequence with Configured Test Boundaries (238 aa)",
        "fixture_kind": "PUBLIC_SEQUENCE_WITH_SYNTHETIC_BOUNDARIES",
        "sequence_provenance": "Canonical 238 aa Superfolder GFP (sfGFP) amino acid sequence.",
        "boundary_provenance": (
            "Configured test positions b=[70, 145] for full-length beta-barrel algorithmic stress testing; "
            "not experimentally measured ribosome pausing or domain boundary coordinates."
        ),
        "protein": (
            "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTL"
            "VTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLV"
            "NRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLAD"
            "HYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        ),
        "boundaries": [70, 145],
        "tunnel_length": 35,
        "half_width": 2,
    },
    {
        "id": "hemopexin_subsequence_206aa",
        "name": "Human Hemopexin N-Terminal 206 aa Subsequence with Configured Test Boundaries",
        "fixture_kind": "PUBLIC_SEQUENCE_WITH_SYNTHETIC_BOUNDARIES",
        "sequence_provenance": (
            "206 aa N-terminal subsequence of mature human Hemopexin (UniProt P02790 residues 24-229); "
            "truncated fixture, not the full 439 aa two-propeller mature protein."
        ),
        "boundary_provenance": (
            "Configured test positions b=[65, 135] within blade repeats for algorithmic comparison; "
            "not the inter-propeller hinge of full-length Hemopexin."
        ),
        "protein": (
            "TPLPPTSAHGNVAEGETKPDPDVTERCSDGWSFDATTLDDNGTMLFFKGEFVWKSHKWDR"
            "ELISERWKNFPSPVDAAFRQGHNSVFLIKGDKVWVYPPEKKEKGYPKLLQDEFPGIPSPL"
            "DAAVECHRGECQAEGVLFFQGDREWFWDLATGTMKERSWPAVGNCSSALRWLGRYYCFQG"
            "NQFLRFDPVRGEVPPRYPRDVRDYFM"
        ),
        "boundaries": [65, 135],
        "tunnel_length": 35,
        "half_width": 2,
    },
    {
        "id": "cterm_boundary_clipping_55aa",
        "name": "C-Terminal Out-of-Bounds Boundary Fixture (55 aa)",
        "fixture_kind": "SYNTHETIC_ALGORITHMIC_FIXTURE",
        "sequence_provenance": "55 aa truncated N-terminal GFP subsequence for C-terminal clipping verification.",
        "boundary_provenance": (
            "Synthetic boundary b_1=30 with L=35, delta=2 so shifted window [63..67] > 55 tests "
            "strict out-of-bounds exclusion (empty 4DP window without C-terminal clamping)."
        ),
        "protein": "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLP",
        "boundaries": [30],
        "tunnel_length": 35,
        "half_width": 2,
    },
)


def run_dp_4dp_benchmark_panel(
    codon_weights: dict[str, float] | None = None,
    panel: Sequence[dict[str, Any]] = BENCHMARK_PROTEIN_PANEL,
    weight: float = 2.5,
    target_tau: float = 0.45,
    verify_invariants_for_all: bool = False,
    independent_state_replay: bool = True,
) -> dict[str, Any]:
    """Runs the matched-condition DP / 3DP / 4DP comparison across the benchmark panel."""
    weights = (
        codon_weights if codon_weights is not None else load_codon_usage_table().codon_weights
    )

    exhaustive_gt = verify_small_instance_exhaustive_dp_4dp(codon_weights=weights)

    entries: list[dict[str, Any]] = []
    invariants_executed_count = 0
    invariants_not_run_count = 0

    for idx, item in enumerate(panel):
        comp = run_matched_dp_comparison(
            protein_sequence=item["protein"],
            boundaries=item["boundaries"],
            codon_weights=weights,
            tunnel_length=item.get("tunnel_length", 35),
            half_width=item.get("half_width", 2),
            weight=weight,
            target_tau=target_tau,
            independent_state_replay=independent_state_replay,
        )
        if verify_invariants_for_all or len(item["protein"]) <= 100 or idx == 0:
            invariants = verify_dp_4dp_collapse_invariants(
                protein_sequence=item["protein"],
                boundaries=item["boundaries"],
                codon_weights=weights,
                tunnel_length=item.get("tunnel_length", 35),
                half_width=item.get("half_width", 2),
                target_tau=target_tau,
            )
            invariants_executed_count += 1
        else:
            invariants = make_not_run_invariants_record(
                "Skipped on >100 aa fixture when verify_invariants_for_all=False; not claimed as executed."
            )
            invariants_not_run_count += 1

        entries.append(
            {
                "id": item["id"],
                "name": item["name"],
                "fixture_kind": item.get("fixture_kind", "SYNTHETIC_ALGORITHMIC_FIXTURE"),
                "sequence_provenance": item.get("sequence_provenance", "Unspecified fixture"),
                "boundary_provenance": item.get("boundary_provenance", "Synthetic test boundary"),
                "comparison": comp,
                "invariants": invariants,
            }
        )

    all_feasible = all(e["comparison"]["all_constraints_satisfied"] for e in entries)
    all_self_preferred = all(
        e["comparison"]["candidate_self_preference_verified"] for e in entries
    )
    executed_invariants = [
        e["invariants"]["all_invariants_verified"]
        for e in entries
        if e["invariants"]["execution_status"] == "EXECUTED"
    ]
    all_executed_invariants_ok = bool(executed_invariants) and all(executed_invariants)
    all_targets_invariants_verified: bool | None = (
        all_executed_invariants_ok if invariants_not_run_count == 0 else None
    )

    return {
        "evidence_classification": "UNVERIFIED_PUBLIC_CORPUS",
        "biological_validation": "NOT_AVAILABLE",
        "benchmark_scope_note": (
            "Matched-condition algorithmic comparison prototype across synthetic and public-subsequence "
            "fixtures with configured test boundaries. Does not constitute biological domain boundary "
            "validation or wet-lab expression/folding evidence."
        ),
        "small_instance_exhaustive_ground_truth": exhaustive_gt,
        "benchmark_target_count": len(entries),
        "all_targets_feasible": all_feasible,
        "all_targets_candidate_self_preference_verified": all_self_preferred,
        "invariants_executed_target_count": invariants_executed_count,
        "invariants_not_run_target_count": invariants_not_run_count,
        "all_executed_invariants_verified": all_executed_invariants_ok,
        "all_targets_invariants_verified": all_targets_invariants_verified,
        "entries": entries,
    }


ANNOTATED_PUBLIC_PROTEIN_PANEL: tuple[dict[str, Any], ...] = (
    {
        "id": "ubiquitin_1ubq_76aa",
        "annotation_kind": 'PFAM_DOMAIN_END',
        "source_snapshot_date": "2026-10-29",
        "source_snapshots": [{'url': 'https://rest.uniprot.org/uniprotkb/P0CG48?format=json', 'sha256': '2a70e3f891cdcbce7c424244846cf87d1677cecdb86c9b39614b3faad4833d81'}, {'url': 'https://data.rcsb.org/rest/v1/core/polymer_entity/1UBQ/1', 'sha256': '2d62cd6a8280b12897a9c2a43f95867165543f0ae81fbf72fad2500903f0e7e0'}, {'url': 'https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00240/protein/uniprot/P0CG48/', 'sha256': 'e663e2840894bae3cbed1d459c78e58954a5c4565393db5d7bcc47f4e5e6bd7e'}],
        "name": "Human Ubiquitin (UniProt P0CG48 / PDB 1UBQ Chain A, 76 aa)",
        "fixture_kind": "ANNOTATED_PUBLIC_PROTEIN",
        "sequence_provenance": (
            "UniProtKB P0CG48 (UBC_HUMAN first cleaved monomer 1-76, sequence version 3, entry version 138) / "
            "RCSB PDB 1UBQ Chain A (76 aa)."
        ),
        "coordinate_system": (
            "1-based mature monomer / PDB 1UBQ Chain A residue numbering [1..76] "
            "(precursor offset = 0 for monomer unit)."
        ),
        "boundary_annotation_source": (
            "InterPro IPR000626 / Pfam PF00240 (Ubiquitin domain, residues 3-74; "
            "remaining mature-chain residues 75-76). Annotated domain C-terminal "
            "boundary b_1 = 74."
        ),
        "boundary_provenance": (
            "Domain annotation with verified source, sequence version, and coordinate system "
            "(Pfam PF00240 span 3..74 on 1-based monomer [1..76], b=[74]). With L=35, delta=2, "
            "static 3DP window is D_1=[72..76] whereas shifted 4DP window [107..111] > 76 is "
            "strictly out-of-bounds (W_1 = empty). Database annotation does not by itself prove "
            "optimal ribosomal pausing or wet-lab expression superiority."
        ),
        "protein": (
            "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYN"
            "IQKESTLHLVLRLRGG"
        ),
        "boundaries": [74],
        "tunnel_length": 35,
        "half_width": 2,
    },
    {
        "id": "adalimumab_lc_4nyl_214aa",
        "annotation_kind": 'CATH_DOMAIN_END',
        "source_snapshot_date": "2026-10-29",
        "source_snapshots": [{'url': 'https://data.rcsb.org/rest/v1/core/polymer_entity/4NYL/1', 'sha256': '5958c0e0d344ab6b76ff1203bda2030cd700b4b6d93a05a345b1f462f7679f4f'}, {'url': 'https://www.ebi.ac.uk/pdbe/api/mappings/4nyl', 'sha256': '8fc603d643e5de96aeaf644cb29fcfa32b865527b86f1fc5b6ebb5937b913ef6'}],
        "name": "Adalimumab Fab Light Chain (RCSB PDB 4NYL Entity 1 Chain L, 214 aa)",
        "fixture_kind": "ANNOTATED_PUBLIC_PROTEIN",
        "sequence_provenance": (
            "RCSB PDB 4NYL polymer Entity 1, 214 aa; author chain L maps to label asym A. "
            "Sequence checked against the source response retrieved 2026-10-29."
        ),
        "coordinate_system": "1-based polymer entity sequence positions [1..214]; use residue_number, not author numbering.",
        "boundary_annotation_source": (
            "PDBe 4nyl CATH 2.60.40.10: chain L/entity 1 domain 4nylL01 spans 2..107, "
            "domain 4nylL02 spans 108..212; selected end b=107. InterPro/HMMER spans differ."
        ),
        "boundary_provenance": (
            "Source-specific CATH domain-end annotation; not a universal immunoglobulin junction. "
            "D=[105..109], W=[140..144] for the configured L=35, half-width=2. "
            "Annotation does not establish optimal pausing or expression effects."
        ),
        "protein": (
            "DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPS"
            "RFSGSGSGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIKRTVAAPSVFIFPP"
            "SDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLT"
            "LSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
        ),
        "boundaries": [107],
        "tunnel_length": 35,
        "half_width": 2,
    },
    {
        "id": "hemopexin_n_propeller_p02790_235aa",
        "annotation_kind": 'UNIPROT_REPEAT_END_NOT_FOLDING_DOMAIN',
        "source_snapshot_date": "2026-10-29",
        "source_snapshots": [{'url': 'https://rest.uniprot.org/uniprotkb/P02790?format=json', 'sha256': '267dc28668399805d9ba2bb2e26d1ffbb0a4b6bd6d9c81f6dbce7b1facfd852d'}],
        "name": "Human Hemopexin Mature N-Propeller + Hinge Construct (UniProt P02790, 235 aa)",
        "fixture_kind": "ANNOTATED_PUBLIC_PROTEIN",
        "sequence_provenance": (
            "UniProtKB P02790 (HEMO_HUMAN, sequence version 2, entry version 223, snapshot 2026-10-29) precursor "
            "residues 24-258 corresponding to mature chain residues 1-235 after cleavage of "
            "the 23-aa N-terminal signal peptide (precursor 1-23)."
        ),
        "coordinate_system": (
            "1-based mature polypeptide residue indexing [1..235], mapped via "
            "mature_pos = uniprot_precursor_pos - 23."
        ),
        "boundary_annotation_source": (
            "UniProtKB P02790 FT REPEAT annotations (PROSITE PS51642 / Pfam PF00045 Hemopexin "
            "repeats 1-3 in the N-terminal 4-blade beta-propeller): Repeat 1 = precursor 53..93 "
            "(mature 30..70, C-terminus b_1 = 70); Repeat 2 = precursor 94..139 (mature 71..116, "
            "C-terminus b_2 = 116); Repeat 3 = precursor 140..184 (mature 117..161, C-terminus "
            "b_3 = 161)."
        ),
        "boundary_provenance": (
            "Repeat-end annotations (NOT independent folding-domain ends) with source and coordinate mapping "
            "(UniProt P02790 entry version 223 Hemopexin repeats 1-3 mapped by -23 signal peptide "
            "offset to mature 1-based b=[70, 116, 161]). Produces static 3DP windows "
            "D=[68..72, 114..118, 159..163] and shifted 4DP windows W=[103..107, 149..153, 194..198]. "
            "Database annotations do not by themselves prove optimal ribosomal pausing or wet-lab "
            "expression superiority."
        ),
        "protein": (
            "TPLPPTSAHGNVAEGETKPDPDVTERCSDGWSFDATTLDDNGTMLFFKGEFVWKSHKWDR"
            "ELISERWKNFPSPVDAAFRQGHNSVFLIKGDKVWVYPPEKKEKGYPKLLQDEFPGIPSPL"
            "DAAVECHRGECQAEGVLFFQGDREWFWDLATGTMKERSWPAVGNCSSALRWLGRYYCFQG"
            "NQFLRFDPVRGEVPPRYPRDVRDYFMPCPGRGHGHRNGTGHGNSTHHGPEYMRCS"
        ),
        "boundaries": [70, 116, 161],
        "tunnel_length": 35,
        "half_width": 2,
    },
)


def run_step2_annotated_public_protein_comparison(
    step1_result: dict[str, Any] | None = None,
    codon_weights: dict[str, float] | None = None,
    panel: Sequence[dict[str, Any]] = ANNOTATED_PUBLIC_PROTEIN_PANEL,
    weight: float = 2.5,
    target_tau: float = 0.45,
    verify_invariants_for_all: bool = False,
    independent_state_replay: bool = True,
) -> dict[str, Any]:
    """STEP 2 Suite: Runs 1DP, 3DP, and 4DP comparison on source-bound public proteins with domain/repeat annotations, gated on STEP 1 passing."""
    weights = (
        codon_weights if codon_weights is not None else load_codon_usage_table().codon_weights
    )
    step1 = (
        step1_result
        if step1_result is not None
        else run_step1_exhaustive_verification_suite(codon_weights=weights)
    )
    if not step1.get("step1_gate_passed", False):
        raise RuntimeError(
            "STEP 2 is gated on STEP 1 computational correctness passing; "
            "step1_gate_passed is False."
        )

    entries: list[dict[str, Any]] = []
    invariants_executed_count = 0
    invariants_not_run_count = 0

    for item in panel:
        comp = run_matched_dp_comparison(
            protein_sequence=item["protein"],
            boundaries=item["boundaries"],
            codon_weights=weights,
            tunnel_length=item.get("tunnel_length", 35),
            half_width=item.get("half_width", 2),
            weight=weight,
            target_tau=target_tau,
            independent_state_replay=independent_state_replay,
        )
        if verify_invariants_for_all or len(item["protein"]) <= 100:
            invariants = verify_dp_4dp_collapse_invariants(
                protein_sequence=item["protein"],
                boundaries=item["boundaries"],
                codon_weights=weights,
                tunnel_length=item.get("tunnel_length", 35),
                half_width=item.get("half_width", 2),
                target_tau=target_tau,
            )
            invariants_executed_count += 1
        else:
            invariants = make_not_run_invariants_record(
                "Skipped on >100 aa annotated protein when verify_invariants_for_all=False; not claimed as executed."
            )
            invariants_not_run_count += 1

        entries.append(
            {
                "id": item["id"],
                "name": item["name"],
                "fixture_kind": item["fixture_kind"],
                "sequence_provenance": item["sequence_provenance"],
                "coordinate_system": item["coordinate_system"],
                "boundary_annotation_source": item["boundary_annotation_source"],
                "boundary_provenance": item["boundary_provenance"],
                "annotation_kind": item["annotation_kind"],
                "source_snapshot_date": item["source_snapshot_date"],
                "source_snapshots": item["source_snapshots"],
                "comparison": comp,
                "invariants": invariants,
            }
        )

    all_feasible = all(e["comparison"]["all_constraints_satisfied"] for e in entries)
    all_self_preferred = all(
        e["comparison"]["candidate_self_preference_verified"] for e in entries
    )
    all_state_parity = all(
        e["comparison"]["state_graph_parity_verified"] for e in entries
    )
    executed_invariants = [
        e["invariants"]["all_invariants_verified"]
        for e in entries
        if e["invariants"]["execution_status"] == "EXECUTED"
    ]
    all_executed_invariants_ok = bool(executed_invariants) and all(executed_invariants)
    all_targets_invariants_verified: bool | None = (
        all_executed_invariants_ok if invariants_not_run_count == 0 else None
    )

    return {
        "step": "STEP_2_ANNOTATED_PUBLIC_PROTEIN_COMPARISON",
        "evidence_classification": "UNVERIFIED_PUBLIC_CORPUS",
        "biological_validation": "NOT_AVAILABLE",
        "scope_note": (
            "Evaluates candidate design differences and common computational metrics across "
            "1DP, 3DP (L=0), and 4DP (L=35) on source-bound public proteins using source-specific "
            "domain ends or repeat ends (not equivalent to independent folding-domain ends). Database annotations "
            "carry uncertainty and do not directly prove optimal pause locations or wet-lab "
            "expression/yield superiority."
        ),
        "step1_prerequisite_passed": True,
        "annotated_target_count": len(entries),
        "all_targets_feasible": all_feasible,
        "all_targets_candidate_self_preference_verified": all_self_preferred,
        "all_targets_state_graph_parity_verified": all_state_parity,
        "invariants_executed_target_count": invariants_executed_count,
        "invariants_not_run_target_count": invariants_not_run_count,
        "all_executed_invariants_verified": all_executed_invariants_ok,
        "all_targets_invariants_verified": all_targets_invariants_verified,
        "entries": entries,
    }


def run_job510_two_step_validation(
    codon_weights: dict[str, float] | None = None,
    verify_invariants_for_all: bool = False,
    independent_state_replay: bool = True,
) -> dict[str, Any]:
    """Executes STEP 1 (exhaustive computational correctness verification) followed sequentially by STEP 2 (annotated public protein comparison)."""
    weights = (
        codon_weights if codon_weights is not None else load_codon_usage_table().codon_weights
    )
    step1 = run_step1_exhaustive_verification_suite(codon_weights=weights)
    step2 = run_step2_annotated_public_protein_comparison(
        step1_result=step1,
        codon_weights=weights,
        verify_invariants_for_all=verify_invariants_for_all,
        independent_state_replay=independent_state_replay,
    )
    return {
        "job_id": "510",
        "evidence_classification": "UNVERIFIED_PUBLIC_CORPUS",
        "biological_validation": "NOT_AVAILABLE",
        "step1_exhaustive_verification": step1,
        "step2_annotated_public_protein_comparison": step2,
    }

