"""Equivalence, parity, and invariant tests for hoisted DP engine loops."""

from __future__ import annotations

import math

import pytest
from Bio.Seq import Seq

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE, load_codon_usage_table
from factorforge.constraints.type_iis import (
    get_canonical_forbidden_motifs,
    resolve_restriction_site_motifs,
)
from factorforge.engines.dp_v2 import DPV2Optimizer
from factorforge.engines.dp_v2_4 import DPV240Optimizer
from factorforge.engines.structural_pause import StructuralPauseDPOptimizer

SYNTHETIC_PROTEIN_SHORT = "MASTWQHLEKDPVNRSGYF"
_AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
# Independently generated synthetic 472 aa sequence (non-biological motif progression)
SYNTHETIC_PROTEIN_LONG_472 = "M" + "".join(
    _AA_ALPHABET[(i * 7 + (i // 11) * 3) % len(_AA_ALPHABET)] for i in range(471)
)


def _reference_dpv2_uncached(protein: str, weights: dict[str, float], motifs, gc_min=0.40, gc_max=0.47):
    """Uncached reference DPV2 DP matching pre-hoist iteration order."""
    aa_to_codons: dict[str, list[str]] = {}
    for codon, aa in STANDARD_GENETIC_CODE.items():
        if aa != "*":
            aa_to_codons.setdefault(aa, []).append(codon)
    automaton = DPV2Optimizer(forbidden_motifs=motifs).automaton
    clean_protein = "".join(protein.upper().split()).rstrip("*")
    n_codons = len(clean_protein)
    total_nt = n_codons * 3
    min_gc = math.ceil(gc_min * total_nt)
    max_gc = math.floor(gc_max * total_nt)
    layers = [{(0, 0): (0.0, None, None, None)}]
    for aa in clean_protein:
        cur = layers[-1]
        nxt = {}
        for (g, s), (score, _, _, _) in cur.items():
            for codon in aa_to_codons[aa]:
                ns, bad = automaton.step_codon(s, codon)
                if bad:
                    continue
                ng = g + codon.count("G") + codon.count("C")
                w = weights.get(codon, 1e-4)
                nsc = score + math.log(max(w, 1e-6))
                k = (ng, ns)
                if k not in nxt:
                    nxt[k] = (nsc, codon, g, s)
                else:
                    es, ec, _, _ = nxt[k]
                    if nsc > es + 1e-12 or (abs(nsc - es) <= 1e-12 and ec is not None and codon < ec):
                        nxt[k] = (nsc, codon, g, s)
        layers.append(nxt)
    mid = (min_gc + max_gc) / 2.0
    in_band = [item for item in layers[-1].items() if min_gc <= item[0][0] <= max_gc]
    best = max(in_band, key=lambda item: (item[1][0], -abs(item[0][0] - mid), -item[0][0], -item[0][1]))
    codons = []
    cg, cs = best[0]
    for idx in range(n_codons, 0, -1):
        _, c, pg, ps = layers[idx][(cg, cs)]
        codons.append(c)
        cg, cs = pg, ps
    codons.reverse()
    return "".join(codons), best[1][0]


def _reference_3dp_uncached(
    protein: str,
    weights: dict[str, float],
    motifs,
    positions: list[int],
    target_adaptiveness: float = 0.40,
    penalty_weight: float = 2.0,
    gc_min: float = 0.35,
    gc_max: float = 0.55,
):
    """Uncached pre-hoist reference solver for StructuralPauseDPOptimizer (3DP)."""
    aa_to_codons: dict[str, list[str]] = {}
    for codon, aa in STANDARD_GENETIC_CODE.items():
        if aa != "*":
            aa_to_codons.setdefault(aa, []).append(codon)
    pos_set = set(positions)
    automaton = StructuralPauseDPOptimizer(
        positions=positions,
        target_adaptiveness=target_adaptiveness,
        penalty_weight=penalty_weight,
        structure_reference="synthetic-3dp",
        forbidden_motifs=motifs,
    ).automaton
    clean_protein = "".join(protein.upper().split()).rstrip("*")
    n_codons = len(clean_protein)
    total_nt = n_codons * 3
    min_gc = math.ceil(gc_min * total_nt)
    max_gc = math.floor(gc_max * total_nt)
    layers = [{(0, 0): (0.0, None, None, None)}]
    for idx, aa in enumerate(clean_protein, start=1):
        cur = layers[-1]
        nxt = {}
        in_pause = idx in pos_set
        for (g, s), (score, _, _, _) in cur.items():
            for codon in aa_to_codons[aa]:
                ns, bad = automaton.step_codon(s, codon)
                if bad:
                    continue
                ng = g + codon.count("G") + codon.count("C")
                w = weights.get(codon, 1e-4)
                penalty = penalty_weight * abs(w - target_adaptiveness) if in_pause else 0.0
                nsc = score + math.log(max(w, 1e-6)) - penalty
                k = (ng, ns)
                if k not in nxt:
                    nxt[k] = (nsc, codon, g, s)
                else:
                    es, ec, _, _ = nxt[k]
                    if nsc > es + 1e-12 or (abs(nsc - es) <= 1e-12 and ec is not None and codon < ec):
                        nxt[k] = (nsc, codon, g, s)
        layers.append(nxt)
    mid = (min_gc + max_gc) / 2.0
    in_band = [item for item in layers[-1].items() if min_gc <= item[0][0] <= max_gc]
    best = max(in_band, key=lambda item: (item[1][0], -abs(item[0][0] - mid), -item[0][0], -item[0][1]))
    codons = []
    cg, cs = best[0]
    for layer_idx in range(n_codons, 0, -1):
        _, c, pg, ps = layers[layer_idx][(cg, cs)]
        codons.append(c)
        cg, cs = pg, ps
    codons.reverse()
    return "".join(codons), best[1][0]


def _reference_4dp_uncached(
    protein: str,
    weights: dict[str, float],
    enzyme_names: tuple[str, ...],
    pause_positions: list[int],
    pause_target_adaptiveness: float = 0.45,
    pause_penalty_weight: float = 2.0,
    gc_min: float = 0.40,
    gc_max: float = 0.47,
):
    """Uncached pre-hoist reference solver for DPV240Optimizer (4DP, exact step_codon)."""
    aa_to_codons: dict[str, list[str]] = {}
    for codon, aa in STANDARD_GENETIC_CODE.items():
        if aa != "*":
            aa_to_codons.setdefault(aa, []).append(codon)
    opt = DPV240Optimizer(
        forbidden_enzymes=enzyme_names,
        pause_positions=pause_positions,
        pause_target_adaptiveness=pause_target_adaptiveness,
        pause_penalty_weight=pause_penalty_weight,
    )
    clean_protein = "".join(protein.upper().split()).rstrip("*")
    n_codons = len(clean_protein)
    total_nt = n_codons * 3
    min_gc = math.ceil(gc_min * total_nt)
    max_gc = math.floor(gc_max * total_nt)
    ramp_k = min(opt.ramp_length_codons, n_codons)
    ramp_nt = ramp_k * 3
    apply_init = n_codons >= opt.ramp_length_codons
    nom_min = math.ceil(opt.initiation_gc_min * ramp_nt)
    nom_max = math.floor(opt.initiation_gc_max * ramp_nt)
    if nom_min > nom_max:
        mid_init = round(opt.nominal_ramp_gc * ramp_nt)
        nom_min = nom_max = mid_init
    ramp_min_poss = sum(min(c.count("G") + c.count("C") for c in aa_to_codons[aa]) for aa in clean_protein[:ramp_k])
    ramp_max_poss = sum(max(c.count("G") + c.count("C") for c in aa_to_codons[aa]) for aa in clean_protein[:ramp_k])
    init_min = max(nom_min, ramp_min_poss)
    init_max = min(nom_max, ramp_max_poss)
    if init_min > init_max:
        if ramp_min_poss > nom_max:
            init_min, init_max = ramp_min_poss, min(ramp_min_poss + 2, ramp_max_poss)
        else:
            init_min, init_max = max(ramp_max_poss - 2, ramp_min_poss), ramp_max_poss

    active_pause_map = opt._resolve_active_pause_map(clean_protein)
    layers = [{(0, 0): (0.0, None, None, None, 0)}]
    for i, aa in enumerate(clean_protein):
        cur = layers[-1]
        nxt = {}
        for (g, s), (score, _, _, _, prev_rank) in cur.items():
            for codon in aa_to_codons[aa]:
                ns, bad = opt.automaton.step_codon(s, codon)
                if bad:
                    continue
                ng = g + codon.count("G") + codon.count("C")
                fit = opt.compute_codon_log_fitness(i, codon, weights, active_pause_map)
                nsc = score + fit
                k = (ng, ns)
                if k not in nxt:
                    nxt[k] = (nsc, codon, g, s, prev_rank)
                else:
                    esc, ec, _, _, epr = nxt[k]
                    if nsc > esc + 1e-12 or (abs(nsc - esc) <= 1e-12 and (prev_rank, codon) < (epr, ec)):
                        nxt[k] = (nsc, codon, g, s, prev_rank)
        if i == ramp_k - 1 and apply_init:
            filt = {k: v for k, v in nxt.items() if init_min <= k[0] <= init_max}
            if filt:
                nxt = filt
            else:
                d_min = min(
                    0 if init_min <= k[0] <= init_max else min(abs(k[0] - init_min), abs(k[0] - init_max))
                    for k in nxt
                )
                nxt = {
                    k: v
                    for k, v in nxt.items()
                    if (0 if init_min <= k[0] <= init_max else min(abs(k[0] - init_min), abs(k[0] - init_max))) == d_min
                }
        order = sorted(nxt.keys(), key=lambda st: (nxt[st][4], nxt[st][1]))
        rmap = {st: idx for idx, st in enumerate(order)}
        layers.append({st: (sc, cd, pg, ps, rmap[st]) for st, (sc, cd, pg, ps, _) in nxt.items()})

    valid_finals = []
    for (g, s), val in layers[-1].items():
        curr = s
        ok = True
        for ch in "TAA":
            curr = opt.automaton.step_nucleotide(curr, ch)
            if curr in opt.automaton.terminal_states:
                ok = False
                break
        if ok:
            valid_finals.append(((g, s), val))
    mid = (min_gc + max_gc) / 2.0
    in_band = [item for item in valid_finals if min_gc <= item[0][0] <= max_gc]
    best_item = max(in_band, key=lambda item: (item[1][0], -abs(item[0][0] - mid), -item[1][4]))
    best_state, (best_score, _, _, _, _) = best_item
    codons = []
    cg, cs = best_state
    for layer_idx in range(n_codons, 0, -1):
        val = layers[layer_idx][(cg, cs)]
        codons.append(val[1])
        cg, cs = val[2], val[3]
    codons.reverse()
    return "".join(codons), best_score


@pytest.mark.parametrize(
    ("gc_min", "gc_max", "enzymes"),
    [
        (0.35, 0.55, ("BsaI", "BsmBI", "NotI")),
        (0.40, 0.47, ("BsaI", "BsmBI", "SapI", "NotI")),
        (0.44, 0.46, ("BsaI", "BsmBI")),
    ],
)
def test_hoisted_dp_v2_matches_uncached_reference_bit_for_bit(gc_min, gc_max, enzymes):
    table = load_codon_usage_table()
    motifs = resolve_restriction_site_motifs(enzymes)
    expected_seq, expected_score = _reference_dpv2_uncached(
        SYNTHETIC_PROTEIN_SHORT + "*", table.codon_weights, motifs, gc_min, gc_max
    )
    actual = DPV2Optimizer(forbidden_motifs=motifs).optimize(
        SYNTHETIC_PROTEIN_SHORT + "*", table.codon_weights, target_gc_min=gc_min, target_gc_max=gc_max
    )
    assert actual["sequence"] == expected_seq
    assert actual["score"] == pytest.approx(expected_score, rel=1e-12, abs=1e-12)


@pytest.mark.parametrize("penalty_weight", [0.0, 0.5, 2.0, 5.0])
@pytest.mark.parametrize("target_adaptiveness", [0.30, 0.40, 0.50])
@pytest.mark.parametrize("gc_band", [(0.35, 0.55), (0.40, 0.47)])
def test_hoisted_structural_pause_3dp_matches_uncached_reference_and_scores(
    penalty_weight, target_adaptiveness, gc_band
):
    table = load_codon_usage_table()
    motifs = get_canonical_forbidden_motifs()
    positions = [4, 9, 14]
    gc_min, gc_max = gc_band
    expected_seq, expected_score = _reference_3dp_uncached(
        SYNTHETIC_PROTEIN_SHORT,
        table.codon_weights,
        motifs,
        positions=positions,
        target_adaptiveness=target_adaptiveness,
        penalty_weight=penalty_weight,
        gc_min=gc_min,
        gc_max=gc_max,
    )
    actual = StructuralPauseDPOptimizer(
        positions=positions,
        target_adaptiveness=target_adaptiveness,
        penalty_weight=penalty_weight,
        structure_reference="synthetic-3dp",
        forbidden_motifs=motifs,
    ).optimize(
        SYNTHETIC_PROTEIN_SHORT, table.codon_weights, target_gc_min=gc_min, target_gc_max=gc_max
    )
    assert actual["sequence"] == expected_seq
    assert actual["score"] == pytest.approx(expected_score, rel=1e-12, abs=1e-12)
    assert str(Seq(actual["sequence"]).translate()) == SYNTHETIC_PROTEIN_SHORT


@pytest.mark.parametrize("penalty_weight", [0.0, 0.5, 2.0, 5.0])
@pytest.mark.parametrize("target_adaptiveness", [0.30, 0.40, 0.50])
@pytest.mark.parametrize(
    ("gc_band", "enzymes"),
    [
        ((0.35, 0.55), ("BsaI", "BsmBI", "NotI")),
        ((0.40, 0.47), ("BsaI", "BsmBI", "SapI", "NotI")),
    ],
)
def test_hoisted_dp_v2_4_matches_uncached_reference_and_scores(
    penalty_weight, target_adaptiveness, gc_band, enzymes
):
    table = load_codon_usage_table()
    protein = SYNTHETIC_PROTEIN_SHORT + "MKTVR*"
    positions = [5, 11, 18]
    gc_min, gc_max = gc_band
    expected_seq, expected_score = _reference_4dp_uncached(
        protein,
        table.codon_weights,
        enzyme_names=enzymes,
        pause_positions=positions,
        pause_target_adaptiveness=target_adaptiveness,
        pause_penalty_weight=penalty_weight,
        gc_min=gc_min,
        gc_max=gc_max,
    )
    actual = DPV240Optimizer(
        forbidden_enzymes=enzymes,
        pause_positions=positions,
        pause_target_adaptiveness=target_adaptiveness,
        pause_penalty_weight=penalty_weight,
    ).optimize(
        protein_sequence=protein,
        codon_weights=table.codon_weights,
        target_gc_min=gc_min,
        target_gc_max=gc_max,
    )
    assert actual["sequence"] == expected_seq
    assert actual["score"] == pytest.approx(expected_score, rel=1e-12, abs=1e-12)
    assert str(Seq(actual["sequence"]).translate()) == protein.rstrip("*")


def test_hoisted_tie_breaking_equal_codon_weights_matches_uncached_reference():
    equal_weights = {c: 0.5 for c, aa in STANDARD_GENETIC_CODE.items() if aa != "*"}
    motifs = resolve_restriction_site_motifs(("BsaI", "BsmBI"))
    protein = "MLLSSPPGGAAK*"
    exp_3dp, exp_3dp_sc = _reference_3dp_uncached(
        protein, equal_weights, motifs, positions=[3, 7], penalty_weight=2.0, gc_min=0.40, gc_max=0.55
    )
    act_3dp = StructuralPauseDPOptimizer(
        positions=[3, 7],
        target_adaptiveness=0.40,
        penalty_weight=2.0,
        structure_reference="synthetic-tie",
        forbidden_motifs=motifs,
    ).optimize(protein, equal_weights, target_gc_min=0.40, target_gc_max=0.55)
    assert act_3dp["sequence"] == exp_3dp
    assert act_3dp["score"] == pytest.approx(exp_3dp_sc, rel=1e-12, abs=1e-12)

    exp_4dp, exp_4dp_sc = _reference_4dp_uncached(
        protein, equal_weights, ("BsaI", "BsmBI"), pause_positions=[3, 7], pause_penalty_weight=2.0, gc_min=0.40, gc_max=0.55
    )
    act_4dp = DPV240Optimizer(
        forbidden_enzymes=("BsaI", "BsmBI"), pause_positions=[3, 7], pause_target_adaptiveness=0.45, pause_penalty_weight=2.0
    ).optimize(protein, equal_weights, target_gc_min=0.40, target_gc_max=0.55)
    assert act_4dp["sequence"] == exp_4dp
    assert act_4dp["score"] == pytest.approx(exp_4dp_sc, rel=1e-12, abs=1e-12)


def test_hoisted_structural_pause_matches_dpv2_at_zero_penalty():
    table = load_codon_usage_table()
    motifs = get_canonical_forbidden_motifs()
    base = DPV2Optimizer(forbidden_motifs=motifs).optimize(
        SYNTHETIC_PROTEIN_LONG_472, table.codon_weights, target_gc_min=0.35, target_gc_max=0.55
    )
    zero_3dp = StructuralPauseDPOptimizer(
        positions=[215, 250],
        target_adaptiveness=0.40,
        penalty_weight=0.0,
        structure_reference="synthetic-472",
        forbidden_motifs=motifs,
    ).optimize(
        SYNTHETIC_PROTEIN_LONG_472, table.codon_weights, target_gc_min=0.35, target_gc_max=0.55
    )
    assert zero_3dp["sequence"] == base["sequence"]
    assert zero_3dp["score"] == pytest.approx(base["score"], rel=1e-12, abs=1e-12)
    assert str(Seq(zero_3dp["sequence"]).translate()) == SYNTHETIC_PROTEIN_LONG_472


def test_hoisted_dp_v2_4_472aa_deterministic_and_valid():
    table = load_codon_usage_table()
    motifs = resolve_restriction_site_motifs(("BsaI", "BsmBI", "BbsI", "SapI", "NotI"))
    pause_kwargs = {
        "pause_positions": [215, 250, 251],
        "pause_target_adaptiveness": 0.35,
        "pause_penalty_weight": 3.5,
    }
    run_a = DPV240Optimizer(forbidden_motifs=motifs, **pause_kwargs).optimize(
        protein_sequence=SYNTHETIC_PROTEIN_LONG_472,
        codon_weights=table.codon_weights,
        target_gc_min=0.35,
        target_gc_max=0.55,
    )
    run_b = DPV240Optimizer(forbidden_motifs=motifs, **pause_kwargs).optimize(
        protein_sequence=SYNTHETIC_PROTEIN_LONG_472,
        codon_weights=table.codon_weights,
        target_gc_min=0.35,
        target_gc_max=0.55,
    )
    assert run_a["sequence"] == run_b["sequence"]
    assert run_a["score"] == pytest.approx(run_b["score"], rel=1e-12, abs=1e-12)
    assert len(run_a["sequence"]) == 472 * 3
    assert str(Seq(run_a["sequence"]).translate()) == SYNTHETIC_PROTEIN_LONG_472
