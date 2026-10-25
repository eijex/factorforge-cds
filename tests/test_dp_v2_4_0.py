"""Comprehensive Verification Suite for DP v2.4.0 Optimizer.

Verifies:
1. Tropical Semiring Max-Plus Algebraic Invariants.
2. Zero Parity Invariant with DP v2.3.0 baseline (bit-identical sequence and score).
3. TASEP Collision Smoothing and Shock Wave Potential Reduction.
4. 1D Wasserstein Optimal Transport Cadence Metric.
5. Transfer Matrix Feasibility & Topological Entropy Guard.
6. Full Construct Assembly & Type IIS Motif Invariants.
"""

from __future__ import annotations
import math

from factorforge.engines.dp_v2_3 import DPV230Optimizer
from factorforge.engines.dp_v2_4 import (
    DPV240Optimizer,
    DPV240EngineAdapter,
    TropicalSemiring,
    TASEPKineticsEngine,
    WassersteinCadenceEvaluator,
    TransferMatrixFeasibilityGuard,
)

BENCHMARK_CODON_WEIGHTS = {
    # Ala
    "GCT": 0.85, "GCC": 0.60, "GCA": 0.50, "GCG": 0.15,
    # Cys
    "TGT": 0.75, "TGC": 0.55,
    # Asp
    "GAT": 0.90, "GAC": 0.45,
    # Glu
    "GAA": 0.80, "GAG": 0.65,
    # Phe
    "TTT": 0.70, "TTC": 0.60,
    # Gly
    "GGT": 0.80, "GGA": 0.50, "GGG": 0.30, "GGC": 0.40,
    # His
    "CAT": 0.75, "CAC": 0.45,
    # Ile
    "ATT": 0.85, "ATC": 0.55, "ATA": 0.20,
    # Lys
    "AAA": 0.80, "AAG": 0.70,
    # Leu
    "TTG": 0.75, "TTA": 0.60, "CTT": 0.50, "CTA": 0.35, "CTG": 0.45, "CTC": 0.40,
    # Met
    "ATG": 1.00,
    # Asn
    "AAT": 0.85, "AAC": 0.50,
    # Pro
    "CCT": 0.80, "CCA": 0.60, "CCG": 0.20, "CCC": 0.40,
    # Gln
    "CAA": 0.85, "CAG": 0.50,
    # Arg
    "AGA": 0.85, "AGG": 0.50, "CGT": 0.40, "CGA": 0.30, "CGG": 0.20, "CGC": 0.25,
    # Ser
    "TCT": 0.85, "TCA": 0.60, "AGT": 0.55, "TCG": 0.25, "AGC": 0.50, "TCC": 0.45,
    # Thr
    "ACT": 0.85, "ACA": 0.60, "ACC": 0.50, "ACG": 0.20,
    # Val
    "GTT": 0.85, "GTG": 0.60, "GTA": 0.35, "GTC": 0.45,
    # Trp
    "TGG": 1.00,
    # Tyr
    "TAT": 0.80, "TAC": 0.50,
}

TEST_PROTEIN = "MEFGLSWLFLVAILKGVQCEVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDSSWYVFDYWGQGTLVTVSS"


def test_tropical_semiring_algebra():
    """Verify Tropical Max-Plus Semiring axioms."""
    assert TropicalSemiring.oplus(3.0, 5.0) == 5.0
    assert TropicalSemiring.oplus(5.0, 3.0) == 5.0
    assert TropicalSemiring.oplus(TropicalSemiring.ZERO, 7.0) == 7.0

    assert TropicalSemiring.otimes(3.0, 5.0) == 8.0
    assert TropicalSemiring.otimes(TropicalSemiring.ONE, 4.5) == 4.5
    assert TropicalSemiring.otimes(TropicalSemiring.ZERO, 4.5) == TropicalSemiring.ZERO


def test_zero_parity_invariant_with_dp_v2_3():
    """Verify bit-identical parity between DP v2.4.0 and DP v2.3.0 when kinetic smoothing is disabled."""
    opt_v23 = DPV230Optimizer(auto_annotate_pauses=False)
    opt_v24 = DPV240Optimizer(auto_annotate_pauses=False, enable_tasep_smoothing=False)

    res_v23 = opt_v23.optimize(
        protein_sequence=TEST_PROTEIN,
        codon_weights=BENCHMARK_CODON_WEIGHTS,
        target_gc_min=0.40,
        target_gc_max=0.47,
    )
    res_v24 = opt_v24.optimize(
        protein_sequence=TEST_PROTEIN,
        codon_weights=BENCHMARK_CODON_WEIGHTS,
        target_gc_min=0.40,
        target_gc_max=0.47,
    )

    assert res_v24["sequence"] == res_v23["sequence"], "DP v2.4.0 sequence diverged from DP v2.3.0!"
    assert math.isclose(res_v24["score"], res_v23["score"], rel_tol=1e-12), "Scores diverged!"
    assert res_v24["gc_percent"] == res_v23["gc_percent"], "GC% diverged!"
    assert res_v24["cai_global"] == res_v23["cai_global"], "Global CAI diverged!"


def test_tasep_collision_smoothing_activation():
    """Verify TASEP upstream ramp smoothing distributes deceleration across a multi-codon window."""
    # Place a deep pause at residue 30
    opt_raw = DPV240Optimizer(
        pause_positions=[30],
        pause_target_adaptiveness=0.30,
        pause_penalty_weight=3.0,
        enable_tasep_smoothing=False,
    )
    opt_smooth = DPV240Optimizer(
        pause_positions=[30],
        pause_target_adaptiveness=0.30,
        pause_penalty_weight=3.0,
        enable_tasep_smoothing=True,
        tasep_smoothing_window=3,
    )

    res_raw = opt_raw.optimize(
        protein_sequence=TEST_PROTEIN,
        codon_weights=BENCHMARK_CODON_WEIGHTS,
        target_gc_min=0.40,
        target_gc_max=0.47,
    )
    res_smooth = opt_smooth.optimize(
        protein_sequence=TEST_PROTEIN,
        codon_weights=BENCHMARK_CODON_WEIGHTS,
        target_gc_min=0.40,
        target_gc_max=0.47,
    )

    # In raw mode, only position 30 is active
    raw_active_pos = {p["position"] for p in res_raw["pause_conditioning"]["details"]}
    assert raw_active_pos == {30}

    # In smoothed mode, positions 28, 29 (30 - 2, 30 - 1) receive upstream buffer ramps
    smooth_active_pos = {p["position"] for p in res_smooth["pause_conditioning"]["details"]}
    assert 29 in smooth_active_pos or 28 in smooth_active_pos, "TASEP buffer ramp was not placed!"

    # Verify shock potential metrics are present
    assert "tasep" in res_smooth["kinetics_and_transport"]
    assert "collision_index" in res_smooth["kinetics_and_transport"]["tasep"]

    # Direct test of TASEPKineticsEngine helper
    test_smoothed = TASEPKineticsEngine.compute_smoothed_pause_map(
        TEST_PROTEIN, {30: [("membrane", 3.0, 0.30)]}, smoothing_window=3
    )
    assert 29 in test_smoothed or 28 in test_smoothed


def test_wasserstein_1d_cadence_metric():
    """Verify 1D Wasserstein cadence metric evaluates translation trajectory differences."""
    opt = DPV240Optimizer(auto_annotate_pauses=True, enable_tasep_smoothing=True)
    res = opt.optimize(
        protein_sequence=TEST_PROTEIN,
        codon_weights=BENCHMARK_CODON_WEIGHTS,
        target_gc_min=0.40,
        target_gc_max=0.47,
    )

    ot_dist = res["kinetics_and_transport"]["wasserstein_1d_cadence_distance"]
    assert isinstance(ot_dist, float)
    assert ot_dist >= 0.0, "Wasserstein distance must be non-negative!"

    # Direct test of WassersteinCadenceEvaluator helper
    direct_w1 = WassersteinCadenceEvaluator.evaluate_1d_wasserstein(
        "GCTGCTGCT", BENCHMARK_CODON_WEIGHTS, target_pause_map={2: [("test", 1.0, 0.5)]}
    )
    assert direct_w1 >= 0.0


def test_transfer_matrix_feasibility_guard():
    """Verify A Priori Transfer Matrix feasibility check."""
    opt = DPV240Optimizer()
    feasibility = TransferMatrixFeasibilityGuard.check_feasibility(
        TEST_PROTEIN, opt.automaton, ramp_length=15
    )
    assert feasibility["feasible"] is True
    assert feasibility["topological_entropy"] > 0.0
    assert feasibility["surviving_states_at_ramp"] > 0


def test_full_construct_invariants_and_adapter():
    """Verify full construct Type IIS site elimination and DPV240EngineAdapter."""
    adapter = DPV240EngineAdapter()
    result = adapter.optimize_cds(
        protein_sequence="MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAEEGPKLVASQAALA",
        codon_weights=BENCHMARK_CODON_WEIGHTS,
        target_gc_min=0.40,
        target_gc_max=0.47,
        upstream_context="TTCTAG",
        stop_codon="TAA",
        downstream_context="AAGCTT",
    )

    dna = result["sequence"]
    full_construct = "TTCTAG" + dna + "TAA" + "AAGCTT"

    for forbidden in ("GGTCTC", "GAGACC", "CGTCTC", "GAGACG", "GCGGCCGC"):
        assert forbidden not in full_construct, f"Forbidden motif {forbidden} violated in construct!"
    assert result["engine_version"] == "2.4.0"
