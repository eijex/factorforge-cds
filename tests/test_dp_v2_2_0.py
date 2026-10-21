"""Unit tests for DP v2.2.0 Unified Optimizer and StemCell Integration."""

import math
import pytest
from factorforge.analysis.metrics import (
    amino_acid_identity,
    load_codon_usage_table,
    translate_dna,
)
from factorforge.engines.dp_v2_1_1 import DPV211Optimizer
from factorforge.engines.dp_v2_2 import DPV220Optimizer
from factorforge.engines.stemcell_3dp import StemCell3DPOptimizer


@pytest.fixture
def codon_weights():
    table = load_codon_usage_table()
    return table.codon_weights


def test_dp_v2_2_zero_penalty_parity(codon_weights):
    """Verify that when lambda=0 or pause_positions is empty, DP v2.2.0 exactly matches DP v2.1.1."""
    protein = "MKWVTFISLLFLFSSAYSRGVIKEDIVMTQ"

    opt_211 = DPV211Optimizer()
    res_211 = opt_211.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=0.40,
        target_gc_max=0.47,
        upstream_context="GCTAGC",
        stop_codon="TAA",
        downstream_context="ACTAGT",
    )

    opt_220 = DPV220Optimizer(
        pause_positions=None,
        pause_penalty_weight=0.0,
    )
    res_220 = opt_220.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=0.40,
        target_gc_max=0.47,
        upstream_context="GCTAGC",
        stop_codon="TAA",
        downstream_context="ACTAGT",
    )

    # Invariant: Bit-identical sequence output and score
    assert res_220["sequence"] == res_211["sequence"], "Sequences must be bit-identical"
    assert math.isclose(res_220["score"], res_211["score"], abs_tol=1e-9), "Scores must match"
    assert res_220["gc_feasible"] == res_211["gc_feasible"]
    assert res_220["initiation_gc_status"] == res_211["initiation_gc_status"]


def test_dp_v2_2_pause_penalty_induction(codon_weights):
    """Verify that a positive penalty weight at a designated residue induces codon pausing."""
    # Test protein with Arginine (R) at position 18:
    # R has 6 synonymous codons with varying adaptiveness:
    # e.g., CGT, CGC, CGA, CGG, AGA, AGG
    protein = "MKWVTFISLLFLFSSAYSRGVIKEDIVMTQ"
    pos = 18  # Residue R

    opt_baseline = DPV220Optimizer(
        pause_positions=None,
        pause_penalty_weight=0.0,
    )
    res_base = opt_baseline.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=0.35,
        target_gc_max=0.55,
    )
    base_codon = res_base["sequence"][(pos - 1) * 3 : pos * 3]
    base_w = codon_weights.get(base_codon, 1e-4)

    # Now enforce pause towards target adaptiveness tau = 0.45 with high penalty lambda = 5.0
    target_tau = 0.45
    opt_paused = DPV220Optimizer(
        pause_positions=[pos],
        pause_target_adaptiveness=target_tau,
        pause_penalty_weight=5.0,
        structure_reference="domain_hinge_reviewed",
    )
    res_paused = opt_paused.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=0.35,
        target_gc_max=0.55,
    )
    paused_codon = res_paused["sequence"][(pos - 1) * 3 : pos * 3]
    paused_w = codon_weights.get(paused_codon, 1e-4)

    assert translate_dna(res_paused["sequence"]) == protein
    assert res_paused["pause_conditioning"]["positions"] == [pos]
    assert pos in res_paused["pause_conditioning"]["details"]

    # Either codon changed to closer to tau, or delta to tau was minimized
    delta_base = abs(base_w - target_tau)
    delta_paused = abs(paused_w - target_tau)
    assert delta_paused <= delta_base + 1e-6


def test_dp_v2_2_synonymous_envelope_clamping(codon_weights):
    """Verify that DP v2.2.0 correctly clamps extreme-GC prefixes without error."""
    # 20 Ws: each W codon is TGG (66.7% GC)
    protein = "WWWWWWWWWWWWWWWWWWWW"
    opt = DPV220Optimizer()
    res = opt.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=0.60,
        target_gc_max=0.70,
    )
    assert res["initiation_gc_status"] == "SYNONYMOUS_ENVELOPE_CLAMPED"
    assert res["sequence"] == "TGG" * 20
    assert translate_dna(res["sequence"]) == protein


def test_dp_v2_2_reachability_relaxation(codon_weights):
    """Verify automaton reachability relaxation when motifs isolate active GC band."""
    opt = DPV220Optimizer(
        forbidden_motifs=["GGT", "GGA"],
        homopolymer_max_run=10,
    )
    res = opt.optimize(
        protein_sequence="GGGGGGGGGGGGGGGGGGGG",
        codon_weights=codon_weights,
        target_gc_min=0.50,
        target_gc_max=1.00,
    )
    assert res["initiation_gc_status"] == "AUTOMATON_REACHABILITY_RELAXED"
    assert res["gc_feasible"] is True


def test_dp_v2_2_stemcell_optimizer_integration():
    """Verify that StemCell3DPOptimizer executes cleanly using DPV220Optimizer backend."""
    protein = "MKWVTFISLLFLFSSAYSRGVIKEDIVMTQ"
    stem_opt = StemCell3DPOptimizer(
        target_gc_min=0.44,
        target_gc_max=0.52,
    )
    result = stem_opt.optimize(
        protein_sequence=protein,
        left_flank="ACCGTT",
        right_flank="TGACCA",
    )
    assert result["translation_identity"] == 1.0
    assert result["dna_length_nt"] == len(protein) * 3
    assert result["gc_in_band"] is True
    assert result["motif_clean"] is True
    assert "thermodynamics_37c" in result
