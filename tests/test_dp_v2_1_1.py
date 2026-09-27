"""Unit tests for DP v2.1.1 Local Composition Guard Optimizer."""

import pytest
from factorforge.analysis.metrics import (
    amino_acid_identity,
    load_codon_usage_table,
    translate_dna,
)
from factorforge.constraints.type_iis import get_canonical_forbidden_motifs
from factorforge.engines.dp_v2_1_1 import DPV211Optimizer


@pytest.fixture
def codon_weights():
    table = load_codon_usage_table()
    return table.codon_weights


def test_dp_v2_1_1_basic_properties(codon_weights):
    """Test standard protein optimization under DP v2.1.1."""
    # Test protein: 30 AA
    protein = "MKWVTFISLLFLFSSAYSRGVIKEDIVMTQ"
    optimizer = DPV211Optimizer()

    result = optimizer.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=0.40,
        target_gc_max=0.47,
    )

    dna = result["sequence"]
    assert len(dna) == len(protein) * 3
    assert translate_dna(dna) == protein
    assert amino_acid_identity(protein, dna) == 1.0
    assert 40.0 <= result["gc_percent"] <= 47.0
    assert result["gc_feasible"] is True

    # Check 5' initiation GC band (first 15 codons = 45 nt)
    ramp_dna = dna[:45]
    ramp_gc_count = ramp_dna.count("G") + ramp_dna.count("C")
    assert 9 <= ramp_gc_count <= 13
    assert 20.0 <= result["gc_5p_45nt_percent"] <= 28.89 + 1e-3


def test_dp_v2_1_1_homopolymer_rejection(codon_weights):
    """Verify that DP v2.1.1 strictly prevents >= 6-mer homopolymers."""
    # Poly-phenylalanine / poly-lysine / poly-proline which might tempt TTTTTT or AAAAAA
    protein = "MFFFFFFFFFFFFFKKKKKKKKKKPPPPPP"
    optimizer = DPV211Optimizer(homopolymer_max_run=5)

    result = optimizer.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=0.35,
        target_gc_max=0.55,
    )

    dna = result["sequence"]
    assert translate_dna(dna) == protein

    # Assert no 6-mers exist
    assert "AAAAAA" not in dna
    assert "TTTTTT" not in dna
    assert "CCCCCC" not in dna
    assert "GGGGGG" not in dna
    assert result["max_homopolymer_run"] <= 5


def test_dp_v2_1_1_type_iis_avoidance(codon_weights):
    """Verify full construct Type IIS site avoidance."""
    protein = "MEVQLVESGGGLVQPGGSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS"

    optimizer = DPV211Optimizer(
        forbidden_enzymes={"BsaI", "BsmBI", "BpiI", "SapI"},
        homopolymer_max_run=5,
    )

    result = optimizer.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=0.40,
        target_gc_max=0.45,
        upstream_context="ACGTACGT",
        stop_codon="TAA",
        downstream_context="GCTAGCTA",
    )

    full_seq = result["full_construct"]
    motifs = get_canonical_forbidden_motifs({"BsaI", "BsmBI", "BpiI", "SapI"})
    for m in motifs:
        assert m not in full_seq


def test_dp_v2_1_1_metrics_emission(codon_weights):
    """Verify all 20-row compatible metrics are emitted."""
    protein = "MKWVTFISLLFLFSSAYSRGVIKEDIVMTQSPDSLAVSLGERATINCKSSQSVLYSSNNKNYLAWYQQKPGQPPKLLIYWASTRESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQYYSTPYTFGQGTKVEIK"
    optimizer = DPV211Optimizer()

    result = optimizer.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=0.40,
        target_gc_max=0.47,
        upstream_context="AACCGAATTC",
    )

    assert "gc_constraint_status" in result
    assert result["gc_constraint_status"] in (
        "LOWER_BOUND_ACTIVE",
        "UPPER_BOUND_ACTIVE",
        "INTERIOR",
    )
    assert "dist_to_gc_min" in result
    assert "dist_to_gc_max" in result
    assert "min_50bp_gc_percent" in result
    assert "outlier_50bp_gc_count" in result
    assert "mfe_5p_window_kcal_mol" in result
    assert result["mfe_5p_status"] == "not_computed"
    assert result["mfe_5p_reason"] == "insufficient_upstream_context"
    assert result["mfe_5p_upstream_context_nt"] == 10
    assert "max_homopolymer_run" in result
    assert result["engine_version"] == "2.1.1"


def test_dp_v2_1_1_short_protein_uses_short_active_ramp(codon_weights):
    """Proteins shorter than 15 codons are constrained over their full length."""
    protein = "MKWVTF"
    result = DPV211Optimizer().optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=0.40,
        target_gc_max=0.60,
    )

    assert translate_dna(result["sequence"]) == protein
    assert result["ramp_length_codons"] == len(protein)
    assert result["initiation_gc_status"] == "NOT_APPLIED_SHORT_SEQUENCE"
    assert result["mfe_5p_fold_window_nt"] == 0


def test_dp_v2_1_1_high_gc_prefix_clamps_to_synonymous_envelope(codon_weights):
    """Intrinsic high-GC prefixes remain solvable and disclose band clamping."""
    protein = "MGGGGGGGGGGGGGG"
    result = DPV211Optimizer().optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=0.60,
        target_gc_max=0.90,
    )

    assert translate_dna(result["sequence"]) == protein
    assert result["initiation_gc_status"] == "SYNONYMOUS_ENVELOPE_CLAMPED"
    assert result["initiation_gc_active_min_count"] >= 29
    assert result["max_homopolymer_run"] <= 5


def test_dp_v2_1_1_automaton_relaxation_is_exact_and_reported(codon_weights):
    """A motif conflict selects only the nearest reachable initiation layer."""
    optimizer = DPV211Optimizer(
        forbidden_motifs=["TTCTTC"],
        include_reverse_complement=False,
        ramp_length_codons=2,
    )
    result = optimizer.optimize(
        protein_sequence="FF",
        codon_weights=codon_weights,
        target_gc_min=0.0,
        target_gc_max=1.0,
    )

    assert result["initiation_gc_status"] == "AUTOMATON_REACHABILITY_RELAXED"
    assert result["initiation_gc_active_min_count"] == 1
    assert result["initiation_gc_active_max_count"] == 1
    assert result["sequence"] != "TTCTTC"
    assert result["max_homopolymer_run"] <= 5
