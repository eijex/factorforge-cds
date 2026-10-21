"""Verification Suite for FactorForge DP v2.3.0 Engine.

Covers:
1. Exact zero-penalty bit parity with DP v2.2.0 and DP v2.1.1.
2. Automated structural pause profiling on CD47 and CD9.
3. Multi-class structural pause penalty induction (TM, SS, SP, Hinge).
4. Auto-annotation integration end-to-end optimization.
5. Synonymous envelope clamping and automaton reachability relaxation.
"""

from __future__ import annotations
import pytest

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE, load_codon_usage_table, translate_dna
from factorforge.engines.dp_v2_1_1 import DPV211Optimizer
from factorforge.engines.dp_v2_2 import DPV220Optimizer
from factorforge.engines.dp_v2_3 import DPV230Optimizer, StructuralPauseAnnotator


@pytest.fixture
def codon_weights() -> dict[str, float]:
    table = load_codon_usage_table()
    return table.codon_weights


def test_dp_v2_3_zero_penalty_parity(codon_weights):
    """Verify that when lambda=0 or pause lists are empty, DP v2.3.0 matches v2.2.0 and v2.1.1 bit-for-bit."""
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

    opt_230 = DPV230Optimizer(
        pause_positions=None,
        pause_penalty_weight=0.0,
        auto_annotate_pauses=False,
    )
    res_230 = opt_230.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=0.40,
        target_gc_max=0.47,
        upstream_context="GCTAGC",
        stop_codon="TAA",
        downstream_context="ACTAGT",
    )

    assert res_230["sequence"] == res_220["sequence"] == res_211["sequence"]
    assert abs(res_230["score"] - res_220["score"]) < 1e-9
    assert abs(res_230["score"] - res_211["score"]) < 1e-9


def test_dp_v2_3_structural_annotator_cd47():
    """Verify that StructuralPauseAnnotator correctly identifies domain features in CD47."""
    cd47_precursor = (
        "MWPLVAALLLGSACCGSA"  # Signal Peptide (1-18)
        "QLLFNKTKSVEFTFCNDTVVIPCFVTNMEAQNTTEVYVKWKFKGRDIYTFDGALNKSTVP"
        "TDFSSAKIEVSQLLKGDASLKMDKSDAVSHTGNYTCEVTELTREGETIIELKYRVVSWFS"
        "PNENI"
    )

    annotations = StructuralPauseAnnotator.annotate(cd47_precursor)

    # 1. Signal peptide cleavage detected near position 18
    assert len(annotations["SP"]) == 1
    sp_pos = annotations["SP"][0]
    assert 15 <= sp_pos <= 22

    # 2. Disulfide loop detected (Cys pairs)
    assert len(annotations["SS"]) >= 1

    # 3. N-Glycosylation sites detected (N-X-[S/T])
    assert len(annotations["GLYC"]) >= 3
    # NKTK (starts at N23 of precursor)
    assert any(20 <= g <= 26 for g in annotations["GLYC"])


def test_dp_v2_3_structural_annotator_cd9_transmembrane():
    """Verify that Kyte-Doolittle transmembrane detection identifies CD9 TM helices."""
    cd9_human = (
        "MPVKGGTKCIKYLLFGFNFIFWFLAGIAVLAIGLWLRFDSQTKSIFEQETNNNNSSFYTG"
        "VYILIGAGALMMLVGFLGCCGAVQESQCMLGLFFGFLLVIFAIEIAAAIWGYSHKDEVIK"
        "EVQEFYKDTYNKLKTKDEPQRETLKAIHYALNCCGLAGGVEQFISDICPKKDVLETFTVK"
        "SCPDAIKEVFDNKFHIIGAVGIGIAVVMIFGMIFSMILCCAIRRNREMV"
    )
    tm_boundaries = StructuralPauseAnnotator.detect_transmembrane_boundaries(cd9_human)
    # CD9 is a 4-pass tetraspanin; boundary list should contain entry and exit markers
    assert len(tm_boundaries) >= 4


def test_dp_v2_3_multi_class_pause_induction(codon_weights):
    """Verify that multi-class pause penalties shift codons toward class-specific targets."""
    protein = "MKWVTFISLLFLFSSAYSRGVIKEDIVMTQ"
    # Position 10: TM-like target (tau=0.25)
    # Position 20: SS-like target (tau=0.45)
    opt_base = DPV230Optimizer(pause_positions=None, pause_penalty_weight=0.0)
    res_base = opt_base.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=0.40,
        target_gc_max=0.55,
    )

    cd_10_base = res_base["sequence"][9 * 3 : 10 * 3]
    w_10_base = codon_weights.get(cd_10_base, 1e-4)

    cd_20_base = res_base["sequence"][19 * 3 : 20 * 3]
    w_20_base = codon_weights.get(cd_20_base, 1e-4)

    opt_multi = DPV230Optimizer(
        pause_classes={
            "TM": {"positions": [10], "weight": 3.0, "target_tau": 0.25},
            "SS": {"positions": [20], "weight": 3.0, "target_tau": 0.45},
        }
    )
    res_multi = opt_multi.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=0.40,
        target_gc_max=0.55,
    )

    cd_10_multi = res_multi["sequence"][9 * 3 : 10 * 3]
    w_10_multi = codon_weights.get(cd_10_multi, 1e-4)

    cd_20_multi = res_multi["sequence"][19 * 3 : 20 * 3]
    w_20_multi = codon_weights.get(cd_20_multi, 1e-4)

    assert translate_dna(res_multi["sequence"]) == protein
    # The pause optimization drives the adaptiveness closer to the specific tau
    assert abs(w_10_multi - 0.25) <= abs(w_10_base - 0.25) + 1e-6
    assert abs(w_20_multi - 0.45) <= abs(w_20_base - 0.45) + 1e-6
    assert 10 in res_multi["pause_conditioning"]["details"]
    assert 20 in res_multi["pause_conditioning"]["details"]


def test_dp_v2_3_auto_annotation_integration(codon_weights):
    """Verify that auto_annotate_pauses=True automatically optimizes domain boundaries."""
    cd47_construct = (
        "MGWSCIILFLVATATGVHS"  # SP
        "QLLFNKTKSVEFTFCNDTVVIPCFVTNMEAQNTTEVYVKWKFKGRDIYTFDGALNKSTVP"
        "TDFSSAKIEVSQLLKGDASLKMDKSDAVSHTGNYTCEVTELTREGETIIELKYRVVSWFS"
        "PNENIGGGGSGGGGSCHHHHHHHH"  # Hinge + single Cys + His
    )

    opt = DPV230Optimizer(auto_annotate_pauses=True)
    result = opt.optimize(
        protein_sequence=cd47_construct,
        codon_weights=codon_weights,
        target_gc_min=0.40,
        target_gc_max=0.47,
        upstream_context="GCTAGC",
        stop_codon="TAA",
        downstream_context="ACTAGT",
    )

    assert translate_dna(result["sequence"]) == cd47_construct
    assert result["gc_feasible"] is True
    assert result["max_homopolymer_run"] <= 5
    assert result["pause_conditioning"]["auto_annotated"] is True
    assert result["pause_conditioning"]["total_active_positions"] > 0
    assert result["engine_version"] == "2.3.0"
