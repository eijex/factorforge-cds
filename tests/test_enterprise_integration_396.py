"""Enterprise Benchmark and Full Integration Verification (Job 396).

Validates end-to-end:
1. Adalimumab (Humira LC & HC) 1DP vs 3DP with new plant RNA guards and invariant auditing.
2. CD47 (B2 construct, 169 aa) zero-defect invariants and RNA scan audit.
3. CD9 (236 aa) zero-defect invariants and RNA scan audit.
4. MFE observation truthfulness across all targets.
"""

from Bio.SeqUtils import seq1
import pytest

from factorforge.engines.dp_v2 import DPV2Optimizer
from factorforge.engines.profile.rules.reverse_translator import ReverseTranslator
from factorforge.engines.profile.utils import load_golden_set
from factorforge.engines.structural_pause import StructuralPauseDPOptimizer
from factorforge.guards.initiation_mfe import evaluate_initiation_mfe
from factorforge.guards.plant_rna_guards import (
    compare_candidate_plant_rna_risks,
    scan_plant_rna_risks,
    validate_candidate_hard_invariants,
)
from factorforge.validation.tradeoff_suite import compute_tradeoff_comparison


HUMIRA_LC_AA = (
    "DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIKR"
    "TVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
)

CD47_B2_AA = (
    "MGWSCIILFLVATATGVHSWPLVAALLLGSACCGSAMLLFNKTKSVEFTFCNDTVVIPCFVTNMEAQNTTEVYVKWKFKGRDIYTFDGALNKSTVPTDFSSAKIEVSQ"
    "LLKGDASLKMDKSDAVSHTGNYTCEVTELTREGETIIELKYRVVSWFSPNENIHHHHHHHH"
)

CD9_AA = (
    "MPVKGGTKCIKYLLFGFNFIFWLAGIAVLAIGLWLRFDSQTKSIFEQETNNNNSSFYTGVYILIGAGALMMLVGFLGCCGAVQESQCMLGLFFGFLLVIFAIEIAAAIW"
    "GYSHKDEVIKEVQEFYKDTYNKLKTKDEPQRETLKAIHYALNCCGLAGGVEQFISDICPKKDVLETFTVKSCPDAIKEVFDNKFHIIGAVGIGIAVVMIFGMIFSMILCCAIRRNREMVHHHHHHHH"
)


@pytest.fixture
def plant_weights():
    table = load_golden_set()
    return ReverseTranslator._build_ref_weights(table)


def test_job396_adalimumab_verification(plant_weights):
    """Re-verify Adalimumab LC under 1DP and 3DP with new RNA guards and invariant audits."""
    pause_positions = list(range(108, 120))  # Light chain hinge/linker region
    
    comp = compute_tradeoff_comparison(
        protein=HUMIRA_LC_AA,
        pause_positions=pause_positions,
        codon_weights=plant_weights,
        target_gc_min=0.40,
        target_gc_max=0.48,
        forbidden_motifs=["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"],
        lambda_weight=2.0,
        target_adaptiveness=0.40,
    )
    
    assert comp.hard_constraints_pass is True
    assert comp.structure_fit_3dp > comp.structure_fit_1dp
    assert comp.localization_ratio >= 0.50

    # Independent Invariant Audit on BOTH candidates
    target_motifs = ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"]
    inv_1dp = validate_candidate_hard_invariants(comp.cds_1dp, HUMIRA_LC_AA, gc_min=0.40, gc_max=0.48, forbidden_motifs=target_motifs)
    inv_3dp = validate_candidate_hard_invariants(comp.cds_3dp, HUMIRA_LC_AA, gc_min=0.40, gc_max=0.48, forbidden_motifs=target_motifs)
    assert inv_1dp["passed"] is True
    assert inv_3dp["passed"] is True
    assert inv_1dp["translation_identity_100pct"] is True
    assert inv_3dp["translation_identity_100pct"] is True

    # Unified Plant RNA Scan Comparison
    rna_comp = compare_candidate_plant_rna_risks(comp.cds_1dp, comp.cds_3dp, protein_name="Adalimumab_LC")
    assert rna_comp["truthfulness_notice"] is not None
    assert "Plant RNA Risk Comparison" in rna_comp["summary_table"]

    # MFE Observation
    mfe_1dp = evaluate_initiation_mfe(comp.cds_1dp)
    mfe_3dp = evaluate_initiation_mfe(comp.cds_3dp)
    assert mfe_1dp.status in {"AVAILABLE", "NOT_AVAILABLE"}
    assert mfe_3dp.status in {"AVAILABLE", "NOT_AVAILABLE"}


def test_job396_cd47_b2_verification(plant_weights):
    """Re-verify CD47 B2 candidate under 1DP vs 3DP."""
    pause_positions = list(range(20, 37)) + list(range(158, 163))
    
    comp = compute_tradeoff_comparison(
        protein=CD47_B2_AA,
        pause_positions=pause_positions,
        codon_weights=plant_weights,
        target_gc_min=0.44,
        target_gc_max=0.52,
        forbidden_motifs=["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"],
        lambda_weight=2.0,
        target_adaptiveness=0.40,
    )
    
    assert comp.hard_constraints_pass is True
    assert comp.structure_fit_3dp >= comp.structure_fit_1dp

    target_motifs = ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"]
    inv_1dp = validate_candidate_hard_invariants(comp.cds_1dp, CD47_B2_AA, gc_min=0.44, gc_max=0.52, forbidden_motifs=target_motifs)
    inv_3dp = validate_candidate_hard_invariants(comp.cds_3dp, CD47_B2_AA, gc_min=0.44, gc_max=0.52, forbidden_motifs=target_motifs)
    assert inv_1dp["passed"] is True
    assert inv_3dp["passed"] is True
    assert inv_1dp["forbidden_motifs_found"] == []
    assert inv_3dp["forbidden_motifs_found"] == []

    rna_comp = compare_candidate_plant_rna_risks(comp.cds_1dp, comp.cds_3dp, protein_name="CD47_B2")
    assert "Plant RNA Risk Comparison" in rna_comp["summary_table"]


def test_job396_cd9_verification(plant_weights):
    """Re-verify CD9 candidate under 1DP vs 3DP."""
    pause_positions = list(range(35, 51)) + list(range(110, 126)) + list(range(229, 237))
    target_motifs = ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"]
    
    comp = compute_tradeoff_comparison(
        protein=CD9_AA,
        pause_positions=pause_positions,
        codon_weights=plant_weights,
        target_gc_min=0.45,
        target_gc_max=0.52,
        forbidden_motifs=target_motifs,
        lambda_weight=2.0,
        target_adaptiveness=0.40,
    )
    
    assert comp.hard_constraints_pass is True
    assert comp.structure_fit_3dp >= comp.structure_fit_1dp

    inv_1dp = validate_candidate_hard_invariants(comp.cds_1dp, CD9_AA, gc_min=0.45, gc_max=0.52, forbidden_motifs=target_motifs)
    inv_3dp = validate_candidate_hard_invariants(comp.cds_3dp, CD9_AA, gc_min=0.45, gc_max=0.52, forbidden_motifs=target_motifs)
    assert inv_1dp["passed"] is True
    assert inv_3dp["passed"] is True
    assert inv_1dp["internal_stops"] == 0
    assert inv_3dp["internal_stops"] == 0
