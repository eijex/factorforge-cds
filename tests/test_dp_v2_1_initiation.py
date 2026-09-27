"""Tests for DP v2.1 Position-Dependent Initiation-Aware Exact DP Optimizer.

Verifies:
1. Hard Invariants (AA identity, frame, stops, Type IIS = 0 across U + CDS + STOP + R, GC bounds).
2. Hard GC UnsatisfiableDesignError fail-fast enforcement (and exploratory fallback).
3. Full Construct STOP & Flank Junction Automaton Veto.
4. True Full-Path Deterministic LexRank Tie-Breaking.
5. Exact Discrete Optimality via Brute-Force Equivalence on small design space.
6. 5' Initiation Window MFE Relaxation & Harmonization vs DP v2.0 Baseline.
7. SharedEvaluator integration (InitiationMetrics, context_digest, and initiation_hairpin_risk check).
"""

import itertools
import math
import pytest
from typing import Dict, List

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE, translate_dna, calculate_cai, calculate_gc
from factorforge.engines.dp_v2 import DPV2Optimizer
from factorforge.engines.dp_v2_1 import DPV21Optimizer, AA_TO_CODONS
from factorforge.engines.profile.utils import load_golden_set
from factorforge.engines.profile.rules.reverse_translator import ReverseTranslator
from factorforge.evaluation.evaluator import SharedEvaluator
from factorforge.evaluation.models import CheckEnforcement, CheckResultValue
from factorforge.utils.exceptions import UnsatisfiableDesignError

# Humira Sequences for Production Benchmark
HUMIRA_LC_AA = (
    "DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIKR"
    "TVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
)

HUMIRA_HC_AA = (
    "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLD"
    "YWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVE"
    "PKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNK"
    "ALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYT"
    "QKSLSLSPGK"
)


@pytest.fixture
def plant_codon_weights() -> Dict[str, float]:
    """Load N. benthamiana codon weights from golden set."""
    table = load_golden_set()
    return ReverseTranslator._build_ref_weights(table)


def test_dp_v2_1_hard_invariants(plant_codon_weights):
    """Verify that DP v2.1 strictly satisfies all hard biological & assembly invariants across full construct."""
    optimizer = DPV21Optimizer()
    res = optimizer.optimize(
        protein_sequence=HUMIRA_LC_AA,
        codon_weights=plant_codon_weights,
        target_gc_min=0.40,
        target_gc_max=0.47,
        upstream_context="GGTC",  # Flank prefix testing BsaI junction exclusion
        stop_codon="TAA",
        downstream_context="GAGA",
    )

    dna = res["sequence"]
    full_construct = res["full_construct"]

    # 1. Translational Identity
    assert translate_dna(dna) == HUMIRA_LC_AA
    assert len(dna) == len(HUMIRA_LC_AA) * 3

    # 2. No Internal Stops
    assert "*" not in translate_dna(dna)

    # 3. Full Construct Type IIS Exclusion (U + CDS + STOP + R)
    assert full_construct == "GGTC" + dna + "TAA" + "GAGA"
    for forbidden in ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG", "GAAGAC", "GTCTTC"]:
        assert forbidden not in full_construct, f"Forbidden motif {forbidden} found in construct!"
        assert forbidden not in dna, f"Forbidden motif {forbidden} found in CDS"

    # 4. GC Envelope Compliance
    assert res["gc_feasible"] is True
    assert 40.0 <= res["gc_percent"] <= 47.0
    assert res["constraint_scope"] == "FULL_CONSTRUCT"


def test_dp_v2_1_hard_gc_unsatisfiable_error(plant_codon_weights):
    """Verify that DP v2.1 raises UnsatisfiableDesignError when target GC band is unreachable."""
    optimizer = DPV21Optimizer()
    
    # Impossible GC band: [0.10, 0.15] for a standard protein
    with pytest.raises(UnsatisfiableDesignError) as exc_info:
        optimizer.optimize(
            protein_sequence="MWKDEC",
            codon_weights=plant_codon_weights,
            target_gc_min=0.10,
            target_gc_max=0.15,
            allow_nearest_infeasible=False,
        )
    assert "No synonymous candidate satisfies the requested GC band" in str(exc_info.value)

    # When exploratory allow_nearest_infeasible=True is requested, fallback is permitted
    fallback_res = optimizer.optimize(
        protein_sequence="MWKDEC",
        codon_weights=plant_codon_weights,
        target_gc_min=0.10,
        target_gc_max=0.15,
        allow_nearest_infeasible=True,
    )
    assert fallback_res["gc_feasible"] is False
    assert fallback_res["sequence"] is not None


def test_dp_v2_1_full_construct_stop_and_flank_veto(plant_codon_weights):
    """Verify that DP v2.1 avoids forbidden motifs across STOP codon and downstream flanks."""
    optimizer = DPV21Optimizer()
    
    # If CDS ends with codon ending in 'GGT' and stop codon is 'CTC' or downstream is 'CTC', BsaI 'GGTCTC' could form.
    # Here stop_codon="TAA", downstream_context="GAGACC" (which is forbidden in flank itself -> should raise ValueError)
    with pytest.raises(ValueError) as exc:
        optimizer.optimize(
            protein_sequence="MD",
            codon_weights=plant_codon_weights,
            upstream_context="GAGACC",  # Contains BsaI
        )
    assert "Upstream context" in str(exc.value)

    # Clean construct junction test
    res = optimizer.optimize(
        protein_sequence="MD",
        codon_weights=plant_codon_weights,
        target_gc_min=0.30,
        target_gc_max=0.60,
        upstream_context="GGTC",
        stop_codon="TAA",
        downstream_context="GAGT",
    )
    full = res["full_construct"]
    assert "GGTCTC" not in full
    assert "GAGACC" not in full


def test_dp_v2_1_full_path_lexrank_tie_breaking(plant_codon_weights):
    """Verify that identical-score paths are broken deterministically by full prefix LexRank."""
    optimizer = DPV21Optimizer(ramp_length_codons=2)
    
    # Run 10 times and assert identical outputs
    results = [
        optimizer.optimize(
            protein_sequence="MMWWKK",
            codon_weights=plant_codon_weights,
            target_gc_min=0.30,
            target_gc_max=0.60,
        )["sequence"]
        for _ in range(10)
    ]
    assert len(set(results)) == 1, "Full-path LexRank failed to produce 100% deterministic output across runs!"


def test_dp_v2_1_bruteforce_equivalence(plant_codon_weights):
    """Verify exact discrete optimality by comparing DP v2.1 against exhaustive brute-force enumeration."""
    short_peptide = "MWKDEC"  # 6 AA -> ~ 1 * 1 * 2 * 2 * 2 * 2 = 16 states, fast enumeration
    optimizer = DPV21Optimizer(ramp_length_codons=3)

    dp_res = optimizer.optimize(
        protein_sequence=short_peptide,
        codon_weights=plant_codon_weights,
        target_gc_min=0.35,
        target_gc_max=0.55,
    )

    # Exhaustive enumeration of all synonymous candidates
    synonymous_lists = [AA_TO_CODONS[aa] for aa in short_peptide]
    best_bf_score = -float("inf")
    best_bf_dna = None

    for candidate_codons in itertools.product(*synonymous_lists):
        dna = "".join(candidate_codons)
        gc_pct = calculate_gc(dna)
        if not (35.0 <= gc_pct <= 55.0):
            continue

        # Score with exact same piecewise objective
        score = 0.0
        for i, codon in enumerate(candidate_codons):
            score += optimizer.compute_codon_log_fitness(i, codon, plant_codon_weights)

        if score > best_bf_score + 1e-9:
            best_bf_score = score
            best_bf_dna = dna

    assert best_bf_dna is not None
    assert math.isclose(dp_res["score"], best_bf_score, rel_tol=1e-6)
    assert dp_res["sequence"] == best_bf_dna


def test_dp_v2_1_initiation_mfe_relaxation_vs_dp_v2_0(plant_codon_weights):
    """Compare DP v2.1 against DP v2.0 baseline on Humira LC and verify 5' structural relaxation."""
    evaluator = SharedEvaluator(codon_weights=plant_codon_weights)

    # 1. Run DP v2.0 Baseline
    opt_v20 = DPV2Optimizer()
    res_v20 = opt_v20.optimize(HUMIRA_LC_AA, plant_codon_weights, target_gc_min=0.40, target_gc_max=0.47)
    eval_v20 = evaluator.evaluate_candidate(res_v20["sequence"], HUMIRA_LC_AA)

    # 2. Run DP v2.1 Initiation-Aware DP
    opt_v21 = DPV21Optimizer(ramp_length_codons=15)
    res_v21 = opt_v21.optimize(HUMIRA_LC_AA, plant_codon_weights, target_gc_min=0.40, target_gc_max=0.47)
    eval_v21 = evaluator.evaluate_candidate(res_v21["sequence"], HUMIRA_LC_AA)

    # Assertions
    # Both must pass all hard integrity & GC checks
    assert eval_v20.passed is True
    assert eval_v21.passed is True

    # DP v2.1 must maintain high body CAI
    assert res_v21["cai_body"] >= 0.980

    # 5' Ramp GC% in DP v2.1 should be relaxed compared to rigid high-GC codons
    assert res_v21["gc_5p_ramp_percent"] <= res_v20["gc_percent"] + 5.0

    # Evaluator metrics verification
    assert eval_v21.metrics.cai_5p_ramp is not None
    assert eval_v21.metrics.cai_body is not None
    assert eval_v21.metrics.cai_body >= 0.980
    assert eval_v21.metrics.context_digest is not None


def test_evaluator_initiation_hairpin_risk_check(plant_codon_weights):
    """Verify SharedEvaluator properly emits PASS, WARNING, or INDETERMINATE for 5' initiation window MFE."""
    evaluator = SharedEvaluator(codon_weights=plant_codon_weights)
    
    optimizer = DPV21Optimizer()
    res = optimizer.optimize(HUMIRA_LC_AA, plant_codon_weights)
    eval_result = evaluator.evaluate_candidate(
        candidate_dna=res["sequence"],
        expected_protein=HUMIRA_LC_AA,
        upstream_context="GGTC",
        stop_codon="TAA",
        downstream_context="GAGA",
        initiation_mfe_threshold=-12.0,
    )

    # Check presence of initiation metrics
    assert eval_result.metrics.cai_5p_ramp is not None
    assert eval_result.metrics.cai_body is not None
    assert eval_result.metrics.context_digest is not None
    
    # Check initiation_hairpin_risk check presence
    check_names = [c.check_name for c in eval_result.checks]
    assert "initiation_hairpin_risk" in check_names
    hairpin_check = next(c for c in eval_result.checks if c.check_name == "initiation_hairpin_risk")
    assert hairpin_check.result in (CheckResultValue.PASS, CheckResultValue.WARNING, CheckResultValue.INDETERMINATE)
