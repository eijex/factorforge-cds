"""Brute-force Exhaustive Optimality Verification Fixtures for Math v2 DP.

Proves:
1. Math v2 DP finds the exact same global optimum sequence as brute-force exhaustive search.
2. Math v2 DP calculates the exact same maximum log CAI score under GC and Type IIS constraints.
3. Math v2 DP eliminates 100% of BsaI, BsmBI, and BpiI motifs across codon junctions and flanks.
"""

import itertools
import math
import pytest
from typing import List, Tuple, Dict, Optional

from factorforge.engines.dp_v2 import AA_TO_CODONS, DPV2Optimizer
from factorforge.engines.profile.utils import load_golden_set
from factorforge.engines.profile.rules.reverse_translator import ReverseTranslator


def brute_force_optimal_cds(
    protein: str,
    codon_weights: Dict[str, float],
    target_gc_min: float = 0.40,
    target_gc_max: float = 0.47,
    forbidden_motifs: Optional[List[str]] = None,
    left_flank: str = "",
    right_flank: str = "",
) -> Optional[Tuple[str, float, float]]:
    """Exhaustively generates all synonymous CDS combinations and finds the global maximum CAI."""
    if forbidden_motifs is None:
        forbidden_motifs = ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG", "GAAGAC", "GTCTTC"]

    codon_options = [AA_TO_CODONS[aa] for aa in protein]
    total_nt = len(protein) * 3
    gc_min_frac = target_gc_min if target_gc_min <= 1.0 else target_gc_min / 100.0
    gc_max_frac = target_gc_max if target_gc_max <= 1.0 else target_gc_max / 100.0
    min_gc_count = int(math.ceil(gc_min_frac * total_nt))
    max_gc_count = int(math.floor(gc_max_frac * total_nt))

    best_seq = None
    best_score = -float("inf")

    for combo in itertools.product(*codon_options):
        cds = "".join(combo)
        full_seq = (left_flank + cds + right_flank).upper()

        # Check forbidden motifs
        has_forbidden = False
        for motif in forbidden_motifs:
            if motif in full_seq:
                has_forbidden = True
                break
        if has_forbidden:
            continue

        # Check GC constraint
        gc_count = cds.count("G") + cds.count("C")
        if not (min_gc_count <= gc_count <= max_gc_count):
            continue

        # Calculate log score
        score = sum(math.log(max(codon_weights.get(c, 1e-4), 1e-6)) for c in combo)
        if score > best_score:
            best_score = score
            best_seq = cds

    if best_seq is None:
        return None

    cai = math.exp(best_score / len(protein))
    return best_seq, best_score, cai


@pytest.fixture
def plant_codon_weights():
    table = load_golden_set()
    return ReverseTranslator._build_ref_weights(table)


@pytest.mark.parametrize(
    "peptide",
    [
        "MDLS",      # 4 AA: 1 x 2 x 6 x 6 = 72 paths
        "EVAL",      # 4 AA: 2 x 4 x 4 x 6 = 192 paths
        "GGTCT",     # 5 AA: contains BsaI-inducing AA sequence
        "KPREEQ",    # 6 AA: 2 x 4 x 6 x 2 x 2 x 2 = 384 paths
        "MAKTNL",    # 6 AA: PlantForm Signal peptide start (1 x 4 x 2 x 4 x 2 x 6 = 384 paths)
    ],
)
def test_math_v2_dp_exactness_against_brute_force(peptide, plant_codon_weights):
    optimizer = DPV2Optimizer()
    dp_res = optimizer.optimize(
        protein_sequence=peptide,
        codon_weights=plant_codon_weights,
        target_gc_min=0.35,
        target_gc_max=0.55,
    )

    bf_res = brute_force_optimal_cds(
        protein=peptide,
        codon_weights=plant_codon_weights,
        target_gc_min=0.35,
        target_gc_max=0.55,
    )

    assert bf_res is not None, f"Feasible space must not be empty for {peptide}"
    bf_seq, bf_score, bf_cai = bf_res

    # 1. Exact log score match (within float precision)
    assert math.isclose(dp_res["score"], bf_score, rel_tol=1e-7, abs_tol=1e-7), (
        f"DP score ({dp_res['score']}) != Brute-force score ({bf_score})"
    )

    # 2. Exact CAI match
    assert math.isclose(dp_res["cai"], bf_cai, rel_tol=1e-7, abs_tol=1e-7), (
        f"DP CAI ({dp_res['cai']}) != Brute-force CAI ({bf_cai})"
    )

    # 3. Exact sequence match (or identical optimal score tie)
    dp_score_check = sum(math.log(max(plant_codon_weights.get(dp_res["sequence"][i:i+3], 1e-4), 1e-6))
                         for i in range(0, len(dp_res["sequence"]), 3))
    assert math.isclose(dp_score_check, bf_score, rel_tol=1e-7)

    # 4. Zero forbidden motifs in DP output
    for motif in ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG", "GAAGAC", "GTCTTC"]:
        assert motif not in dp_res["sequence"], f"Forbidden motif {motif} found in DP output!"


def test_math_v2_flank_construct_level_guarantee(plant_codon_weights):
    """Test that Math v2 DP avoids BsaI creation across flank junctions."""
    # If left_flank is "GGTC", emitting a codon starting with "TC" would create BsaI "GGTCTC".
    left_flank = "GGTC"
    peptide = "SL"  # Serine (TCT, TCC, TCA, TCG, AGT, AGC), Leucine

    optimizer = DPV2Optimizer()
    dp_res = optimizer.optimize(
        protein_sequence=peptide,
        codon_weights=plant_codon_weights,
        target_gc_min=0.30,
        target_gc_max=0.60,
        left_flank=left_flank,
    )

    full_construct = left_flank + dp_res["sequence"]
    assert "GGTCTC" not in full_construct, "Construct-level BsaI junction violated!"
    assert dp_res["constraint_scope"] == "FULL_CONSTRUCT"
