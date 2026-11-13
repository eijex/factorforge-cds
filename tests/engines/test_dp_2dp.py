from __future__ import annotations

from Bio.Seq import Seq

from factorforge.analysis.metrics import load_codon_usage_table
from factorforge.engines.dp_2dp import DP2DPOptimizer
from factorforge.engines.dp_v2_1_1 import DPV211Optimizer


def test_dp_2dp_preserves_protein_and_reduces_dinucleotide_junctions() -> None:
    protein = "MAKTNLFLFLIFSLLLSLSSAEVQLVESGGGLVQPGGSLRLSCAASGFTFD"
    table = load_codon_usage_table()

    res_1dp = DPV211Optimizer().optimize(
        protein_sequence=protein,
        codon_weights=table.codon_weights,
        target_gc_min=0.42,
        target_gc_max=0.55,
    )
    res_2dp = DP2DPOptimizer().optimize(
        protein_sequence=protein,
        codon_weights=table.codon_weights,
        target_gc_min=0.42,
        target_gc_max=0.55,
    )

    assert str(Seq(res_2dp["cds_only"]).translate()) == protein
    assert res_2dp["homopolymers"] == 0
    assert res_2dp["design_tier"] == "2DP"
    assert res_2dp["biological_validation"] == "NOT_AVAILABLE"
    # 2DP actively penalizes CpG and TpA junctions relative to 1DP
    cpg_1dp = sum(1 for i in range(len(res_1dp["sequence"]) - 1) if res_1dp["sequence"][i : i + 2] == "CG")
    tpa_1dp = sum(1 for i in range(len(res_1dp["sequence"]) - 1) if res_1dp["sequence"][i : i + 2] == "TA")
    assert res_2dp["cpg_count"] + res_2dp["tpa_count"] <= cpg_1dp + tpa_1dp
