# factorforge/tests/benchmarks/test_discovery_smoke.py
"""
CI Smoke Test for FactorForge Discovery Pilot Benchmark (Phase 283B).
Executes a fast subset of 3 targets across novelty tiers (Class A, B, C) in < 5 seconds.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

# Ensure repository root is in sys.path for benchmarks package import
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pytest

from benchmarks.run_discovery_benchmark import (
    compute_codon_distance,
    compute_pairwise_hamming,
    compute_slate_diversity_metrics,
    evaluate_slate_condition,
    run_benchmark,
)
from factorforge.discovery.filter import HardConstraintFilter
from factorforge.discovery.schemas import SlateCandidate, TraitVector
from factorforge.discovery.traits import TraitVectorExtractor


def test_diversity_metric_invariants():
    """Verify mathematical properties of diversity metric computations."""
    seq1 = "ATGGCGAAATTTTAA"
    seq2 = "ATGGCGAAATTTTAA"
    seq3 = "ATGGCCAAGTTCTAG"
    
    assert compute_pairwise_hamming(seq1, seq2) == 0.0
    assert compute_codon_distance(seq1, seq2) == 0.0
    
    h_dist = compute_pairwise_hamming(seq1, seq3)
    assert 0.0 < h_dist < 1.0
    
    c_dist = compute_codon_distance(seq1, seq3)
    assert 0.0 < c_dist < 1.0
    
    # Test candidate slate metrics
    tv = TraitVector(
        cai_golden_set=0.8,
        global_gc_percent=0.5,
        initiation_mfe_kcal_mol=-15.0,
        local_50bp_gc_min=0.4,
        local_50bp_gc_max=0.6,
        homopolymer_max_run=3,
        synthesis_penalty=0.0,
        codon_context_prior=0.85,
    )
    c1 = SlateCandidate(
        rank=1,
        candidate_id="c1",
        generation_contract="dp_v2_1_1_initiation_favored",
        strategy_cluster="Cluster-Initiation: Open 5' Structure",
        utility_score=0.9,
        trait_vector=tv,
        sequence_dna=seq1,
        sequence_digest="sha256_1",
        rationale="test",
        is_pareto_optimal=True,
    )
    c2 = SlateCandidate(
        rank=2,
        candidate_id="c2",
        generation_contract="profile_local_guarded",
        strategy_cluster="Cluster-Composition: Local Guarded",
        utility_score=0.85,
        trait_vector=tv,
        sequence_dna=seq3,
        sequence_digest="sha256_2",
        rationale="test",
        is_pareto_optimal=True,
    )
    
    d_mean, d_min, c_mean, strat_cov = compute_slate_diversity_metrics([c1, c2])
    assert d_mean == d_min == h_dist
    assert c_mean == c_dist
    assert strat_cov == 2


def test_discovery_smoke_fast_execution(tmp_path: Path):
    """
    CI smoke test running on a 3-target subset fixture (RbcS, sfGFP, VP28).
    Ensures end-to-end execution completes fast and produces valid output schemas.
    """
    smoke_targets = [
        {
            "benchmark_id": "SMK-A01",
            "class": "Class A",
            "protein_name": "RbcS (Rubisco Small Subunit)",
            "sequence": "MASSVLSSAAVATRSNVAQANMVAPFTGLKSAASFPVSRKQNLDITSIASNGGRVQCMRVWPPIGKKKFETLSYLPDLTDSELAKEVDYLIRNKWIPCVEFELEHGFVYREHGNSPGYYDGRYWTMWKLPMFGCTDATQVLAEVEEAKKAYPQAWIRIIGFDNVRQVQCISFIAYKPEGY",
            "construct_length": 180,
            "max_identity_nb_ref": 1.0,
        },
        {
            "benchmark_id": "SMK-B01",
            "class": "Class B",
            "protein_name": "sfGFP",
            "sequence": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
            "construct_length": 238,
            "max_identity_nb_ref": 0.32,
        },
        {
            "benchmark_id": "SMK-C01",
            "class": "Class C",
            "protein_name": "VP28 Viral Antigen",
            "sequence": "MDLSFTLSVVSAILAITAVIAVFIVIFRYHNTVTKTIETHTDNIETNMDENLRIPVTAEVGSGYFKMTDVSFDSDTLGKIKIRNGKSDAQMKEEDADLVITPVEGRALEVTVGQNLTFEGTFKVWNNTSRKINITGMQMVPKINPSKAFVGSSNTSSFTPVSIDEDEVGTFVCGTTFGAPIAATAGGNLFDMYVHVTYSGTETE",
            "construct_length": 204,
            "max_identity_nb_ref": 0.12,
        },
    ]
    
    smoke_fixture_path = tmp_path / "smoke_panel.json"
    with open(smoke_fixture_path, "w", encoding="utf-8") as f:
        json.dump({"manifest_version": "1.0.0", "reference_corpus_id": "NB-EXPRESSION-REF-v1", "targets": smoke_targets}, f)
        
    out_dir = tmp_path / "smoke_results"
    summary, rows = run_benchmark(
        panel_path=smoke_fixture_path,
        out_dir=out_dir,
        host="nbenthamiana",
        seed=42,
    )
    
    assert summary["target_count"] == 3
    assert (out_dir / "benchmark_summary.json").exists()
    assert (out_dir / "benchmark_details.csv").exists()
    assert (out_dir / "discovery_pilot_benchmark_report.md").exists()
    
    # Check that Top-3 and Top-5 slates strictly pass SlateValid for smoke targets
    for target_res in summary["detailed_results"]:
        assert target_res["conditions"]["single_dp_v2_1_1"]["slate_valid"] is True
        assert target_res["conditions"]["slate_top3"]["slate_valid"] is True
        assert target_res["conditions"]["slate_top5"]["slate_valid"] is True
        assert target_res["conditions"]["slate_top3"]["strategy_coverage"] >= 2


def test_full_panel_benchmark_integrity(tmp_path: Path):
    """
    Full 9-target benchmark verification test.
    Validates complete novelty panel (Class A, B, C) for 100% emitted validity,
    strategy coverage, and Top-K completeness metrics.
    """
    fixture_path = REPO_ROOT / "benchmarks" / "fixtures" / "discovery_novelty_panel.json"
    assert fixture_path.exists(), f"Fixture missing at {fixture_path}"

    out_dir = tmp_path / "full_panel_results"
    summary, rows = run_benchmark(
        panel_path=fixture_path,
        out_dir=out_dir,
        host="nbenthamiana",
        seed=42,
    )

    assert summary["target_count"] == 9
    assert (out_dir / "benchmark_summary.json").exists()
    assert (out_dir / "benchmark_details.csv").exists()
    assert (out_dir / "discovery_pilot_benchmark_report.md").exists()

    # Verify 100% emitted-candidate hard constraint validity across all conditions and targets
    for target_res in summary["detailed_results"]:
        for cond_name in ["single_dp_v2_1_1", "slate_top1", "slate_top3", "slate_top5"]:
            c_info = target_res["conditions"][cond_name]
            assert c_info["emitted_candidate_valid_rate"] == 1.0, f"Target {target_res['protein_name']} failed validity in {cond_name}"
            for cand in c_info["candidates"]:
                assert cand["is_valid"] is True
                assert len(cand["validation_errors"]) == 0

    # Verify Top-1 and Top-3 100% completeness, and Top-5 88.9% (8/9) completeness
    top1_complete = sum(1 for t in summary["detailed_results"] if t["conditions"]["slate_top1"]["slate_complete"])
    top3_complete = sum(1 for t in summary["detailed_results"] if t["conditions"]["slate_top3"]["slate_complete"])
    top5_complete = sum(1 for t in summary["detailed_results"] if t["conditions"]["slate_top5"]["slate_complete"])

    assert top1_complete == 9
    assert top3_complete == 9
    assert top5_complete >= 8  # Full panel achieves >= 88.9% (now 100% with enhanced diversity trajectories)

