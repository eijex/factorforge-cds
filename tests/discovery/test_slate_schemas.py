"""Unit tests for Candidate Slate schemas and JSON serialization (Job 283A)."""

import json
from factorforge.discovery.schemas import (
    CandidateSlate,
    SlateCandidate,
    SlateSummary,
    TargetMetadata,
    TraitVector,
)


def test_slate_schema_serialization():
    target_meta = TargetMetadata(
        target_name="Target-mAb-A-LC",
        mature_protein_aa_length=214,
        construct_aa_length=235,
        signal_peptide_included=True,
        novelty_class="Class_B_RemoteHomolog",
        uncertainty_status="not_calibrated",
        evidence_tier="COMPUTATIONAL_HEURISTIC",
    )

    trait = TraitVector(
        cai_golden_set=0.988,
        global_gc_percent=40.2,
        initiation_mfe_kcal_mol=-8.4,
        local_50bp_gc_min=24.0,
        local_50bp_gc_max=54.0,
        homopolymer_max_run=5,
        synthesis_penalty=0.0,
        codon_context_prior=0.88,
    )

    candidate = SlateCandidate(
        rank=1,
        candidate_id="cand_01",
        generation_contract="dp_v2_1_1_initiation_favored",
        strategy_cluster="Cluster-A: Initiation-Optimized DP",
        utility_score=0.942,
        trait_vector=trait,
        sequence_dna="ATGGCTTGGTAA",
        sequence_digest="sha256:4a8b7c",
        rationale="Test rationale",
        is_pareto_optimal=True,
    )

    summary = SlateSummary(
        generated_pool_size=10,
        feasible_pool_size=8,
        pareto_front_size=4,
        top_k_count=1,
        diversity_index=0.75,
    )

    slate = CandidateSlate(
        run_id="FF-SLATE-TEST",
        target_metadata=target_meta,
        slate_summary=summary,
        candidates=[candidate],
    )

    data = slate.to_dict()
    assert data["$schema"] == "https://eijex.com/schemas/factorforge/candidate-slate-v0.1.json"
    assert data["run_id"] == "FF-SLATE-TEST"
    assert data["target_metadata"]["mature_protein_aa_length"] == 214
    assert data["target_metadata"]["construct_aa_length"] == 235
    assert data["target_metadata"]["signal_peptide_included"] is True
    assert data["target_metadata"]["evidence_tier"] == "COMPUTATIONAL_HEURISTIC"
    assert len(data["candidates"]) == 1
    assert data["candidates"][0]["trait_vector"]["cai_golden_set"] == 0.988

    # Verify JSON serializability
    json_str = json.dumps(data)
    assert "FF-SLATE-TEST" in json_str
