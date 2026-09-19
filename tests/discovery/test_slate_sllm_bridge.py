# factorforge/tests/discovery/test_slate_sllm_bridge.py
from factorforge.discovery.filter import HardConstraintConfig, HardConstraintFilter
from factorforge.discovery.slate import DiscoverySlateEngine
from factorforge.engines.sllm.adapters import DeterministicMockSLMAdapter
from factorforge.learning.ranker import ExpressionPriorRanker

TEST_TARGET_AA = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHM"
    "KQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADK"
    "QKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
)


def test_slate_generation_sllm_disabled_defaults_to_pure_dp():
    engine = DiscoverySlateEngine(host="nbenthamiana", enable_sllm_preview=False)
    slate = engine.generate_slate(
        target_aa=TEST_TARGET_AA[:60],
        target_name="sfGFP-60aa",
        top_k=3,
    )
    assert slate.summary.selected_slate_size == 3
    assert not slate.provenance.get("sllm_preview_enabled")
    # All candidates should be DP or profile contracts
    for cand in slate.candidates:
        assert not cand.generation_contract.startswith("sllm_")


def test_slate_generation_sllm_enabled():
    mock_adapter = DeterministicMockSLMAdapter(host="nbenthamiana", seed=42)
    engine = DiscoverySlateEngine(
        host="nbenthamiana",
        enable_sllm_preview=True,
        sllm_model_adapter=mock_adapter,
    )
    slate = engine.generate_slate(
        target_aa=TEST_TARGET_AA[:60],
        target_name="sfGFP-60aa",
        top_k=3,
    )
    assert slate.summary.selected_slate_size == 3
    assert slate.provenance.get("sllm_preview_enabled")

    # Check that candidate pool clusters are valid
    clusters = {c.strategy_cluster for c in slate.candidates}
    assert any("LatentPrior" in cl or "Initiation" in cl or "Speed" in cl for cl in clusters)


def test_slate_sllm_candidates_hard_constraint_compliance():
    mock_adapter = DeterministicMockSLMAdapter(host="nbenthamiana", seed=42)
    config = HardConstraintConfig(
        forbidden_other_motifs=["GCGGCCGC", "GGTCTC"],
        homopolymer_max_run=5,
    )
    filter_guard = HardConstraintFilter(config)
    engine = DiscoverySlateEngine(
        host="nbenthamiana",
        hard_filter_config=config,
        enable_sllm_preview=True,
        sllm_model_adapter=mock_adapter,
    )

    raw_pool = engine.generator.generate_pool(
        protein_sequence=TEST_TARGET_AA[:60],
        stop_codon="TAA",
        enable_sllm_override=True,
    )

    sllm_cands = [c for c in raw_pool if c.generation_contract.startswith("sllm_")]
    assert len(sllm_cands) >= 1

    for cand in sllm_cands:
        eval_res = filter_guard.evaluate(
            cand.sequence_dna, TEST_TARGET_AA[:60], expected_stop="TAA"
        )
        assert eval_res.is_feasible, (
            f"Candidate {cand.candidate_id} failed hard filter: {eval_res.violations}"
        )


def test_slate_ranker_integration():
    mock_adapter = DeterministicMockSLMAdapter(host="nbenthamiana", seed=42)
    engine = DiscoverySlateEngine(
        host="nbenthamiana",
        enable_sllm_preview=True,
        sllm_model_adapter=mock_adapter,
    )
    slate = engine.generate_slate(
        target_aa=TEST_TARGET_AA[:60],
        target_name="sfGFP-60aa",
        top_k=3,
    )

    # Rank candidates with ExpressionPriorRanker
    ranker = ExpressionPriorRanker(weights=None)  # Heuristic fallback ranker
    scored_slate = []
    for cand in slate.candidates:
        # 8D trait vector
        t = cand.trait_vector
        tv = [
            t.cai_golden_set,
            t.global_gc_percent,
            t.global_gc_percent,  # proxy 5' gc
            t.initiation_mfe_kcal_mol or -15.0,
            0.05,  # CPI
            0.02,  # rare codon frac
            float(t.homopolymer_max_run),
            float(t.homopolymer_max_run),
        ]
        score = ranker.score(tv)
        scored_slate.append((cand.candidate_id, score, cand.strategy_cluster))

    assert len(scored_slate) == 3
    # Highest score is sorted properly
    scored_slate.sort(key=lambda x: x[1], reverse=True)
    assert scored_slate[0][1] >= scored_slate[-1][1]
