# factorforge/tests/learning/test_expression_prior_ranker.py
from factorforge.learning.ranker import (
    ExpressionPriorRanker,
    RankerConfig,
)
from factorforge.registry.model_registry import ModelRegistry


def test_ranker_fit_and_inference():
    # Construct synthetic pairs where higher 5' MFE (index 3) and higher CAI (index 0) is better
    train_pairs = [
        {
            "pair_id": "P1",
            "trait_vector_a": [0.82, 44.0, 37.5, -6.2, 0.02, 0.04, 3.0, 2.0],  # H1: Open 5' MFE
            "trait_vector_b": [0.88, 45.2, 40.1, -12.4, 0.01, 0.03, 4.0, 3.0],  # H0: Stable MFE
            "preference_label": 1,
            "tier_weight": 1.0,
        },
        {
            "pair_id": "P2",
            "trait_vector_a": [0.88, 45.2, 40.1, -12.4, 0.01, 0.03, 4.0, 3.0],  # H0
            "trait_vector_b": [0.85, 48.0, 44.0, -18.5, 0.01, 0.02, 3.0, 3.0],  # H2
            "preference_label": 1,
            "tier_weight": 0.8,
        },
    ]

    ranker = ExpressionPriorRanker(config=RankerConfig(max_epochs=150, learning_rate=0.1))
    weights = ranker.fit(train_pairs=train_pairs)

    assert ranker.is_trained()
    assert weights.train_pairwise_acc == 1.0

    # Score candidates: H1 should score highest
    c_h1 = {"id": "H1", "trait_vector": [0.82, 44.0, 37.5, -6.2, 0.02, 0.04, 3.0, 2.0]}
    c_h0 = {"id": "H0", "trait_vector": [0.88, 45.2, 40.1, -12.4, 0.01, 0.03, 4.0, 3.0]}
    c_h2 = {"id": "H2", "trait_vector": [0.85, 48.0, 44.0, -18.5, 0.01, 0.02, 3.0, 3.0]}

    ranked = ranker.rank_candidates([c_h2, c_h0, c_h1])
    assert ranked[0]["id"] == "H1"
    assert ranked[0]["expression_prior_score"] > ranked[1]["expression_prior_score"]


def test_model_registry_lifecycle(tmp_path):
    registry = ModelRegistry(registry_dir=str(tmp_path))

    # Stage candidate
    rec = registry.register_candidate(
        model_id="M-001",
        model_version="0.1.0",
        weights_sha256="sha_test_weights",
        training_snapshot_id="SNAP-01",
        validation_metrics={"pairwise_acc": 0.85},
    )
    assert rec.status == "STAGED_CANDIDATE"
    assert registry.get_active_model() is None

    # Promote candidate (simulated vs actual)
    sim_promoted = registry.promote_model(
        model_id="M-001",
        promoted_by="Dr. Chemist",
        rationale="Simulation benchmark verified",
        validation_bundle_sha256="sha_bundle_01",
        is_simulated=True,
    )
    assert sim_promoted.status == "SIMULATION_ONLY"
    assert registry.get_active_model() is None

    # Actual production promotion
    prod_promoted = registry.promote_model(
        model_id="M-001",
        promoted_by="Lead Bioprocess Engineer",
        rationale="Empirical PlantForm validation complete",
        validation_bundle_sha256="sha_bundle_01",
        is_simulated=False,
    )
    assert prod_promoted.status == "ACTIVE_PRODUCTION"
    assert registry.get_active_model().model_id == "M-001"
