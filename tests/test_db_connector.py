"""PyTest suite for upgraded Context-Aware FactorForgeDBConnector."""

import tempfile
from factorforge.db.connector import FactorForgeDBConnector


def test_db_connector_save_candidate_and_evaluations():
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tmp:
        db_path = tmp.name

    db = FactorForgeDBConnector(db_path=db_path)
    res = db.save_candidate_with_evaluations(
        candidate_data={
            "engine": "lm",
            "model_version": "v3.5.0-SynCodonLM-V2",
            "optimized_sequence": "ATGGCTAAATGGTAA",
            "cai": 0.84,
            "gc_percent": 41.2,
            "type2is_clean": True,
        },
        evaluations=[
            {"constraint_code": "AA_IDENTITY", "status": "PASS", "observed_value": 1.0},
            {"constraint_code": "TYPE_IIS_CLEAN", "status": "PASS", "observed_value": 0},
        ],
        construct_id="CD9-MOCLO-TEST",
    )

    assert res["package_id"] >= 1
    assert res["candidate_id"] >= 1
    assert len(res["frozen_manifest_hash"]) == 64
