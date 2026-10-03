"""Unit tests for Job 351: Codon Closed-Loop Evidence Foundation."""

import os
import shutil
import tempfile
from pathlib import Path
import pytest

from factorforge.closed_loop import (
    PolicyRecord,
    ConstructRecord,
    ConstructSetRecord,
    ExperimentRecord,
    MeasurementRecord,
    ReadinessStatus,
    ReadinessEvaluator,
    ModelRecommendationGate,
    EvidenceLedger,
)


@pytest.fixture
def temp_ledger(tmp_path):
    ledger = EvidenceLedger(ledger_dir=tmp_path / "test_ledger")
    return ledger


class TestJob351ClosedLoopFoundation:
    def test_initial_readiness_is_no_data(self, temp_ledger):
        report = temp_ledger.get_readiness()
        assert report.status == ReadinessStatus.NO_DATA
        assert report.real_empirical_records == 0
        assert report.model_training_permitted is False

        # Model Recommendation Gate must strictly reject recommendation requests
        rec = ModelRecommendationGate.query_next_ratio(report)
        assert rec["recommendation_status"] == "UNAVAILABLE"
        assert rec["suggested_profile"] is None
        assert "NO_DATA" in rec["reason"]

    def test_synthetic_measurement_is_strictly_isolated(self, temp_ledger):
        # Ingest a synthetic/mock measurement
        mock_record = MeasurementRecord(
            measurement_id="MOCK-EXP-001-YIELD",
            experiment_id="EXP-MOCK-ROUND1",
            construct_set_id="SET-BALANCED-W0",
            biological_replicate_id="BIO-REP-01",
            technical_replicate_id="TECH-REP-01",
            measurement_type="protein_yield",
            value=1.20,
            unit="mg/L",
            qc_status="PASS",
            is_synthetic=True,  # Synthetic test flag
            tag="test_only",
            notes="Pipeline plumbing verification mock fixture",
        )
        temp_ledger.record_measurement(mock_record)

        # Readiness must still be NO_DATA for real training
        report = temp_ledger.get_readiness()
        assert report.total_records == 1
        assert report.synthetic_test_records == 1
        assert report.real_empirical_records == 0
        assert report.status == ReadinessStatus.NO_DATA
        assert report.model_training_permitted is False

        # Recommender must detect mock data and refuse active learning
        rec = ModelRecommendationGate.query_next_ratio(report)
        assert rec["recommendation_status"] == "UNAVAILABLE"
        assert rec["mock_data_detected"] is True

    def test_end_to_end_lineage_trace(self, temp_ledger):
        # 1. Register Policy
        pol = PolicyRecord(
            policy_id="doug_balanced_v1",
            canonical_sha256="abc123canonicalhash",
            target_distribution={"R": {"AGA": 0.5, "AGG": 0.5}},
            region_definition="mature_chain",
            feedback_source_id="doug_2026_10_02",
        )
        temp_ledger.register_policy(pol)

        # 2. Register LC & HC Constructs
        c_lc = ConstructRecord(
            construct_id="PF-HUMIRA-LC-BALANCED-W0",
            cds_sha256="hash_cds_lc",
            protein_sha256="hash_prot_lc",
            chain="LC",
            watermark_pair_id="HUMIRA-BALANCED-PAIR",
            watermark_applied=False,
            achieved_distribution={"R": {"AGA": 0.556, "AGG": 0.444}},
            policy_id="doug_balanced_v1",
        )
        c_hc = ConstructRecord(
            construct_id="PF-HUMIRA-HC-BALANCED-W0",
            cds_sha256="hash_cds_hc",
            protein_sha256="hash_prot_hc",
            chain="HC",
            watermark_pair_id="HUMIRA-BALANCED-PAIR",
            watermark_applied=False,
            achieved_distribution={"R": {"AGA": 0.50, "AGG": 0.50}},
            policy_id="doug_balanced_v1",
        )
        temp_ledger.register_construct(c_lc)
        temp_ledger.register_construct(c_hc)

        # 3. Register Co-expressed Construct Set
        cset = ConstructSetRecord(
            construct_set_id="SET-HUMIRA-BALANCED-W0",
            light_construct_id="PF-HUMIRA-LC-BALANCED-W0",
            heavy_construct_id="PF-HUMIRA-HC-BALANCED-W0",
            watermark_applied=False,
            description="Co-expressed Adalimumab balanced W0 twin pair",
        )
        temp_ledger.register_construct_set(cset)

        # 4. Register Experiment
        exp = ExperimentRecord(
            experiment_id="EXP-PLANTFORM-ROUND1",
            construct_set_id="SET-HUMIRA-BALANCED-W0",
            host="Nicotiana benthamiana",
            batch_id="BATCH-2026-10-OCT",
            experiment_date="2026-10-15",
        )
        temp_ledger.register_experiment(exp)

        # 5. Record Measurement
        meas = MeasurementRecord(
            measurement_id="MEAS-EXP1-YIELD-001",
            experiment_id="EXP-PLANTFORM-ROUND1",
            construct_set_id="SET-HUMIRA-BALANCED-W0",
            biological_replicate_id="BIO-REP-01",
            technical_replicate_id="TECH-REP-01",
            measurement_type="protein_yield",
            value=14.5,
            unit="mg/L",
            qc_status="PASS",
            is_synthetic=False,
        )
        temp_ledger.record_measurement(meas)

        # Verify complete unbroken lineage trace
        lineage = temp_ledger.verify_lineage("MEAS-EXP1-YIELD-001")
        assert lineage["lineage_intact"] is True
        assert lineage["experiment"]["experiment_id"] == "EXP-PLANTFORM-ROUND1"
        assert lineage["construct_set"]["construct_set_id"] == "SET-HUMIRA-BALANCED-W0"
        assert lineage["light_construct"]["construct_id"] == "PF-HUMIRA-LC-BALANCED-W0"
        assert lineage["heavy_construct"]["construct_id"] == "PF-HUMIRA-HC-BALANCED-W0"
        assert lineage["policy"]["policy_id"] == "doug_balanced_v1"

    def test_gcs_cloud_sync_plumbing(self, temp_ledger):
        # Register a test policy and check GCS sync capability
        pol = PolicyRecord(
            policy_id="test_gcs_policy",
            canonical_sha256="test_hash",
            target_distribution={"R": {"AGA": 0.5, "AGG": 0.5}},
            feedback_source_id="test",
        )
        temp_ledger.register_policy(pol)

        key_path = r"C:\Users\munky\.gcp\factorforge-runner-key.json"
        if os.path.exists(key_path):
            synced = temp_ledger.sync_to_gcs(
                bucket_name="factorforge-fair-conduit",
                destination_folder="validationhub/_test_sync",
                gcp_key_path=key_path,
            )
            assert "policies.jsonl" in synced
            assert synced["policies.jsonl"].startswith("gs://factorforge-fair-conduit/")
            assert "ledger_manifest.json" in synced

    def test_compound_biological_replicate_counting(self):
        # Two different constructs using the same biological replicate label (e.g. BIO-REP-01)
        # must be counted as 2 distinct biological observations, not 1!
        m1 = MeasurementRecord(
            measurement_id="M1",
            experiment_id="EXP1",
            construct_set_id="SET-A",
            biological_replicate_id="BIO-REP-01",
            technical_replicate_id="T1",
            measurement_type="protein_yield",
            value=10.0,
            unit="mg/L",
            qc_status="PASS",
            is_synthetic=False,
        )
        m2 = MeasurementRecord(
            measurement_id="M2",
            experiment_id="EXP1",
            construct_set_id="SET-B",
            biological_replicate_id="BIO-REP-01",
            technical_replicate_id="T1",
            measurement_type="protein_yield",
            value=12.0,
            unit="mg/L",
            qc_status="PASS",
            is_synthetic=False,
        )
        report = ReadinessEvaluator.evaluate([m1, m2])
        assert report.distinct_biological_replicates == 2
        assert report.distinct_construct_sets == 2

    def test_referential_integrity_violation_in_lineage(self, temp_ledger):
        # Register Experiment pointing to SET-A
        exp = ExperimentRecord(
            experiment_id="EXP-MISMATCH",
            construct_set_id="SET-A",
            host="Nicotiana benthamiana",
            batch_id="BATCH-01",
            experiment_date="2026-10-03",
        )
        temp_ledger.register_experiment(exp)

        # Register measurement pointing to SET-B (Mismatch!)
        meas = MeasurementRecord(
            measurement_id="MEAS-MISMATCH",
            experiment_id="EXP-MISMATCH",
            construct_set_id="SET-B",
            biological_replicate_id="BIO-01",
            technical_replicate_id="T1",
            measurement_type="protein_yield",
            value=5.0,
            unit="mg/L",
            qc_status="PASS",
            is_synthetic=False,
        )
        temp_ledger.record_measurement(meas)

        lineage = temp_ledger.verify_lineage("MEAS-MISMATCH")
        assert lineage["lineage_intact"] is False
        assert any("Referential mismatch" in v or "does not exist" in v for v in lineage["violations"])

    def test_model_active_requires_explicit_human_approval(self):
        # Create a report with MODEL_ACTIVE status
        report = ReadinessEvaluator.evaluate([])
        # Artificially test MODEL_ACTIVE gate
        report.status = ReadinessStatus.MODEL_ACTIVE
        report.model_training_permitted = True

        # Without human approval token: must return APPROVAL_REQUIRED
        gate_res = ModelRecommendationGate.query_next_ratio(report)
        assert gate_res["recommendation_status"] == "APPROVAL_REQUIRED"

        # With explicit approval token: activates
        gate_active = ModelRecommendationGate.query_next_ratio(
            report, explicit_human_approval_token="HUMAN-SIGNOFF-2026-OCT-PLANTFORM"
        )
        assert gate_active["recommendation_status"] == "ACTIVE"
        assert gate_active["approved_token"] == "HUMAN-SIGNOFF-2026-OCT-PLANTFORM"
