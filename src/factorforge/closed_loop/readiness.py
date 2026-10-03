"""Readiness State Machine & Active Learner Safety Gate."""

from __future__ import annotations

from enum import Enum
from typing import Dict, List, Optional
from pydantic import BaseModel
from factorforge.closed_loop.contracts import MeasurementRecord, ExperimentRecord


class ReadinessStatus(str, Enum):
    """5-State Closed-Loop Readiness State Machine."""

    NO_DATA = "NO_DATA"
    DATA_COLLECTION = "DATA_COLLECTION"
    EXPLORATORY = "EXPLORATORY"
    MODEL_ELIGIBLE = "MODEL_ELIGIBLE"
    MODEL_ACTIVE = "MODEL_ACTIVE"


# Provisional software thresholds (Not statistical guarantees; subject to future power analysis)
PROVISIONAL_READINESS_THRESHOLDS: Dict[str, int] = {
    "min_bioreps_collection": 3,
    "min_sets_collection": 2,
    "min_bioreps_exploratory": 6,
    "min_sets_exploratory": 4,
}


class ReadinessReport(BaseModel):
    """Evaluation summary of dataset readiness for machine learning."""

    status: ReadinessStatus
    total_records: int
    synthetic_test_records: int
    real_empirical_records: int
    distinct_construct_sets: int
    distinct_biological_replicates: int
    distinct_batches: int
    model_training_permitted: bool
    status_explanation: str
    thresholds_applied: Dict[str, int] = PROVISIONAL_READINESS_THRESHOLDS


class ReadinessEvaluator:
    """Evaluates readiness of experimental records strictly guarding against premature ML claims."""

    @staticmethod
    def evaluate(
        measurements: List[MeasurementRecord],
        experiments: Optional[Dict[str, object]] = None,
        thresholds: Optional[Dict[str, int]] = None,
    ) -> ReadinessReport:
        active_thresholds = thresholds or PROVISIONAL_READINESS_THRESHOLDS
        synthetic_count = sum(1 for m in measurements if m.is_synthetic)
        real_records = [m for m in measurements if not m.is_synthetic and m.qc_status == "PASS"]
        real_count = len(real_records)

        distinct_sets = len(set(m.construct_set_id for m in real_records))

        # Hardening 1: Compound biological replicate key to avoid ID collision across constructs/experiments
        distinct_bioreps = len(set(
            (m.experiment_id, m.construct_set_id, m.biological_replicate_id)
            for m in real_records
        ))

        # Hardening 2: Resolve true batch_id from ExperimentRecord if provided
        if experiments:
            batch_ids = set()
            for m in real_records:
                exp = experiments.get(m.experiment_id)
                if exp:
                    b_id = getattr(exp, "batch_id", None) if hasattr(exp, "batch_id") else exp.get("batch_id")
                    batch_ids.add(b_id or m.experiment_id)
                else:
                    batch_ids.add(m.experiment_id)
            distinct_batches = len(batch_ids)
        else:
            distinct_batches = len(set(m.experiment_id for m in real_records))

        if real_count == 0:
            status = ReadinessStatus.NO_DATA
            explanation = "No empirical wet-lab measurements recorded yet. Learning models remain offline."
            permitted = False
        elif (
            distinct_bioreps < active_thresholds["min_bioreps_collection"]
            or distinct_sets < active_thresholds["min_sets_collection"]
        ):
            status = ReadinessStatus.DATA_COLLECTION
            explanation = (
                f"Empirical data collection in progress ({real_count} records). "
                "Insufficient biological replicates or construct sets for causal modeling."
            )
            permitted = False
        elif (
            distinct_bioreps < active_thresholds["min_bioreps_exploratory"]
            or distinct_sets < active_thresholds["min_sets_exploratory"]
        ):
            status = ReadinessStatus.EXPLORATORY
            explanation = "Exploratory evidence available for descriptive review, but insufficient for active learning optimization."
            permitted = False
        else:
            status = ReadinessStatus.MODEL_ELIGIBLE
            explanation = (
                "Dataset meets provisional software coverage thresholds for formal model training approval review."
            )
            permitted = True

        return ReadinessReport(
            status=status,
            total_records=len(measurements),
            synthetic_test_records=synthetic_count,
            real_empirical_records=real_count,
            distinct_construct_sets=distinct_sets,
            distinct_biological_replicates=distinct_bioreps,
            distinct_batches=distinct_batches,
            model_training_permitted=permitted,
            status_explanation=explanation,
            thresholds_applied=active_thresholds,
        )


class ModelRecommendationGate:
    """Safety gate preventing premature or hallucinated model ratio recommendations."""

    @staticmethod
    def query_next_ratio(
        report: ReadinessReport,
        explicit_human_approval_token: Optional[str] = None,
    ) -> Dict[str, object]:
        """Query for model-recommended next-generation codon distribution."""
        if not report.model_training_permitted or report.status != ReadinessStatus.MODEL_ACTIVE:
            return {
                "recommendation_status": "UNAVAILABLE",
                "current_readiness": report.status.value,
                "reason": (
                    f"Model-guided ratio optimization is unavailable (Current state: {report.status.value}). "
                    "Empirical wet-lab measurements must be accumulated, statistically reviewed, "
                    "and explicitly approved by human sign-off before activating the Active Learner."
                ),
                "suggested_profile": None,
                "mock_data_detected": report.synthetic_test_records > 0,
            }

        # Hardening 5: Explicit human sign-off token required to transition from MODEL_ELIGIBLE to ACTIVE
        if not explicit_human_approval_token:
            return {
                "recommendation_status": "APPROVAL_REQUIRED",
                "current_readiness": ReadinessStatus.MODEL_ACTIVE.value,
                "reason": "Dataset is MODEL_ACTIVE but requires verified human sign-off token to trigger candidate generation.",
                "suggested_profile": None,
                "mock_data_detected": False,
            }

        return {
            "recommendation_status": "ACTIVE",
            "current_readiness": report.status.value,
            "suggested_profile": "learned_active_learner_v1",
            "approved_token": explicit_human_approval_token,
            "reason": "Empirically validated and human-approved closed-loop recommendation available.",
        }
