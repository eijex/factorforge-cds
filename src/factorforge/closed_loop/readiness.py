"""Readiness State Machine & Active Learner Safety Gate."""

from __future__ import annotations

from enum import Enum
from typing import Dict, List, Optional
from pydantic import BaseModel
from factorforge.closed_loop.contracts import MeasurementRecord


class ReadinessStatus(str, Enum):
    """5-State Closed-Loop Readiness State Machine."""

    NO_DATA = "NO_DATA"
    DATA_COLLECTION = "DATA_COLLECTION"
    EXPLORATORY = "EXPLORATORY"
    MODEL_ELIGIBLE = "MODEL_ELIGIBLE"
    MODEL_ACTIVE = "MODEL_ACTIVE"


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


class ReadinessEvaluator:
    """Evaluates readiness of experimental records strictly guarding against premature ML claims."""

    @staticmethod
    def evaluate(measurements: List[MeasurementRecord]) -> ReadinessReport:
        synthetic_count = sum(1 for m in measurements if m.is_synthetic)
        real_records = [m for m in measurements if not m.is_synthetic and m.qc_status == "PASS"]
        real_count = len(real_records)

        distinct_sets = len(set(m.construct_set_id for m in real_records))
        distinct_bioreps = len(set(m.biological_replicate_id for m in real_records))
        distinct_batches = len(set(m.experiment_id for m in real_records))

        if real_count == 0:
            status = ReadinessStatus.NO_DATA
            explanation = "No empirical wet-lab measurements recorded yet. Learning models remain offline."
            permitted = False
        elif distinct_bioreps < 3 or distinct_sets < 2:
            status = ReadinessStatus.DATA_COLLECTION
            explanation = f"Empirical data collection in progress ({real_count} records). Insufficient biological replicates for causal modeling."
            permitted = False
        elif distinct_bioreps < 6 or distinct_sets < 4:
            status = ReadinessStatus.EXPLORATORY
            explanation = "Exploratory evidence available for descriptive review, but insufficient for active learning optimization."
            permitted = False
        else:
            status = ReadinessStatus.MODEL_ELIGIBLE
            explanation = "Dataset meets minimum replicate and construct coverage for formal model training approval."
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
        )


class ModelRecommendationGate:
    """Safety gate preventing premature or hallucinated model ratio recommendations."""

    @staticmethod
    def query_next_ratio(report: ReadinessReport) -> Dict[str, object]:
        """Query for model-recommended next-generation codon distribution."""
        if not report.model_training_permitted or report.status != ReadinessStatus.MODEL_ACTIVE:
            return {
                "recommendation_status": "UNAVAILABLE",
                "current_readiness": report.status.value,
                "reason": (
                    f"Model-guided ratio optimization is unavailable (Current state: {report.status.value}). "
                    "Empirical wet-lab measurements must be accumulated and verified before activating the Active Learner."
                ),
                "suggested_profile": None,
                "mock_data_detected": report.synthetic_test_records > 0,
            }

        # In production, once MODEL_ACTIVE is explicitly approved:
        return {
            "recommendation_status": "ACTIVE",
            "current_readiness": report.status.value,
            "suggested_profile": "learned_active_learner_v1",
        }
