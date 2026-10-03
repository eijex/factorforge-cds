"""Synthetic regression fixtures; no collaborator data or private targets."""

import pytest

from factorforge.closed_loop.contracts import MeasurementRecord
from factorforge.closed_loop.readiness import (
    ModelRecommendationGate,
    ReadinessEvaluator,
    ReadinessStatus,
)
from factorforge.codon.policy import FACTORFORGE_RECOMMENDED_PRESET, CodonDistributionPolicy
from factorforge.codon.qa import InSilicoQA


@pytest.mark.parametrize(
    "weights", [{"AGA": float("nan"), "AGG": 0.5}, {"AGA": float("inf")}, {"aga": 1.0}]
)
def test_invalid_policy_values(weights):
    with pytest.raises(ValueError):
        CodonDistributionPolicy(policy_id="synthetic", distributions={"R": weights})


@pytest.mark.parametrize("cds,protein", [("ATGTAAA", "M"), ("", ""), ("ATGN", "M")])
def test_partial_frame_never_passes(cds, protein):
    assert not InSilicoQA.evaluate(cds, protein).overall_passed


def test_exploratory_preset_not_a_learned_model():
    assert FACTORFORGE_RECOMMENDED_PRESET.source_type == "custom"


def measurement(**updates):
    values = {
        "measurement_id": "synthetic",
        "experiment_id": "e",
        "construct_set_id": "s",
        "biological_replicate_id": "b",
        "technical_replicate_id": "t",
        "measurement_type": "protein_yield",
        "value": 0.0,
        "unit": "mg/L",
    }
    return MeasurementRecord(**(values | updates))


def test_missing_and_test_only_records_are_not_empirical():
    for m in [
        measurement(value=None),
        measurement(tag="test_only"),
        measurement(is_synthetic=True),
    ]:
        assert ReadinessEvaluator.evaluate([m]).real_empirical_records == 0
    assert ReadinessEvaluator.evaluate([measurement(value=0.0)]).real_empirical_records == 1
    with pytest.raises(ValueError):
        measurement(value=float("nan"))


def test_arbitrary_approval_token_cannot_invent_model():
    report = ReadinessEvaluator.evaluate([])
    report.status = ReadinessStatus.MODEL_ACTIVE
    report.model_training_permitted = True
    result = ModelRecommendationGate.query_next_ratio(report, "synthetic-token")
    assert result["recommendation_status"] == "UNAVAILABLE"
    assert result["suggested_profile"] is None
    assert "approved_token" not in result


def test_conflicting_identity_and_unlinked_measurements_do_not_count(tmp_path):
    from factorforge.closed_loop.ledger import EvidenceLedger

    ledger = EvidenceLedger(tmp_path)
    record = measurement()
    ledger.record_measurement(record)
    ledger.record_measurement(record)
    assert len(ledger.list_measurements()) == 1
    assert ledger.get_readiness().real_empirical_records == 0
    with pytest.raises(ValueError, match="Conflicting identity"):
        ledger.record_measurement(record.model_copy(update={"value": 1.0}))


def test_short_or_mismatched_baseline_rejected():
    from factorforge.codon.policy import DOUG_BALANCED_PRESET
    from factorforge.codon.twin_generator import TwinCandidateGenerator

    generator = TwinCandidateGenerator()
    with pytest.raises(ValueError, match="Baseline"):
        generator.generate("LC", "M", "ATGAAATAA", DOUG_BALANCED_PRESET)
    with pytest.raises(ValueError, match="frozen region"):
        generator.generate("LC", "M", "ATGTAA", DOUG_BALANCED_PRESET)


def test_allocator_rejects_invalid_direct_calls():
    from factorforge.codon.allocator import DeterministicAllocator

    with pytest.raises(ValueError):
        DeterministicAllocator.allocate("R", -1, {"AGA": 1.0})
    with pytest.raises(ValueError):
        DeterministicAllocator.allocate("R", 9, {"AGA": float("nan")})
