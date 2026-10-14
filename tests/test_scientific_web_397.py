"""Synthetic regression checks for observation truthfulness and web integration."""

import math
import sys
from types import SimpleNamespace

import pytest

from factorforge.guards.initiation_mfe import evaluate_initiation_mfe
from factorforge.rules.registry import RuleRegistry
from factorforge.sops.config import SopProfile
from factorforge.structure.auto_annotator import recommend_structure_pause_regions
from tests.test_auto_annotator import _generate_synthetic_pdb
from tests.api.test_structure_web import PROTEIN, post, request
from api.optimize import handler


def test_flanked_construct_is_not_mislabeled_as_sense_cds():
    result = handler.attach_scientific_observations(object.__new__(handler),
        {"optimized_sequence":"GCCACCATGGCTTGA"}, expected_protein="MA")
    review = result["scientific_review"]
    assert review["scope"] == "not_established"
    assert review["rna_risks"]["summary_status"] == "NOT_AVAILABLE"
    assert review["initiation_mfe"]["status"] == "NOT_AVAILABLE"
    assert review["initiation_mfe"]["mfe_kcal_mol"] is None


def required_policy():
    return {
        "$schema": "factorforge-sop-v1",
        "profile_id": "synthetic_required_mfe",
        "sop_name": "Synthetic required observation",
        "version": "1.0.0",
        "status": "EXAMPLE_TEMPLATE",
        "default_enforcement": "IGNORE",
        "unknown_rule_policy": "ERROR",
        "rules": {"policy.initiation_mfe.v1": "HARD_FAIL"},
    }


@pytest.mark.parametrize("value", [None, float("nan"), float("inf"), True])
def test_required_unavailable_blocks(value):
    registry = RuleRegistry(sop_profile=SopProfile.from_dict(required_policy()))
    result = registry.evaluate_sequence("ATGGCTTGA", context={"initiation_mfe": value})
    assert not result["all_passed"]
    assert result["evaluation_status"] == "INCOMPLETE"
    assert result["required_unavailable_count"] == 1
    assert result["hard_fails"][0]["result"]["passed"] is None


@pytest.mark.parametrize(
    "energy,structure",
    [(float("nan"), "." * 48), (float("inf"), "." * 48), (-2, "(" * 48), (-2, ".")],
)
def test_invalid_fold_never_becomes_available(monkeypatch, energy, structure):
    monkeypatch.setitem(
        sys.modules,
        "RNA",
        SimpleNamespace(
            md=lambda: SimpleNamespace(),
            fold_compound=lambda *a: SimpleNamespace(mfe=lambda: (structure, energy)),
        ),
    )
    result = evaluate_initiation_mfe("ATG" + "GCA" * 20)
    assert result.status == "ERROR"
    assert result.mfe_kcal_mol is None


def test_real_per_call_temperature_and_start_anchor():
    RNA = pytest.importorskip("RNA")
    before = RNA.cvar.temperature
    result = evaluate_initiation_mfe(
        "ATG" + "GCA" * 20, upstream_context="C" * 80, temperature_c=42
    )
    assert result.status == "AVAILABLE"
    assert math.isfinite(result.mfe_kcal_mol)
    assert (result.window_start, result.window_end) == (-4, 44)
    assert result.model_details["temperature"] == 42
    assert RNA.cvar.temperature == before


def test_absent_tool_and_legacy_api_are_explicit(monkeypatch):
    monkeypatch.setitem(sys.modules, "RNA", None)
    assert evaluate_initiation_mfe("ATGGCTTGA").status == "NOT_AVAILABLE"
    monkeypatch.setitem(sys.modules, "RNA", SimpleNamespace(fold=lambda x: ("." * len(x), -1)))
    assert evaluate_initiation_mfe("ATGGCTTGA").status == "NOT_AVAILABLE"


@pytest.mark.parametrize(
    "source,missing_ca,status",
    [
        (None, False, "NOT_AVAILABLE"),
        (False, False, "NO_CANDIDATES"),
        (True, True, "NOT_AVAILABLE"),
    ],
)
def test_no_invented_recommendation(source, missing_ca, status):
    content = _generate_synthetic_pdb(PROTEIN)
    if missing_ca:
        content = content.replace("  CA ", "  N  ")
    result = recommend_structure_pause_regions(PROTEIN, content, is_alphafold=source)
    assert result.status == status
    assert not result.recommended_positions
    assert not result.reviewed


def test_recommend_action_does_not_approve_or_generate():
    data = request(action="recommend_structure")
    data["structure_context"].update(
        content=_generate_synthetic_pdb(PROTEIN, [90] * 4 + [40] * 6 + [90] * 9), is_alphafold=True
    )
    status, result = post(data)
    assert status == 200, result
    assert not result["recommendation"]["reviewed"]
    assert result["recommendation"]["recommended_positions"]
    assert "candidates" not in result
    assert result["approval_status"] == "not_collected"


@pytest.mark.parametrize("paired", [False, True])
def test_actual_observations_attached_to_web_result(paired):
    status, result = post(request() if paired else {"sequence": PROTEIN})
    assert status == 200, result
    assert result["scientific_review"]["rna_risks"]["evidence_class"] == "SEQUENCE_SCAN"
    mfe = result["scientific_review"]["initiation_mfe"]
    assert mfe["status"] == "AVAILABLE"
    assert math.isfinite(mfe["mfe_kcal_mol"])
    if paired:
        assert all(c["scientific_review"] for c in result["candidates"])


@pytest.mark.parametrize("paired", [False, True])
def test_required_absent_tool_is_http_422(monkeypatch, paired):
    monkeypatch.setitem(sys.modules, "RNA", None)
    data = request() if paired else {"sequence": PROTEIN}
    data["sop_profile"] = required_policy()
    status, result = post(data)
    assert status == 422, result
    assert result["error_code"] == "REQUIRED_CHECK_UNAVAILABLE"
    assert not result["success"] and "candidates" not in result


@pytest.mark.parametrize("paired", [False, True])
def test_measured_required_failure_does_not_display_pass(monkeypatch, paired):
    monkeypatch.setitem(
        sys.modules,
        "RNA",
        SimpleNamespace(
            md=lambda: SimpleNamespace(),
            fold_compound=lambda seq, md: SimpleNamespace(mfe=lambda: ("." * len(seq), -99)),
        ),
    )
    data = request() if paired else {"sequence": PROTEIN}
    data["sop_profile"] = required_policy()
    status, result = post(data)
    assert status == 200, result
    assert result["policy_decision"] == "BLOCK"
    assert result["automated_decision"] == "FAIL"
    assert result["decision_summary"]["required_failure_count"] >= 1


@pytest.mark.parametrize(
    "field,value", [("format", None), ("model", True), ("is_alphafold", "yes")]
)
def test_recommend_invalid_parameters_are_400(field, value):
    data = request(action="recommend_structure")
    data["structure_context"][field] = value
    assert post(data)[0] == 400
