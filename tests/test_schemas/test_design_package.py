"""DesignPackage schema tests."""

from api.optimize import handler
from factorforge.schemas import DesignPackage


def minimal_payload() -> dict:
    return {
        "construct_id": "CF-20260525-000000",
        "created_at": "2026-05-25T00:00:00Z",
        "cds_design": {
            "engine": "factorforge_cds",
            "product_version": "3.1.3",
            "host_profile": "nbenthamiana",
            "objective": "feasibility_best",
            "cai": 0.91,
            "gc_percent": 52.4,
        },
        "provenance": {
            "input_sequence_hash": "sha256:abc",
            "output_cds_hash": "sha256:def",
            "parameter_hash": "sha256:ghi",
        },
    }


def test_minimal_construction():
    """Can create a package with the minimum required fields."""
    pkg = DesignPackage(**minimal_payload())

    assert pkg.construct_id == "CF-20260525-000000"
    assert pkg.design_package_version == "1.0"


def test_optional_target_null():
    """Target can be omitted or left null."""
    payload = minimal_payload()
    payload["construct_id"] = "CF-20260525-000001"

    pkg = DesignPackage(**payload)

    assert pkg.target is None or pkg.target.name is None


def test_schema_export():
    """JSON Schema export works."""
    schema = DesignPackage.model_json_schema()

    assert "properties" in schema
    assert "construct_id" in schema["properties"]


def test_wet_lab_feedback_default():
    """Wet-lab feedback defaults to pending with no submissions."""
    payload = minimal_payload()
    payload["construct_id"] = "CF-20260525-000002"

    pkg = DesignPackage(**payload)

    assert pkg.wet_lab_feedback.status == "pending"
    assert pkg.wet_lab_feedback.submissions == []


def test_api_response_validates_as_design_package():
    """The optimize API success response can validate as a DesignPackage."""
    request_handler = handler.__new__(handler)

    response = request_handler.optimize_feasibility_best(
        sequence="MKT",
        profile="balanced",
        host_profile="nbenthamiana",
        constraints={"gc_min": 40.0, "gc_max": 55.0},
        kozak=False,
        dinuc=False,
        return_candidates=True,
    )

    pkg = DesignPackage.model_validate(response)

    assert pkg.construct_id.startswith("CF-")
    assert pkg.cds_design.host_profile == "nbenthamiana"
    assert pkg.cds_design.profile == "feasibility_best"
