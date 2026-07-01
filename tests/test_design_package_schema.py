"""JSON Schema coverage for the public Open Bio Design Package contract."""

import json
from copy import deepcopy
from pathlib import Path

import jsonschema
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = json.loads(
    (ROOT / "src/factorforge/schemas/design_package.schema.json").read_text(encoding="utf-8")
)


def test_schema_is_draft_2020_12() -> None:
    assert SCHEMA["$schema"] == "https://json-schema.org/draft/2020-12/schema"
    jsonschema.Draft202012Validator.check_schema(SCHEMA)


def test_valid_design_package(design_package: dict) -> None:
    jsonschema.validate(design_package, SCHEMA)


def test_mfe_status_reason_is_optional_and_enum_checked(design_package: dict) -> None:
    without_reason = deepcopy(design_package)
    without_reason["metrics"].pop("mfe_status_reason", None)
    jsonschema.validate(without_reason, SCHEMA)

    with_reason = deepcopy(design_package)
    with_reason["metrics"]["mfe_status_reason"] = "missing_dependency"
    jsonschema.validate(with_reason, SCHEMA)

    invalid = deepcopy(design_package)
    invalid["metrics"]["mfe_status_reason"] = "dependency_missing"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(invalid, SCHEMA)


@pytest.mark.parametrize(
    "field",
    [
        "design_package_version",
        "design_id",
        "sequence_type",
        "input_type",
        "host_profile",
        "optimization",
        "metrics",
        "validation",
        "evidence",
        "claim_boundary",
    ],
)
def test_required_top_level_fields(design_package: dict, field: str) -> None:
    invalid = deepcopy(design_package)
    invalid.pop(field)
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(invalid, SCHEMA)
