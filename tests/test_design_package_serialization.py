"""Serialization and public-sequence privacy rules for Design Packages."""

import json
from copy import deepcopy

import jsonschema
import pytest

from test_design_package_schema import SCHEMA


def test_none_round_trips_as_json_null(design_package: dict) -> None:
    encoded = json.dumps(design_package)
    assert '"mfe_kcal_mol": null' in encoded
    assert json.loads(encoded)["metrics"]["mfe_kcal_mol"] is None


@pytest.mark.parametrize("placeholder", [0.0, "0.0", "N/A", "null", ""])
def test_missing_mfe_placeholders_are_rejected(design_package: dict, placeholder) -> None:
    invalid = deepcopy(design_package)
    invalid["metrics"]["mfe_kcal_mol"] = placeholder
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(invalid, SCHEMA)


@pytest.mark.parametrize("field", ["sequence", "raw_sequence", "input_sequence", "output_cds"])
def test_public_package_rejects_raw_sequence_fields(design_package: dict, field: str) -> None:
    invalid = deepcopy(design_package)
    invalid[field] = "ATGAAACCC"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(invalid, SCHEMA)
