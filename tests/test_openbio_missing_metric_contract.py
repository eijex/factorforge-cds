"""MFE missing-value contract for Open Bio public outputs."""

from copy import deepcopy

import jsonschema
import pytest

from test_design_package_schema import SCHEMA


@pytest.mark.parametrize("status", ["not_computed", "unavailable", "failed"])
def test_noncomputed_mfe_requires_null_and_unused(design_package: dict, status: str) -> None:
    package = deepcopy(design_package)
    package["metrics"].update(mfe_status=status, mfe_kcal_mol=None, mfe_used=False)
    jsonschema.validate(package, SCHEMA)


@pytest.mark.parametrize("status", ["not_computed", "unavailable", "failed"])
def test_noncomputed_mfe_rejects_numeric_value(design_package: dict, status: str) -> None:
    package = deepcopy(design_package)
    package["metrics"].update(mfe_status=status, mfe_kcal_mol=-1.0, mfe_used=False)
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(package, SCHEMA)


def test_computed_mfe_accepts_real_float(design_package: dict) -> None:
    package = deepcopy(design_package)
    package["metrics"].update(mfe_status="computed", mfe_kcal_mol=-12.3, mfe_used=True)
    jsonschema.validate(package, SCHEMA)


def test_computed_mfe_rejects_zero_placeholder(design_package: dict) -> None:
    package = deepcopy(design_package)
    package["metrics"].update(mfe_status="computed", mfe_kcal_mol=0.0, mfe_used=True)
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(package, SCHEMA)
