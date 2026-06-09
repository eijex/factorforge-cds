"""Cross-field semantic rules for the public Design Package."""

from copy import deepcopy

import jsonschema
import pytest

from test_design_package_schema import SCHEMA


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("validation", "aa_identity"), 0.99),
        (("validation", "forbidden_type_iis_site_count"), 1),
        (("metrics", "mfe_kcal_mol"), -12.3),
        (("metrics", "mfe_used"), True),
        (("claim_boundary", "no_yield_claim"), False),
        (("host_profile", "ncbi_taxonomy_id"), 4197),
    ],
)
def test_invalid_cross_field_semantics(design_package: dict, path: tuple[str, str], value) -> None:
    invalid = deepcopy(design_package)
    invalid[path[0]][path[1]] = value
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(invalid, SCHEMA)


def test_by2_registry_metadata_is_supported(design_package: dict) -> None:
    package = deepcopy(design_package)
    package["host_profile"] = {
        "id": "by2",
        "display_name": "BY-2 / N. tabacum",
        "scientific_name": "Nicotiana tabacum",
        "ncbi_taxonomy_id": 4097,
        "status": "experimental",
    }
    jsonschema.validate(package, SCHEMA)
