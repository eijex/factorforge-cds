# factorforge/tests/test_validation_contract_compat.py
"""Cross-repo compatibility test: the runtime registry must still cover every
check the frozen validation-check taxonomy snapshot (validation_contract_v1.yaml)
claimed to exist. The registry is the ongoing source of truth; the snapshot is a
pinned historical reference and may be a subset of the registry, never a
superset."""

from __future__ import annotations

from pathlib import Path

import yaml

from factorforge.validation_registry import REVIEW_CONTRACT_CHECK_IDS, VALIDATION_CHECKS

_FIXTURE_PATH = (
    Path(__file__).parent / "fixtures" / "contracts" / "validation_contract_v1.yaml"
)


def _load_contract() -> dict:
    with open(_FIXTURE_PATH, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def test_fixture_file_exists() -> None:
    assert _FIXTURE_PATH.exists()


def test_frozen_review_contract_ids_are_a_subset_of_runtime_review_ids() -> None:
    contract = _load_contract()
    contract_review_ids = {
        check["check_id"]
        for check in contract["checks"]
        if check["primary_domain"] in {"advisory_scan", "assembly_review"}
    }
    assert contract_review_ids <= set(REVIEW_CONTRACT_CHECK_IDS)


def test_frozen_contract_check_count_matches_registry_expectation() -> None:
    contract = _load_contract()
    assert len(contract["checks"]) == 17


def test_review_contract_check_ids_length_is_11() -> None:
    assert len(REVIEW_CONTRACT_CHECK_IDS) == 11


def test_review_contract_ids_are_all_present_in_the_full_registry() -> None:
    registry_ids = {check["check_id"] for check in VALIDATION_CHECKS}
    assert set(REVIEW_CONTRACT_CHECK_IDS) <= registry_ids
