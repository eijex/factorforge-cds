import json
from pathlib import Path

import pytest

from factorforge.rules.models import EnforcementLevel
from factorforge.rules.registry import RuleRegistry
from factorforge.sops.config import SopProfile


SOP_PATH = Path(__file__).parents[1] / "configs" / "sops" / "default_conservative.json"


def test_default_conservative_sop_is_valid_and_hash_stable():
    known = [rule.rule_id for rule in RuleRegistry().list_rules()]
    first = SopProfile(SOP_PATH, known_rule_ids=known)
    second = SopProfile.from_dict(json.loads(SOP_PATH.read_text()), known_rule_ids=known)
    assert first.profile_id == "default_conservative_plant_expression"
    assert first.digest == second.digest
    assert first.get_enforcement("rna.cryptic_splice.v1") == EnforcementLevel.WARNING


def test_sop_unknown_rule_fails_closed():
    data = json.loads(SOP_PATH.read_text())
    data["rules"]["unknown.rule"] = "WARNING"
    with pytest.raises(ValueError, match="Unknown SOP rule IDs"):
        SopProfile.from_dict(data, known_rule_ids=[rule.rule_id for rule in RuleRegistry().list_rules()])


def test_sop_policy_overrides_registry_without_changing_detector():
    data = json.loads(SOP_PATH.read_text())
    data["rules"]["rna.polya_motifs.v1"] = "HARD_FAIL"
    profile = SopProfile.from_dict(data, known_rule_ids=[rule.rule_id for rule in RuleRegistry().list_rules()])
    result = RuleRegistry(sop_profile=profile).evaluate_sequence("ATGGCTAATAAAGCTGCTGCT")
    assert "rna.polya_motifs.v1" in [item["rule_id"] for item in result["hard_fails"]]
