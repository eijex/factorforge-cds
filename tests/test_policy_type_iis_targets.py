import re
from pathlib import Path

import pytest

from factorforge.rules.registry import RuleRegistry


@pytest.mark.parametrize(
    "rule_id,forward,reverse",
    [
        ("assembly.type_iis.bpii.v1", "GAAGAC", "GTCTTC"),
        ("assembly.type_iis.sapi.v1", "GAAGAGC", "GCTCTTC"),
    ],
)
def test_policy_type_iis_targets_registered_and_scan_both_strands(rule_id, forward, reverse):
    rule = RuleRegistry().get_rule(rule_id)
    assert rule is not None
    evaluator = rule.evaluator_fn
    assert evaluator("ATG" + forward + "AAA", {})["forward_positions"] == [3]
    assert evaluator("ATG" + reverse.lower() + "AAA", {})["reverse_positions"] == [3]
    assert evaluator("ATGACCAACTAA", {})["passed"] is True


def test_web_presets_only_reference_registered_rules():
    source = (Path(__file__).parents[1] / "web/js/app.js").read_text(encoding="utf-8")
    rule_ids = set(re.findall(r"'([a-z_]+\.[a-z_.]+\.v\d+)':", source))
    known = {rule.rule_id for rule in RuleRegistry().list_rules()}
    assert rule_ids
    assert rule_ids <= known
