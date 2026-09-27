# factorforge/tests/test_rule_registry.py
"""Comprehensive Tests for Rule Registry, Scopes, Authority Attribution, and Enforcement Levels."""

import pytest
from factorforge.rules.models import (
    AuthorityType,
    EnforcementLevel,
    RuleAuthority,
    RuleCategory,
    RuleDefinition,
    RuleScope,
)
from factorforge.rules.registry import RuleRegistry


def test_rule_registry_initialization() -> None:
    registry = RuleRegistry()
    rules = registry.list_rules()
    assert len(rules) >= 8
    
    bsai_rule = registry.get_rule("assembly.type_iis.bsai.v1")
    assert bsai_rule is not None
    assert bsai_rule.enforcement == EnforcementLevel.HARD_FAIL
    assert bsai_rule.authority.authority_type == AuthorityType.EIJEX_INTERNAL_POLICY


def test_rule_evaluation_hard_fail_on_bsai_and_bsmbi() -> None:
    registry = RuleRegistry()
    
    # Sequence with BsaI site (GGTCTC) and BsmBI site (CGTCTC)
    dirty_seq = "ATGGCAAAAGGTCTCAAACTTCGTCTCTGA"
    res = registry.evaluate_sequence(dirty_seq, assembly_method="golden_gate")
    
    assert res["all_passed"] is False
    assert res["hard_fail_count"] == 2
    failed_rule_ids = [f["rule_id"] for f in res["hard_fails"]]
    assert "assembly.type_iis.bsai.v1" in failed_rule_ids
    assert "assembly.type_iis.bsmbi.v1" in failed_rule_ids


def test_rule_evaluation_hard_fail_on_reading_frame() -> None:
    registry = RuleRegistry()
    
    # Sequence not multiple of 3
    bad_frame_seq = "ATGGCATCAAC" # length 11
    res = registry.evaluate_sequence(bad_frame_seq)
    
    assert res["all_passed"] is False
    failed_rule_ids = [f["rule_id"] for f in res["hard_fails"]]
    assert "biological.reading_frame.v1" in failed_rule_ids


def test_rule_evaluation_warnings_polya_and_are() -> None:
    registry = RuleRegistry()

    # Clean assembly (len 33, starts with ATG) but contains premature polyA (AATAAA) and ARE element (ATTTA)
    seq = "ATGGCAAAAGCTAATAAAGCTATTTACACTGAA"  # length 33
    res = registry.evaluate_sequence(seq)

    # Hard fails must be 0, but warnings must be present
    assert res["hard_fail_count"] == 0
    assert res["warning_count"] >= 2
    warning_ids = [w["rule_id"] for w in res["warnings"]]
    assert "rna.polya_motifs.v1" in warning_ids
    assert "rna.au_rich_elements.v1" in warning_ids


def test_rule_evaluation_clean_sequence() -> None:
    registry = RuleRegistry()
    
    # Clean sequence with valid ATG, multiple of 3, no BsaI/BsmBI, no polyA, no ARE, no 6-homopolymers
    clean_seq = "ATGGCATCAACTCAATCTTCTACTGCATGA"
    res = registry.evaluate_sequence(clean_seq, assembly_method="golden_braid")
    
    assert res["all_passed"] is True
    assert res["hard_fail_count"] == 0


def test_rule_scope_filtering() -> None:
    registry = RuleRegistry()
    
    # Add custom rule scoped only to yeast
    yeast_rule = RuleDefinition(
        rule_id="host.yeast.specific.v1",
        name="Yeast Specific Rule",
        description="Only applies to S. cerevisiae",
        category=RuleCategory.INTERNAL_POLICY,
        enforcement=EnforcementLevel.WARNING,
        authority=RuleAuthority(
            authority_type=AuthorityType.EIJEX_INTERNAL_POLICY,
            source_name="yeast_guidelines",
        ),
        scope=RuleScope(target_hosts={"Saccharomyces cerevisiae"}),
    )
    registry.register(yeast_rule)
    
    plant_rules = registry.list_rules(host="Nicotiana benthamiana")
    assert "host.yeast.specific.v1" not in [r.rule_id for r in plant_rules]
    
    yeast_rules = registry.list_rules(host="Saccharomyces cerevisiae")
    assert "host.yeast.specific.v1" in [r.rule_id for r in yeast_rules]


def test_ruleset_digest_stability() -> None:
    registry1 = RuleRegistry()
    registry2 = RuleRegistry()
    
    digest1 = registry1.compute_digest()
    digest2 = registry2.compute_digest()
    
    assert digest1.startswith("sha256:")
    assert digest1 == digest2
