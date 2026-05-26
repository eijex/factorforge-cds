"""Tests for profile engine registration and v2 compatibility aliases."""

from __future__ import annotations

from factorforge.engines import EngineRegistry, register_builtin_engines
from factorforge.engines.profile import (
    ConstructBuilder,
    OptimizationPipeline,
    RuleBasedOptimizer,
)
from factorforge.engines.v2 import (
    ConstructBuilder as V2ConstructBuilder,
    OptimizationPipeline as V2OptimizationPipeline,
    RuleBasedOptimizer as V2RuleBasedOptimizer,
)


def test_profile_exports_stable_engine_classes() -> None:
    assert RuleBasedOptimizer is V2RuleBasedOptimizer
    assert OptimizationPipeline is V2OptimizationPipeline
    assert ConstructBuilder is V2ConstructBuilder


def test_registry_registers_profile_and_v2_alias() -> None:
    EngineRegistry.clear()
    register_builtin_engines()

    profile_optimizer = EngineRegistry.get("profile")
    v2_optimizer = EngineRegistry.get("v2")
    metadata = EngineRegistry.list_with_metadata()

    assert isinstance(profile_optimizer, RuleBasedOptimizer)
    assert isinstance(v2_optimizer, RuleBasedOptimizer)
    assert metadata["profile"]["role"] == "stable_profile_engine"
    assert metadata["v2"]["role"] == "compatibility_alias"
    assert metadata["v2"]["alias_for"] == "profile"
    assert metadata["v2"]["deprecated"] is True
