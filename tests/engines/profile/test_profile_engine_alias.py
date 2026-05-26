"""Tests for profile engine registration."""

from __future__ import annotations

from factorforge.engines import EngineRegistry, register_builtin_engines
from factorforge.engines.profile import (
    ConstructBuilder,
    OptimizationPipeline,
    RuleBasedOptimizer,
)


def test_profile_exports_stable_engine_classes() -> None:
    assert RuleBasedOptimizer.__name__ == "RuleBasedOptimizer"
    assert OptimizationPipeline.__name__ == "OptimizationPipeline"
    assert ConstructBuilder.__name__ == "ConstructBuilder"


def test_registry_registers_profile_engine() -> None:
    EngineRegistry.clear()
    register_builtin_engines()

    profile_optimizer = EngineRegistry.get("profile")
    metadata = EngineRegistry.list_with_metadata()

    assert isinstance(profile_optimizer, RuleBasedOptimizer)
    assert metadata["profile"]["role"] == "stable_profile_engine"
