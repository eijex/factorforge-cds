"""Optimization Engines"""

from __future__ import annotations

from typing import Type

from factorforge.core.interfaces import OptimizerEngine

from .registry import EngineRegistry


def _load_v1() -> Type[OptimizerEngine]:
    raise ImportError(
        "FactorForge v1 is archived. Install with: pip install factorforge[v1]"
    )


def _load_v3() -> Type[OptimizerEngine]:
    from .v3 import V3Optimizer

    return V3Optimizer  # type: ignore[return-value]


def register_builtin_engines() -> None:
    """Register bundled engines with lazy loaders for archived/ML engines."""
    from .v2 import RuleBasedOptimizer

    EngineRegistry.register(
        "v2",
        RuleBasedOptimizer,
        metadata={
            "version": "2.5.3",
            "engine_type": "rule_based",
            "role": "legacy_fallback",
            "stable": True,
        },
    )
    EngineRegistry.register_lazy(
        "v1",
        _load_v1,
        metadata={
            "version": "archived",
            "engine_type": "rule_based",
            "role": "archived",
            "stable": False,
        },
    )
    EngineRegistry.register_lazy(
        "v3",
        _load_v3,
        metadata={
            "version": "alpha",
            "engine_type": "ml",
            "role": "experimental",
            "stable": False,
        },
    )


__all__ = ["EngineRegistry", "register_builtin_engines"]
