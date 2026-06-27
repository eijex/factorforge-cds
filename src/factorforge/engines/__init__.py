"""Optimization Engines"""

from __future__ import annotations

from .registry import EngineRegistry


def register_builtin_engines() -> None:
    """Register bundled engines."""
    from .profile import RuleBasedOptimizer

    EngineRegistry.register(
        "profile",
        RuleBasedOptimizer,
        metadata={
            "version": "3.2.6",
            "engine_type": "profile_rule_based",
            "role": "stable_profile_engine",
            "stable": True,
        },
    )


__all__ = ["EngineRegistry", "register_builtin_engines"]
