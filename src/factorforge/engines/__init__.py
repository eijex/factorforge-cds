"""Optimization Engines"""

from __future__ import annotations

from .registry import EngineRegistry


def register_builtin_engines() -> None:
    """Register bundled engines."""
    from .profile import RuleBasedOptimizer
    from .lm.inference import LMEngineAdapter

    EngineRegistry.register(
        "profile",
        RuleBasedOptimizer,
        metadata={
            "version": "3.4.6",
            "engine_type": "profile_rule_based",
            "role": "stable_profile_engine",
            "stable": True,
        },
    )
    EngineRegistry.register(
        "lm",
        LMEngineAdapter,
        metadata={
            "version": "3.5.0",
            "engine_type": "constrained_beam_search_lm",
            "role": "experimental_lm_engine",
            "stable": False,
            "status": "work_in_progress",
        },
    )
    EngineRegistry.register(
        "slm",
        LMEngineAdapter,
        metadata={
            "version": "3.5.0",
            "engine_type": "constrained_beam_search_lm",
            "role": "experimental_slm_engine",
            "stable": False,
            "status": "work_in_progress",
        },
    )


__all__ = ["EngineRegistry", "register_builtin_engines"]
