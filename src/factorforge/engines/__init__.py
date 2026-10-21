"""Optimization Engines"""

from __future__ import annotations

from .registry import EngineRegistry
from factorforge.registry.versioning import engine_metadata


def _load_dp_v2_2():
    """Load the opt-in core only when explicitly requested."""
    from .dp_v2_2 import DPV220EngineAdapter

    return DPV220EngineAdapter


def _load_dp_v2_3():
    """Load DP v2.3.0 engine when requested."""
    from .dp_v2_3 import DPV230EngineAdapter

    return DPV230EngineAdapter


def register_builtin_engines() -> None:
    """Register bundled engines."""
    from .profile.optimizer import RuleBasedOptimizer
    from .dp_adapter import DPEngineAdapter

    EngineRegistry.register_lazy(
        "dp_v2_2",
        _load_dp_v2_2,
        metadata={**engine_metadata("dp_v2_2"), "stable": False},
    )

    EngineRegistry.register_lazy(
        "dp_v2_3",
        _load_dp_v2_3,
        metadata={**engine_metadata("dp_v2_3"), "stable": False},
    )

    EngineRegistry.register(
        "profile",
        RuleBasedOptimizer,
        metadata={
            **engine_metadata("profile"),
            "engine_type": "profile_rule_based",
            "role": "stable_profile_engine",
            "stable": True,
        },
    )
    EngineRegistry.register(
        "dp",
        DPEngineAdapter,
        metadata={
            **engine_metadata("dp"),
            "engine_type": "deterministic_constrained_optimizer",
            "role": "stable_dp_engine",
            "stable": True,
        },
    )

    try:
        from .lm.inference import LMEngineAdapter

        EngineRegistry.register(
            "lm",
            LMEngineAdapter,
            metadata={
                **engine_metadata("slm"),
                "engine_type": "constrained_beam_search_lm",
                "role": "experimental_lm_engine",
                "stable": False,
            },
        )
        EngineRegistry.register(
            "slm",
            LMEngineAdapter,
            metadata={
                **engine_metadata("slm"),
                "engine_type": "constrained_beam_search_lm",
                "role": "experimental_slm_engine",
                "stable": False,
            },
        )
    except Exception:
        # LM engines require PyTorch/heavy dependencies which may not be present in lightweight/serverless bundles
        pass


__all__ = ["EngineRegistry", "register_builtin_engines"]
