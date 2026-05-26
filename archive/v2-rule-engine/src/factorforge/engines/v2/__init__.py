"""Compatibility alias for the FactorForge profile engine.

New code should import from `factorforge.engines.profile`. This package remains
available so existing `factorforge.engines.v2` imports continue to work.
"""

from factorforge.engines.profile import OptimizationPipeline, RuleBasedOptimizer
from factorforge.engines.profile.construct_builder import ConstructBuilder

__version__ = "3.1.3"

__all__ = ["ConstructBuilder", "OptimizationPipeline", "RuleBasedOptimizer"]
