"""
FactorForge profile engine - Rule-based Engine

Production system (2026)
Plant-specific rule-based optimization
"""

from factorforge.registry.versioning import engine_version

__version__ = engine_version("profile")

from .optimizer import RuleBasedOptimizer
from .pipeline import OptimizationPipeline
from .construct_builder import ConstructBuilder

__all__ = ["ConstructBuilder", "OptimizationPipeline", "RuleBasedOptimizer"]
