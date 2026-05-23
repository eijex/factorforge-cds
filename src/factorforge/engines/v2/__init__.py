"""
FactorForge v2 - Rule-based Engine

Production system (2026)
Plant-specific rule-based optimization
"""

__version__ = "3.0.0"

from .optimizer import RuleBasedOptimizer
from .pipeline import OptimizationPipeline

__all__ = ["OptimizationPipeline", "RuleBasedOptimizer"]
