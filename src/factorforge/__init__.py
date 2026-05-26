"""
FactorForge - Codon Optimization Platform

profile: constraint-aware rule/profile engine
v2: compatibility alias for profile
v3: experimental engine namespace
"""

__version__ = "3.1.3"
__author__ = "Eijex"

# Auto-register engines (safe when running from source tree)
try:
    from .engines import EngineRegistry, register_builtin_engines

    register_builtin_engines()
    __all__ = ["EngineRegistry"]
except Exception:
    __all__ = ["__version__"]
