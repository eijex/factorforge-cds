"""
FactorForge - Codon Optimization Platform

v1_archived: Rule-based v1 (Archived)
v2: Rule-based (Production) — engine version 2.5.3
v3: ML engine / v3-alpha (ESM2 + BART, in development)
"""

__version__ = "1.0.0"
__author__ = "Eijex"

# Auto-register engines (safe when running from source tree)
try:
    from .engines import EngineRegistry, register_builtin_engines

    register_builtin_engines()
    __all__ = ["EngineRegistry"]
except Exception:
    __all__ = ["__version__"]
