"""
FactorForge - Codon Optimization Platform

v2: Rule-based engine (Production)
v3: ML engine (ESM2 + BART)
"""

__version__ = "3.1.1"
__author__ = "Eijex"

# Auto-register engines (safe when running from source tree)
try:
    from .engines import EngineRegistry, register_builtin_engines

    register_builtin_engines()
    __all__ = ["EngineRegistry"]
except Exception:
    __all__ = ["__version__"]
