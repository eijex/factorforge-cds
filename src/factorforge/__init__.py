"""
FactorForge - Codon Optimization Platform

v2: Rule-based (Production)
v1: BPE Tokenizer (Archived)
"""

__version__ = "2.5.2"
__author__ = "Mun-Kyu Kim"

# Auto-register engines (safe when running from source tree)
try:
    from .engines import EngineRegistry, register_builtin_engines

    register_builtin_engines()
    __all__ = ["EngineRegistry"]
except Exception:
    __all__ = ["__version__"]
