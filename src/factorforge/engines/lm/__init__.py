"""FactorForge experimental LM scaffold.

Inactive research helpers for tokenizer/model-shape and deterministic decoding-mask
experiments. This package is not a supported public optimization engine.
"""

from factorforge.engines.lm.tokenizer.control_tokenizer import FactorForgeControlTokenizer

__all__ = ["FactorForgeControlTokenizer"]
