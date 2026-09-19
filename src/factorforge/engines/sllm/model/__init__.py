# factorforge/src/factorforge/engines/sllm/model/__init__.py
"""Model architectures, baselines, and ONNX export for FactorForge sLLM (Job 285B)."""

from factorforge.engines.sllm.model.baselines import (
    CodonBigramMarkovBaseline,
    CodonFrequencyBaseline,
)

__all__ = [
    "CodonFrequencyBaseline",
    "CodonBigramMarkovBaseline",
]

try:
    from factorforge.engines.sllm.model.arch import (  # noqa: F401
        CompactCodonTransformer,
        compute_synonymous_masked_loss,
    )
    from factorforge.engines.sllm.model.export_onnx import (  # noqa: F401
        export_to_onnx,
        quantize_onnx_model,
        verify_onnx_parity,
    )

    __all__.extend(
        [
            "CompactCodonTransformer",
            "compute_synonymous_masked_loss",
            "export_to_onnx",
            "quantize_onnx_model",
            "verify_onnx_parity",
        ]
    )
except ImportError:
    pass
