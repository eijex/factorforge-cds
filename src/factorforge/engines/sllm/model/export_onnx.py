# factorforge/src/factorforge/engines/sllm/model/export_onnx.py
"""Two-Stage ONNX Export & INT8 Quantization Pipeline (Job 285B).

Stage 1: Export PyTorch CompactCodonTransformer to dynamic-axes ONNX (FP32).
Stage 2: Apply dynamic INT8 quantization and verify behavioral equivalence.
"""

from __future__ import annotations
import os
from typing import Tuple
import numpy as np
import onnx
import onnxruntime as ort
import torch

from factorforge.engines.sllm.model.arch import CompactCodonTransformer


def export_to_onnx(
    model: CompactCodonTransformer,
    output_onnx_path: str,
    max_len: int = 512,
) -> str:
    """Stage 1: Exports PyTorch model to dynamic-axes FP32 ONNX."""
    model.eval()
    os.makedirs(os.path.dirname(os.path.abspath(output_onnx_path)), exist_ok=True)

    # Dummy inputs
    dummy_codon = torch.zeros((1, 8), dtype=torch.long)
    dummy_aa = torch.zeros((1, 8), dtype=torch.long)

    try:
        torch.onnx.export(
            model,
            (dummy_codon, dummy_aa),
            output_onnx_path,
            input_names=["codon_ids", "aa_ids"],
            output_names=["logits"],
            dynamic_axes={
                "codon_ids": {0: "batch_size", 1: "seq_len"},
                "aa_ids": {0: "batch_size", 1: "seq_len"},
                "logits": {0: "batch_size", 1: "seq_len"},
            },
            opset_version=14,
            do_constant_folding=True,
            dynamo=False,
        )
    except TypeError:
        # Fallback for older torch versions without dynamo flag
        torch.onnx.export(
            model,
            (dummy_codon, dummy_aa),
            output_onnx_path,
            input_names=["codon_ids", "aa_ids"],
            output_names=["logits"],
            dynamic_axes={
                "codon_ids": {0: "batch_size", 1: "seq_len"},
                "aa_ids": {0: "batch_size", 1: "seq_len"},
                "logits": {0: "batch_size", 1: "seq_len"},
            },
            opset_version=14,
            do_constant_folding=True,
        )

    # Verify ONNX model validity
    onnx_model = onnx.load(output_onnx_path)
    onnx.checker.check_model(onnx_model)
    return output_onnx_path


def quantize_onnx_model(
    input_onnx_path: str,
    output_quant_path: str,
) -> str:
    """Stage 2: Quantizes FP32 ONNX model to dynamic INT8."""
    from onnxruntime.quantization import QuantType, quantize_dynamic

    os.makedirs(os.path.dirname(os.path.abspath(output_quant_path)), exist_ok=True)
    quantize_dynamic(
        model_input=input_onnx_path,
        model_output=output_quant_path,
        weight_type=QuantType.QInt8,
    )
    return output_quant_path


def verify_onnx_parity(
    torch_model: CompactCodonTransformer,
    onnx_path: str,
    test_length: int = 16,
) -> Tuple[float, bool]:
    """Verifies numerical parity between PyTorch and ONNX runtime (MSE < 1e-4)."""
    torch_model.eval()

    dummy_codon = torch.randint(0, 64, (1, test_length), dtype=torch.long)
    dummy_aa = torch.randint(0, 20, (1, test_length), dtype=torch.long)

    with torch.no_grad():
        torch_logits = torch_model(dummy_codon, dummy_aa).cpu().numpy()

    session = ort.InferenceSession(onnx_path, providers=["CPUExecutionProvider"])
    ort_inputs = {
        "codon_ids": dummy_codon.numpy().astype(np.int64),
        "aa_ids": dummy_aa.numpy().astype(np.int64),
    }
    ort_logits = session.run(None, ort_inputs)[0]

    mse = float(np.mean((torch_logits - ort_logits) ** 2))
    is_valid = mse < 1e-4
    return mse, is_valid
