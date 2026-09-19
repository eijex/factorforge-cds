# factorforge/src/factorforge/engines/sllm/adapters.py
"""Model Adapter Interfaces and Implementations for sLLM Generative Runtime (Job 284A).

Decouples core inference from external model formats and network transport.
"""

from __future__ import annotations
from abc import ABC, abstractmethod
import math
import os
from typing import List, Optional
import numpy as np

from factorforge.analysis.metrics import load_codon_usage_table
from factorforge.engines.sllm.interfaces import (
    CODON_TO_INDEX,
    CodonLogits,
)


class SLMModelAdapter(ABC):
    """Abstract interface for genomic language model logit generation."""

    @property
    @abstractmethod
    def model_name(self) -> str:
        """Name or identifier of the underlying model."""
        raise NotImplementedError

    @abstractmethod
    def predict_logits(
        self,
        prefix_codons: List[str],
        target_aa: str,
        position: int,
    ) -> CodonLogits:
        """Predicts raw unconstrained 64-dimensional logits for position t given context."""
        raise NotImplementedError


class DeterministicMockSLMAdapter(SLMModelAdapter):
    """Deterministic, biophysically realistic model adapter for CI/CD and testing.

    Uses host codon usage preferences with positional pseudo-context modulation
    to produce stable, reproducible logit distributions.
    """

    def __init__(self, host: str = "nbenthamiana", seed: int = 42) -> None:
        self.host = host
        self.seed = seed
        table = load_codon_usage_table()
        self.codon_weights = table.codon_weights

    @property
    def model_name(self) -> str:
        return f"mock_sllm_v1_{self.host}"

    def predict_logits(
        self,
        prefix_codons: List[str],
        target_aa: str,
        position: int,
    ) -> CodonLogits:
        # Base logit values from codon weights
        scores = np.zeros(64, dtype=np.float32)
        for codon, idx in CODON_TO_INDEX.items():
            w = self.codon_weights.get(codon, 0.05)
            # Log-scale weighting
            base_score = math.log(max(1e-4, w))

            # Subtle context modulation based on position & last codon
            if prefix_codons:
                last_codon = prefix_codons[-1]
                # Small deterministic pseudo-co-occurrence bonus
                bonus = (hash(last_codon + codon + str(position)) % 100) / 1000.0
            else:
                bonus = 0.0

            scores[idx] = base_score + bonus

        return CodonLogits(scores)


class OnnxSLMAdapter(SLMModelAdapter):
    """Lightweight ONNX Runtime inference adapter for local CPU/GPU execution (Job 285B/C)."""

    def __init__(
        self,
        onnx_model_path: Optional[str] = None,
        full_protein_aa: Optional[str] = None,
    ) -> None:
        if onnx_model_path is None:
            # Default auto-discovery in data/models
            possible_paths = [
                os.path.join(
                    os.path.dirname(__file__),
                    "..",
                    "..",
                    "..",
                    "data",
                    "models",
                    "sllm_nbent_v1.quant.onnx",
                ),
                os.path.join(
                    os.path.dirname(__file__),
                    "..",
                    "..",
                    "..",
                    "data",
                    "models",
                    "sllm_nbent_v1.onnx",
                ),
                r"c:\Work\eijex\factorforge\data\models\sllm_nbent_v1.quant.onnx",
                r"c:\Work\eijex\factorforge\data\models\sllm_nbent_v1.onnx",
            ]
            for p in possible_paths:
                norm_p = os.path.abspath(p)
                if os.path.exists(norm_p):
                    onnx_model_path = norm_p
                    break

        if onnx_model_path is None or not os.path.exists(onnx_model_path):
            raise FileNotFoundError(f"ONNX model file not found at {onnx_model_path}")

        self.model_path = onnx_model_path
        self.full_protein_aa = full_protein_aa

        try:
            import onnxruntime as ort

            # Configure efficient CPU execution
            opts = ort.SessionOptions()
            opts.intra_op_num_threads = 1
            opts.inter_op_num_threads = 1
            opts.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
            self.session = ort.InferenceSession(
                onnx_model_path, opts, providers=["CPUExecutionProvider"]
            )
        except ImportError as e:
            raise ImportError("onnxruntime package is required to use OnnxSLMAdapter.") from e

        from factorforge.engines.sllm.data.tokenizer import CodonSequenceTokenizer

        self.tokenizer = CodonSequenceTokenizer()

    @property
    def model_name(self) -> str:
        return f"onnx_{os.path.basename(self.model_path)}"

    def set_target_protein(self, protein_aa: str) -> None:
        """Sets full target protein sequence for contextual positional decoding."""
        self.full_protein_aa = protein_aa

    def predict_logits(
        self,
        prefix_codons: List[str],
        target_aa: str,
        position: int,
    ) -> CodonLogits:
        # Build autoregressive codon input: [BOS_CODON_ID] + prefix_codons
        from factorforge.engines.sllm.data.tokenizer import BOS_CODON_ID, AA_TO_INDEX

        prefix_ids = [CODON_TO_INDEX.get(c, 0) for c in prefix_codons]
        input_codons = [BOS_CODON_ID] + prefix_ids

        # Build amino acid sequence context up to current position
        if self.full_protein_aa is not None and position < len(self.full_protein_aa):
            aa_context = self.full_protein_aa[: position + 1]
        else:
            # Fallback when full protein context is not attached: synthesize from prefix + target
            aa_context = "".join(["A"] * len(prefix_codons)) + target_aa

        aa_ids = [AA_TO_INDEX.get(a, 0) for a in aa_context]
        seq_len = min(len(input_codons), len(aa_ids))
        input_codons = input_codons[:seq_len]
        aa_ids = aa_ids[:seq_len]

        # Context window truncation (max 128 codons) for fast, unbounded sequence length decoding
        MAX_CONTEXT = 128
        if len(input_codons) > MAX_CONTEXT:
            input_codons = input_codons[-MAX_CONTEXT:]
            aa_ids = aa_ids[-MAX_CONTEXT:]

        codon_tensor = np.array([input_codons], dtype=np.int64)
        aa_tensor = np.array([aa_ids], dtype=np.int64)

        inputs = {
            "codon_ids": codon_tensor,
            "aa_ids": aa_tensor,
        }
        outputs = self.session.run(None, inputs)
        raw_last_step = outputs[0][0, -1, :64].astype(np.float32)
        return CodonLogits(raw_last_step)
