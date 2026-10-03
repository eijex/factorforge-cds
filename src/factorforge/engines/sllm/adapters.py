# factorforge/src/factorforge/engines/sllm/adapters.py
"""Model Adapter Interfaces and Implementations for sLLM Generative Runtime (Job 284A).

Decouples core inference from external model formats and network transport.
"""

from __future__ import annotations

import hashlib
import math
from abc import ABC, abstractmethod
from pathlib import Path

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
        prefix_codons: list[str],
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
        prefix_codons: list[str],
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


class DeadEndError(ValueError):
    """No legal synonymous codon is available."""


class OnnxSLMAdapter(SLMModelAdapter):
    """Optional ONNX adapter; invalid artifacts never silently become mock models."""

    def __init__(
        self,
        model_path=None,
        expected_sha256=None,
        vocab=None,
        use_int8=True,
        *,
        onnx_model_path=None,
        full_protein_aa=None,
    ):
        if expected_sha256 is not None and vocab is None and len(expected_sha256) != 64:
            if full_protein_aa is not None:
                raise ValueError("Conflicting protein context arguments")
            full_protein_aa, expected_sha256 = expected_sha256, None
        if model_path is not None and onnx_model_path is not None:
            raise ValueError("Supply one model path")
        path = model_path or onnx_model_path
        if path is None:
            root = Path(__file__).resolve().parents[3] / "data" / "models"
            path = next(
                (
                    str(p)
                    for p in (root / "sllm_nbent_v1.quant.onnx", root / "sllm_nbent_v1.onnx")
                    if p.is_file()
                ),
                None,
            )
        if path is None or not Path(path).is_file():
            raise FileNotFoundError("Optional ONNX model is not installed")
        self.model_path = str(path)
        self.full_protein_aa = full_protein_aa
        self.vocab = dict(CODON_TO_INDEX if vocab is None else vocab)
        if (
            set(self.vocab) != set(CODON_TO_INDEX)
            or any(type(i) is not int for i in self.vocab.values())
            or set(self.vocab.values()) != set(range(64))
        ):
            raise ValueError("Vocabulary must bijectively map 64 codons to IDs 0..63")
        self.inv_vocab = {i: c for c, i in self.vocab.items()}
        # Load exactly the checked bytes, not a path that can change after verification.
        artifact = Path(path).read_bytes()
        if expected_sha256 is not None:
            if len(expected_sha256) != 64 or any(
                c not in "0123456789abcdefABCDEF" for c in expected_sha256
            ):
                raise ValueError("Expected SHA-256 must be 64 hexadecimal characters")
            if hashlib.sha256(artifact).hexdigest() != expected_sha256.lower():
                raise ValueError("Model SHA-256 mismatch")
        try:
            import onnxruntime as ort
        except ImportError as exc:
            raise ImportError("Install factorforge-cds[ml] to use the ONNX adapter") from exc
        opts = ort.SessionOptions()
        opts.intra_op_num_threads = 1
        opts.inter_op_num_threads = 1
        # use_int8 is retained for compatibility; quantization belongs to the artifact.
        self.session = ort.InferenceSession(artifact, opts, providers=["CPUExecutionProvider"])
        names = {i.name for i in self.session.get_inputs()}
        if names == {"input_ids"}:
            self._input_contract = "input_ids"
        elif names == {"codon_ids", "aa_ids"}:
            if self.vocab != CODON_TO_INDEX:
                raise ValueError("Paired-input models require the canonical vocabulary")
            self._input_contract = "paired"
        else:
            raise ValueError("Unsupported ONNX input contract")

    @property
    def model_name(self):
        return f"onnx_{Path(self.model_path).name}"

    def set_target_protein(self, protein_aa):
        self.full_protein_aa = protein_aa

    def predict_logits(self, prefix_codons, target_aa, position):
        from factorforge.engines.sllm.interfaces import STANDARD_GENETIC_CODE

        if target_aa not in STANDARD_GENETIC_CODE or position != len(prefix_codons):
            raise ValueError("Target amino acid and position must match context")
        if any(c not in self.vocab for c in prefix_codons):
            raise ValueError("Unknown context codon")
        if self._input_contract == "input_ids":
            ids = [self.vocab[c] for c in prefix_codons] or [0]
            inputs = {"input_ids": np.asarray([ids[-128:]], dtype=np.int64)}
        else:
            from factorforge.engines.sllm.data.tokenizer import AA_TO_INDEX, BOS_CODON_ID

            reverse = {c: a for a, codons in STANDARD_GENETIC_CODE.items() for c in codons}
            observed = "".join(reverse[c] for c in prefix_codons) + target_aa
            context = (
                observed if self.full_protein_aa is None else self.full_protein_aa[: position + 1]
            )
            if context != observed or any(a not in AA_TO_INDEX for a in context):
                raise ValueError("Protein context must match prefix and target")
            ids = [BOS_CODON_ID] + [self.vocab[c] for c in prefix_codons]
            inputs = {
                "codon_ids": np.asarray([ids[-128:]], dtype=np.int64),
                "aa_ids": np.asarray([[AA_TO_INDEX[a] for a in context[-128:]]], dtype=np.int64),
            }
        outputs = self.session.run(None, inputs)
        if not outputs:
            raise ValueError("ONNX model returned no logits")
        raw = np.asarray(outputs[0])
        if (
            raw.ndim != 3
            or raw.shape[:2] != next(iter(inputs.values())).shape
            or raw.shape[2] != 64
        ):
            raise ValueError("Expected logits shape (1, context_length, 64)")
        scores = raw[0, -1, :]
        if not np.all(np.isfinite(scores)):
            raise ValueError("ONNX logits must be finite")
        canonical = np.asarray([scores[self.vocab[c]] for c in CODON_TO_INDEX], dtype=np.float32)
        if not np.all(np.isfinite(canonical)):
            raise ValueError("Logits exceed the supported float32 range")
        return CodonLogits(canonical)

    def predict_next_codon(
        self,
        context_codons,
        target_aa,
        allowed_codons,
        decoding_mode="sampling",
        temperature=1.0,
        seed=42,
    ):
        from factorforge.engines.sllm.interfaces import STANDARD_GENETIC_CODE

        if not allowed_codons:
            raise DeadEndError("No allowed synonymous codons")
        if target_aa not in STANDARD_GENETIC_CODE or any(
            c not in STANDARD_GENETIC_CODE[target_aa] for c in allowed_codons
        ):
            raise ValueError("Allowed codons must preserve the target amino acid")
        if decoding_mode not in {"sampling", "greedy"}:
            raise ValueError("Unknown decoding mode")
        if not np.isfinite(temperature) or temperature <= 0:
            raise ValueError("Temperature must be finite and positive")
        scores = self.predict_logits(context_codons, target_aa, len(context_codons)).scores
        ids = sorted({CODON_TO_INDEX[c] for c in allowed_codons})
        if decoding_mode == "greedy":
            chosen = ids[int(np.argmax(scores[ids]))]
        else:
            weights = scores[ids].astype(np.float64)
            weights = (weights - weights.max()) / temperature
            probs = np.exp(weights)
            probs /= probs.sum()
            chosen = int(np.random.default_rng(seed).choice(ids, p=probs))
        return next(c for c, i in CODON_TO_INDEX.items() if i == chosen)
