"""Synthetic artifact contracts only; not a trained-model evaluation."""

import hashlib
import subprocess
import sys

import numpy as np
import pytest

from factorforge.engines.sllm.adapters import DeadEndError, OnnxSLMAdapter, SLMModelAdapter
from factorforge.engines.sllm.interfaces import CODON_TO_INDEX


def fixture_model(tmp_path, *, paired=False, width=64, nonfinite=False):
    onnx = pytest.importorskip("onnx")
    pytest.importorskip("onnxruntime")
    helper, tensor = onnx.helper, onnx.TensorProto
    names = ["codon_ids", "aa_ids"] if paired else ["input_ids"]
    inputs = [helper.make_tensor_value_info(n, tensor.INT64, [1, "length"]) for n in names]
    scores = np.arange(width, dtype=np.float32)
    if nonfinite:
        scores[0] = np.nan
    nodes = [
        helper.make_node("Shape", [names[0]], ["context_shape"]),
        helper.make_node("Concat", ["context_shape", "width"], ["output_shape"], axis=0),
        helper.make_node("Expand", ["scores", "output_shape"], ["logits"]),
    ]
    graph = helper.make_graph(
        nodes,
        "synthetic-contract-only",
        inputs,
        [helper.make_tensor_value_info("logits", tensor.FLOAT, [1, "length", width])],
        initializer=[
            onnx.numpy_helper.from_array(scores, "scores"),
            onnx.numpy_helper.from_array(np.array([width], dtype=np.int64), "width"),
        ],
    )
    model = helper.make_model(graph, opset_imports=[helper.make_opsetid("", 13)])
    model.ir_version = 8
    path = tmp_path / "synthetic.onnx"
    onnx.save(model, path)
    return str(path), hashlib.sha256(path.read_bytes()).hexdigest()


def test_import_without_onnx_runtime():
    script = """
import sys
sys.modules['onnxruntime'] = None
from factorforge.engines.sllm.adapters import DeterministicMockSLMAdapter, SLMModelAdapter
from factorforge.engines.sllm.decoder import ConstrainedBeamDecoder
assert isinstance(DeterministicMockSLMAdapter(), SLMModelAdapter)
"""
    subprocess.run([sys.executable, "-c", script], check=True, capture_output=True)


def test_input_ids_contract_and_synonymous_sampling(tmp_path):
    path, digest = fixture_model(tmp_path)
    adapter = OnnxSLMAdapter(path, digest, CODON_TO_INDEX)
    assert isinstance(adapter, SLMModelAdapter)
    assert adapter.predict_logits(["ATG", "GCT"], "R", 2).scores.shape == (64,)
    allowed = ["AGA", "AGG"]
    np.random.seed(7)
    before = np.random.get_state()
    choices = [adapter.predict_next_codon(["ATG"], "R", allowed, seed=12) for _ in range(3)]
    assert len(set(choices)) == 1 and choices[0] in allowed
    after = np.random.get_state()
    assert before[0] == after[0] and np.array_equal(before[1], after[1])
    assert before[2:] == after[2:]
    assert adapter.predict_next_codon([], "R", allowed, decoding_mode="greedy") == "AGG"
    with pytest.raises(DeadEndError):
        adapter.predict_next_codon([], "R", [])
    with pytest.raises(ValueError, match="preserve"):
        adapter.predict_next_codon([], "R", ["ATG"])
    for temperature in [0, -1, float("nan"), float("inf")]:
        with pytest.raises(ValueError, match="Temperature"):
            adapter.predict_next_codon([], "R", allowed, temperature=temperature)


def test_legacy_paired_contract(tmp_path):
    path, _ = fixture_model(tmp_path, paired=True)
    adapter = OnnxSLMAdapter(onnx_model_path=path, full_protein_aa="MAK")
    assert adapter.predict_logits(["ATG", "GCT"], "K", 2).scores.shape == (64,)
    positional = OnnxSLMAdapter(path, "MAK")
    assert positional.model_name.startswith("onnx_")
    positional.set_target_protein("MAG")
    with pytest.raises(ValueError, match="Protein context"):
        positional.predict_logits(["ATG", "GCT"], "K", 2)


def test_checksum_and_loading_fail_closed(tmp_path):
    path, digest = fixture_model(tmp_path)
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        OnnxSLMAdapter(path, "0" * 64, CODON_TO_INDEX)
    with pytest.raises(ValueError, match="hexadecimal"):
        OnnxSLMAdapter(path, "invalid", CODON_TO_INDEX)
    with pytest.raises(FileNotFoundError):
        OnnxSLMAdapter(str(tmp_path / "absent.onnx"), digest, CODON_TO_INDEX)
    # Correct checksum does not turn an invalid artifact into successful mock inference.
    invalid = tmp_path / "invalid.onnx"
    invalid.write_bytes(b"synthetic invalid model")
    with pytest.raises(Exception) as caught:
        OnnxSLMAdapter(
            str(invalid), hashlib.sha256(invalid.read_bytes()).hexdigest(), CODON_TO_INDEX
        )
    assert not isinstance(caught.value, FileNotFoundError)


def test_vocab_remapping_and_validation(tmp_path):
    path, digest = fixture_model(tmp_path)
    shuffled = {c: 63 - i for c, i in CODON_TO_INDEX.items()}
    adapter = OnnxSLMAdapter(path, digest, shuffled)
    scores = adapter.predict_logits([], "M", 0).scores
    assert scores[CODON_TO_INDEX["ATG"]] == shuffled["ATG"]
    for vocab in [{"ATG": 0}, dict.fromkeys(CODON_TO_INDEX, 0), {**CODON_TO_INDEX, "ATG": True}]:
        with pytest.raises(ValueError, match="Vocabulary"):
            OnnxSLMAdapter(path, digest, vocab)
    with pytest.raises(ValueError, match="Unknown context"):
        adapter.predict_logits(["NNN"], "M", 1)


@pytest.mark.parametrize("width,nonfinite,match", [(63, False, "shape"), (64, True, "finite")])
def test_invalid_model_outputs(tmp_path, width, nonfinite, match):
    path, digest = fixture_model(tmp_path, width=width, nonfinite=nonfinite)
    adapter = OnnxSLMAdapter(path, digest, CODON_TO_INDEX)
    with pytest.raises(ValueError, match=match):
        adapter.predict_logits([], "M", 0)
