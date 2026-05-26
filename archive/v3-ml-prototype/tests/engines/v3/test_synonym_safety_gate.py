"""Pre-Run-4 synonym-mask and constrained-decoding safety tests."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

torch = pytest.importorskip("torch")

ROOT = Path(__file__).resolve().parents[3]
TRAINING = ROOT / "scripts" / "3_training"
if str(TRAINING) not in sys.path:
    sys.path.insert(0, str(TRAINING))

from dataset import CodonDataset, collate_fn  # noqa: E402
from factorforge.engines.v3.inference.constrained_decoder import (  # noqa: E402
    constrained_greedy_decode,
    validate_candidate_or_fallback,
)
from factorforge.engines.v3.synonym_mask import (  # noqa: E402
    build_synonym_token_mask,
    synonymous_codons_for_aa,
)
from factorforge.engines.v3.tokenizer import CodonTokenizer  # noqa: E402
from factorforge.ml.metrics import translate_dna  # noqa: E402
from train_v3_esm2_bart import (  # noqa: E402
    _codon_gc_vector,
    _codon_log_cai_vector,
    _codon_token_mask,
    _compute_expected_gc,
    _compute_expected_log_cai,
)


class FixedLogitDecoder:
    """Decoder that always favors an incompatible token unless masked."""

    def __init__(self, tokenizer: CodonTokenizer) -> None:
        self.tokenizer = tokenizer

    def __call__(self, encoder_hidden_states, decoder_input_ids):
        batch, length = decoder_input_ids.shape
        logits = torch.full((batch, length, len(self.tokenizer.token_to_id)), -10.0)
        logits[:, -1, self.tokenizer.token_to_id["TAA"]] = 100.0
        logits[:, -1, self.tokenizer.token_to_id["ATG"]] = 5.0
        logits[:, -1, self.tokenizer.token_to_id["GCC"]] = 4.0
        return logits


def test_synonym_mask_allows_only_synonymous_codons() -> None:
    tokenizer = CodonTokenizer.default()
    mask = build_synonym_token_mask("MA", tokenizer.token_to_id)

    allowed_m = {
        tokenizer.id_to_token[index]
        for index, allowed in enumerate(mask[0].tolist())
        if allowed
    }
    allowed_a = {
        tokenizer.id_to_token[index]
        for index, allowed in enumerate(mask[1].tolist())
        if allowed
    }

    assert allowed_m == {"ATG"}
    assert allowed_a == set(synonymous_codons_for_aa("A"))


def test_synonym_mask_excludes_special_tokens_and_standard_stops() -> None:
    tokenizer = CodonTokenizer.default()
    mask = build_synonym_token_mask("A", tokenizer.token_to_id)

    for token in tokenizer.special_tokens:
        assert not bool(mask[0, tokenizer.token_to_id[token]])
    for stop in ("TAA", "TAG", "TGA"):
        assert not bool(mask[0, tokenizer.token_to_id[stop]])


def test_dataset_collate_carries_synonym_mask(tmp_path: Path) -> None:
    tokenizer = CodonTokenizer.default()
    embeddings_dir = tmp_path / "emb"
    embeddings_dir.mkdir()
    torch.save({"embeddings": torch.ones((2, 320))}, embeddings_dir / "p1.pt")
    training_jsonl = tmp_path / "train.jsonl"
    training_jsonl.write_text(
        json.dumps(
            {
                "protein_id": "p1",
                "sequence": "MA",
                "codon_sequence": "ATG GCC",
            }
        )
        + "\n",
        encoding="utf-8",
    )

    dataset = CodonDataset(
        str(training_jsonl),
        str(embeddings_dir),
        tokenizer.token_to_id,
        tokenizer.bos_token_id,
        tokenizer.eos_token_id,
        tokenizer.unk_token_id,
    )
    batch = collate_fn([dataset[0]])

    assert batch["synonym_mask"].shape == (1, 3, len(tokenizer.token_to_id))
    assert bool(batch["synonym_mask"][0, 0, tokenizer.token_to_id["ATG"]])
    assert not bool(batch["synonym_mask"][0, 0, tokenizer.token_to_id["GCC"]])
    assert not bool(batch["synonym_mask"][0, 2].any())


def test_expected_metrics_with_synonym_mask_ignore_incompatible_codons() -> None:
    tokenizer = CodonTokenizer.default()
    logits = torch.full((1, 1, len(tokenizer.token_to_id)), -10.0)
    logits[0, 0, tokenizer.token_to_id["GCC"]] = 100.0
    logits[0, 0, tokenizer.token_to_id["ATG"]] = 1.0
    labels = torch.tensor([[tokenizer.token_to_id["ATG"]]])
    synonym_mask = build_synonym_token_mask("M", tokenizer.token_to_id).unsqueeze(0)
    token_mask = _codon_token_mask(tokenizer, torch.device("cpu"))

    expected_gc = _compute_expected_gc(
        logits,
        labels,
        _codon_gc_vector(tokenizer, torch.device("cpu")),
        token_mask,
        synonym_mask,
    )
    expected_log_cai = _compute_expected_log_cai(
        logits,
        labels,
        _codon_log_cai_vector(tokenizer, torch.device("cpu")),
        token_mask,
        synonym_mask,
    )

    assert expected_gc.item() == pytest.approx(1 / 3)
    assert torch.isfinite(expected_log_cai)


def test_constrained_decoding_preserves_amino_acid_sequence() -> None:
    tokenizer = CodonTokenizer.default()
    decoder = FixedLogitDecoder(tokenizer)
    ids = constrained_greedy_decode(decoder, torch.zeros((1, 2, 320)), "MA", tokenizer)
    dna = tokenizer.decode(ids[0].tolist())

    assert len(dna) == 6
    assert translate_dna(dna) == "MA"


def test_invalid_candidate_falls_back_to_v2() -> None:
    result = validate_candidate_or_fallback("MA", "GCCGCC")

    assert result["engine"] == "v2"
    assert result["fallback_used"] is True
    assert result["validator"]["passed"] is True
    assert translate_dna(result["dna_sequence"]) == "MA"

