"""Tests for differentiable v3 training loss components."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
torch = pytest.importorskip("torch", reason="torch not installed — skipping ML tests")

ROOT = Path(__file__).resolve().parents[3]
TRAINING = ROOT / "scripts" / "3_training"
if str(TRAINING) not in sys.path:
    sys.path.insert(0, str(TRAINING))

from factorforge.engines.v3.tokenizer import CodonTokenizer
from train_v3_esm2_bart import (  # noqa: E402
    _codon_gc_vector,
    _codon_log_cai_vector,
    _codon_token_mask,
    _compute_bounded_gc_penalty,
    _compute_expected_gc,
    _compute_expected_log_cai,
)


def _logits_for(tokenizer: CodonTokenizer, tokens: list[str], high: float = 80.0) -> torch.Tensor:
    logits = torch.full((1, len(tokens), len(tokenizer.token_to_id)), -80.0)
    for index, token in enumerate(tokens):
        logits[0, index, tokenizer.token_to_id[token]] = high
    return logits


def test_expected_gc_deterministic_logits() -> None:
    tokenizer = CodonTokenizer.default()
    logits = _logits_for(tokenizer, ["ATG", "GCC"])
    labels = torch.tensor([[tokenizer.token_to_id["ATG"], tokenizer.token_to_id["GCC"]]])

    expected_gc = _compute_expected_gc(
        logits,
        labels,
        _codon_gc_vector(tokenizer, logits.device),
        _codon_token_mask(tokenizer, logits.device),
    )

    assert expected_gc.item() == pytest.approx((1 / 3 + 1.0) / 2)


def test_expected_gc_changes_when_logits_favor_high_gc_codons() -> None:
    tokenizer = CodonTokenizer.default()
    labels = torch.tensor([[tokenizer.token_to_id["GCC"]]])
    codon_gc = _codon_gc_vector(tokenizer, torch.device("cpu"))
    mask = _codon_token_mask(tokenizer, torch.device("cpu"))

    low_gc = _compute_expected_gc(_logits_for(tokenizer, ["AAA"]), labels, codon_gc, mask)
    high_gc = _compute_expected_gc(_logits_for(tokenizer, ["GCC"]), labels, codon_gc, mask)

    assert high_gc > low_gc


def test_bounded_gc_penalty() -> None:
    inside = _compute_bounded_gc_penalty(torch.tensor(0.45), gc_low=0.40, gc_high=0.55)
    below = _compute_bounded_gc_penalty(torch.tensor(0.30), gc_low=0.40, gc_high=0.55)
    above = _compute_bounded_gc_penalty(torch.tensor(0.70), gc_low=0.40, gc_high=0.55)

    assert inside.item() == 0.0
    assert below.item() > 0.0
    assert above.item() > 0.0


def test_expected_log_cai_increases_for_preferred_codons() -> None:
    tokenizer = CodonTokenizer.default()
    labels = torch.tensor([[tokenizer.token_to_id["GCC"]]])
    mask = _codon_token_mask(tokenizer, torch.device("cpu"))
    log_cai = torch.full((len(tokenizer.token_to_id),), -20.0)
    log_cai[tokenizer.token_to_id["GCC"]] = 0.0
    log_cai[tokenizer.token_to_id["GCT"]] = torch.log(torch.tensor(0.25))

    low = _compute_expected_log_cai(_logits_for(tokenizer, ["GCT"]), labels, log_cai, mask)
    high = _compute_expected_log_cai(_logits_for(tokenizer, ["GCC"]), labels, log_cai, mask)

    assert high > low


def test_gradients_flow_back_to_logits() -> None:
    tokenizer = CodonTokenizer.default()
    logits = _logits_for(tokenizer, ["AAA"], high=8.0).requires_grad_(True)
    labels = torch.tensor([[tokenizer.token_to_id["AAA"]]])
    expected_gc = _compute_expected_gc(
        logits,
        labels,
        _codon_gc_vector(tokenizer, logits.device),
        _codon_token_mask(tokenizer, logits.device),
    )
    loss = _compute_bounded_gc_penalty(expected_gc, gc_low=0.40, gc_high=0.55)
    loss.backward()

    assert logits.grad is not None
    assert torch.any(logits.grad != 0)


def test_special_tokens_and_padding_do_not_affect_expected_gc_or_cai() -> None:
    tokenizer = CodonTokenizer.default()
    logits = _logits_for(tokenizer, ["GCC", "AAA"], high=20.0)
    logits[0, 0, tokenizer.pad_token_id] = 100.0
    logits[0, 1, tokenizer.token_to_id["GCC"]] = 100.0
    labels = torch.tensor([[tokenizer.token_to_id["GCC"], -100]])
    mask = _codon_token_mask(tokenizer, logits.device)

    expected_gc = _compute_expected_gc(logits, labels, _codon_gc_vector(tokenizer, logits.device), mask)
    expected_log_cai = _compute_expected_log_cai(
        logits,
        labels,
        _codon_log_cai_vector(tokenizer, logits.device),
        mask,
    )

    assert expected_gc.item() == pytest.approx(1.0)
    assert torch.isfinite(expected_log_cai)

