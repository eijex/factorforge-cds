"""BART decoder skeleton tests for v3."""

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("transformers")
pytest.importorskip("transformers.models.bart.modeling_bart")

from factorforge.engines.v3.modeling_bart_decoder import BartDecoderSkeleton


def test_decoder_output_shapes():
    vocab_size = 97
    model = BartDecoderSkeleton(vocab_size=vocab_size)
    batch_size = 2
    encoder_len = 5
    target_len = 7

    encoder_hidden_states = torch.randn(batch_size, encoder_len, 320)
    decoder_input_ids = torch.randint(0, vocab_size, (batch_size, target_len))

    logits = model(
        encoder_hidden_states=encoder_hidden_states,
        decoder_input_ids=decoder_input_ids,
    )

    assert logits.shape == (batch_size, target_len, vocab_size)


def test_generate_greedy():
    vocab_size = 97
    model = BartDecoderSkeleton(vocab_size=vocab_size)
    batch_size = 2
    encoder_len = 4

    encoder_hidden_states = torch.randn(batch_size, encoder_len, 320)
    generated = model.generate(
        encoder_hidden_states=encoder_hidden_states,
        max_new_tokens=6,
        bos_token_id=1,
        eos_token_id=2,
        pad_token_id=0,
    )

    assert generated.shape[0] == batch_size
    assert generated.shape[1] >= 2
