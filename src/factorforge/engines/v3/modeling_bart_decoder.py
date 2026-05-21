"""BART decoder skeleton for FactorForge v3."""

from __future__ import annotations

from typing import TYPE_CHECKING

try:
    import torch
    from torch import nn
    from transformers import BartConfig
    from transformers.models.bart.modeling_bart import BartDecoder
except ImportError as exc:  # pragma: no cover - exercised in ML installs
    raise ImportError(
        "ML dependencies not installed. Install with: pip install -e \".[ml]\""
    ) from exc

if TYPE_CHECKING:
    from torch import Tensor


class BartDecoderSkeleton(nn.Module):
    """Lightweight BART decoder that consumes unified encoder embeddings."""

    def __init__(
        self,
        vocab_size: int,
        d_model: int = 256,
        encoder_dim: int = 320,  # ESM2 t6_8M per-token embedding dim
        decoder_layers: int = 4,
        decoder_attention_heads: int = 4,
        ffn_dim: int = 1024,
        max_position_embeddings: int = 512,
        dropout: float = 0.1,
        pad_token_id: int = 0,
        bos_token_id: int = 1,
        eos_token_id: int = 2,
    ) -> None:
        super().__init__()
        self.config = BartConfig(
            vocab_size=vocab_size,
            d_model=d_model,
            decoder_layers=decoder_layers,
            decoder_attention_heads=decoder_attention_heads,
            decoder_ffn_dim=ffn_dim,
            max_position_embeddings=max_position_embeddings,
            encoder_layers=1,
            encoder_attention_heads=decoder_attention_heads,
            dropout=dropout,
            activation_function="gelu",
            pad_token_id=pad_token_id,
            bos_token_id=bos_token_id,
            eos_token_id=eos_token_id,
            decoder_start_token_id=bos_token_id,
            is_encoder_decoder=True,
        )

        self.encoder_projection = nn.Linear(encoder_dim, d_model)
        self.decoder = BartDecoder(self.config)
        self.lm_head = nn.Linear(d_model, vocab_size, bias=False)
        self.register_buffer("final_logits_bias", torch.zeros((1, vocab_size)))

    def forward(
        self,
        encoder_hidden_states: "Tensor",
        decoder_input_ids: "Tensor",
        attention_mask: "Tensor | None" = None,
        encoder_attention_mask: "Tensor | None" = None,
    ) -> "Tensor":
        """Teacher-forcing forward pass returning logits [B, T, vocab_size]."""
        encoder_states = self.encoder_projection(encoder_hidden_states)
        outputs = self.decoder(
            input_ids=decoder_input_ids,
            attention_mask=attention_mask,
            encoder_hidden_states=encoder_states,
            encoder_attention_mask=encoder_attention_mask,
            use_cache=False,
            return_dict=True,
        )
        logits = self.lm_head(outputs.last_hidden_state) + self.final_logits_bias.to(
            outputs.last_hidden_state.device
        )
        return logits

    def generate(
        self,
        encoder_hidden_states: "Tensor",
        max_new_tokens: int,
        bos_token_id: int,
        eos_token_id: int,
        pad_token_id: int | None = None,
    ) -> "Tensor":
        """Greedy decoding for a few steps (no beam search)."""
        batch_size = encoder_hidden_states.shape[0]
        device = encoder_hidden_states.device

        decoder_input_ids = torch.full(
            (batch_size, 1),
            bos_token_id,
            dtype=torch.long,
            device=device,
        )

        for _ in range(max_new_tokens):
            logits = self.forward(
                encoder_hidden_states=encoder_hidden_states,
                decoder_input_ids=decoder_input_ids,
                attention_mask=None,
                encoder_attention_mask=None,
            )
            next_token = torch.argmax(logits[:, -1, :], dim=-1)
            decoder_input_ids = torch.cat(
                [decoder_input_ids, next_token.unsqueeze(1)],
                dim=1,
            )
            if torch.all(next_token == eos_token_id):
                break

        if pad_token_id is not None:
            decoder_input_ids = decoder_input_ids.clone()
            for idx in range(decoder_input_ids.shape[0]):
                row = decoder_input_ids[idx]
                eos_positions = (row == eos_token_id).nonzero(as_tuple=False)
                if eos_positions.numel() > 0:
                    eos_position = int(eos_positions[0])
                    if eos_position + 1 < row.shape[0]:
                        row[eos_position + 1 :] = pad_token_id
        return decoder_input_ids
