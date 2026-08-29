"""FactorForge-LM mBART Model Architecture (v3.5.0 ML Track).

Implements a Seq2Seq Encoder-Decoder Transformer architecture tailored for
host-aware codon generation and infilling, supporting control-token conditioning
and FactorForge logit masking.
"""

from typing import Optional, Tuple, Dict, Any, List
import math
try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
    class DummyTorch:
        Tensor = Any
    torch = DummyTorch()
    nn = type("nn", (), {"Module": object})()



from factorforge.engines.lm.tokenizer.control_tokenizer import FactorForgeControlTokenizer


class SynCodonLMConfig:
    """Configuration class for SynCodonLM / FactorForge-LM."""

    def __init__(
        self,
        vocab_size: int = 98,
        d_model: int = 256,
        encoder_layers: int = 6,
        decoder_layers: int = 6,
        encoder_attention_heads: int = 8,
        decoder_attention_heads: int = 8,
        d_ff: int = 1024,
        dropout: float = 0.1,
        max_position_embeddings: int = 1024,
        pad_token_id: int = 0,
        bos_token_id: int = 2,
        eos_token_id: int = 3,
    ):
        self.vocab_size = vocab_size
        self.d_model = d_model
        self.encoder_layers = encoder_layers
        self.decoder_layers = decoder_layers
        self.encoder_attention_heads = encoder_attention_heads
        self.decoder_attention_heads = decoder_attention_heads
        self.d_ff = d_ff
        self.dropout = dropout
        self.max_position_embeddings = max_position_embeddings
        self.pad_token_id = pad_token_id
        self.bos_token_id = bos_token_id
        self.eos_token_id = eos_token_id


class PositionalEncoding(nn.Module):
    """Sinusoidal positional encoding."""

    def __init__(self, d_model: int, max_len: int = 1024):
        super().__init__()
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer("pe", pe.unsqueeze(0))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.pe[:, :x.size(1)]


class SynCodonLMModel(nn.Module):
    """Seq2Seq Transformer Encoder-Decoder for Codon Optimization (SynCodonLM-V2)."""

    def __init__(self, config: SynCodonLMConfig):
        super().__init__()
        self.config = config

        self.token_embedding = nn.Embedding(config.vocab_size, config.d_model, padding_idx=config.pad_token_id)
        self.pos_encoding = PositionalEncoding(config.d_model, config.max_position_embeddings)
        self.dropout = nn.Dropout(config.dropout)

        # PyTorch Transformer
        self.transformer = nn.Transformer(
            d_model=config.d_model,
            nhead=config.encoder_attention_heads,
            num_encoder_layers=config.encoder_layers,
            num_decoder_layers=config.decoder_layers,
            dim_feedforward=config.d_ff,
            dropout=config.dropout,
            batch_first=True,
        )

        # Output LM Head mapping hidden states to vocab logits
        self.lm_head = nn.Linear(config.d_model, config.vocab_size)

    def forward(
        self,
        input_ids: torch.Tensor,
        decoder_input_ids: torch.Tensor,
        src_key_padding_mask: Optional[torch.Tensor] = None,
        tgt_key_padding_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """Forward pass generating vocabulary logits."""
        # Embeddings + Positional Encodings
        src_emb = self.dropout(self.pos_encoding(self.token_embedding(input_ids)))
        tgt_emb = self.dropout(self.pos_encoding(self.token_embedding(decoder_input_ids)))

        # Target Causal Mask for Decoder
        tgt_len = decoder_input_ids.size(1)
        tgt_mask = torch.triu(torch.full((tgt_len, tgt_len), float("-inf"), device=input_ids.device), diagonal=1)

        # Transformer Pass
        out = self.transformer(
            src=src_emb,
            tgt=tgt_emb,
            tgt_mask=tgt_mask,
            src_key_padding_mask=src_key_padding_mask,
            tgt_key_padding_mask=tgt_key_padding_mask,
        )

        logits = self.lm_head(out)
        return logits


if __name__ == "__main__":
    tokenizer = FactorForgeControlTokenizer()
    config = SynCodonLMConfig(vocab_size=tokenizer.vocab_size)
    print(f"SynCodonLM Config initialized. Vocab size: {config.vocab_size}, Model dim: {config.d_model}")

    if TORCH_AVAILABLE:
        model = SynCodonLMModel(config)
        src = torch.tensor([[2, 5, 8, 11, 23, 13, 21, 31, 3]])
        tgt = torch.tensor([[2, 45, 62, 18]])
        logits = model(src, tgt)
        print(f"SynCodonLM PyTorch Model initialized successfully!")
        print(f"Input batch shape: {src.shape}, Output logits shape: {logits.shape}")
        assert logits.shape == (1, 4, tokenizer.vocab_size), "Logits shape mismatch!"
        print("SynCodonLM Forward Pass Verified Successfully!")
    else:
        print("PyTorch not installed in environment - SynCodonLMConfig structure verified successfully!")


