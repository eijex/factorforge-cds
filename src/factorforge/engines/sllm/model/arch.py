# factorforge/src/factorforge/engines/sllm/model/arch.py
"""Compact PyTorch Causal Codon Transformer Architecture (Job 285B).

Lightweight autoregressive decoder predicting 64 codon logits conditioned on
protein amino acid tokens and prior codon history.
"""

from __future__ import annotations
import math
import torch
import torch.nn as nn
import torch.nn.functional as F

from factorforge.engines.sllm.data.tokenizer import (
    AA_VOCAB_SIZE,
    CODON_VOCAB_SIZE,
    PAD_AA_ID,
    PAD_CODON_ID,
)


class MultiHeadAttentionBlock(nn.Module):
    """Clean multi-head attention module with guaranteed dynamic shape ONNX exportability."""

    def __init__(self, d_model: int, nhead: int, dropout: float = 0.1) -> None:
        super().__init__()
        self.d_model = d_model
        self.nhead = nhead
        self.head_dim = d_model // nhead
        assert self.head_dim * nhead == d_model, (
            f"d_model ({d_model}) must be divisible by nhead ({nhead})"
        )

        self.q_proj = nn.Linear(d_model, d_model)
        self.k_proj = nn.Linear(d_model, d_model)
        self.v_proj = nn.Linear(d_model, d_model)
        self.out_proj = nn.Linear(d_model, d_model)
        self.dropout = nn.Dropout(dropout)

    def forward(
        self, q_x: torch.Tensor, kv_x: torch.Tensor, is_causal: bool = False
    ) -> torch.Tensor:
        B, L_q, _ = q_x.shape
        _, L_kv, _ = kv_x.shape

        q = self.q_proj(q_x).view(B, L_q, self.nhead, self.head_dim).transpose(1, 2)
        k = self.k_proj(kv_x).view(B, L_kv, self.nhead, self.head_dim).transpose(1, 2)
        v = self.v_proj(kv_x).view(B, L_kv, self.nhead, self.head_dim).transpose(1, 2)

        scores = torch.matmul(q, k.transpose(-2, -1)) / math.sqrt(self.head_dim)

        if is_causal:
            mask = torch.triu(
                torch.full((L_q, L_kv), float("-10000.0"), device=q_x.device), diagonal=1
            )
            scores = scores + mask

        weights = F.softmax(scores, dim=-1)
        weights = self.dropout(weights)

        attn = torch.matmul(weights, v)
        attn = attn.transpose(1, 2).contiguous().view(B, L_q, self.d_model)
        return self.out_proj(attn)


class TransformerDecoderLayer(nn.Module):
    """Transformer Decoder Block with Causal Self-Attention and Cross-Attention."""

    def __init__(
        self, d_model: int, nhead: int, dim_feedforward: int, dropout: float = 0.1
    ) -> None:
        super().__init__()
        self.self_attn = MultiHeadAttentionBlock(d_model, nhead, dropout=dropout)
        self.cross_attn = MultiHeadAttentionBlock(d_model, nhead, dropout=dropout)
        self.linear1 = nn.Linear(d_model, dim_feedforward)
        self.linear2 = nn.Linear(dim_feedforward, d_model)
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)
        self.norm3 = nn.LayerNorm(d_model)
        self.dropout = nn.Dropout(dropout)
        self.dropout1 = nn.Dropout(dropout)
        self.dropout2 = nn.Dropout(dropout)
        self.dropout3 = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor, memory: torch.Tensor) -> torch.Tensor:
        sa_out = self.self_attn(x, x, is_causal=True)
        x = self.norm1(x + self.dropout1(sa_out))
        ca_out = self.cross_attn(x, memory, is_causal=False)
        x = self.norm2(x + self.dropout2(ca_out))
        ff_out = self.linear2(self.dropout(F.gelu(self.linear1(x))))
        x = self.norm3(x + self.dropout3(ff_out))
        return x


class CompactCodonTransformer(nn.Module):
    """Compact Causal Transformer for Host Codon-Context Prior."""

    def __init__(
        self,
        d_model: int = 128,
        nhead: int = 4,
        num_layers: int = 3,
        dim_feedforward: int = 256,
        max_len: int = 1024,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.d_model = d_model
        self.max_len = max_len

        self.codon_embedding = nn.Embedding(CODON_VOCAB_SIZE, d_model, padding_idx=PAD_CODON_ID)
        self.aa_embedding = nn.Embedding(AA_VOCAB_SIZE, d_model, padding_idx=PAD_AA_ID)
        self.pos_embedding = nn.Embedding(max_len, d_model)

        self.layers = nn.ModuleList(
            [
                TransformerDecoderLayer(d_model, nhead, dim_feedforward, dropout)
                for _ in range(num_layers)
            ]
        )

        self.output_head = nn.Linear(d_model, 64)
        self.layer_norm = nn.LayerNorm(d_model)
        self._init_weights()

    def _init_weights(self) -> None:
        for p in self.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)

    def forward(
        self,
        codon_ids: torch.Tensor,  # Shape: (B, L)
        aa_ids: torch.Tensor,  # Shape: (B, L)
    ) -> torch.Tensor:
        """Forward pass generating 64 codon logits per position."""
        B, L = codon_ids.shape
        device = codon_ids.device

        positions = torch.arange(0, L, device=device).unsqueeze(0).expand(B, L)
        tgt = self.codon_embedding(codon_ids) + self.pos_embedding(positions)
        memory = self.aa_embedding(aa_ids) + self.pos_embedding(positions)

        x = tgt
        for layer in self.layers:
            x = layer(x, memory)

        x = self.layer_norm(x)
        logits = self.output_head(x)  # (B, L, 64)
        return logits


def compute_synonymous_masked_loss(
    logits: torch.Tensor,  # (B, L, 64)
    target_codons: torch.Tensor,  # (B, L)
    synonymous_masks: torch.Tensor,  # (B, L, 64) boolean mask
    pad_mask: torch.Tensor,  # (B, L) boolean mask (True for valid tokens)
) -> torch.Tensor:
    """Computes Cross-Entropy loss strictly over valid synonymous codons.

    Invalid non-synonymous codons receive -1e4 logit penalty so they do not
    compete in the denominator of Softmax.
    """
    mask_penalty = torch.where(synonymous_masks, 0.0, -10000.0)
    masked_logits = logits + mask_penalty

    B, L, C = masked_logits.shape
    flat_logits = masked_logits.view(-1, C)
    flat_targets = target_codons.view(-1)
    flat_pad = pad_mask.view(-1)

    active_logits = flat_logits[flat_pad]
    active_targets = flat_targets[flat_pad]

    if active_targets.numel() == 0:
        return torch.tensor(0.0, device=logits.device, requires_grad=True)

    loss = F.cross_entropy(active_logits, active_targets)
    return loss
