"""Constrained v3 decoding and fallback helpers."""

from __future__ import annotations

from typing import Any

import torch

from factorforge.engines.v3.inference.v2_adapter import optimize_with_v2
from factorforge.engines.v3.synonym_mask import build_synonym_token_mask, normalize_protein_sequence
from factorforge.engines.v3.tokenizer import CodonTokenizer
from factorforge.utils.validation import validate_candidate_sequence


def constrained_greedy_decode(
    decoder: Any,
    encoder_hidden_states: torch.Tensor,
    protein_sequence: str,
    codon_tokenizer: CodonTokenizer,
) -> torch.Tensor:
    """Greedy decode one codon per amino acid under a synonym mask."""
    protein = normalize_protein_sequence(protein_sequence)
    mask = build_synonym_token_mask(
        protein,
        codon_tokenizer.token_to_id,
        device=encoder_hidden_states.device,
    )
    batch_size = int(encoder_hidden_states.shape[0])
    if batch_size != 1:
        raise ValueError("constrained_greedy_decode currently supports batch_size=1")

    decoder_input_ids = torch.full(
        (1, 1),
        codon_tokenizer.bos_token_id,
        dtype=torch.long,
        device=encoder_hidden_states.device,
    )
    generated: list[torch.Tensor] = []
    for position in range(len(protein)):
        logits = decoder(
            encoder_hidden_states=encoder_hidden_states,
            decoder_input_ids=decoder_input_ids,
        )
        next_logits = logits[:, -1, :].masked_fill(~mask[position].unsqueeze(0), -1.0e9)
        next_token = torch.argmax(next_logits, dim=-1)
        generated.append(next_token)
        decoder_input_ids = torch.cat([decoder_input_ids, next_token.unsqueeze(1)], dim=1)

    eos = torch.tensor(
        [[codon_tokenizer.eos_token_id]],
        dtype=torch.long,
        device=encoder_hidden_states.device,
    )
    return torch.cat([decoder_input_ids, eos], dim=1)


def validate_candidate_or_fallback(
    protein_sequence: str,
    dna_sequence: str,
    fallback_options: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return a valid v3 candidate or deterministic v2 fallback metadata."""
    protein = normalize_protein_sequence(protein_sequence)
    validator = validate_candidate_sequence(protein, dna_sequence)
    if validator["passed"]:
        return {
            "engine": "v3",
            "protein_sequence": protein,
            "dna_sequence": dna_sequence,
            "validator": validator,
            "fallback_used": False,
            "warnings": list(validator["warnings"]),
            "errors": [],
        }

    fallback = optimize_with_v2(protein, options=fallback_options)
    fallback["fallback_used"] = True
    fallback.setdefault("metadata", {})["fallback_reason"] = validator["errors"]
    return fallback

