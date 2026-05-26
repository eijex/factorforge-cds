"""Synonymous codon masks for v3 training and constrained decoding."""

from __future__ import annotations

from typing import Mapping

import torch

from factorforge.ml.metrics import STANDARD_GENETIC_CODE


AA_TO_CODONS: dict[str, tuple[str, ...]] = {}
for _codon, _aa in STANDARD_GENETIC_CODE.items():
    AA_TO_CODONS.setdefault(_aa, tuple())
for _aa in list(AA_TO_CODONS):
    AA_TO_CODONS[_aa] = tuple(
        codon for codon, codon_aa in STANDARD_GENETIC_CODE.items() if codon_aa == _aa
    )


def normalize_protein_sequence(protein_sequence: str) -> str:
    """Normalize a protein sequence for codon-mask construction."""
    protein = "".join(protein_sequence.upper().split())
    if protein.endswith("*"):
        protein = protein[:-1]
    if not protein:
        raise ValueError("protein_sequence must not be empty")
    invalid = [aa for aa in protein if aa not in AA_TO_CODONS or aa == "*"]
    if invalid:
        raise ValueError(f"No synonymous codons for amino acid(s): {''.join(sorted(set(invalid)))}")
    return protein


def synonymous_codons_for_aa(amino_acid: str) -> tuple[str, ...]:
    """Return synonymous non-stop codons for one standard amino acid."""
    aa = amino_acid.upper()
    codons = AA_TO_CODONS.get(aa, tuple())
    if not codons or aa == "*":
        raise ValueError(f"No synonymous codons for amino acid: {amino_acid}")
    return codons


def build_synonym_token_mask(
    protein_sequence: str,
    token_to_id: Mapping[str, int],
    *,
    device: torch.device | None = None,
) -> torch.Tensor:
    """Build a boolean mask [protein_length, vocab_size] for synonymous codons only."""
    protein = normalize_protein_sequence(protein_sequence)
    vocab_size = max(token_to_id.values()) + 1
    mask = torch.zeros((len(protein), vocab_size), dtype=torch.bool, device=device)
    for index, aa in enumerate(protein):
        for codon in synonymous_codons_for_aa(aa):
            token_id = token_to_id.get(codon)
            if token_id is not None:
                mask[index, token_id] = True
        if not bool(mask[index].any()):
            raise ValueError(f"No tokenizer codons available for amino acid: {aa}")
    return mask

