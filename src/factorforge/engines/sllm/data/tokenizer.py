# factorforge/src/factorforge/engines/sllm/data/tokenizer.py
"""Codon and Amino Acid sequence tokenizer for FactorForge sLLM (Job 285A).

Provides:
1. Exact mapping of 64 standard codons to vocabulary IDs [0..63].
2. Standard 20 amino acids mapping to IDs [0..19].
3. Special tokens: PAD, BOS, EOS.
4. Synonymous codon index masks for constrained loss and inference.
"""

from __future__ import annotations
from typing import Dict, List, Tuple
import numpy as np

from factorforge.engines.sllm.interfaces import (
    AA_TO_INDICES,
    CODON_TO_INDEX,
)

# Amino acid vocabulary (20 standard AAs)
AA_VOCAB: Tuple[str, ...] = (
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
)
AA_TO_INDEX: Dict[str, int] = {aa: idx for idx, aa in enumerate(AA_VOCAB)}

# Special tokens
PAD_CODON_ID: int = 64
BOS_CODON_ID: int = 65
EOS_CODON_ID: int = 66

PAD_AA_ID: int = 20
BOS_AA_ID: int = 21
EOS_AA_ID: int = 22

CODON_VOCAB_SIZE: int = 67
AA_VOCAB_SIZE: int = 23


class CodonSequenceTokenizer:
    """Bidirectional tokenizer between DNA CDS, protein sequences, and discrete token IDs."""

    def __init__(self) -> None:
        self.codon_to_id = dict(CODON_TO_INDEX)
        self.id_to_codon = {idx: c for c, idx in self.codon_to_id.items()}

        self.aa_to_id = dict(AA_TO_INDEX)
        self.id_to_aa = {idx: a for a, idx in self.aa_to_id.items()}

        # Build lookup table: AA ID -> binary mask array of shape (64,)
        self._synonymous_masks: Dict[int, np.ndarray] = {}
        for aa, aa_id in self.aa_to_id.items():
            mask = np.zeros(64, dtype=bool)
            valid_indices = AA_TO_INDICES.get(aa, [])
            mask[valid_indices] = True
            self._synonymous_masks[aa_id] = mask

    def encode_codons(
        self,
        cds_sequence: str,
        add_bos: bool = False,
        add_eos: bool = False,
    ) -> List[int]:
        """Encodes DNA CDS string into list of codon token IDs (each 3-nt)."""
        clean_cds = "".join(cds_sequence.upper().split())
        if len(clean_cds) % 3 != 0:
            raise ValueError(f"CDS length ({len(clean_cds)}) must be divisible by 3.")

        tokens: List[int] = []
        if add_bos:
            tokens.append(BOS_CODON_ID)

        for i in range(0, len(clean_cds), 3):
            codon = clean_cds[i : i + 3]
            if codon not in self.codon_to_id:
                raise ValueError(
                    f"Invalid non-standard or ambiguous codon '{codon}' at offset {i}."
                )
            tokens.append(self.codon_to_id[codon])

        if add_eos:
            tokens.append(EOS_CODON_ID)

        return tokens

    def decode_codons(
        self,
        token_ids: List[int],
        skip_special: bool = True,
    ) -> str:
        """Decodes list of codon token IDs into DNA sequence string."""
        codons: List[str] = []
        for tid in token_ids:
            if tid in self.id_to_codon:
                codons.append(self.id_to_codon[tid])
            elif not skip_special:
                if tid == BOS_CODON_ID:
                    codons.append("<BOS>")
                elif tid == EOS_CODON_ID:
                    codons.append("<EOS>")
                elif tid == PAD_CODON_ID:
                    codons.append("<PAD>")
        return "".join(codons)

    def encode_amino_acids(
        self,
        aa_sequence: str,
        add_bos: bool = False,
        add_eos: bool = False,
    ) -> List[int]:
        """Encodes single-letter amino acid string into list of AA token IDs."""
        clean_aa = "".join(aa_sequence.upper().split()).rstrip("*")
        tokens: List[int] = []
        if add_bos:
            tokens.append(BOS_AA_ID)

        for i, aa in enumerate(clean_aa):
            if aa not in self.aa_to_id:
                raise ValueError(f"Invalid amino acid '{aa}' at position {i}.")
            tokens.append(self.aa_to_id[aa])

        if add_eos:
            tokens.append(EOS_AA_ID)

        return tokens

    def decode_amino_acids(
        self,
        token_ids: List[int],
        skip_special: bool = True,
    ) -> str:
        """Decodes list of AA token IDs into amino acid string."""
        aas: List[str] = []
        for tid in token_ids:
            if tid in self.id_to_aa:
                aas.append(self.id_to_aa[tid])
            elif not skip_special:
                if tid == BOS_AA_ID:
                    aas.append("<BOS>")
                elif tid == EOS_AA_ID:
                    aas.append("<EOS>")
                elif tid == PAD_AA_ID:
                    aas.append("<PAD>")
        return "".join(aas)

    def get_synonymous_mask(self, aa_id: int) -> np.ndarray:
        """Returns boolean mask of shape (64,) where True indicates valid synonymous codon."""
        return self._synonymous_masks.get(aa_id, np.zeros(64, dtype=bool))
