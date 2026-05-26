"""FactorForge v3 training dataset for per-token ESM2 embeddings."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import torch
from torch.nn.utils.rnn import pad_sequence
from torch.utils.data import Dataset

from factorforge.engines.v3.synonym_mask import build_synonym_token_mask


class CodonDataset(Dataset):
    """Load paired ESM2 per-token embeddings and codon token IDs."""

    def __init__(
        self,
        training_jsonl: str,
        embeddings_dir: str,
        codon_to_id: dict[str, int],
        bos_token_id: int,
        eos_token_id: int,
        unk_token_id: int,
        max_length: int = 512,
    ) -> None:
        self.embeddings_dir = Path(embeddings_dir)
        self.codon_to_id = codon_to_id
        self.bos_token_id = bos_token_id
        self.eos_token_id = eos_token_id
        self.unk_token_id = unk_token_id
        self.max_length = max_length
        self.samples: list[dict[str, Any]] = []

        with open(training_jsonl, encoding="utf-8") as handle:
            for line in handle:
                item = json.loads(line.strip())
                emb_path = self.embeddings_dir / f"{item['protein_id']}.pt"
                if emb_path.exists():
                    self.samples.append(item)

        print(f"Dataset: {len(self.samples)} samples loaded")

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, idx: int) -> dict[str, torch.Tensor]:
        sample = self.samples[idx]
        emb_data = torch.load(
            self.embeddings_dir / f"{sample['protein_id']}.pt",
            map_location="cpu",
            weights_only=True,
        )
        embeddings = emb_data["embeddings"]

        codon_sequence = sample.get("codon_sequence") or sample.get("dna_seq")
        if not isinstance(codon_sequence, str):
            raise ValueError("Sample must include 'codon_sequence' or 'dna_seq'")

        codon_ids = self._encode_codons(codon_sequence)[: self.max_length]
        embeddings = embeddings[: len(codon_ids)]
        protein_sequence = self._extract_protein_sequence(sample)
        protein_for_mask = protein_sequence[: len(codon_ids)]
        if len(protein_for_mask) != len(codon_ids):
            raise ValueError("Protein and codon sequence lengths do not match for synonym mask")
        codon_mask = build_synonym_token_mask(protein_for_mask, self.codon_to_id)
        eos_mask = torch.zeros((1, len(self.codon_to_id)), dtype=torch.bool)
        synonym_mask = torch.cat([codon_mask, eos_mask], dim=0)

        decoder_input_ids = torch.tensor(
            [self.bos_token_id, *codon_ids],
            dtype=torch.long,
        )
        labels = torch.tensor([*codon_ids, self.eos_token_id], dtype=torch.long)

        return {
            "encoder_hidden_states": embeddings,
            "decoder_input_ids": decoder_input_ids,
            "labels": labels,
            "synonym_mask": synonym_mask,
        }

    def _encode_codons(self, sequence: str) -> list[int]:
        seq = sequence.strip().upper().replace("U", "T").replace(" ", "")
        codons = [seq[i : i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]
        return [self.codon_to_id.get(codon, self.unk_token_id) for codon in codons]

    def _extract_protein_sequence(self, sample: dict[str, Any]) -> str:
        sequence = (
            sample.get("protein_sequence")
            or sample.get("protein_seq")
            or sample.get("sequence")
        )
        if not isinstance(sequence, str) or not sequence.strip():
            raise ValueError(
                "Sample must include 'protein_sequence', 'protein_seq', or 'sequence' "
                "to build synonym masks"
            )
        return "".join(sequence.upper().split()).rstrip("*")


def collate_fn(batch: list[dict[str, torch.Tensor]]) -> dict[str, torch.Tensor]:
    """Pad variable-length ESM2 and codon sequences in a batch."""
    encoder_states = pad_sequence(
        [item["encoder_hidden_states"] for item in batch],
        batch_first=True,
    )
    decoder_ids = pad_sequence(
        [item["decoder_input_ids"] for item in batch],
        batch_first=True,
        padding_value=0,
    )
    labels = pad_sequence(
        [item["labels"] for item in batch],
        batch_first=True,
        padding_value=-100,
    )
    synonym_masks = pad_sequence(
        [item["synonym_mask"] for item in batch],
        batch_first=True,
        padding_value=False,
    )

    return {
        "encoder_hidden_states": encoder_states,
        "decoder_input_ids": decoder_ids,
        "labels": labels,
        "synonym_mask": synonym_masks,
    }
