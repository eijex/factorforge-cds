"""
Tokenizer utilities for FactorForge v3.

Provides deterministic AA and codon tokenizers with stable token-id mappings.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

SPECIAL_TOKENS: tuple[str, ...] = ("[PAD]", "[BOS]", "[EOS]", "[MASK]", "[UNK]")
AA_TOKENS: tuple[str, ...] = (
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
CODON_BASE_ORDER: tuple[str, ...] = ("T", "C", "A", "G")


def _build_codons() -> tuple[str, ...]:
    return tuple(
        f"{first}{second}{third}"
        for first in CODON_BASE_ORDER
        for second in CODON_BASE_ORDER
        for third in CODON_BASE_ORDER
    )


CODON_TOKENS: tuple[str, ...] = _build_codons()


@dataclass(frozen=True)
class TokenizerMetadata:
    kind: str
    token_to_id: dict[str, int]
    special_tokens: tuple[str, ...]


class BaseTokenizer:
    def __init__(self, metadata: TokenizerMetadata) -> None:
        self.kind = metadata.kind
        self.token_to_id = dict(metadata.token_to_id)
        self.special_tokens = metadata.special_tokens
        self.id_to_token = {idx: token for token, idx in self.token_to_id.items()}

        self.pad_token_id = self._require_token_id("[PAD]")
        self.bos_token_id = self._require_token_id("[BOS]")
        self.eos_token_id = self._require_token_id("[EOS]")
        self.mask_token_id = self._require_token_id("[MASK]")
        self.unk_token_id = self._require_token_id("[UNK]")

    def _require_token_id(self, token: str) -> int:
        if token not in self.token_to_id:
            raise ValueError(f"Tokenizer missing required token: {token}")
        return self.token_to_id[token]

    def encode(self, sequence: str, add_special_tokens: bool = True) -> list[int]:
        tokens = self._tokenize(sequence)
        ids = [self.token_to_id.get(token, self.unk_token_id) for token in tokens]
        if add_special_tokens:
            return [self.bos_token_id, *ids, self.eos_token_id]
        return ids

    def decode(self, ids: Sequence[int], skip_special_tokens: bool = True) -> str:
        tokens: list[str] = []
        for idx in ids:
            token = self.id_to_token.get(int(idx))
            if token is None:
                continue
            if skip_special_tokens and token in self.special_tokens:
                continue
            tokens.append(token)
        return "".join(tokens)

    def mapping_hash(self) -> str:
        payload = json.dumps(self.token_to_id, sort_keys=True).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "kind": self.kind,
            "special_tokens": list(self.special_tokens),
            "token_to_id": self.token_to_id,
        }
        path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    def _tokenize(self, sequence: str) -> list[str]:
        raise NotImplementedError


class AATokenizer(BaseTokenizer):
    @classmethod
    def default(cls) -> "AATokenizer":
        token_to_id = _build_token_map((*SPECIAL_TOKENS, *AA_TOKENS))
        return cls(TokenizerMetadata(kind="aa", token_to_id=token_to_id, special_tokens=SPECIAL_TOKENS))

    @classmethod
    def load(cls, path: Path) -> "AATokenizer":
        metadata = _load_metadata(path)
        if metadata.kind != "aa":
            raise ValueError(f"Tokenizer kind mismatch: expected 'aa', got {metadata.kind!r}")
        return cls(metadata)

    def _tokenize(self, sequence: str) -> list[str]:
        return [token for token in sequence.strip().upper() if token]


class CodonTokenizer(BaseTokenizer):
    @classmethod
    def default(cls) -> "CodonTokenizer":
        token_to_id = _build_token_map((*SPECIAL_TOKENS, *CODON_TOKENS))
        return cls(
            TokenizerMetadata(kind="codon", token_to_id=token_to_id, special_tokens=SPECIAL_TOKENS)
        )

    @classmethod
    def load(cls, path: Path) -> "CodonTokenizer":
        metadata = _load_metadata(path)
        if metadata.kind != "codon":
            raise ValueError(f"Tokenizer kind mismatch: expected 'codon', got {metadata.kind!r}")
        return cls(metadata)

    def encode(
        self, sequence: str, add_special_tokens: bool = True, strict: bool = True
    ) -> list[int]:
        seq = sequence.strip().upper().replace("U", "T")
        if strict and len(seq) % 3 != 0:
            raise ValueError("Codon sequence length must be a multiple of 3")
        tokens = [seq[i : i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]
        ids = [self.token_to_id.get(token, self.unk_token_id) for token in tokens]
        if add_special_tokens:
            return [self.bos_token_id, *ids, self.eos_token_id]
        return ids

    def _tokenize(self, sequence: str) -> list[str]:
        return [sequence]


def _build_token_map(tokens: Iterable[str]) -> dict[str, int]:
    token_to_id: dict[str, int] = {}
    for idx, token in enumerate(tokens):
        token_to_id[token] = idx
    return token_to_id


def _load_metadata(path: Path) -> TokenizerMetadata:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("Tokenizer JSON must be an object")

    kind = raw.get("kind")
    if not isinstance(kind, str):
        raise ValueError("Tokenizer JSON missing 'kind'")

    raw_special = raw.get("special_tokens")
    if not isinstance(raw_special, list) or not all(isinstance(token, str) for token in raw_special):
        raise ValueError("Tokenizer JSON missing valid 'special_tokens'")
    special_tokens = tuple(raw_special)

    raw_map = raw.get("token_to_id")
    if not isinstance(raw_map, dict):
        raise ValueError("Tokenizer JSON missing 'token_to_id'")

    token_to_id: dict[str, int] = {}
    for token, idx in raw_map.items():
        if not isinstance(token, str) or not isinstance(idx, int):
            raise ValueError("Tokenizer JSON has invalid token mapping")
        token_to_id[token] = idx

    return TokenizerMetadata(kind=kind, token_to_id=token_to_id, special_tokens=special_tokens)
