"""Metrics for FactorForge v3 (CAI, perplexity, GC%)."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - for type checkers only
    import torch


def _default_codon_table_path() -> Path:
    return Path(__file__).resolve().parents[4] / "data" / "nbenthamiana_codons.json"


@dataclass(frozen=True)
class CodonUsageTable:
    codon_to_aa: dict[str, str]
    codon_weights: dict[str, float]
    best_codon_for_aa: dict[str, str]
    source: str | None = None


def load_codon_usage_table(path: Path | None = None) -> CodonUsageTable:
    table_path = path or _default_codon_table_path()
    raw = json.loads(table_path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("Codon table JSON must be an object")

    codons = raw.get("codons")
    if not isinstance(codons, dict):
        raise ValueError("Codon table JSON missing 'codons'")

    codon_to_aa: dict[str, str] = {}
    codon_freq: dict[str, float] = {}
    for codon, entry in codons.items():
        if not isinstance(entry, dict) or not isinstance(codon, str):
            continue
        aa = entry.get("aa")
        freq = entry.get("frequency")
        if not isinstance(aa, str) or not isinstance(freq, (int, float)):
            continue
        codon_to_aa[codon] = aa
        codon_freq[codon] = float(freq)

    codon_weights = _build_codon_weights(codon_to_aa, codon_freq)
    best_codon_for_aa = _best_codon_map(codon_to_aa, codon_freq)
    source = raw.get("source") if isinstance(raw.get("source"), str) else None

    return CodonUsageTable(
        codon_to_aa=codon_to_aa,
        codon_weights=codon_weights,
        best_codon_for_aa=best_codon_for_aa,
        source=source,
    )


def compute_cai(dna_sequence: str, table: CodonUsageTable) -> float:
    seq = dna_sequence.upper().replace("U", "T")
    codon_count = len(seq) // 3
    if codon_count == 0:
        return 0.0

    weights: list[float] = []
    for i in range(codon_count):
        codon = seq[i * 3 : i * 3 + 3]
        weight = table.codon_weights.get(codon)
        if weight is None or weight <= 0:
            return 0.0
        weights.append(weight)

    log_sum = sum(math.log(w) for w in weights)
    return math.exp(log_sum / len(weights))


def compute_gc(dna_sequence: str) -> float:
    seq = dna_sequence.upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return (gc / len(seq)) * 100.0


def compute_perplexity(
    logits: "torch.Tensor",
    labels: "torch.Tensor",
    ignore_index: int = -100,
) -> float:
    torch, functional = _require_torch()
    vocab_size = int(logits.shape[-1])
    loss = functional.cross_entropy(
        logits.view(-1, vocab_size),
        labels.view(-1),
        ignore_index=ignore_index,
    )
    return float(torch.exp(loss).item())


def _build_codon_weights(
    codon_to_aa: dict[str, str],
    codon_freq: dict[str, float],
) -> dict[str, float]:
    by_aa: dict[str, list[float]] = {}
    for codon, aa in codon_to_aa.items():
        if aa == "*":
            continue
        by_aa.setdefault(aa, []).append(codon_freq.get(codon, 0.0))

    weights: dict[str, float] = {}
    for codon, aa in codon_to_aa.items():
        if aa == "*":
            continue
        max_freq = max(by_aa.get(aa, [0.0]))
        freq = codon_freq.get(codon, 0.0)
        weights[codon] = freq / max_freq if max_freq > 0 else 0.0
    return weights


def _best_codon_map(
    codon_to_aa: dict[str, str],
    codon_freq: dict[str, float],
) -> dict[str, str]:
    best: dict[str, tuple[str, float]] = {}
    for codon, aa in codon_to_aa.items():
        if aa == "*":
            continue
        current = best.get(aa)
        freq = codon_freq.get(codon, 0.0)
        if current is None or freq > current[1]:
            best[aa] = (codon, freq)
    return {aa: codon for aa, (codon, _) in best.items()}


def _require_torch():
    try:
        import torch
        from torch.nn import functional
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "ML dependencies not installed. Install with: pip install -e \".[ml]\""
        ) from exc
    return torch, functional
