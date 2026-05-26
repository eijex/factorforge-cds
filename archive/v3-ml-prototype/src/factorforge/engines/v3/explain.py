"""FDA-style explainability report helpers for FactorForge v3."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from typing import Any

DEFAULT_RESTRICTION_SITES: dict[str, str] = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC",
    "NotI": "GCGGCCGC",
}


@dataclass(frozen=True)
class ExplainabilityInputs:
    aa_sequence: str
    dna_sequence: str
    metrics: dict[str, float]
    model_id: str
    tokenizer_hash: str
    seed: int
    config: dict[str, Any]
    post_guard: dict[str, Any] | None = None


def build_fda_report(inputs: ExplainabilityInputs) -> dict[str, Any]:
    aa_seq = inputs.aa_sequence
    dna_seq = inputs.dna_sequence

    return {
        "model": {
            "id": inputs.model_id,
            "hash": _hash_string(inputs.model_id),
        },
        "tokenizer": {
            "hash": inputs.tokenizer_hash,
        },
        "inputs": {
            "aa_length": len(aa_seq),
        },
        "outputs": {
            "dna_length": len(dna_seq),
            "cai": float(inputs.metrics.get("cai", 0.0)),
            "gc_percent": float(inputs.metrics.get("gc_content", 0.0)),
        },
        "constraint_checks": _constraint_checks(dna_seq),
        "post_guard": inputs.post_guard or {},
        "determinism": {
            "seed": inputs.seed,
            "config": inputs.config,
        },
    }


def write_fda_report(path: str, report: dict[str, Any]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)


def _constraint_checks(sequence: str) -> dict[str, Any]:
    return {
        "poly_a_runs": _scan_poly_a(sequence),
        "restriction_sites": _scan_restriction_sites(sequence, DEFAULT_RESTRICTION_SITES),
        "splice_like_motifs": _scan_splice_like(sequence),
        "homopolymers": _scan_homopolymers(sequence),
        "repeats": _scan_repeats(sequence),
    }


def _scan_poly_a(sequence: str, min_len: int = 6) -> list[dict[str, int]]:
    matches = re.finditer(rf"A{{{min_len},}}", sequence)
    return [{"start": m.start(), "end": m.end()} for m in matches]


def _scan_restriction_sites(
    sequence: str, sites: dict[str, str]
) -> list[dict[str, Any]]:
    hits: list[dict[str, Any]] = []
    for enzyme, motif in sites.items():
        positions = [m.start() for m in re.finditer(motif, sequence)]
        if positions:
            hits.append({"enzyme": enzyme, "motif": motif, "positions": positions})
    return hits


def _scan_splice_like(sequence: str) -> list[dict[str, int]]:
    matches = re.finditer(r"GT[ACGT]{2,20}AG", sequence)
    return [{"start": m.start(), "end": m.end()} for m in matches]


def _scan_homopolymers(sequence: str, min_len: int = 6) -> list[dict[str, Any]]:
    runs: list[dict[str, Any]] = []
    for base in "ACGT":
        for match in re.finditer(rf"{base}{{{min_len},}}", sequence):
            runs.append({"base": base, "start": match.start(), "end": match.end()})
    return runs


def _scan_repeats(sequence: str, kmer: int = 6) -> list[dict[str, Any]]:
    seen: dict[str, list[int]] = {}
    for idx in range(0, max(0, len(sequence) - kmer + 1)):
        token = sequence[idx : idx + kmer]
        seen.setdefault(token, []).append(idx)

    repeats: list[dict[str, Any]] = []
    for token, positions in seen.items():
        if len(positions) > 1:
            repeats.append({"kmer": token, "positions": positions})
    return repeats


def _hash_string(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()
