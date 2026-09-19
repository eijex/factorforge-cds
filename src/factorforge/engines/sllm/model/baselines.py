# factorforge/src/factorforge/engines/sllm/model/baselines.py
"""Statistical baseline models for codon sequence prediction (Job 285B).

Provides:
- Baseline 0: Host Marginal Codon Frequency
- Baseline 1: Host Amino-Acid Conditional Codon Frequency
- Baseline 2: Codon Bigram / Trigram Markov Model with Laplace smoothing
"""

from __future__ import annotations
import math
from typing import Dict, List, Optional, Tuple
import numpy as np

from factorforge.engines.sllm.interfaces import AA_TO_INDICES, CODON_TO_INDEX, STANDARD_GENETIC_CODE


class CodonFrequencyBaseline:
    """Baseline 0 & 1: Evaluates marginal and AA-conditional codon frequencies."""

    def __init__(self) -> None:
        self.marginal_counts = np.zeros(64, dtype=np.float64)
        self.aa_counts = {aa: np.zeros(64, dtype=np.float64) for aa in STANDARD_GENETIC_CODE.keys()}
        self.total_codons = 0

    def fit(self, train_records: List[Dict[str, str]]) -> None:
        """Counts codon occurrences from training records."""
        for rec in train_records:
            cds = rec["coding_cds"]
            aa_seq = rec["protein_aa"]
            for i in range(0, len(cds), 3):
                codon = cds[i : i + 3]
                pos = i // 3
                if pos >= len(aa_seq):
                    break
                aa = aa_seq[pos]
                if codon in CODON_TO_INDEX and aa in self.aa_counts:
                    idx = CODON_TO_INDEX[codon]
                    self.marginal_counts[idx] += 1.0
                    self.aa_counts[aa][idx] += 1.0
                    self.total_codons += 1

    def evaluate_nll(
        self, val_records: List[Dict[str, str]], mode: str = "conditional"
    ) -> Tuple[float, float]:
        """Calculates synonymous-masked NLL and perplexity on held-out records."""
        total_nll = 0.0
        total_tokens = 0

        for rec in val_records:
            cds = rec["coding_cds"]
            aa_seq = rec["protein_aa"]
            for i in range(0, len(cds), 3):
                codon = cds[i : i + 3]
                pos = i // 3
                if pos >= len(aa_seq) or codon not in CODON_TO_INDEX:
                    continue
                aa = aa_seq[pos]
                valid_indices = AA_TO_INDICES.get(aa, [])
                if not valid_indices:
                    continue

                actual_idx = CODON_TO_INDEX[codon]

                if mode == "marginal":
                    freqs = self.marginal_counts[valid_indices] + 1e-6
                else:  # conditional
                    freqs = self.aa_counts[aa][valid_indices] + 1e-6

                probs = freqs / np.sum(freqs)
                valid_idx_map = {idx: p for idx, p in zip(valid_indices, probs)}
                prob = valid_idx_map.get(actual_idx, 1e-6)
                total_nll += -math.log(max(1e-9, prob))
                total_tokens += 1

        mean_nll = total_nll / max(1, total_tokens)
        perplexity = math.exp(mean_nll)
        return float(mean_nll), float(perplexity)


class CodonBigramMarkovBaseline:
    """Baseline 2: Codon Bigram Markov Model: P(c_t | c_{t-1}, a_t)."""

    def __init__(self, smoothing: float = 0.1) -> None:
        self.smoothing = smoothing
        # Shape: (64, 64) -> transition counts from previous_codon to current_codon
        self.transition_counts = np.zeros((64, 64), dtype=np.float64)

    def fit(self, train_records: List[Dict[str, str]]) -> None:
        """Fits bigram transition matrix."""
        for rec in train_records:
            cds = rec["coding_cds"]
            prev_idx: Optional[int] = None
            for i in range(0, len(cds), 3):
                codon = cds[i : i + 3]
                if codon in CODON_TO_INDEX:
                    curr_idx = CODON_TO_INDEX[codon]
                    if prev_idx is not None:
                        self.transition_counts[prev_idx, curr_idx] += 1.0
                    prev_idx = curr_idx

    def evaluate_nll(self, val_records: List[Dict[str, str]]) -> Tuple[float, float]:
        """Calculates synonymous-masked NLL and perplexity on held-out records."""
        total_nll = 0.0
        total_tokens = 0

        for rec in val_records:
            cds = rec["coding_cds"]
            aa_seq = rec["protein_aa"]
            prev_idx: Optional[int] = None

            for i in range(0, len(cds), 3):
                codon = cds[i : i + 3]
                pos = i // 3
                if pos >= len(aa_seq) or codon not in CODON_TO_INDEX:
                    continue
                aa = aa_seq[pos]
                valid_indices = AA_TO_INDICES.get(aa, [])
                if not valid_indices:
                    continue

                actual_idx = CODON_TO_INDEX[codon]

                if prev_idx is not None:
                    counts = self.transition_counts[prev_idx, valid_indices] + self.smoothing
                else:
                    counts = np.ones(len(valid_indices), dtype=np.float64)

                probs = counts / np.sum(counts)
                valid_idx_map = {idx: p for idx, p in zip(valid_indices, probs)}
                prob = valid_idx_map.get(actual_idx, 1e-6)
                total_nll += -math.log(max(1e-9, prob))
                total_tokens += 1
                prev_idx = actual_idx

        mean_nll = total_nll / max(1, total_tokens)
        perplexity = math.exp(mean_nll)
        return float(mean_nll), float(perplexity)
