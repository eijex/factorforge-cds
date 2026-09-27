# factorforge/src/factorforge/scoring/quantization.py
"""Canonical Fixed-Point Quantization and Deterministic Tie-Breaking for FactorForge."""

from __future__ import annotations

from decimal import Decimal, ROUND_HALF_EVEN
import hashlib
from typing import List, Tuple, Union


class CanonicalQuantizer:
    """Fixed-point integer quantizer for DP additive objective terms.
    
    Guarantees cross-platform, byte-identical integer arithmetic:
        Q_i(c) = RoundCanonical( l_i(c) * Scale )
        J_int(X) = sum_i Q_i(c_i)
    """

    def __init__(self, scale: int = 1_000_000, rounding: str = "ROUND_HALF_EVEN"):
        self.scale = scale
        self.rounding = rounding
        self._scale_dec = Decimal(str(scale))

    def quantize(self, value: Union[float, Decimal, str]) -> int:
        """Convert a floating-point score term into a canonical fixed-point integer.
        
        Uses Decimal with ROUND_HALF_EVEN to eliminate platform-dependent floating-point rounding biases.
        """
        if isinstance(value, float):
            # Convert float through canonical string format to avoid float representation artifacts
            dec = Decimal(f"{value:.8f}")
        elif isinstance(value, str):
            dec = Decimal(value)
        else:
            dec = value

        scaled = dec * self._scale_dec
        rounded = scaled.quantize(Decimal("1"), rounding=ROUND_HALF_EVEN)
        return int(rounded)

    def quantize_additive_vector(self, values: List[float]) -> List[int]:
        """Quantize an array of additive score terms."""
        return [self.quantize(v) for v in values]

    def compute_digest(self) -> str:
        """Fingerprint the quantizer configuration."""
        raw = f"scale:{self.scale};rounding:{self.rounding}".encode("utf-8")
        return f"sha256:{hashlib.sha256(raw).hexdigest()}"


class PathRank:
    """Full-prefix deterministic lexicographic rank calculator for DP tie-breaking."""

    CODON_ALPHABET = sorted([
        "TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",
        "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG",
        "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG",
        "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG",
        "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG",
        "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG",
        "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
        "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"
    ])
    CODON_PRIORITY = {codon: idx for idx, codon in enumerate(CODON_ALPHABET)}

    @classmethod
    def compare_prefix_paths(cls, path_a: List[str], path_b: List[str]) -> int:
        """Compare two full prefix paths lexicographically.
        
        Returns:
            -1 if path_a < path_b (path_a has higher priority / lower index)
             1 if path_a > path_b
             0 if path_a == path_b
        """
        for c_a, c_b in zip(path_a, path_b):
            rank_a = cls.CODON_PRIORITY.get(c_a, 999)
            rank_b = cls.CODON_PRIORITY.get(c_b, 999)
            if rank_a < rank_b:
                return -1
            if rank_a > rank_b:
                return 1
        if len(path_a) < len(path_b):
            return -1
        if len(path_a) > len(path_b):
            return 1
        return 0

    @classmethod
    def best_path(cls, candidate_paths: List[List[str]]) -> List[str]:
        """Select the uniquely best path among equal-score candidates using full-prefix ordering."""
        if not candidate_paths:
            raise ValueError("Candidate paths list is empty")
        
        best = candidate_paths[0]
        for candidate in candidate_paths[1:]:
            if cls.compare_prefix_paths(candidate, best) < 0:
                best = candidate
        return best
