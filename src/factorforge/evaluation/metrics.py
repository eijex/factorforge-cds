"""Metric calculators for Shared Evaluation.

Includes PolyA, Homopolymer, and Repeat detectors.
"""
import re
from typing import List, Optional

def calculate_polya_motifs(dna_sequence: str, threshold: int = 6) -> List[str]:
    """Find PolyA-like motifs (A runs) of length >= threshold."""
    seq = dna_sequence.upper()
    return re.findall(f"A{{{threshold},}}", seq)

def calculate_homopolymers(dna_sequence: str, threshold: int = 6) -> List[str]:
    """Find homopolymers (runs of any same nucleotide) of length >= threshold."""
    seq = dna_sequence.upper()
    matches = []
    for nt in "ATGC":
        matches.extend(re.findall(f"{nt}{{{threshold},}}", seq))
    return matches

def calculate_direct_repeats(dna_sequence: str, min_length: int = 10) -> List[str]:
    """Find simple direct repeats of min_length. (Naïve approach for evaluation)."""
    seq = dna_sequence.upper()
    repeats = set()
    for i in range(len(seq) - min_length * 2 + 1):
        kmer = seq[i:i+min_length]
        if seq.count(kmer) > 1:
            repeats.add(kmer)
    return list(repeats)
