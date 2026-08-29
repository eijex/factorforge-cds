"""Sequence risk scanners (PolyA, Repeats, Homopolymers, etc.)

Extracted from legacy RuleEngine to serve as shared evaluation metrics.
"""
import re
from typing import Dict, Any, List, Set, Tuple

# Polyadenylation signals (true biological motifs)
POLYA_PATTERNS = {
    # Tier 1 (canonical, highest efficiency in plants/eukaryotes)
    "AATAAA": "canonical",
    "ATTAAA": "variant_1",
    "AGTAAA": "variant_2",
    # Tier 2 (lower-frequency but functional in plants)
    "AATACA": "tier_2",
    "AAGAAA": "tier_2",
    "AATGAA": "tier_2",
}

UNSTABLE_MOTIFS = {
    "ATTTA": "ARE (AU-rich element)", 
    "WWWWWW": "W=A/T, 6+ in a row"
}

def scan_polya(seq: str) -> List[Dict[str, Any]]:
    """Detect biological PolyA signal motifs (e.g. AATAAA)."""
    seq_upper = seq.upper()
    violations: List[Dict[str, Any]] = []
    
    for pattern, pattern_type in POLYA_PATTERNS.items():
        pos = 0
        while True:
            idx = seq_upper.find(pattern, pos)
            if idx == -1:
                break
            violations.append({
                "type": "polya_signal",
                "pattern": pattern,
                "pattern_type": pattern_type,
                "position": idx,
                "context": seq_upper[max(0, idx - 10):min(len(seq_upper), idx + len(pattern) + 10)]
            })
            pos = idx + 1
            
    return violations

def scan_homopolymers(seq: str, threshold: int = 6) -> List[Dict[str, Any]]:
    """Find homopolymers (runs of any single nucleotide)."""
    seq_upper = seq.upper()
    violations: List[Dict[str, Any]] = []
    for nt in "ATGC":
        pattern = nt * threshold
        pos = 0
        while True:
            idx = seq_upper.find(pattern, pos)
            if idx == -1:
                break
            match = re.match(f"{nt}+", seq_upper[idx:])
            length = match.end() if match else threshold
            violations.append({
                "type": "homopolymer",
                "base": nt,
                "length": length,
                "position": idx,
                "context": seq_upper[max(0, idx - 10):min(len(seq_upper), idx + length + 10)]
            })
            pos = idx + length
    return violations

def scan_repeats(seq: str, min_length: int = 20, max_distance: int = 1000) -> List[Dict[str, Any]]:
    """Detect direct or inverted repeats causing assembly/stability issues.
    
    Direct repeats (exact k-mer matches) >= 20 bp are known to cause issues with 
    DNA synthesis (e.g. Twist/IDT) and homologous recombination in host organisms.
    """
    seq_upper = seq.upper()
    violations: List[Dict[str, Any]] = []
    
    # Naive k-mer search for direct repeats
    # For a real implementation, a suffix tree or kmer hash map is better.
    # We will use a fast localized search.
    seen_kmers = {}
    for i in range(len(seq_upper) - min_length + 1):
        kmer = seq_upper[i:i + min_length]
        if kmer in seen_kmers:
            dist = i - seen_kmers[kmer]
            if dist <= max_distance:
                violations.append({
                    "type": "direct_repeat",
                    "pattern": kmer,
                    "position1": seen_kmers[kmer],
                    "position2": i,
                    "distance": dist
                })
        seen_kmers[kmer] = i
        
    return violations
