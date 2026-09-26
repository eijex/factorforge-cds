
import itertools
from typing import Dict, Any, List, Optional
from factorforge.engines.balanced_optimizer import STANDARD_GENETIC_CODE, AA_TO_CODONS
from factorforge.analysis.metrics import translate_dna, _codons

class TargetedRemediator:
    """
    Layer 1 Engine component for minimal-edit synonymous targeted remediation.
    It takes a DNA sequence and a target window, and explores the local synonymous
    codon space to find a replacement that satisfies a validation callback
    while strictly minimizing the Hamming distance to the original sequence.
    """
    
    def __init__(self, include_reverse_complement: bool = False):
        self.include_reverse_complement = include_reverse_complement
        
    def remediate_window(
        self,
        dna_sequence: str,
        target_start: int,
        target_end: int,
        validation_callback,
        window_expansion: int = 6  # Expand search space by N nt on each side
    ) -> Dict[str, Any]:
        """
        Attempt to remediate a target region by synonymous edits.
        
        Args:
            dna_sequence: The full original DNA sequence.
            target_start: 1-based start coordinate of the finding to break.
            target_end: 1-based end coordinate of the finding to break.
            validation_callback: A function f(mutated_seq: str) -> bool that returns True if the finding is resolved.
            window_expansion: How many nucleotides to expand the mutational search space to catch the full codon boundaries.
        """
        
        # 0-based coordinates
        s = max(0, target_start - 1 - window_expansion)
        e = min(len(dna_sequence), target_end + window_expansion)
        
        # Snap to codon boundaries
        codon_start_idx = s // 3
        codon_end_idx = (e + 2) // 3
        
        nt_start = codon_start_idx * 3
        nt_end = codon_end_idx * 3
        
        prefix = dna_sequence[:nt_start]
        suffix = dna_sequence[nt_end:]
        target_window = dna_sequence[nt_start:nt_end]
        
        target_aa = translate_dna(target_window)
        
        # Generate synonymous codon combinations for the window
        possible_codons_per_aa = [AA_TO_CODONS.get(aa, [target_window[i*3:i*3+3]]) for i, aa in enumerate(target_aa)]
        
        # To avoid combinatorial explosion if the window is huge, we cap the search
        # or sort combinations by hamming distance from the original window early.
        
        best_seq = dna_sequence
        min_hamming = float('inf')
        found = False
        
        # For small windows (e.g. 5-7 AAs), itertools.product is fine.
        # 6 AAs with avg 3 codons each = 729 combinations.
        for codon_combo in itertools.product(*possible_codons_per_aa):
            candidate_window = "".join(codon_combo)
            
            # Calculate local hamming distance
            hamming = sum(1 for a, b in zip(target_window, candidate_window) if a != b)
            
            if hamming == 0:
                continue # Original sequence
                
            if hamming >= min_hamming:
                continue # Optimization: skip if we already found a closer one
                
            candidate_seq = prefix + candidate_window + suffix
            
            # Check if this resolves the finding (and doesn't introduce bad things if callback checks it)
            if validation_callback(candidate_seq):
                min_hamming = hamming
                best_seq = candidate_seq
                found = True
                
        return {
            "resolved": found,
            "original_sequence": dna_sequence,
            "remediated_sequence": best_seq,
            "hamming_distance": min_hamming if found else 0,
            "window_nt_start": nt_start,
            "window_nt_end": nt_end
        }
