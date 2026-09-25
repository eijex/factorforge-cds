import math
from typing import Dict, Any, Optional, Set, List
from collections import defaultdict
from factorforge.engines.sllm.automaton import AutomatonCompiler, CompiledAutomaton

STANDARD_GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

AA_TO_CODONS: Dict[str, List[str]] = {}
for codon, aa in STANDARD_GENETIC_CODE.items():
    if aa != "*":
        AA_TO_CODONS.setdefault(aa, []).append(codon)

def get_canonical_forbidden_motifs(forbidden_enzymes=None):
    # Fallback mock for canonical motifs
    return ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"]

class BalancedOptimizer:
    """Deterministic Distribution-Aware Optimizer (Greedy with GC and Automaton Support)."""

    def __init__(
        self,
        forbidden_motifs: Optional[List[str]] = None,
        forbidden_enzymes: Optional[Set[str]] = None,
        include_reverse_complement: bool = True,
    ) -> None:
        if forbidden_motifs is not None:
            motifs = forbidden_motifs
        else:
            motifs = get_canonical_forbidden_motifs(forbidden_enzymes)

        self.automaton: CompiledAutomaton = AutomatonCompiler.compile(
            motifs, include_rc=include_reverse_complement
        )

    def optimize(
        self,
        protein_sequence: str,
        codon_weights: Dict[str, float],
        target_gc_min: float = 0.40,
        target_gc_max: float = 0.47,
        left_flank: str = "",
        right_flank: str = "",
    ) -> Dict[str, Any]:
        """Runs a deterministic greedy search optimizing distribution divergence and GC constraints."""
        protein = "".join(protein_sequence.upper().split()).rstrip("*")
        if not protein:
            raise ValueError("protein_sequence must not be empty")

        for aa in protein:
            if aa not in AA_TO_CODONS:
                raise ValueError(f"Unsupported amino acid residue: {aa}")

        # Basic setup
        s0 = 0
        for char in left_flank.upper():
            s0 = self.automaton.step_nucleotide(s0, char)
            if s0 in self.automaton.terminal_states:
                raise ValueError("Left flank contains forbidden motif")

        n_codons = len(protein)
        total_nt = n_codons * 3

        # State initialization
        current_counts = defaultdict(int)
        aa_totals = defaultdict(int)
        
        target_mid_gc_count = (target_gc_min + target_gc_max) / 2.0 * total_nt
        
        seq_codons = []
        curr_s = s0
        curr_g = 0
        score_sum = 0.0

        for i, aa in enumerate(protein):
            synonymous_codons = AA_TO_CODONS[aa]
            best_codon = None
            best_score = -float('inf')
            best_next_s = None
            
            for codon in synonymous_codons:
                # 1. Automaton check
                next_s, is_forbidden = self.automaton.step_codon(curr_s, codon)
                if is_forbidden:
                    continue
                
                # 2. Distribution score (Minimize divergence)
                # Score = TargetFreq - CurrentFreq
                target_freq = codon_weights.get(codon, 1e-4)
                new_aa_total = aa_totals[aa] + 1
                curr_freq = current_counts[codon] / new_aa_total
                dist_score = target_freq - curr_freq
                
                # 3. GC Score (drive toward mid GC)
                codon_gc = codon.count('G') + codon.count('C')
                next_g = curr_g + codon_gc
                expected_gc_after_this = next_g + (n_codons - i - 1) * 1.5 # estimate remaining
                gc_diff = abs(expected_gc_after_this - target_mid_gc_count)
                
                # 4. Total heuristic score (Scale dist_score * 100 to dominate GC slightly)
                total_score = dist_score * 100 - gc_diff
                
                if total_score > best_score:
                    best_score = total_score
                    best_codon = codon
                    best_next_s = next_s
                elif abs(total_score - best_score) < 1e-9:
                    if best_codon is not None and codon < best_codon:
                        best_score = total_score
                        best_codon = codon
                        best_next_s = next_s

            if best_codon is None:
                raise RuntimeError(f"Dead end at position {i} for AA {aa}")
                
            seq_codons.append(best_codon)
            current_counts[best_codon] += 1
            aa_totals[aa] += 1
            curr_s = best_next_s
            curr_g += best_codon.count('G') + best_codon.count('C')
            
            w = codon_weights.get(best_codon, 1e-4)
            score_sum += math.log(max(w, 1e-6))

        optimized_dna = "".join(seq_codons)
        final_gc_percent = curr_g / total_nt * 100.0
        final_cai = math.exp(score_sum / n_codons)

        return {
            "sequence": optimized_dna,
            "cai": final_cai,
            "gc_percent": final_gc_percent,
            "gc_feasible": (target_gc_min <= (final_gc_percent/100.0) <= target_gc_max),
            "constraint_scope": "CDS_ONLY",
            "score": score_sum,
            "length_nt": total_nt,
            "length_aa": n_codons,
        }
