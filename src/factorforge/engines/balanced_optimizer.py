import math
import hashlib
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
    return ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"]

class BalancedOptimizer:
    def __init__(
        self,
        forbidden_motifs: Optional[List[str]] = None,
        forbidden_enzymes: Optional[Set[str]] = None,
        include_reverse_complement: bool = True,
    ) -> None:
        motifs = forbidden_motifs if forbidden_motifs is not None else get_canonical_forbidden_motifs(forbidden_enzymes)
        self.automaton: CompiledAutomaton = AutomatonCompiler.compile(motifs, include_rc=include_reverse_complement)

    def optimize(
        self,
        protein_sequence: str,
        codon_weights: Dict[str, float],
        target_gc_min: float = 0.40,
        target_gc_max: float = 0.47,
        left_flank: str = "",
        right_flank: str = "",
        seed: str = "default_seed_307",
        lambda_cai: float = 0.70
    ) -> Dict[str, Any]:
        protein = "".join(protein_sequence.upper().split()).rstrip("*")
        if not protein:
            raise ValueError("protein_sequence must not be empty")

        for aa in protein:
            if aa not in AA_TO_CODONS:
                raise ValueError(f"Unsupported amino acid residue: {aa}")

        s0 = 0
        for char in left_flank.upper():
            s0 = self.automaton.step_nucleotide(s0, char)
            if s0 in self.automaton.terminal_states:
                raise ValueError("Left flank contains forbidden motif")

        n_codons = len(protein)
        total_nt = n_codons * 3

        true_probs = {}
        for aa, codons in AA_TO_CODONS.items():
            total_w = sum(codon_weights.get(c, 0.0) for c in codons)
            if total_w > 0:
                for c in codons:
                    true_probs[c] = codon_weights.get(c, 0.0) / total_w
            else:
                for c in codons:
                    true_probs[c] = 1.0 / len(codons)

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
                next_s, is_forbidden = self.automaton.step_codon(curr_s, codon)
                if is_forbidden:
                    continue
                
                target_prob = true_probs.get(codon, 1e-4)
                w = codon_weights.get(codon, 1e-4)
                cai_score = math.log(max(w, 1e-6))
                
                new_aa_total = aa_totals[aa] + 1
                curr_freq = (current_counts[codon] + 1) / new_aa_total
                dist_score = -abs(curr_freq - target_prob)
                
                concentration_penalty = 0.0
                if curr_freq > 0.8 and new_aa_total > 5:
                    concentration_penalty = (curr_freq - 0.8) * 2.0
                
                codon_gc = codon.count('G') + codon.count('C')
                expected_gc = curr_g + codon_gc + (n_codons - i - 1) * 1.5 
                gc_diff = abs(expected_gc - target_mid_gc_count) / total_nt
                
                hash_val = int(hashlib.sha256(f"{seed}:{i}:{codon}".encode('utf-8')).hexdigest()[:8], 16)
                noise = (hash_val / 0xFFFFFFFF) * 0.001
                
                total_score = (lambda_cai * cai_score) + ((1.0 - lambda_cai) * dist_score * 5.0) - concentration_penalty - (gc_diff * 2.0) + noise
                
                if total_score > best_score:
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
