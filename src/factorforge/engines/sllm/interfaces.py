import math
import numpy as np
from dataclasses import dataclass

# 방금 만든 Aho-Corasick 오토마톤 컴파일러 임포트
from factorforge.engines.sllm.automaton import CompiledAutomaton

# =====================================================================
# Canonical Codon Vocabulary Contract
# =====================================================================
CODON_VOCAB = (
    'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT',
    'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
    'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT',
    'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
    'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
    'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
    'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
    'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT'
)
CODON_TO_INDEX = {codon: idx for idx, codon in enumerate(CODON_VOCAB)}

STANDARD_GENETIC_CODE = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTT', 'TTC'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'], 'K': ['AAA', 'AAG'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M': ['ATG'],
    'N': ['AAT', 'AAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'W': ['TGG'], 'Y': ['TAT', 'TAC'], '*': ['TAA', 'TAG', 'TGA']
}

AA_TO_INDICES = {
    aa: [CODON_TO_INDEX[c] for c in codons]
    for aa, codons in STANDARD_GENETIC_CODE.items()
}

# =====================================================================
# State & Logits Contracts
# =====================================================================

@dataclass(frozen=True)
class ConstraintState:
    """Immutable state for safe Beam Search tracking."""
    automaton_node: int
    position: int

class CodonLogits:
    """Adapter Contract enforcing exactly 64 elements and true immutability."""
    def __init__(self, raw_scores: np.ndarray):
        if raw_scores.shape != (64,):
            raise ValueError("CodonLogits must strictly have 64 elements.")
        self.scores = np.array(raw_scores, dtype=np.float32, copy=True)
        self.scores.setflags(write=False)

# =====================================================================
# Processors
# =====================================================================

class BaseConstraintProcessor:
    def __call__(self, position: int, target_aa: str, state: ConstraintState, logits: CodonLogits) -> CodonLogits:
        raise NotImplementedError
        
    def get_next_state(self, position: int, state: ConstraintState, chosen_codon: str) -> ConstraintState:
        """AI가 최종 선택한 코돈을 기반으로 다음 상태(State)를 반환합니다."""
        return state

class SynonymousMaskProcessor(BaseConstraintProcessor):
    def __call__(self, position: int, target_aa: str, state: ConstraintState, logits: CodonLogits) -> CodonLogits:
        valid_indices = AA_TO_INDICES.get(target_aa, [])
        if not valid_indices:
            raise ValueError(f"Invalid amino acid '{target_aa}' at pos {position}.")
            
        masked = np.copy(logits.scores)
        for i in range(64):
            if i not in valid_indices:
                masked[i] = -math.inf
        return CodonLogits(masked)

class AutomatonConstraintProcessor(BaseConstraintProcessor):
    def __init__(self, automaton: CompiledAutomaton):
        self.automaton = automaton

    def __call__(self, position: int, target_aa: str, state: ConstraintState, logits: CodonLogits) -> CodonLogits:
        masked = np.copy(logits.scores)
        
        for i in range(64):
            if masked[i] == -math.inf: 
                continue
            
            # Aho-Corasick 상태 머신을 통과시켜 Veto 판정 시뮬레이션
            candidate_codon = CODON_VOCAB[i]
            _, is_forbidden = self.automaton.step_codon(state.automaton_node, candidate_codon)
            
            if is_forbidden:
                masked[i] = -math.inf
            
        return CodonLogits(masked)
        
    def get_next_state(self, position: int, state: ConstraintState, chosen_codon: str) -> ConstraintState:
        """선택된 코돈을 입력받아 진짜로 상태(Node)를 다음 단계로 이동시킵니다."""
        next_node, is_forbidden = self.automaton.step_codon(state.automaton_node, chosen_codon)
        if is_forbidden:
            raise ValueError(f"Fatal: Chosen codon '{chosen_codon}' results in a forbidden state!")
            
        return ConstraintState(
            automaton_node=next_node,
            position=position + 1
        )

class ConstraintProcessorPipeline:
    """Enforces explicit ordering of constraint processors and global state tracking."""
    def __init__(self, processors: list[BaseConstraintProcessor]):
        self.processors = processors

    def __call__(self, position: int, target_aa: str, state: ConstraintState, logits: CodonLogits) -> CodonLogits:
        current_logits = logits
        for processor in self.processors:
            current_logits = processor(position, target_aa, state, current_logits)
        
        if np.all(current_logits.scores == -math.inf):
            raise RuntimeError("BEAM_DEAD_END")
            
        return current_logits
        
    def get_next_state(self, position: int, state: ConstraintState, chosen_codon: str) -> ConstraintState:
        """모든 프로세서의 상태를 순차적으로 전이시킵니다."""
        current_state = state
        for processor in self.processors:
            current_state = processor.get_next_state(position, current_state, chosen_codon)
        return current_state
