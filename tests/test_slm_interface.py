import math
import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from factorforge.engines.sllm.interfaces import (
    CodonLogits, 
    SynonymousMaskProcessor, 
    AutomatonConstraintProcessor,
    ConstraintProcessorPipeline,
    ConstraintState,
    CODON_VOCAB
)
from factorforge.engines.sllm.automaton import AutomatonCompiler

class MockSLMDecoder:
    def generate(self, protein_seq: str, mock_raw_logits: list, pipeline: ConstraintProcessorPipeline) -> list:
        selected_indices = []
        # 초기 상태(State 0) 설정
        state = ConstraintState(automaton_node=0, position=0)
        
        for pos, aa in enumerate(protein_seq):
            current_logits = CodonLogits(mock_raw_logits[pos])
            
            # 1. Pipeline을 통과하며 위험 경로 차단 (Predictive Masking)
            masked_logits = pipeline(pos, aa, state, current_logits)
            
            # 2. 가장 높은 확률의 코돈 선택 (Argmax)
            best_idx = int(np.argmax(masked_logits.scores))
            if masked_logits.scores[best_idx] == -np.inf:
                raise RuntimeError("BEAM_DEAD_END")
                
            selected_codon = CODON_VOCAB[best_idx]
            selected_indices.append(best_idx)
            
            # 3. 방금 찾아낸 Missing Link! 상태 업데이트 (State Transition)
            state = pipeline.get_next_state(pos, state, selected_codon)
            
        return selected_indices

def test_translation_invariant_masking():
    raw_scores = np.random.rand(64)
    logits = CodonLogits(raw_scores)
    processor = SynonymousMaskProcessor()
    state = ConstraintState(0, 0)
    
    masked_logits = processor(0, 'F', state, logits)
    
    # 'F'는 TTC, TTT 단 2개만 존재해야 함
    valid_count = np.sum(masked_logits.scores > -math.inf)
    assert valid_count == 2

def test_automaton_cross_codon_veto():
    """코돈 경계를 넘나드는 BsaI(GGTCTC) 모티프가 완벽히 차단(Veto)되는지 테스트"""
    # GGTCTC를 금지 모티프로 등록
    automaton = AutomatonCompiler.compile(["GGTCTC"])
    pipeline = ConstraintProcessorPipeline([
        SynonymousMaskProcessor(),
        AutomatonConstraintProcessor(automaton)
    ])
    decoder = MockSLMDecoder()
    
    # G(Glycine) -> L(Leucine) 순서로 생성
    # 1. G(Glycine) 위치의 Logits: GGT를 압도적 1위로 설정
    mock_logits_0 = np.zeros(64)
    mock_logits_0[CODON_VOCAB.index("GGT")] = 99.0 
    
    # 2. L(Leucine) 위치의 Logits: CTC를 1위로 설정하지만, 결합 시 BsaI(GGT+CTC)가 완성됨
    mock_logits_1 = np.zeros(64)
    mock_logits_1[CODON_VOCAB.index("CTC")] = 99.0 # Veto 대상
    mock_logits_1[CODON_VOCAB.index("CTT")] = 50.0 # 안전한 2위 코돈
    
    selected = decoder.generate("GL", [mock_logits_0, mock_logits_1], pipeline)
    
    # 결과 검증: AI가 99.0점을 준 CTC를 오토마톤이 차단하고 CTT를 선택했는가?
    assert CODON_VOCAB[selected[0]] == "GGT"
    assert CODON_VOCAB[selected[1]] == "CTT" 
