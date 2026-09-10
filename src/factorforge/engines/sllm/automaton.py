from typing import Iterable, Set, Dict, Tuple
from collections import deque

class CompiledAutomaton:
    """
    Aho-Corasick 알고리즘이 적용된 상태 머신 (결정론적 유한 오토마타 - DFA).
    특정 상태(Node)에서 염기(A,T,C,G)나 코돈(3bp)을 입력받아 다음 상태로 빠르고 안전하게 전이합니다.
    """
    def __init__(self):
        self.transitions: Dict[int, Dict[str, int]] = {0: {}}
        self.fail_links: Dict[int, int] = {0: 0}
        self.terminal_states: Set[int] = set()
        self.state_counter = 0

    def add_state(self) -> int:
        self.state_counter += 1
        self.transitions[self.state_counter] = {}
        self.fail_links[self.state_counter] = 0
        return self.state_counter

    def step_nucleotide(self, state: int, char: str) -> int:
        """단일 염기에 대한 상태 전이 (실패 링크 자동 추적)"""
        curr = state
        while curr != 0 and char not in self.transitions[curr]:
            curr = self.fail_links[curr]
        return self.transitions[curr].get(char, 0)

    def step_codon(self, state: int, codon: str) -> Tuple[int, bool]:
        """
        코돈(3bp) 입력에 대한 상태 전이.
        반환값: (next_state, is_forbidden)
        """
        curr = state
        for char in codon:
            curr = self.step_nucleotide(curr, char)
            if curr in self.terminal_states:
                return curr, True  # 금지 모티프 도달 (Veto)
        return curr, False

class AutomatonCompiler:
    """
    제약 조건 모티프(BsaI 등) 리스트를 입력받아 CompiledAutomaton 인스턴스로 컴파일하는 빌더.
    역방향 서열(Reverse Complement)을 자동으로 생성하여 양방향 검사를 100% 보장합니다.
    """
    COMPLEMENT = str.maketrans('ATCG', 'TAGC')

    @classmethod
    def reverse_complement(cls, seq: str) -> str:
        return seq.translate(cls.COMPLEMENT)[::-1]

    @classmethod
    def compile(cls, motifs: Iterable[str], include_rc: bool = True) -> CompiledAutomaton:
        automaton = CompiledAutomaton()
        
        # 1. 정방향 및 역방향(RC) 모티프 수집
        all_motifs = set(motifs)
        if include_rc:
            all_motifs.update({cls.reverse_complement(m) for m in motifs})

        # 2. Trie(접두사 트리) 구축
        for motif in all_motifs:
            curr = 0
            for char in motif.upper():
                if char not in automaton.transitions[curr]:
                    new_state = automaton.add_state()
                    automaton.transitions[curr][char] = new_state
                curr = automaton.transitions[curr][char]
            automaton.terminal_states.add(curr)

        # 3. Aho-Corasick 실패 링크(Fail Links) 구축 (BFS)
        queue = deque()
        for char, child in automaton.transitions[0].items():
            automaton.fail_links[child] = 0
            queue.append(child)

        while queue:
            curr = queue.popleft()
            
            # 실패 링크가 가리키는 곳이 터미널(금지)이면, 현재 노드도 금지 노드로 취급 (부분 모티프 포함 방지)
            if automaton.fail_links[curr] in automaton.terminal_states:
                automaton.terminal_states.add(curr)

            for char, child in automaton.transitions[curr].items():
                fail_state = automaton.fail_links[curr]
                while fail_state != 0 and char not in automaton.transitions[fail_state]:
                    fail_state = automaton.fail_links[fail_state]
                
                automaton.fail_links[child] = automaton.transitions[fail_state].get(char, 0)
                queue.append(child)

        return automaton
