"""Math v2 DP Diagnostic Bridge for Hybrid Engine Failures."""

from factorforge.engines.sllm.interfaces import STANDARD_GENETIC_CODE
from factorforge.engines.sllm.automaton import AutomatonCompiler

def diagnose_global_feasibility(protein_seq: str, motifs_to_ban: list[str]) -> dict:
    """
    수학적 동적계획법(Exact Math DP)을 사용하여, 주어진 단백질 서열이
    금지 모티프를 완벽히 피하는 것이 '수학적으로 가능한지(Feasible)' 전수 검사합니다.
    """
    automaton = AutomatonCompiler.compile(motifs_to_ban, include_rc=True)
    
    # DP State: 도달 가능한 모든 오토마톤 상태(Node)의 집합
    active_states = {0}
    protein = "".join(protein_seq.upper().split()).rstrip("*")
    
    for pos, aa in enumerate(protein):
        codons = STANDARD_GENETIC_CODE.get(aa, [])
        next_states = set()
        
        for state in active_states:
            for codon in codons:
                next_node, is_forbidden = automaton.step_codon(state, codon)
                if not is_forbidden:
                    next_states.add(next_node)
                    
        if not next_states:
            # 어떤 경로를 택해도 금지 모티프를 피할 수 없음 (수학적 불가능)
            return {
                "global_feasible": False,
                "failure_position": pos + 1,
                "failure_aa": aa,
                "message": f"GLOBAL_INFEASIBLE: Mathematically impossible to avoid motifs at residue {pos+1} ({aa})."
            }
            
        active_states = next_states

    # 완벽히 통과 가능한 경로가 수학적으로 최소 1개 이상 존재함
    return {
        "global_feasible": True,
        "message": "FEASIBLE: Valid paths exist mathematically. Failure is due to SLM_SEARCH_LIMITATION (Beam width too small)."
    }
