# factorforge/src/factorforge/engines/sllm/interfaces.py
"""Interfaces and Processors for Neuro-Symbolic sLLM Constraint Pipeline (Job 284A / 284B).

Enforces:
1. Exact amino acid translation invariant (SynonymousMaskProcessor).
2. Cross-codon restriction site dynamic veto via Aho-Corasick automaton (AutomatonConstraintProcessor).
3. Homopolymer run restriction (HomopolymerConstraintProcessor).
4. 5' Initiation Ramp GC Envelope (InitiationGCConstraintProcessor).
5. Canonical immutable state and logit adapters.
"""

from __future__ import annotations
import math
from dataclasses import dataclass
from typing import List
import numpy as np

from factorforge.engines.sllm.automaton import CompiledAutomaton

# =====================================================================
# Canonical Codon Vocabulary Contract
# =====================================================================
CODON_VOCAB = (
    "AAA",
    "AAC",
    "AAG",
    "AAT",
    "ACA",
    "ACC",
    "ACG",
    "ACT",
    "AGA",
    "AGC",
    "AGG",
    "AGT",
    "ATA",
    "ATC",
    "ATG",
    "ATT",
    "CAA",
    "CAC",
    "CAG",
    "CAT",
    "CCA",
    "CCC",
    "CCG",
    "CCT",
    "CGA",
    "CGC",
    "CGG",
    "CGT",
    "CTA",
    "CTC",
    "CTG",
    "CTT",
    "GAA",
    "GAC",
    "GAG",
    "GAT",
    "GCA",
    "GCC",
    "GCG",
    "GCT",
    "GGA",
    "GGC",
    "GGG",
    "GGT",
    "GTA",
    "GTC",
    "GTG",
    "GTT",
    "TAA",
    "TAC",
    "TAG",
    "TAT",
    "TCA",
    "TCC",
    "TCG",
    "TCT",
    "TGA",
    "TGC",
    "TGG",
    "TGT",
    "TTA",
    "TTC",
    "TTG",
    "TTT",
)
CODON_TO_INDEX = {codon: idx for idx, codon in enumerate(CODON_VOCAB)}

STANDARD_GENETIC_CODE = {
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "C": ["TGT", "TGC"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "F": ["TTT", "TTC"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"],
    "I": ["ATT", "ATC", "ATA"],
    "K": ["AAA", "AAG"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "M": ["ATG"],
    "N": ["AAT", "AAC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "Q": ["CAA", "CAG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "W": ["TGG"],
    "Y": ["TAT", "TAC"],
    "*": ["TAA", "TAG", "TGA"],
}

AA_TO_INDICES = {
    aa: [CODON_TO_INDEX[c] for c in codons] for aa, codons in STANDARD_GENETIC_CODE.items()
}

# =====================================================================
# State & Logits Contracts
# =====================================================================


@dataclass(frozen=True)
class ConstraintState:
    """Immutable state for safe Beam Search and lookahead tracking."""

    automaton_node: int = 0
    position: int = 0
    trailing_nt: str = ""  # Last few nucleotides to track homopolymer runs across codons
    gc_count: int = 0  # Cumulative GC count for initiation ramp tracking


class CodonLogits:
    """Adapter Contract enforcing exactly 64 elements and true immutability."""

    def __init__(self, raw_scores: np.ndarray):
        if raw_scores.shape != (64,):
            raise ValueError(f"CodonLogits must strictly have 64 elements, got {raw_scores.shape}.")
        self.scores = np.array(raw_scores, dtype=np.float32, copy=True)
        self.scores.setflags(write=False)


# =====================================================================
# Constraint Processors
# =====================================================================


class BaseConstraintProcessor:
    """Base interface for neuro-symbolic logit processors and state transformers."""

    def __call__(
        self,
        position: int,
        target_aa: str,
        state: ConstraintState,
        logits: CodonLogits,
    ) -> CodonLogits:
        raise NotImplementedError

    def get_next_state(
        self,
        position: int,
        state: ConstraintState,
        chosen_codon: str,
    ) -> ConstraintState:
        """Computes the next state transition given chosen codon."""
        return state


class SynonymousMaskProcessor(BaseConstraintProcessor):
    """Hard Veto: masks all non-synonymous codons to -inf (100% AA fidelity guarantee)."""

    def __call__(
        self,
        position: int,
        target_aa: str,
        state: ConstraintState,
        logits: CodonLogits,
    ) -> CodonLogits:
        valid_indices = AA_TO_INDICES.get(target_aa, [])
        if not valid_indices:
            raise ValueError(f"Invalid amino acid '{target_aa}' at pos {position}.")

        masked = np.copy(logits.scores)
        for i in range(64):
            if i not in valid_indices:
                masked[i] = -math.inf
        return CodonLogits(masked)


class AutomatonConstraintProcessor(BaseConstraintProcessor):
    """Hard Veto: evaluates Aho-Corasick automaton transitions to block forbidden motifs across codon boundaries."""

    def __init__(self, automaton: CompiledAutomaton):
        self.automaton = automaton

    def __call__(
        self,
        position: int,
        target_aa: str,
        state: ConstraintState,
        logits: CodonLogits,
    ) -> CodonLogits:
        masked = np.copy(logits.scores)

        for i in range(64):
            if masked[i] == -math.inf:
                continue

            candidate_codon = CODON_VOCAB[i]
            _, is_forbidden = self.automaton.step_codon(state.automaton_node, candidate_codon)

            if is_forbidden:
                masked[i] = -math.inf

        return CodonLogits(masked)

    def get_next_state(
        self,
        position: int,
        state: ConstraintState,
        chosen_codon: str,
    ) -> ConstraintState:
        next_node, is_forbidden = self.automaton.step_codon(state.automaton_node, chosen_codon)
        if is_forbidden:
            raise ValueError(f"Fatal: Chosen codon '{chosen_codon}' results in a forbidden state!")

        return ConstraintState(
            automaton_node=next_node,
            position=state.position,
            trailing_nt=state.trailing_nt,
            gc_count=state.gc_count,
        )


class HomopolymerConstraintProcessor(BaseConstraintProcessor):
    """Hard Veto: prevents homopolymer runs exceeding configured limit (e.g. > 5-nt)."""

    def __init__(self, max_run: int = 5):
        self.max_run = max_run

    def _has_excess_homopolymer(self, prefix: str, codon: str) -> bool:
        full = prefix + codon
        if not full:
            return False
        current_char = full[0]
        current_len = 1
        for char in full[1:]:
            if char == current_char:
                current_len += 1
                if current_len > self.max_run:
                    return True
            else:
                current_char = char
                current_len = 1
        return False

    def __call__(
        self,
        position: int,
        target_aa: str,
        state: ConstraintState,
        logits: CodonLogits,
    ) -> CodonLogits:
        masked = np.copy(logits.scores)
        trailing = state.trailing_nt

        for i in range(64):
            if masked[i] == -math.inf:
                continue
            codon = CODON_VOCAB[i]
            if self._has_excess_homopolymer(trailing, codon):
                masked[i] = -math.inf

        return CodonLogits(masked)

    def get_next_state(
        self,
        position: int,
        state: ConstraintState,
        chosen_codon: str,
    ) -> ConstraintState:
        new_trailing = (state.trailing_nt + chosen_codon)[-self.max_run :]
        return ConstraintState(
            automaton_node=state.automaton_node,
            position=state.position,
            trailing_nt=new_trailing,
            gc_count=state.gc_count,
        )


class InitiationGCConstraintProcessor(BaseConstraintProcessor):
    """Hard Veto: ensures 5' 45-nt (first 15 codons) GC% stays within allowed bounds."""

    def __init__(
        self,
        ramp_length_codons: int = 15,
        gc_min_frac: float = 0.15,
        gc_max_frac: float = 0.45,
    ):
        self.ramp_length_codons = ramp_length_codons
        self.max_gc_count = int(math.floor(gc_max_frac * (ramp_length_codons * 3)))
        self.min_gc_count = int(math.ceil(gc_min_frac * (ramp_length_codons * 3)))

    def __call__(
        self,
        position: int,
        target_aa: str,
        state: ConstraintState,
        logits: CodonLogits,
    ) -> CodonLogits:
        if position >= self.ramp_length_codons:
            return logits

        masked = np.copy(logits.scores)
        cur_gc = state.gc_count

        for i in range(64):
            if masked[i] == -math.inf:
                continue
            codon = CODON_VOCAB[i]
            codon_gc = codon.count("G") + codon.count("C")
            cand_total_gc = cur_gc + codon_gc

            # Do not exceed max GC count during ramp
            if cand_total_gc > self.max_gc_count:
                masked[i] = -math.inf

        return CodonLogits(masked)

    def get_next_state(
        self,
        position: int,
        state: ConstraintState,
        chosen_codon: str,
    ) -> ConstraintState:
        codon_gc = chosen_codon.count("G") + chosen_codon.count("C")
        return ConstraintState(
            automaton_node=state.automaton_node,
            position=state.position,
            trailing_nt=state.trailing_nt,
            gc_count=state.gc_count + codon_gc,
        )


class ConstraintProcessorPipeline:
    """Enforces ordered sequential constraint processing and state transitions."""

    def __init__(self, processors: List[BaseConstraintProcessor]):
        self.processors = processors

    def __call__(
        self,
        position: int,
        target_aa: str,
        state: ConstraintState,
        logits: CodonLogits,
    ) -> CodonLogits:
        current_logits = logits
        for processor in self.processors:
            current_logits = processor(position, target_aa, state, current_logits)

        return current_logits

    def get_next_state(
        self,
        position: int,
        state: ConstraintState,
        chosen_codon: str,
    ) -> ConstraintState:
        current_state = state
        for processor in self.processors:
            current_state = processor.get_next_state(position, current_state, chosen_codon)
        return current_state
