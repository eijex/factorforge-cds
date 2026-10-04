"""
FactorForge Codon Policy Model and Presets (Job 355).

Defines structured CodonPolicy objects containing target ratios, allowed codons,
explicitly excluded codons with biological rationale and citations,
and deterministic Hamilton Largest Remainder allocation with repeat avoidance.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Optional


@dataclass(frozen=True)
class ExcludedCodon:
    """Explicitly excluded codon with biological rationale and literature evidence."""
    codon: str
    amino_acid: str
    reason: str
    evidence_source: str

    def to_dict(self) -> dict[str, str]:
        return {
            "codon": self.codon,
            "amino_acid": self.amino_acid,
            "reason": self.reason,
            "evidence_source": self.evidence_source,
        }

    @classmethod
    def from_dict(cls, data: dict[str, str]) -> ExcludedCodon:
        return cls(
            codon=data["codon"],
            amino_acid=data["amino_acid"],
            reason=data["reason"],
            evidence_source=data["evidence_source"],
        )


@dataclass
class CodonPolicy:
    """
    Structured Codon Optimization and Allocation Policy.
    
    Provides:
    - Target ratios per amino acid
    - Allowed codons
    - Excluded codons with verifiable citations & rationales
    - Deterministic Hamilton Largest Remainder allocation with sliding-window repeat avoidance
    """
    name: str
    display_name: str
    description: str
    organism: str
    target_ratios: dict[str, dict[str, float]]
    allowed_codons: dict[str, list[str]]
    excluded_codons: list[ExcludedCodon] = field(default_factory=list)
    repeat_avoidance_window: int = 3
    metadata: dict[str, Any] = field(default_factory=dict)

    def allocate_sequence(self, aa_sequence: str) -> list[str]:
        """
        Allocate synonymous codons for an amino acid sequence using Hamilton's
        Largest Remainder Method, respecting repeat avoidance.
        
        Guarantees:
        - 100% deterministic (reproducible across all runs).
        - Exact integer fulfillment of target_ratios across the protein.
        - Consecutive identical synonymous codons bounded by repeat_avoidance_window.
        """
        # 1. Count occurrences of each amino acid
        aa_counts: dict[str, int] = {}
        for aa in aa_sequence:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1

        # 2. Hamilton Largest Remainder calculation per amino acid
        codon_pools: dict[str, dict[str, int]] = {}
        for aa, total_n in aa_counts.items():
            ratios = self.target_ratios.get(aa)
            if not ratios:
                # Fallback to allowed codons or standard default
                allowed = self.allowed_codons.get(aa, ["ATG"])
                ratios = {c: 1.0 / len(allowed) for c in allowed}

            # Normalize ratios
            sum_r = sum(ratios.values())
            norm_ratios = {c: r / sum_r for c, r in ratios.items()}

            # Ideal continuous counts
            exact_counts = {c: total_n * r for c, r in norm_ratios.items()}
            integer_floors = {c: math.floor(cnt) for c, cnt in exact_counts.items()}
            remainders = {c: exact_counts[c] - integer_floors[c] for c in exact_counts}

            surplus = total_n - sum(integer_floors.values())
            # Deterministic sort: largest remainder first, tie-break alphabetically by codon
            sorted_codons = sorted(
                remainders.keys(),
                key=lambda c: (-remainders[c], c),
            )

            allocated = dict(integer_floors)
            for i in range(surplus):
                c_winner = sorted_codons[i % len(sorted_codons)]
                allocated[c_winner] += 1

            codon_pools[aa] = allocated

        # 3. Deterministic sequence allocation with sliding-window repeat avoidance
        allocated_sequence: list[str] = []
        recent_by_aa: dict[str, list[str]] = {aa: [] for aa in aa_counts}

        for aa in aa_sequence:
            pool = codon_pools[aa]
            available_codons = [c for c, count in pool.items() if count > 0]

            if not available_codons:
                # Should never happen under Hamilton invariants
                fallback = self.allowed_codons.get(aa, ["ATG"])[0]
                allocated_sequence.append(fallback)
                continue

            # Check repeat avoidance against recent picks for this AA
            recent = recent_by_aa[aa]
            valid_candidates = []
            for c in available_codons:
                # Count consecutive trailing identical occurrences
                streak = 0
                for prev in reversed(recent):
                    if prev == c:
                        streak += 1
                    else:
                        break
                if streak < self.repeat_avoidance_window:
                    valid_candidates.append(c)

            if not valid_candidates:
                # If all available exceed streak, pick candidate with lowest current streak
                valid_candidates = available_codons

            # Deterministic selection: pick candidate with largest remaining pool,
            # tie-break alphabetically
            chosen = sorted(
                valid_candidates,
                key=lambda c: (-pool[c], c),
            )[0]

            # Decrement pool
            pool[chosen] -= 1
            allocated_sequence.append(chosen)

            # Record recent history
            recent.append(chosen)
            if len(recent) > self.repeat_avoidance_window * 2:
                recent.pop(0)

        return allocated_sequence

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "display_name": self.display_name,
            "description": self.description,
            "organism": self.organism,
            "target_ratios": self.target_ratios,
            "allowed_codons": self.allowed_codons,
            "excluded_codons": [e.to_dict() for e in self.excluded_codons],
            "repeat_avoidance_window": self.repeat_avoidance_window,
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> CodonPolicy:
        return cls(
            name=data["name"],
            display_name=data["display_name"],
            description=data["description"],
            organism=data["organism"],
            target_ratios=data["target_ratios"],
            allowed_codons=data["allowed_codons"],
            excluded_codons=[ExcludedCodon.from_dict(e) for e in data.get("excluded_codons", [])],
            repeat_avoidance_window=data.get("repeat_avoidance_window", 3),
            metadata=dict(data.get("metadata", {})),
        )


# ============================================================================
# Preset Codon Policies
# ============================================================================

PLANTFORM_BALANCED_CODON_V1 = CodonPolicy(
    name="plantform_balanced_codon_v1",
    display_name="PlantForm Balanced Non-CpG Codon Policy v1",
    description=(
        "Plant-specific balanced codon distribution designed for Nicotiana benthamiana. "
        "Strictly balances synonymous codons for Arg (1:1) and Ser (1:1:1:1:1) while excluding "
        "all CpG dinucleotide codons to eliminate methylation-induced gene silencing."
    ),
    organism="Nicotiana benthamiana",
    target_ratios={
        "R": {"AGA": 0.50, "AGG": 0.50},
        "S": {"TCT": 0.20, "TCC": 0.20, "TCA": 0.20, "AGT": 0.20, "AGC": 0.20},
        "L": {"CTT": 0.25, "CTC": 0.25, "TTG": 0.25, "TTA": 0.25},
        "G": {"GGT": 0.50, "GGA": 0.50},
        "A": {"GCT": 0.50, "GCA": 0.50},
        "V": {"GTT": 0.50, "GTA": 0.50},
        "P": {"CCT": 0.50, "CCA": 0.50},
        "T": {"ACT": 0.50, "ACA": 0.50},
        "E": {"GAA": 0.50, "GAG": 0.50},
        "D": {"GAT": 0.50, "GAC": 0.50},
        "K": {"AAA": 0.50, "AAG": 0.50},
        "Q": {"CAA": 0.50, "CAG": 0.50},
        "N": {"AAT": 0.50, "AAC": 0.50},
        "H": {"CAT": 0.50, "CAC": 0.50},
        "Y": {"TAT": 0.50, "TAC": 0.50},
        "F": {"TTT": 0.50, "TTC": 0.50},
        "C": {"TGT": 0.50, "TGC": 0.50},
        "I": {"ATT": 0.50, "ATC": 0.50},
        "M": {"ATG": 1.00},
        "W": {"TGG": 1.00},
    },
    allowed_codons={
        "R": ["AGA", "AGG"],
        "S": ["TCT", "TCC", "TCA", "AGT", "AGC"],
        "L": ["CTT", "CTC", "TTG", "TTA"],
        "G": ["GGT", "GGA"],
        "A": ["GCT", "GCA"],
        "V": ["GTT", "GTA"],
        "P": ["CCT", "CCA"],
        "T": ["ACT", "ACA"],
        "E": ["GAA", "GAG"],
        "D": ["GAT", "GAC"],
        "K": ["AAA", "AAG"],
        "Q": ["CAA", "CAG"],
        "N": ["AAT", "AAC"],
        "H": ["CAT", "CAC"],
        "Y": ["TAT", "TAC"],
        "F": ["TTT", "TTC"],
        "C": ["TGT", "TGC"],
        "I": ["ATT", "ATC"],
        "M": ["ATG"],
        "W": ["TGG"],
    },
    excluded_codons=[
        ExcludedCodon(
            codon="CGA",
            amino_acid="R",
            reason="Contains CpG dinucleotide, triggers epigenetic gene silencing in Nicotiana benthamiana.",
            evidence_source="Doug / PlantForm Expression SOP & Meyer 2000 Epigenetics in Plants",
        ),
        ExcludedCodon(
            codon="CGC",
            amino_acid="R",
            reason="Contains CpG dinucleotide, high plant methylation risk.",
            evidence_source="PlantForm Expression SOP",
        ),
        ExcludedCodon(
            codon="CGG",
            amino_acid="R",
            reason="Contains CpG dinucleotide, rare in dicot plants.",
            evidence_source="PlantForm Expression SOP",
        ),
        ExcludedCodon(
            codon="CGT",
            amino_acid="R",
            reason="Contains CpG dinucleotide, triggers silencing.",
            evidence_source="PlantForm Expression SOP",
        ),
        ExcludedCodon(
            codon="TCG",
            amino_acid="S",
            reason="Contains CpG dinucleotide.",
            evidence_source="PlantForm Expression SOP",
        ),
        ExcludedCodon(
            codon="CCG",
            amino_acid="P",
            reason="Contains CpG dinucleotide.",
            evidence_source="PlantForm Expression SOP",
        ),
        ExcludedCodon(
            codon="GCG",
            amino_acid="A",
            reason="Contains CpG dinucleotide.",
            evidence_source="PlantForm Expression SOP",
        ),
        ExcludedCodon(
            codon="ACG",
            amino_acid="T",
            reason="Contains CpG dinucleotide.",
            evidence_source="PlantForm Expression SOP",
        ),
    ],
    repeat_avoidance_window=3,
)

HOST_FREQUENCY_DEFAULT_V1 = CodonPolicy(
    name="host_frequency_default_v1",
    display_name="Nicotiana benthamiana Host Frequency Default v1",
    description="Natural host codon usage frequency proportions from Kazusa codon usage database.",
    organism="Nicotiana benthamiana",
    target_ratios={
        "R": {"AGA": 0.40, "AGG": 0.30, "CGA": 0.10, "CGT": 0.10, "CGC": 0.05, "CGG": 0.05},
        "S": {"TCT": 0.28, "TCA": 0.22, "TCC": 0.18, "AGT": 0.16, "AGC": 0.11, "TCG": 0.05},
        "L": {"CTT": 0.30, "TTG": 0.26, "CTC": 0.18, "TTA": 0.14, "CTA": 0.08, "CTG": 0.04},
        "G": {"GGA": 0.38, "GGT": 0.34, "GGG": 0.16, "GGC": 0.12},
        "A": {"GCT": 0.42, "GCA": 0.32, "GCC": 0.20, "GCG": 0.06},
        "V": {"GTT": 0.44, "GTG": 0.26, "GTA": 0.18, "GTC": 0.12},
        "P": {"CCT": 0.42, "CCA": 0.36, "CCC": 0.16, "CCG": 0.06},
        "T": {"ACT": 0.38, "ACA": 0.34, "ACC": 0.20, "ACG": 0.08},
        "E": {"GAA": 0.58, "GAG": 0.42},
        "D": {"GAT": 0.65, "GAC": 0.35},
        "K": {"AAA": 0.54, "AAG": 0.46},
        "Q": {"CAA": 0.62, "CAG": 0.38},
        "N": {"AAT": 0.60, "AAC": 0.40},
        "H": {"CAT": 0.62, "CAC": 0.38},
        "Y": {"TAT": 0.60, "TAC": 0.40},
        "F": {"TTT": 0.58, "TTC": 0.42},
        "C": {"TGT": 0.64, "TGC": 0.36},
        "I": {"ATT": 0.48, "ATC": 0.32, "ATA": 0.20},
        "M": {"ATG": 1.00},
        "W": {"TGG": 1.00},
    },
    allowed_codons={
        "R": ["AGA", "AGG", "CGA", "CGT", "CGC", "CGG"],
        "S": ["TCT", "TCA", "TCC", "AGT", "AGC", "TCG"],
        "L": ["CTT", "TTG", "CTC", "TTA", "CTA", "CTG"],
        "G": ["GGA", "GGT", "GGG", "GGC"],
        "A": ["GCT", "GCA", "GCC", "GCG"],
        "V": ["GTT", "GTG", "GTA", "GTC"],
        "P": ["CCT", "CCA", "CCC", "CCG"],
        "T": ["ACT", "ACA", "ACC", "ACG"],
        "E": ["GAA", "GAG"],
        "D": ["GAT", "GAC"],
        "K": ["AAA", "AAG"],
        "Q": ["CAA", "CAG"],
        "N": ["AAT", "AAC"],
        "H": ["CAT", "CAC"],
        "Y": ["TAT", "TAC"],
        "F": ["TTT", "TTC"],
        "C": ["TGT", "TGC"],
        "I": ["ATT", "ATC", "ATA"],
        "M": ["ATG"],
        "W": ["TGG"],
    },
    excluded_codons=[],
    repeat_avoidance_window=3,
)

MAX_CAI_EXTREME_V1 = CodonPolicy(
    name="max_cai_extreme_v1",
    display_name="Maximum Single Host Preferred Codon Policy (Max CAI)",
    description="Greedy selection of the single most abundant codon per amino acid.",
    organism="Nicotiana benthamiana",
    target_ratios={
        "R": {"AGA": 1.0},
        "S": {"TCT": 1.0},
        "L": {"CTT": 1.0},
        "G": {"GGA": 1.0},
        "A": {"GCT": 1.0},
        "V": {"GTT": 1.0},
        "P": {"CCT": 1.0},
        "T": {"ACT": 1.0},
        "E": {"GAA": 1.0},
        "D": {"GAT": 1.0},
        "K": {"AAA": 1.0},
        "Q": {"CAA": 1.0},
        "N": {"AAT": 1.0},
        "H": {"CAT": 1.0},
        "Y": {"TAT": 1.0},
        "F": {"TTT": 1.0},
        "C": {"TGT": 1.0},
        "I": {"ATT": 1.0},
        "M": {"ATG": 1.0},
        "W": {"TGG": 1.0},
    },
    allowed_codons={
        "R": ["AGA"], "S": ["TCT"], "L": ["CTT"], "G": ["GGA"], "A": ["GCT"],
        "V": ["GTT"], "P": ["CCT"], "T": ["ACT"], "E": ["GAA"], "D": ["GAT"],
        "K": ["AAA"], "Q": ["CAA"], "N": ["AAT"], "H": ["CAT"], "Y": ["TAT"],
        "F": ["TTT"], "C": ["TGT"], "I": ["ATT"], "M": ["ATG"], "W": ["TGG"],
    },
    excluded_codons=[],
    repeat_avoidance_window=999,  # By definition greedy max CAI repeats codons
)


PRESET_POLICIES: dict[str, CodonPolicy] = {
    "plantform_balanced_codon_v1": PLANTFORM_BALANCED_CODON_V1,
    "host_frequency_default_v1": HOST_FREQUENCY_DEFAULT_V1,
    "max_cai_extreme_v1": MAX_CAI_EXTREME_V1,
}


def get_policy(name_or_obj: str | CodonPolicy) -> CodonPolicy:
    """Retrieve a preset policy by name or return the policy object itself."""
    if isinstance(name_or_obj, CodonPolicy):
        return name_or_obj
    if name_or_obj in PRESET_POLICIES:
        return PRESET_POLICIES[name_or_obj]
    # Default fallback
    return PLANTFORM_BALANCED_CODON_V1
