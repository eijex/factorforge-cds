"""Deterministic Integer Codon Allocator using Largest Remainder Method."""

from __future__ import annotations

import math

from pydantic import BaseModel

from factorforge.codon.policy import CodonDistributionPolicy


class AllocationResult(BaseModel):
    """Result of integer codon allocation for a specific amino acid."""

    aa: str
    total_residues: int
    target_fractions: dict[str, float]
    achieved_counts: dict[str, int]
    achieved_fractions: dict[str, float]


class DeterministicAllocator:
    """Allocates integer counts of synonymous codons deterministically."""

    @staticmethod
    def allocate(
        aa: str, total_count: int, target_distribution: dict[str, float]
    ) -> AllocationResult:
        """Allocate total_count among codons according to largest remainder method.

        Args:
            aa: Amino acid 1-letter code (e.g. 'R', 'S').
            total_count: Total number of residues of this amino acid in target region.
            target_distribution: Target fractions, summing to 1.0.

        Returns:
            AllocationResult with exact integer counts and achieved fractions.
        """
        if type(total_count) is not int or total_count < 0:
            raise ValueError("Residue count must be a nonnegative integer")
        CodonDistributionPolicy(
            policy_id="allocation_validation", distributions={aa: target_distribution}
        )
        if total_count == 0:
            return AllocationResult(
                aa=aa,
                total_residues=0,
                target_fractions=dict(target_distribution),
                achieved_counts={c: 0 for c in target_distribution},
                achieved_fractions={c: 0.0 for c in target_distribution},
            )

        codons = sorted(target_distribution.keys())
        floors: dict[str, int] = {}
        remainders: dict[str, float] = {}

        total_floor = 0
        for c in codons:
            exact_quota = total_count * target_distribution[c]
            floor_val = math.floor(exact_quota)
            rem = exact_quota - floor_val
            floors[c] = floor_val
            remainders[c] = round(rem, 8)
            total_floor += floor_val

        unassigned = total_count - total_floor
        if not 0 <= unassigned <= len(codons):
            raise ValueError("Distribution precision cannot realize the requested integer quota")

        # Sort codons by remainder descending, tie-breaking by lexical codon order
        # Key: (-remainder, codon_string)
        sorted_by_rem = sorted(codons, key=lambda c: (-remainders[c], c))

        achieved_counts = dict(floors)
        for i in range(unassigned):
            c_recipient = sorted_by_rem[i % len(sorted_by_rem)]
            achieved_counts[c_recipient] += 1

        achieved_fractions = {c: round(achieved_counts[c] / total_count, 6) for c in codons}

        return AllocationResult(
            aa=aa,
            total_residues=total_count,
            target_fractions=dict(target_distribution),
            achieved_counts=achieved_counts,
            achieved_fractions=achieved_fractions,
        )
