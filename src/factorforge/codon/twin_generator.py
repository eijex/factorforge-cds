"""Twin Candidate Generator producing W0 (unwatermarked) and W1 (watermarked) pairs.

Enforces:
1. Exact integer quota realization for configured amino acids (Arg, Ser).
2. Complete 5' Signal Peptide preservation (first 60 nt / 20 aa).
3. 100% codon count invariance for Arg & Ser between W0 and W1 (Zero Drift).
4. Watermark embedding restricted exclusively to non-target residues.
5. In-silico QA clearance (Type IIS clean, correct translation, no internal stops).
"""

from __future__ import annotations

import hashlib

from pydantic import BaseModel

from factorforge.codon.allocator import AllocationResult, DeterministicAllocator
from factorforge.codon.policy import (
    AA_TO_CODONS,
    CodonDistributionPolicy,
)
from factorforge.codon.qa import ENZYME_RECOGNITION, InSilicoQA, InSilicoQAReport, translate


def get_preferred_codon(key: str, aa_index: int, available_codons: list[str]) -> str:
    """Deterministic cryptographic pseudo-random codon preference selection."""
    sorted_codons = sorted(available_codons)
    if not sorted_codons:
        return ""
    hash_input = f"{key}_{aa_index}".encode()
    hash_int = int(hashlib.sha256(hash_input).hexdigest(), 16)
    return sorted_codons[hash_int % len(sorted_codons)]


class WatermarkDetector:
    """Evaluates watermark signal presence and statistical detectability (p-value)."""

    def __init__(self, target_protein: str):
        self.target_protein = target_protein.rstrip("*")

    def score_sequence(self, cds: str, key: str) -> float:
        codons = [cds[i : i + 3] for i in range(0, len(cds) - 2, 3)]
        matches = 0
        eligible_positions = 0

        for i, (codon, aa) in enumerate(zip(codons, self.target_protein)):
            available = AA_TO_CODONS.get(aa, [])
            if len(available) > 1:
                eligible_positions += 1
                preferred = get_preferred_codon(key, i, available)
                if codon == preferred:
                    matches += 1

        if eligible_positions == 0:
            return 0.0
        return matches / eligible_positions

    def compute_detectability(
        self, cds: str, correct_key: str, num_null_samples: int = 200
    ) -> dict:
        actual_score = self.score_sequence(cds, correct_key)
        null_scores = []
        for i in range(num_null_samples):
            wrong_key = f"null_key_{i}_{correct_key}"
            null_scores.append(self.score_sequence(cds, wrong_key))

        b = sum(1 for s in null_scores if s >= actual_score)
        p_value = (b + 1) / (num_null_samples + 1)

        return {
            "score": round(actual_score, 4),
            "p_value": round(p_value, 4),
            "null_mean": round(sum(null_scores) / len(null_scores), 4),
            "null_max": round(max(null_scores), 4),
            "is_detected": p_value < 0.05,
        }


class TwinCandidatePair(BaseModel):
    """Container for paired W0 (unwatermarked) and W1 (watermarked) candidates."""

    chain: str
    construct_id_w0: str
    construct_id_w1: str
    cds_w0: str
    cds_w1: str
    protein: str
    policy_id: str
    watermark_key: str
    delta_nt: int
    delta_codons: int
    gc_percent_w0: float
    gc_percent_w1: float
    gc_delta_percent: float
    arg_counts_w0: dict[str, int]
    arg_counts_w1: dict[str, int]
    ser_counts_w0: dict[str, int]
    ser_counts_w1: dict[str, int]
    quotas_identical: bool
    w0_qa: InSilicoQAReport
    w1_qa: InSilicoQAReport
    w0_detection: dict
    w1_detection: dict
    allocations: dict[str, AllocationResult]


class TwinCandidateGenerator:
    """Generates rigorous, matched W0/W1 candidate pairs."""

    def __init__(self, watermark_key: str = "Eijex_FactorForge_Watermark_2026"):
        self.watermark_key = watermark_key

    @staticmethod
    def _has_restriction_site(seq: str) -> bool:
        seq_upper = seq.upper()
        for patterns in ENZYME_RECOGNITION.values():
            for pat in patterns:
                if pat in seq_upper:
                    return True
        return False

    def generate(
        self,
        chain: str,
        protein: str,
        baseline_cds: str,
        policy: CodonDistributionPolicy,
        construct_prefix: str = "FF-HUMIRA",
    ) -> TwinCandidatePair:
        """Generate twin candidate pair W0 and W1 adhering to all invariants."""
        if policy.scope_region != "mature_chain":
            raise ValueError("Only mature_chain generation is currently implemented")
        if (
            not baseline_cds
            or len(baseline_cds) % 3
            or set(baseline_cds) - set("ACGT")
            or translate(baseline_cds) not in {protein.rstrip("*"), protein.rstrip("*") + "*"}
        ):
            raise ValueError("Baseline must match the complete protein and terminal-stop contract")
        clean_protein = protein.rstrip("*")
        sp_aa_len = policy.signal_peptide_aa_len if policy.preserve_signal_peptide else 0
        if sp_aa_len > len(clean_protein):
            raise ValueError("Configured frozen region exceeds the protein length")
        sp_nt_len = sp_aa_len * 3

        sp_cds = baseline_cds[:sp_nt_len]
        mature_prot = clean_protein[sp_aa_len:]

        # Step 1: Run deterministic allocation on mature chain for configured amino acids
        allocations: dict[str, AllocationResult] = {}
        for aa, dist in policy.distributions.items():
            count = mature_prot.count(aa)
            allocations[aa] = DeterministicAllocator.allocate(aa, count, dist)

        # Step 2: Build multiset queues for mature chain Arg and Ser
        codon_pools: dict[str, list[str]] = {}
        for aa, alloc in allocations.items():
            pool = []
            for codon, cnt in alloc.achieved_counts.items():
                pool.extend([codon] * cnt)
            codon_pools[aa] = sorted(pool)

        # Step 3: Construct W0 (Unwatermarked candidate)
        # Distribute Arg and Ser evenly to avoid consecutive repeats
        mature_codons_w0: list[str] = []
        last_codon_for_aa: dict[str, str] = {}

        # Default host preferred codons for N. benthamiana
        host_default = {
            "A": "GCT",
            "C": "TGT",
            "D": "GAT",
            "E": "GAA",
            "F": "TTT",
            "G": "GGA",
            "H": "CAT",
            "I": "ATT",
            "K": "AAG",
            "L": "CTT",
            "M": "ATG",
            "N": "AAT",
            "P": "CCT",
            "Q": "CAA",
            "R": "AGA",
            "S": "TCT",
            "T": "ACT",
            "V": "GTT",
            "W": "TGG",
            "Y": "TAT",
        }

        # Backup synonymous codons for resolving Type IIS sites
        alt_codons = {
            "A": ["GCA", "GCC"],
            "C": ["TGC"],
            "D": ["GAC"],
            "E": ["GAG"],
            "F": ["TTC"],
            "G": ["GGT", "GGC"],
            "H": ["CAC"],
            "I": ["ATC"],
            "K": ["AAA"],
            "L": ["CTC", "TTG"],
            "P": ["CCA", "CCC"],
            "Q": ["CAG"],
            "T": ["ACC", "ACA"],
            "V": ["GTG", "GTC"],
            "Y": ["TAC"],
        }

        for aa in mature_prot:
            if codon_pools.get(aa):
                # Pick a codon different from last pick if possible
                pool = codon_pools[aa]
                picked = None
                for candidate in pool:
                    if candidate != last_codon_for_aa.get(aa):
                        picked = candidate
                        break
                if picked is None:
                    picked = pool[0]
                pool.remove(picked)
                last_codon_for_aa[aa] = picked
                mature_codons_w0.append(picked)
            else:
                chosen = host_default.get(aa, AA_TO_CODONS[aa][0])
                mature_codons_w0.append(chosen)

        # Combine SP + mature + stop codon
        terminal_stop = baseline_cds[-3:] if translate(baseline_cds).endswith("*") else "TAA"
        cds_w0 = sp_cds + "".join(mature_codons_w0) + terminal_stop

        # Resolve any accidental Type IIS restriction sites in W0 by mutating non-target codons
        cds_w0_codons = [cds_w0[i : i + 3] for i in range(0, len(cds_w0), 3)]
        for i in range(sp_aa_len, len(cds_w0_codons) - 1):
            aa = clean_protein[i]
            if aa not in policy.distributions and aa in alt_codons:
                window_start = max(0, (i - 2) * 3)
                window_end = min(len(cds_w0), (i + 3) * 3)
                local_sub = "".join(cds_w0_codons)[window_start:window_end]
                if self._has_restriction_site(local_sub):
                    # Try an alternative codon
                    for alt in alt_codons[aa]:
                        cds_w0_codons[i] = alt
                        new_sub = "".join(cds_w0_codons)[window_start:window_end]
                        if not self._has_restriction_site(new_sub):
                            break

        cds_w0 = "".join(cds_w0_codons)

        # Step 4: Construct W1 (Watermarked counterpart)
        # INVARIANT: Target residues (Arg, Ser) in mature region MUST NOT BE MODIFIED!
        cds_w1_codons = list(cds_w0_codons)
        for i in range(sp_aa_len, len(clean_protein)):
            aa = clean_protein[i]
            # ONLY modify non-target amino acids
            if aa not in policy.distributions:
                available = AA_TO_CODONS.get(aa, [])
                if len(available) > 1:
                    pref = get_preferred_codon(self.watermark_key, i, available)
                    if pref != cds_w1_codons[i]:
                        # Speculatively apply swap
                        old_codon = cds_w1_codons[i]
                        cds_w1_codons[i] = pref

                        # Check if swap introduced a restriction site
                        window_start = max(0, (i - 2) * 3)
                        window_end = min(len(cds_w1_codons) * 3, (i + 3) * 3)
                        local_sub = "".join(cds_w1_codons)[window_start:window_end]
                        if self._has_restriction_site(local_sub):
                            cds_w1_codons[i] = (
                                old_codon  # Revert to maintain restriction site cleanliness
                            )

        cds_w1 = "".join(cds_w1_codons)

        # Step 5: Verify Invariant & QA Checks
        w0_codons_mature = [cds_w0[i : i + 3] for i in range(sp_nt_len, len(cds_w0) - 3, 3)]
        w1_codons_mature = [cds_w1[i : i + 3] for i in range(sp_nt_len, len(cds_w1) - 3, 3)]

        arg_counts_w0 = {c: w0_codons_mature.count(c) for c in AA_TO_CODONS["R"]}
        arg_counts_w1 = {c: w1_codons_mature.count(c) for c in AA_TO_CODONS["R"]}
        ser_counts_w0 = {c: w0_codons_mature.count(c) for c in AA_TO_CODONS["S"]}
        ser_counts_w1 = {c: w1_codons_mature.count(c) for c in AA_TO_CODONS["S"]}

        quotas_identical = (arg_counts_w0 == arg_counts_w1) and (ser_counts_w0 == ser_counts_w1)

        # Calculate deltas between W0 and W1
        delta_nt = sum(1 for c0, c1 in zip(cds_w0, cds_w1) if c0 != c1)
        w0_all_codons = [cds_w0[i : i + 3] for i in range(0, len(cds_w0), 3)]
        w1_all_codons = [cds_w1[i : i + 3] for i in range(0, len(cds_w1), 3)]
        delta_codons = sum(1 for c0, c1 in zip(w0_all_codons, w1_all_codons) if c0 != c1)

        # QA reports
        construct_id_w0 = f"{construct_prefix}-{chain}-W0"
        construct_id_w1 = f"{construct_prefix}-{chain}-W1"

        w0_qa = InSilicoQA.evaluate(cds_w0, clean_protein, construct_id_w0, sp_cds)
        w1_qa = InSilicoQA.evaluate(cds_w1, clean_protein, construct_id_w1, sp_cds)

        # Watermark detectability
        detector = WatermarkDetector(clean_protein)
        w0_detection = detector.compute_detectability(cds_w0, self.watermark_key)
        w1_detection = detector.compute_detectability(cds_w1, self.watermark_key)
        if not w0_qa.overall_passed or not w1_qa.overall_passed or not quotas_identical:
            raise ValueError("Twin candidate failed computational QA or quota checks")
        if not w1_detection["is_detected"]:
            raise ValueError("insufficient_capacity: watermark detector contract not met")

        gc_percent_w0 = w0_qa.gc_percent
        gc_percent_w1 = w1_qa.gc_percent
        gc_delta_percent = round(w1_qa.gc_percent - w0_qa.gc_percent, 2)

        return TwinCandidatePair(
            chain=chain,
            construct_id_w0=construct_id_w0,
            construct_id_w1=construct_id_w1,
            cds_w0=cds_w0,
            cds_w1=cds_w1,
            protein=clean_protein,
            policy_id=policy.policy_id,
            watermark_key=self.watermark_key,
            delta_nt=delta_nt,
            delta_codons=delta_codons,
            gc_percent_w0=gc_percent_w0,
            gc_percent_w1=gc_percent_w1,
            gc_delta_percent=gc_delta_percent,
            arg_counts_w0=arg_counts_w0,
            arg_counts_w1=arg_counts_w1,
            ser_counts_w0=ser_counts_w0,
            ser_counts_w1=ser_counts_w1,
            quotas_identical=quotas_identical,
            w0_qa=w0_qa,
            w1_qa=w1_qa,
            w0_detection=w0_detection,
            w1_detection=w1_detection,
            allocations=allocations,
        )
