"""Diverse Candidate Generator & Hard Invariant Gate for FactorForge Slate v2 (Job 293A).

Implements Stage 0 (Diverse Heuristic Beam Generation v1) and Stage 1 (Hard Invariant Gate)
ensuring that amino acid identity is strictly guaranteed by construction and verified post-generation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import random
import re
from typing import Any, Dict, List, Optional, Set, Tuple

from factorforge.analysis.metrics import (
    STANDARD_GENETIC_CODE,
    calculate_cai,
    calculate_gc,
    calculate_gc_windows,
)
from factorforge.core.host_model import HostModel, STOP_CODONS
from factorforge.utils.restriction_sites import detect_restriction_sites

DEFAULT_FORBIDDEN_SITES = [
    {"name": "BsaI", "sequence": "GGTCTC", "scan_rc": True},
    {"name": "BsmBI", "sequence": "CGTCTC", "scan_rc": True},
    {"name": "NotI", "sequence": "GCGGCCGC", "scan_rc": True},
    {"name": "XhoI", "sequence": "CTCGAG", "scan_rc": True},
]

SUPPORTED_GENERATION_MODES = {"heuristic_beam_v1", "deterministic_beam"}


@dataclass
class CandidateDesign:
    """A generated candidate CDS design with generation metadata."""

    candidate_id: str
    protein_sequence: str
    dna_sequence: str
    coding_sequence: str
    terminal_stop: str
    profile_bias: str
    sha256: str
    cai: float = 0.0
    gc_global: float = 0.0
    gc_5p_ramp: float = 0.0
    rare_codon_count: int = 0
    raw_score: float = 0.0
    metadata: Dict[str, Any] = field(default_factory=dict)


def normalize_protein_sequence(raw_sequence: str) -> str:
    """Normalize input protein sequence once at engine entry.

    Rules:
    - Uppercase & whitespace stripped
    - At most one trailing '*' allowed and stripped from sense sequence
    - Internal '*' strictly rejected with ValueError
    """
    clean = "".join(raw_sequence.upper().split())
    if not clean:
        raise ValueError("Protein sequence is required and cannot be empty.")

    # Check for internal stop codons
    if "*" in clean[:-1]:
        raise ValueError("Internal stop codon '*' is not permitted in protein sequence.")

    # Strip single terminal stop if present
    return clean.rstrip("*")


def translate_dna(dna: str) -> str:
    """Translate DNA coding sequence into amino acid string strictly."""
    dna_clean = "".join(dna.upper().split())
    if len(dna_clean) % 3 != 0:
        raise ValueError(f"DNA sequence length ({len(dna_clean)}) must be divisible by 3.")
    aa_list = []
    for i in range(0, len(dna_clean), 3):
        codon = dna_clean[i : i + 3]
        aa = STANDARD_GENETIC_CODE.get(codon, "X")
        if aa == "*":
            # stop codon
            break
        aa_list.append(aa)
    return "".join(aa_list)


class DiverseCandidateGenerator:
    """Generates a rich, deduplicated pool of candidates across multiple objective biases using HostModel."""

    def __init__(
        self,
        host_model: HostModel,
        target_gc: Optional[float] = None,
    ) -> None:
        self.host_model = host_model
        self.target_gc = target_gc if target_gc is not None else host_model.default_target_gc
        self.aa_to_codons = host_model.aa_to_codons
        self.codon_frequencies = host_model.codon_frequencies
        self.codon_weights = host_model.codon_weights

    def generate_pool(
        self,
        protein_sequence: str,
        target_pool_size: int = 512,
        source_cds: Optional[str] = None,
        stop_policy: str = "append_preferred",
        seed: int = 42,
        generation_mode: str = "heuristic_beam_v1",
    ) -> List[CandidateDesign]:
        """Generate deduplicated candidates across multiple objective biases."""
        # Canonical normalization
        canonical_protein = normalize_protein_sequence(protein_sequence)

        # Validate generation mode
        if generation_mode not in SUPPORTED_GENERATION_MODES:
            raise ValueError(
                f"Unsupported generation_mode '{generation_mode}'. "
                f"Supported modes: {sorted(SUPPORTED_GENERATION_MODES)}"
            )

        # Determine terminal stop codon
        terminal_stop = self._resolve_stop_codon(stop_policy, source_cds)

        candidates: List[CandidateDesign] = []
        seen_hashes: Set[str] = set()

        # Generate candidates across 4 distinct objective biases
        biases = [
            ("balanced", self.target_gc, 0.40, 0.25, 0.20),
            ("high_cai", self.target_gc, 0.80, 0.10, 0.05),
            ("open_5p_ramp", max(0.35, self.target_gc - 0.05), 0.30, 0.15, 0.45),
            ("strict_gc", self.target_gc, 0.25, 0.60, 0.10),
        ]

        per_bias_count = max(64, target_pool_size // len(biases))

        for bias_idx, (bias_name, target_gc_bias, w_cai, w_gc, w_ramp) in enumerate(biases):
            # Pass deterministic seed offset per bias
            bias_seed = (seed + (bias_idx * 10007)) & 0xFFFFFFFF
            bias_candidates = self._beam_search_bias(
                protein=canonical_protein,
                terminal_stop=terminal_stop,
                bias_name=bias_name,
                target_gc=target_gc_bias,
                w_cai=w_cai,
                w_gc=w_gc,
                w_ramp=w_ramp,
                beam_width=per_bias_count,
                seed=bias_seed,
            )
            for cand in bias_candidates:
                if cand.sha256 not in seen_hashes:
                    seen_hashes.add(cand.sha256)
                    candidates.append(cand)

        # Ensure deterministic candidate ordering
        candidates.sort(key=lambda c: (-c.raw_score, c.candidate_id))
        return candidates

    def _resolve_stop_codon(self, stop_policy: str, source_cds: Optional[str]) -> str:
        if stop_policy == "exclude":
            return ""
        if stop_policy.startswith("specified:"):
            specified = stop_policy.split(":", 1)[1].upper().strip()
            if specified in STOP_CODONS:
                return specified
            raise ValueError(f"Invalid stop codon '{specified}'. Allowed: {sorted(STOP_CODONS)}")
        if stop_policy == "preserve_source":
            if not source_cds:
                raise ValueError("stop_policy 'preserve_source' requires non-empty source_cds.")
            clean_source = "".join(source_cds.upper().split())
            if len(clean_source) >= 3 and clean_source[-3:] in STOP_CODONS:
                return clean_source[-3:]
            raise ValueError(f"source_cds does not end with a recognized stop codon: {source_cds[-3:] if len(source_cds)>=3 else source_cds}")
        if stop_policy == "append_preferred":
            return self.host_model.preferred_stop
        raise ValueError(
            f"Unsupported stop_policy '{stop_policy}'. "
            f"Allowed: append_preferred, preserve_source, exclude, specified:<TAA|TAG|TGA>"
        )

    def _beam_search_bias(
        self,
        protein: str,
        terminal_stop: str,
        bias_name: str,
        target_gc: float,
        w_cai: float,
        w_gc: float,
        w_ramp: float,
        beam_width: int,
        seed: int,
    ) -> List[CandidateDesign]:
        """Deterministic heuristic beam search v1 constructing candidates by sampling synonymous codons with seed perturbation."""
        rng = random.Random(seed)
        beam: List[Tuple[float, int, List[str]]] = [(0.0, 0, [])]
        num_residues = len(protein)

        # Forbidden site sequences to avoid at junctions
        forbidden_patterns = ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG", "GCGGCCGC", "CTCGAG"]

        for pos, aa in enumerate(protein):
            possible_codons = sorted(self.aa_to_codons.get(aa, []))
            if not possible_codons:
                raise ValueError(f"Unknown amino acid: {aa} at position {pos + 1}")

            next_beam: List[Tuple[float, int, List[str]]] = []

            for acc_score, gc_count, prefix in beam:
                for codon in possible_codons:
                    codon_weight = self.codon_weights.get(codon, 0.05)
                    codon_freq = self.codon_frequencies.get(codon, 0.05)
                    codon_gc = codon.count("G") + codon.count("C")

                    # Check forbidden site avoidance at junction
                    tail = ("".join(prefix[-3:]) + codon) if len(prefix) >= 3 else ("".join(prefix) + codon)
                    has_forbidden = any(pat in tail for pat in forbidden_patterns)
                    if has_forbidden:
                        continue

                    # Positional 5' ramp bonus (first 10 codons)
                    ramp_bonus = 0.0
                    if pos < 10:
                        au_count = codon.count("A") + codon.count("T")
                        ramp_bonus = (au_count / 3.0) * w_ramp

                    # Running GC fidelity
                    total_bases = (pos + 1) * 3
                    running_gc_ratio = (gc_count + codon_gc) / total_bases
                    gc_fidelity = 1.0 - abs(running_gc_ratio - target_gc)

                    # Bounded seed perturbation for exploration (preserving biological objective dominance)
                    seed_perturb = rng.uniform(-0.015, 0.015)

                    # Step score
                    step_score = (
                        w_cai * codon_weight
                        + w_gc * gc_fidelity
                        + ramp_bonus
                        + (0.01 * codon_freq)
                        + seed_perturb
                    )

                    new_acc_score = acc_score + step_score
                    next_beam.append(
                        (new_acc_score, gc_count + codon_gc, prefix + [codon])
                    )

            # Sort deterministically by descending score, then tie-break by codon string
            next_beam.sort(
                key=lambda item: (-item[0], "".join(item[2]))
            )
            beam = next_beam[:beam_width]

        # Convert top beam paths into CandidateDesign objects
        candidates: List[CandidateDesign] = []
        for idx, (score, _, codons) in enumerate(beam):
            coding_seq = "".join(codons)
            full_dna = coding_seq + terminal_stop
            sha = hashlib.sha256(full_dna.encode("ascii")).hexdigest()
            cai_val = calculate_cai(coding_seq, self.codon_weights)
            gc_val = calculate_gc(coding_seq) / 100.0
            ramp_gc = (calculate_gc(coding_seq[:48]) / 100.0) if len(coding_seq) >= 48 else gc_val

            rare_count = sum(
                1 for c in codons if self.codon_frequencies.get(c, 0.0) < 0.10
            )

            cand = CandidateDesign(
                candidate_id=f"cand_{bias_name}_{idx+1:03d}",
                protein_sequence=protein,
                dna_sequence=full_dna,
                coding_sequence=coding_seq,
                terminal_stop=terminal_stop,
                profile_bias=bias_name,
                sha256=sha,
                cai=round(float(cai_val), 4),
                gc_global=round(float(gc_val), 4),
                gc_5p_ramp=round(float(ramp_gc), 4),
                rare_codon_count=rare_count,
                raw_score=round(float(score / max(1, num_residues)), 5),
                metadata={
                    "bias_name": bias_name,
                    "target_gc": target_gc,
                    "beam_rank": idx + 1,
                    "seed": seed,
                },
            )
            candidates.append(cand)

        return candidates


class HardInvariantGate:
    """Filters candidate designs against strict biological and synthesis constraints."""

    def __init__(
        self,
        forbidden_sites: Optional[List[Dict[str, Any]]] = None,
        extreme_gc_guard: Tuple[float, float] = (20.0, 80.0),
    ) -> None:
        self.forbidden_sites = forbidden_sites or DEFAULT_FORBIDDEN_SITES
        self.min_extreme_gc, self.max_extreme_gc = extreme_gc_guard

    def evaluate(self, candidate: CandidateDesign) -> Tuple[bool, List[str]]:
        """Evaluate a candidate against all hard invariants."""
        violations: List[str] = []

        # 1. Postcondition AA translation verification
        try:
            translated = translate_dna(candidate.coding_sequence)
            if translated != candidate.protein_sequence:
                violations.append(
                    f"AA translation mismatch: expected {len(candidate.protein_sequence)} AA, got {len(translated)} AA"
                )
        except Exception as e:
            violations.append(f"Translation failure: {str(e)}")

        # 2. Forbidden restriction sites screening (Type IIS + Additional)
        hits = detect_restriction_sites(candidate.dna_sequence, self.forbidden_sites, scan_rc=True)
        if hits:
            site_names = sorted({hit.get("name", "Unknown") for hit in hits})
            violations.append(f"Forbidden restriction site detected: {', '.join(site_names)}")

        # 3. Extreme local GC window guard (50-nt sliding window)
        if len(candidate.coding_sequence) >= 50:
            windows = calculate_gc_windows(candidate.coding_sequence, window_size=50, step=10)
            for w in windows:
                w_gc = float(w.get("gc", w.get("gc_percent", 50.0)))
                if w_gc < self.min_extreme_gc or w_gc > self.max_extreme_gc:
                    violations.append(
                        f"Local GC out of extreme bound ({w_gc:.1f}% at nt {w.get('start', 0)})"
                    )
                    break

        return (len(violations) == 0, violations)

    def filter_pool(
        self, candidates: List[CandidateDesign]
    ) -> Tuple[List[CandidateDesign], List[Dict[str, Any]]]:
        """Filter a list of candidates, returning (feasible_candidates, rejected_logs)."""
        feasible: List[CandidateDesign] = []
        rejected: List[Dict[str, Any]] = []

        for cand in candidates:
            passed, violations = self.evaluate(cand)
            if passed:
                feasible.append(cand)
            else:
                rejected.append({
                    "candidate_id": cand.candidate_id,
                    "sha256": cand.sha256,
                    "violations": violations,
                })

        return feasible, rejected

    def filter_candidates(self, candidates: List[CandidateDesign]) -> List[CandidateDesign]:
        """Filter candidates and return list of feasible candidates."""
        feasible, _ = self.filter_pool(candidates)
        return feasible
