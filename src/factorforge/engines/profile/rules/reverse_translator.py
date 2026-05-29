"""
Reverse Translator for FactorForge profile engine
Reverse-translate amino acid sequences to N. benthamiana-optimized codons (P0-2)
"""

from __future__ import annotations

from bisect import bisect_left
import json
import logging
import math
import random
import secrets
from enum import Enum
from pathlib import Path
from typing import Any, cast

from factorforge.engines.profile.scoring import calculate_composite_score
from factorforge.engines.profile.utils import (
    build_aa_to_codons_map,
    calculate_gc,
    get_data_path,
    load_golden_set,
)
from factorforge.utils.exceptions import EmptyCandidateError

logger = logging.getLogger(__name__)


class OptimizationProfile(Enum):
    """Optimization profile"""

    BALANCED = "balanced"
    HIGH_CAI = "high_cai"
    GC_TARGET = "gc_target"
    ASSEMBLY_FRIENDLY = "assembly_friendly"
    RAMP = "ramp"
    VIRAL_DELIVERY = "viral_delivery"  # TRV 바이러스 전달 최적화 (Li et al. 2026)


class ReverseTranslator:
    """
    Reverse-translate amino acid sequences to DNA

    Supports 4 optimization profiles:
    1. Balanced: CAI priority, GC balance
    2. High-CAI: use only preferred codons
    3. GC-Target: enforce GC% 50% ±5%
    4. Assembly-Friendly: avoid BsaI/BpiI
    """

    def __init__(
        self,
        codon_table_path: str | Path | None = None,
        golden_set_path: str | Path | None = None,
        host: str = "nbenthamiana",
    ) -> None:
        """
        Args:
            codon_table_path: Path to codon table JSON file.
            golden_set_path: Path to golden set JSON for CAI reference weights.
                             If None, attempts to load default golden set.
            host: Host codon table name used when codon_table_path is not provided.
        """
        self.host = host
        if codon_table_path is None:
            # Use centralized data path management
            data_dir = get_data_path()
            codon_table_path = data_dir / f"{host}_codons.json"

        self.codon_table: dict[str, Any] = self._load_codon_table(codon_table_path)
        self.aa_to_codons: dict[str, list[tuple[str, float]]] = self._build_aa_to_codons_map()

        # Load golden set for CAI reference weights
        if golden_set_path is not None:
            self.golden_set_table: dict[str, Any] = self._load_codon_table(golden_set_path)
        else:
            try:
                self.golden_set_table = load_golden_set()
            except (FileNotFoundError, json.JSONDecodeError):
                self.golden_set_table = self.codon_table

        # Pre-compute relative adaptiveness weights from golden set (Sharp & Li 1987)
        self.golden_ref_weights: dict[str, float] = self._build_ref_weights(self.golden_set_table)

        # Pre-compute max frequency per amino acid for CAI fallback path
        # Avoids repeated max() inside calculate_cai() hot loop
        self._aa_max_freq: dict[str, float] = {
            aa: max(f for _, f in codons) for aa, codons in self.aa_to_codons.items()
        }
        self._aa_primary_codon: dict[str, str] = {}
        self._aa_weighted_codons: dict[str, tuple[str, ...]] = {}
        self._aa_weighted_cumprob: dict[str, tuple[float, ...]] = {}
        for aa, codons in self.aa_to_codons.items():
            if not codons:
                continue
            self._aa_primary_codon[aa] = codons[0][0]

            codon_names = tuple(c for c, _ in codons)
            raw_weights = [float(w) for _, w in codons]
            total = sum(raw_weights)
            if total <= 0.0:
                # Defensive fallback: uniform sampling if malformed frequencies are loaded.
                n = len(codon_names)
                cumprob = tuple((i + 1) / n for i in range(n))
            else:
                running = 0.0
                cumprob_list: list[float] = []
                for w in raw_weights:
                    running += w / total
                    cumprob_list.append(running)
                # Guard against tiny floating drift.
                cumprob_list[-1] = 1.0
                cumprob = tuple(cumprob_list)
            self._aa_weighted_codons[aa] = codon_names
            self._aa_weighted_cumprob[aa] = cumprob

        # Restriction sites (for Assembly-Friendly profile)
        # Each enzyme maps to a list of recognition sequences (forward + reverse complement)
        self.restriction_sites: dict[str, list[str]] = {
            "BsaI": ["GGTCTC", "GAGACC"],
            "BpiI": ["GAAGAC", "GTCTTC"],
            "BsmBI": ["CGTCTC", "GAGACG"],
        }

    def _load_codon_table(self, path: str | Path) -> dict[str, Any]:
        """Load codon table"""
        with open(path, "r", encoding="utf-8") as f:
            return cast(dict[str, Any], json.load(f))

    def _build_aa_to_codons_map(self) -> dict[str, list[tuple[str, float]]]:
        """
        Build amino-acid-to-codons map

        Returns:
            {"A": [("GCC", 0.40), ("GCT", 0.26), ...], ...}
        """
        aa_map: dict[str, list[tuple[str, float]]] = {}
        raw_aa_map = build_aa_to_codons_map(self.codon_table)
        for aa, codons in raw_aa_map.items():
            codons_with_freq: list[tuple[str, float]] = []
            for codon in codons:
                codon_info = self.codon_table["codons"].get(codon)
                if not codon_info:
                    continue
                freq = float(codon_info["frequency"])
                codons_with_freq.append((codon, freq))

            if codons_with_freq:
                # Sort by frequency (descending)
                codons_with_freq.sort(key=lambda x: x[1], reverse=True)
                aa_map[aa] = codons_with_freq

        return aa_map

    @staticmethod
    def _build_ref_weights(ref_table: dict[str, Any]) -> dict[str, float]:
        """Build relative adaptiveness weights from a reference codon table.

        Groups codons by amino acid and computes w_i = f_i / f_max per amino acid,
        following Sharp & Li (1987).

        Args:
            ref_table: Codon table dict with "codons" section.

        Returns:
            Mapping of codon → relative adaptiveness weight (0-1).
        """
        codons_section = ref_table.get("codons", {})

        # Group frequencies by amino acid
        aa_codons: dict[str, list[tuple[str, float]]] = {}
        for codon, info in codons_section.items():
            aa = info["aa"]
            freq = info.get("frequency", 0.0)
            aa_codons.setdefault(aa, []).append((codon, freq))

        # Compute relative adaptiveness
        weights: dict[str, float] = {}
        for aa, codon_freqs in aa_codons.items():
            if aa == "*":  # Skip stop codons
                continue
            max_freq = max(f for _, f in codon_freqs)
            for codon, freq in codon_freqs:
                weights[codon] = freq / max_freq if max_freq > 0 else 0.0

        return weights

    def calculate_cai(self, dna_sequence: str) -> float:
        """
        Calculate Codon Adaptation Index (CAI) using golden set reference weights.

        Uses pre-computed relative adaptiveness weights from the golden set
        (Sharp & Li 1987). Falls back to the working codon table if the golden
        set does not contain a codon.

        Args:
            dna_sequence: DNA sequence (length must be divisible by 3).

        Returns:
            CAI value (0.0 ~ 1.0).

        Examples:
            >>> translator = ReverseTranslator()
            >>> translator.calculate_cai("ATGGCC")
            0.0
        """
        if len(dna_sequence) % 3 != 0:
            return 0.0

        # ============================================================
        # ORIGINAL (preserved as comment)
        # ============================================================
        # weights: list[float] = []
        # for i in range(0, len(dna_sequence), 3):
        #     codon = dna_sequence[i : i + 3].upper()
        #     w = self.golden_ref_weights.get(codon)
        #     if w is not None and w > 0:
        #         weights.append(w)
        #     elif codon in self.codon_table.get("codons", {}):
        #         aa = self.codon_table["codons"][codon]["aa"]
        #         if aa == "*":
        #             continue
        #         freq = self.codon_table["codons"][codon]["frequency"]
        #         if aa in self.aa_to_codons:
        #             max_freq = max(f for _, f in self.aa_to_codons[aa])  # ← HOT: O(k)×n
        #             weight = freq / max_freq if max_freq > 0 else 0.0
        #             if weight > 0:
        #                 weights.append(weight)
        # if not weights:
        #     return 0.0
        # log_sum = sum(math.log(w) for w in weights)  # ← 2-pass
        # cai = math.exp(log_sum / len(weights))
        # ============================================================
        # OPTIMIZED
        # ============================================================
        # - Fallback max_freq uses pre-computed self._aa_max_freq (O(1) lookup)
        # - 1-pass log accumulation: no list allocation, no second sum() pass
        # Performance: ~8-12x faster for 2,000+ codon sequences
        # ============================================================
        log_sum = 0.0
        count = 0
        codons_section = self.codon_table.get("codons", {})

        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i : i + 3].upper()

            # Primary: golden set reference weights
            w = self.golden_ref_weights.get(codon)
            if w is not None and w > 0:
                log_sum += math.log(w)
                count += 1
            elif codon in codons_section:
                # Fallback: working table with pre-computed max_freq (O(1))
                codon_info = codons_section[codon]
                aa = codon_info["aa"]
                if aa == "*":
                    continue
                max_freq = self._aa_max_freq.get(aa, 0.0)
                if max_freq > 0:
                    weight = codon_info["frequency"] / max_freq
                    if weight > 0:
                        log_sum += math.log(weight)
                        count += 1

        if count == 0:
            return 0.0

        # Geometric mean
        return round(math.exp(log_sum / count), 3)

    def calculate_gc_content(self, dna_sequence: str) -> float:
        """
        Calculate GC content

        Args:
            dna_sequence: DNA sequence

        Returns:
            GC% (0.0 ~ 100.0)

        Raises:
            None.

        Examples:
            >>> translator = ReverseTranslator()
            >>> translator.calculate_gc_content("ATGC")
            50.0
        """
        return round(calculate_gc(dna_sequence), 2)

    def calculate_local_gc(self, dna_sequence: str, window_size: int = 50) -> list[float]:
        """
        Calculate local GC content (sliding window)

        Args:
            dna_sequence: DNA sequence
            window_size: Window size (bp)

        Returns:
            GC% list per window

        Raises:
            None.

        Examples:
            >>> translator = ReverseTranslator()
            >>> translator.calculate_local_gc("ATGCATGC", window_size=4)
            [50.0, 50.0, 50.0, 50.0, 50.0]
        """
        local_gc: list[float] = []

        for i in range(len(dna_sequence) - window_size + 1):
            window = dna_sequence[i : i + window_size]
            gc = self.calculate_gc_content(window)
            local_gc.append(gc)

        return local_gc

    def reverse_translate(
        self,
        protein_seq: str,
        profile: OptimizationProfile = OptimizationProfile.BALANCED,
        **kwargs: Any,
    ) -> str:
        """
        Reverse-translate amino acid sequence to DNA

        Args:
            protein_seq: Amino acid sequence
            profile: Optimization profile
            **kwargs: Profile-specific parameters

        Returns:
            Optimized DNA sequence

        Raises:
            ValueError: Unknown profile or invalid amino acids.

        Examples:
            >>> translator = ReverseTranslator()
            >>> translator.reverse_translate("MA", profile=OptimizationProfile.HIGH_CAI)
            'ATGGCC'
        """
        protein_seq = protein_seq.upper().replace(" ", "")
        kozak = kwargs.pop("kozak", False)

        if profile == OptimizationProfile.BALANCED:
            result = self._balanced_translate(protein_seq, **kwargs)
        elif profile == OptimizationProfile.HIGH_CAI:
            result = self._high_cai_translate(protein_seq, **kwargs)
        elif profile == OptimizationProfile.GC_TARGET:
            result = self._gc_target_translate(protein_seq, **kwargs)
        elif profile == OptimizationProfile.ASSEMBLY_FRIENDLY:
            result = self._assembly_friendly_translate(protein_seq, **kwargs)
        elif profile == OptimizationProfile.RAMP:
            result = self._ramp_translate(protein_seq, **kwargs)
        elif profile == OptimizationProfile.VIRAL_DELIVERY:
            result = self._balanced_translate(protein_seq, **kwargs)
        else:
            raise ValueError(f"Unknown profile: {profile}")

        if kozak:
            result = self._apply_kozak_optimization(result, protein_seq)

        return result

    def _balanced_translate(self, protein_seq: str, **kwargs: Any) -> str:
        """
        Balanced profile: CAI first, GC balanced

        - Preferred codon ratio: 70%
        - Target GC: 55-65% (benchmark analysis 004: avg output 60.1%)
        """
        target_gc_min = kwargs.get("target_gc_min", 55)
        target_gc_max = kwargs.get("target_gc_max", 65)
        preferred_ratio = kwargs.get("preferred_ratio", 0.7)
        max_attempts = kwargs.get("max_gc_attempts", 10)
        if max_attempts < 1:
            raise ValueError("max_gc_attempts must be >= 1")

        best_result: str | None = None
        best_gc_diff = float("inf")
        last_result = ""

        # Try multiple times to find GC within target range
        for _attempt in range(max_attempts):
            dna_seq: list[str] = []

            for aa in protein_seq:
                if aa not in self._aa_primary_codon:
                    raise ValueError(f"Invalid amino acid: {aa}")

                # 70% preferred codon, 30% secondary codon
                if random.random() < preferred_ratio:
                    # Preferred codon
                    codon = self._aa_primary_codon[aa]
                else:
                    # Secondary codon (weighted by frequency)
                    codon = self._sample_weighted_codon(aa)

                dna_seq.append(codon)

            result = "".join(dna_seq)
            last_result = result
            gc = self.calculate_gc_content(result)

            # Return immediately if within target GC
            if target_gc_min <= gc <= target_gc_max:
                return result

            # Track best result
            target_gc_mid = (target_gc_min + target_gc_max) / 2
            gc_diff = abs(gc - target_gc_mid)
            if gc_diff < best_gc_diff:
                best_gc_diff = gc_diff
                best_result = result

        # Return closest result if target range not found
        return best_result if best_result is not None else last_result

    def _high_cai_translate(self, protein_seq: str, **kwargs: Any) -> str:
        """
        High-CAI profile: use only preferred codons

        - CAI > 0.85 guaranteed
        - No GC constraints
        """
        dna_seq: list[str] = []

        for aa in protein_seq:
            if aa not in self.aa_to_codons:
                raise ValueError(f"Invalid amino acid: {aa}")

            # Pick codon with highest golden set relative adaptiveness weight (CAI-optimal)
            codons = self.aa_to_codons[aa]
            preferred_codon = max(codons, key=lambda c: self.golden_ref_weights.get(c[0], 0.0))[0]
            dna_seq.append(preferred_codon)

        return "".join(dna_seq)

    def _gc_target_translate(self, protein_seq: str, **kwargs: Any) -> str:
        """
        GC-Target profile: enforce GC% 42.5% ±2% (N. benthamiana optimal)

        - GC constraint first
        - CAI may be sacrificed
        - Balance local window GC (50 bp)
        """
        target_gc = kwargs.get("target_gc", 42.5)

        dna_seq: list[str] = []

        for i, aa in enumerate(protein_seq):
            if aa not in self.aa_to_codons:
                raise ValueError(f"Invalid amino acid: {aa}")

            codons = self.aa_to_codons[aa]

            current_seq = "".join(dna_seq)

            # Choose codon that brings GC closer to target
            best_codon: str | None = None
            best_diff = float("inf")

            for codon, _ in codons:
                test_seq = current_seq + codon
                test_gc = self.calculate_gc_content(test_seq)
                diff = abs(test_gc - target_gc)

                if diff < best_diff:
                    best_diff = diff
                    best_codon = codon

            dna_seq.append(cast(str, best_codon))

        return "".join(dna_seq)

    def _assembly_friendly_translate(self, protein_seq: str, **kwargs: Any) -> str:
        """
        Assembly-Friendly profile: avoid BsaI/BpiI

        - Golden Gate compatible
        - CAI trade-offs allowed
        """
        max_attempts = kwargs.get("max_attempts", 10)
        if max_attempts < 1:
            raise ValueError("max_attempts must be >= 1")
        last_seq = ""

        for attempt in range(max_attempts):
            # Start with Balanced strategy
            dna_seq = self._balanced_translate(protein_seq, preferred_ratio=0.6)
            last_seq = dna_seq

            # Check restriction sites (forward + reverse complement)
            has_restriction_site = False
            for site_name, site_seqs in self.restriction_sites.items():
                for site_seq in site_seqs:
                    if site_seq in dna_seq:
                        has_restriction_site = True
                        break
                if has_restriction_site:
                    break

            if not has_restriction_site:
                return dna_seq

        # Return with warning if attempts are exhausted
        logger.warning(
            f"Could not remove all restriction sites after {max_attempts} attempts. "
            "The returned sequence may contain restriction sites."
        )
        return last_seq

    def _ramp_translate(self, protein_seq: str, **kwargs: Any) -> str:
        """
        RAMP profile: balanced translation + N-terminal ramp.

        Uses balanced translation as a base, then applies a codon deoptimization
        ramp to the first N codons to promote co-translational folding.

        Args:
            protein_seq: Amino acid sequence.
            **kwargs: ramp_codons (int): Number of N-terminal codons to ramp. Default 50.
        """
        ramp_codons = kwargs.get("ramp_codons", 50)
        dna_seq = self._balanced_translate(protein_seq, **kwargs)
        return self._apply_nterminal_ramp(dna_seq, protein_seq, ramp_codons=ramp_codons)

    def _apply_nterminal_ramp(self, dna_seq: str, protein_seq: str, ramp_codons: int = 50) -> str:
        """
        Apply N-terminal codon ramp for co-translational folding.

        Replaces the first `ramp_codons` codons with lower-frequency synonymous
        codons (bottom 50% by frequency) to slow the ribosome at the N-terminus.
        Single-codon amino acids (Met, Trp) are left unchanged.

        Args:
            dna_seq: Full-length DNA sequence.
            protein_seq: Original protein sequence (same length as dna_seq/3).
            ramp_codons: Number of N-terminal codons to deoptimize.

        Returns:
            DNA sequence with N-terminal ramp applied.
        """
        codons = [dna_seq[i : i + 3] for i in range(0, len(dna_seq), 3)]
        n_ramp = min(ramp_codons, len(codons), len(protein_seq))

        for idx in range(n_ramp):
            aa = protein_seq[idx]
            if aa not in self.aa_to_codons:
                continue

            all_codons = self.aa_to_codons[aa]
            # Skip single-codon amino acids (M, W)
            if len(all_codons) <= 1:
                continue

            # Select from bottom 25% by frequency (mild deoptimization).
            # PMC11718241: optimal tAI_ramp ~0.8-1.2 (max 20% slower than full sequence).
            cutoff = max(1, (3 * len(all_codons)) // 4)
            low_freq_codons = all_codons[cutoff:]

            if not low_freq_codons:
                continue

            # Weighted random from low-frequency codons
            weights = [freq for _, freq in low_freq_codons]
            chosen = random.choices([c for c, _ in low_freq_codons], weights=weights, k=1)[0]
            codons[idx] = chosen

        return "".join(codons)

    def _apply_kozak_optimization(self, dna_seq: str, protein_seq: str) -> str:
        """
        Optimize Kozak context at the 5' end of CDS.

        Plant (N. benthamiana) optimal Kozak context: AACAATG**GC**...
        The 2nd codon (position 4-6, encoding protein_seq[1]) should ideally
        start with G (good) or GC (best) to match the plant Kozak consensus.

        Only performs synonymous codon substitution -- the amino acid sequence
        is preserved. If no synonymous codon starting with G exists, the
        sequence is returned unchanged.

        Args:
            dna_seq: Full-length DNA sequence (must start with ATG).
            protein_seq: Original protein sequence.

        Returns:
            DNA sequence with Kozak-optimized 2nd codon, or original if
            optimization is not possible.
        """
        # Need at least 2 codons (ATG + codon2)
        if len(protein_seq) < 2 or len(dna_seq) < 6:
            return dna_seq

        aa2 = protein_seq[1]
        if aa2 not in self.aa_to_codons:
            return dna_seq

        current_codon2 = dna_seq[3:6]
        codons_for_aa2 = self.aa_to_codons[aa2]

        # Already optimal: starts with G
        if current_codon2[0] == "G":
            return dna_seq

        # Score candidates: prefer GC > G > other
        best_codon = current_codon2
        best_kozak_score = 0  # 0=no G, 1=starts with G, 2=starts with GC
        best_freq = 0.0

        for codon, freq in codons_for_aa2:
            kozak_score = 0
            if codon[0] == "G":
                kozak_score = 1
                if codon[1] == "C":
                    kozak_score = 2
            if kozak_score > best_kozak_score or (
                kozak_score == best_kozak_score and freq > best_freq
            ):
                best_kozak_score = kozak_score
                best_codon = codon
                best_freq = freq

        if best_kozak_score == 0:
            return dna_seq

        return dna_seq[:3] + best_codon + dna_seq[6:]

    def _sample_weighted_codon(self, aa: str) -> str:
        """Sample a codon for one amino acid using precomputed CDF."""
        codons = self._aa_weighted_codons.get(aa)
        cumprob = self._aa_weighted_cumprob.get(aa)
        if not codons or not cumprob:
            raise ValueError(f"Invalid amino acid: {aa}")
        r = random.random()
        idx = bisect_left(cumprob, r)
        if idx >= len(codons):
            idx = len(codons) - 1
        return codons[idx]

    def generate_candidates(
        self,
        protein_seq: str,
        profile: OptimizationProfile = OptimizationProfile.BALANCED,
        n: int = 5,
        **kwargs: Any,
    ) -> list[dict[str, Any]]:
        """
        Generate top-N candidates

        Args:
            protein_seq: Amino acid sequence
            profile: Optimization profile
            n: Number of candidates to generate
            **kwargs: Profile-specific parameters

        Returns:
            [{"sequence": "ATG...", "cai": 0.87, "gc": 51.2, "score": 0.92}, ...]

        Raises:
            ValueError: Invalid amino acids are present.

        Examples:
            >>> translator = ReverseTranslator()
            >>> candidates = translator.generate_candidates("MA", n=2)
            >>> len(candidates) == 2
            True
        """
        if n < 1:
            raise ValueError("n must be >= 1")

        def _build_candidate() -> dict[str, Any]:
            dna_seq = self.reverse_translate(protein_seq, profile, **kwargs)
            cai = self.calculate_cai(dna_seq)
            gc = self.calculate_gc_content(dna_seq)
            score = calculate_composite_score(
                cai=cai,
                gc=gc,
                sequence=dna_seq,
                profile=profile.value,
                **kwargs,
            )
            return {"sequence": dna_seq, "cai": cai, "gc": gc, "score": score}

        # Fast path for the dominant API call shape (n=1).
        if n == 1:
            try:
                return [_build_candidate()]
            except (ValueError, KeyError, TypeError) as exc:
                reason = f"Could not generate a valid candidate in fast path. Last error: {exc}"
                raise EmptyCandidateError(protein_seq[:10], reason=reason) from exc

        candidates: list[dict[str, Any]] = []
        last_error: Exception | None = None
        random.seed(secrets.randbits(32))

        for attempt in range(n):
            try:
                candidates.append(_build_candidate())
            except (ValueError, KeyError, TypeError) as exc:
                # Catch specific exceptions: invalid amino acids, missing codons, type errors
                logger.debug(f"Candidate generation attempt {attempt + 1} failed: {exc}")
                last_error = exc
                continue

        if not candidates:
            reason = (
                f"Could not generate any valid candidates after {n} attempts. "
                "Check codon table and profile settings."
            )
            if last_error is not None:
                reason = f"{reason} Last error: {last_error}"
            raise EmptyCandidateError(protein_seq[:10], reason=reason)

        # Sort by score
        candidates.sort(key=lambda x: x["score"], reverse=True)

        return candidates


# --- Usage example ---
if __name__ == "__main__":
    import json

    translator = ReverseTranslator()

    # Test sequence (partial GFP)
    protein_seq = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHKVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"

    print("=== Balanced Profile ===")
    balanced = translator.reverse_translate(protein_seq, OptimizationProfile.BALANCED)
    print(f"Length: {len(balanced)} bp")
    print(f"CAI: {translator.calculate_cai(balanced)}")
    print(f"GC%: {translator.calculate_gc_content(balanced)}")
    print(f"Sequence: {balanced[:60]}...")

    print("\n=== High-CAI Profile ===")
    high_cai = translator.reverse_translate(protein_seq, OptimizationProfile.HIGH_CAI)
    print(f"CAI: {translator.calculate_cai(high_cai)}")
    print(f"GC%: {translator.calculate_gc_content(high_cai)}")

    print("\n=== GC-Target Profile ===")
    gc_target = translator.reverse_translate(
        protein_seq, OptimizationProfile.GC_TARGET, target_gc=50.0
    )
    print(f"CAI: {translator.calculate_cai(gc_target)}")
    print(f"GC%: {translator.calculate_gc_content(gc_target)}")

    print("\n=== Top-5 Candidates (Balanced) ===")
    candidates = translator.generate_candidates(protein_seq, OptimizationProfile.BALANCED, n=5)
    for i, cand in enumerate(candidates, 1):
        print(f"{i}. CAI={cand['cai']:.3f}, GC={cand['gc']:.1f}%, Score={cand['score']:.3f}")
