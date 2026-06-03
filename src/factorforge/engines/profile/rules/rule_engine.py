"""
Rule Engine for FactorForge profile engine
Plant-aware rule engine - scanning + auto-fix (P0-3)
"""

from __future__ import annotations

import json
import logging
import math
import re
from typing import Any

from factorforge.analysis.metrics import HOMOPOLYMER_SYNTHESIS_WARN_NT
from factorforge.engines.profile.rules.reverse_translator import ReverseTranslator
from factorforge.engines.profile.utils import (
    build_aa_to_codons_map,
    count_dinucleotides,
    get_data_path,
    load_codon_table,
    load_golden_set,
)

logger = logging.getLogger(__name__)


class RuleEngine:
    """
    Plant-aware rule engine

    Features:
    - Detect and remove PolyA signals
    - Detect ARE (AU-rich elements)
    - Detect repeats/homopolymer runs
    - Detect extreme GC regions
    - Detect potential splice sites
    - Auto-fix via synonymous substitutions
    """

    # Pattern definitions
    # PolyA signal patterns: Tier 1 (canonical & high-frequency) + Tier 2 (plant-functional)
    POLYA_PATTERNS = {
        # Tier 1 (canonical & high-frequency variants)
        "AATAAA": "canonical",
        "ATTAAA": "variant_1",
        "AGTAAA": "variant_2",
        # Tier 2 (lower-frequency but functional in plants)
        "AATACA": "variant_3",
        "AAGAAA": "variant_4",
        "AATGAA": "variant_5",
    }
    POLYA_TIER1_PATTERNS = {"AATAAA", "ATTAAA", "AGTAAA"}
    POLYA_TIER2_PATTERNS = {"AATACA", "AAGAAA", "AATGAA"}

    UNSTABLE_MOTIFS = {"ATTTA": "ARE (AU-rich element)", "WWWWWW": "W=A/T, 6+ in a row"}

    def __init__(
        self,
        codon_table: dict[str, Any] | None = None,
        host: str = "nbenthamiana",
    ) -> None:
        """
        Args:
            codon_table: Codon table (loads default if None)
            host: Host codon table name used when codon_table is not provided.
        """
        self.host = host
        if codon_table is None:
            # Use centralized data path management
            codon_table = load_codon_table(host, get_data_path())

        self.codon_table: dict[str, Any] = codon_table
        self.aa_to_codons: dict[str, list[str]] = self._build_aa_to_codons_map()
        self._rare_codon_weights: dict[str, float] = ReverseTranslator._build_ref_weights(
            self.codon_table
        )
        try:
            golden_table = load_golden_set()
            self._golden_cai_weights: dict[str, float] = ReverseTranslator._build_ref_weights(
                golden_table
            )
        except (FileNotFoundError, json.JSONDecodeError):
            self._golden_cai_weights = self._rare_codon_weights

    def _build_aa_to_codons_map(self) -> dict[str, list[str]]:
        """Build amino-acid-to-codons map"""
        return build_aa_to_codons_map(self.codon_table)

    def scan_polya(self, seq: str, window: int = 30) -> list[dict[str, Any]]:
        """
        Detect PolyA signal family

        Args:
            seq: DNA sequence
            window: Window size (bp)

        Returns:
            List of violations

        Raises:
            None.

        Examples:
            >>> engine = RuleEngine()
            >>> engine.scan_polya("AATAAA")
            [{'type': 'polya_signal', 'pattern': 'AATAAA', ...}]
        """
        violations: list[dict[str, Any]] = []
        seq_len = len(seq)
        pattern_hits: dict[str, list[int]] = {}

        # Detect individual patterns
        for pattern, pattern_type in self.POLYA_PATTERNS.items():
            hits: list[int] = []
            pos = 0
            while True:
                idx = seq.find(pattern, pos)
                if idx == -1:
                    break

                hits.append(idx)
                violations.append(
                    {
                        "type": "polya_signal",
                        "pattern": pattern,
                        "pattern_type": pattern_type,
                        "position": idx,
                        "context": seq[max(0, idx - 10) : min(len(seq), idx + len(pattern) + 10)],
                    }
                )
                pos = idx + 1
            pattern_hits[pattern] = hits

        if window < 1 or seq_len < window:
            return violations

        # Precompute per-pattern prefix arrays for fast "pattern exists in window".
        # Semantics match `pattern in window_seq`: count each pattern at most once/window.
        pattern_prefix: dict[str, tuple[int, list[int]]] = {}
        for pattern, hits in pattern_hits.items():
            plen = len(pattern)
            if plen > window or not hits:
                continue
            prefix = [0] * (seq_len + 1)
            for idx in hits:
                prefix[idx + 1] = 1
            for i in range(1, seq_len + 1):
                prefix[i] += prefix[i - 1]
            pattern_prefix[pattern] = (plen, prefix)

        # Add warning if 2+ patterns in 30 bp window
        for i in range(seq_len - window + 1):
            count = 0
            for _pattern, (plen, prefix) in pattern_prefix.items():
                max_start = i + window - plen
                if max_start >= i and (prefix[max_start + 1] - prefix[i]) > 0:
                    count += 1

            if count >= 2:
                window_seq = seq[i : i + window]
                violations.append(
                    {
                        "type": "multiple_polya",
                        "position": i,
                        "window_size": window,
                        "count": count,
                        "context": window_seq,
                        "severity": "high",
                    }
                )

        return violations

    def scan_are(self, seq: str) -> list[dict[str, Any]]:
        """
        Detect ARE (AU-rich element) pattern

        Args:
            seq: DNA sequence

        Returns:
            List of violations

        Raises:
            None.

        Examples:
            >>> engine = RuleEngine()
            >>> engine.scan_are("ATTTA")
            [{'type': 'are_element', ...}]
        """
        violations: list[dict[str, Any]] = []

        # ATTTA pattern
        pos = 0
        while True:
            idx = seq.find("ATTTA", pos)
            if idx == -1:
                break

            violations.append(
                {
                    "type": "are_element",
                    "pattern": "ATTTA",
                    "position": idx,
                    "context": seq[max(0, idx - 10) : min(len(seq), idx + 15)],
                    "severity": "medium",
                }
            )
            pos = idx + 1

        return violations

    def scan_at_runs(self, seq: str, min_length: int = 6) -> list[dict[str, Any]]:
        """
        Detect A/T runs

        Args:
            seq: DNA sequence
            min_length: Minimum length

        Returns:
            List of violations

        Raises:
            None.

        Examples:
            >>> engine = RuleEngine()
            >>> engine.scan_at_runs("AAAAAATTT", min_length=6)
            [{'type': 'at_run', ...}]
        """
        violations: list[dict[str, Any]] = []
        pattern = r"[AT]{" + str(min_length) + r",}"

        for match in re.finditer(pattern, seq):
            violations.append(
                {
                    "type": "at_run",
                    "position": match.start(),
                    "length": len(match.group()),
                    "sequence": match.group(),
                    "context": seq[max(0, match.start() - 5) : min(len(seq), match.end() + 5)],
                    "severity": "medium" if len(match.group()) < 8 else "high",
                }
            )

        return violations

    def scan_homopolymers(
        self, seq: str, min_length: int = HOMOPOLYMER_SYNTHESIS_WARN_NT
    ) -> list[dict[str, Any]]:
        """Detect homopolymer runs for synthesis/manufacturing risk evaluation.

        Uses HOMOPOLYMER_SYNTHESIS_WARN_NT (default 8 nt) — the threshold at
        which gene synthesis vendors flag homopolymers as difficult to synthesize
        with high fidelity.

        For expression stability risk (≥6 nt), see
        factorforge.analysis.metrics.detect_homopolymers() which uses
        HOMOPOLYMER_EXPRESSION_WARN_NT.

        Args:
            seq: DNA sequence
            min_length: Minimum run length to flag (default: HOMOPOLYMER_SYNTHESIS_WARN_NT)

        Examples:
            >>> engine = RuleEngine()
            >>> engine.scan_homopolymers("AAAAAAAA", min_length=8)
            [{'type': 'homopolymer', 'context': 'synthesis', ...}]
        """
        violations: list[dict[str, Any]] = []

        for base in "ATGC":
            pattern = base * min_length
            pos = 0
            while True:
                idx = seq.find(pattern, pos)
                if idx == -1:
                    break

                # Compute actual run length
                actual_length = min_length
                while idx + actual_length < len(seq) and seq[idx + actual_length] == base:
                    actual_length += 1

                violations.append({
                    "type": "homopolymer",
                    "context": "synthesis",
                    "threshold_nt": min_length,
                    "base": base,
                    "position": idx,
                    "length": actual_length,
                    "sequence": base * actual_length,
                    "severity": "high" if actual_length >= 10 else "medium",
                })
                pos = idx + actual_length

        return violations

    def scan_repeats(self, seq: str, min_length: int = 15) -> list[dict[str, Any]]:
        """
        Detect perfect repeats >= 15 bp (recombination risk)

        Args:
            seq: DNA sequence
            min_length: Minimum repeat length

        Returns:
            List of violations

        Raises:
            None.

        Examples:
            >>> engine = RuleEngine()
            >>> engine.scan_repeats("ATGATGATGATGATG", min_length=3)
            [{'type': 'repeat', ...}]
        """
        violations: list[dict[str, Any]] = []
        seen_fragments: dict[str, list[int]] = {}

        for i in range(len(seq) - min_length + 1):
            fragment = seq[i : i + min_length]

            if fragment in seen_fragments:
                # Already found repeat
                seen_fragments[fragment].append(i)
            else:
                # First occurrence
                seen_fragments[fragment] = [i]

        # Report only fragments that appear 2+ times
        for fragment, positions in seen_fragments.items():
            if len(positions) > 1:
                violations.append(
                    {
                        "type": "repeat",
                        "fragment": fragment,
                        "length": len(fragment),
                        "positions": positions,
                        "count": len(positions),
                        "severity": "high" if len(positions) > 2 else "medium",
                    }
                )

        return violations

    def scan_gc_extremes(
        self,
        seq: str,
        window: int = 50,
        min_gc: float = 25,
        max_gc: float = 75,
    ) -> list[dict[str, Any]]:
        """
        Detect extreme GC regions in a sliding local window.

        This is a LOCAL synthesis/extreme-window guard (default 25-75% over a
        50 bp window), NOT the global GC target. Global GC is governed separately
        by the scoring band (GC_OPT_MIN/MAX, ~55-65%) and the API/DP gc_min/gc_max
        constraints. The wide 25-75% band intentionally flags only synthesis-hostile
        local windows; narrowing it toward the global optimum would raise false
        positives against the engine's own output distribution (analysis 004: 55-71%).

        Args:
            seq: DNA sequence
            window: Window size (bp)
            min_gc: Minimum local GC% threshold (synthesis guard, not global target)
            max_gc: Maximum local GC% threshold (synthesis guard, not global target)

        Returns:
            List of violations

        Raises:
            None.

        Examples:
            >>> engine = RuleEngine()
            >>> engine.scan_gc_extremes("GGGGGG", window=3, max_gc=80)
            [{'type': 'gc_extreme', ...}]
        """
        violations: list[dict[str, Any]] = []
        seq_len = len(seq)
        if window < 1 or seq_len < window:
            return violations

        seq_upper = seq.upper()
        gc_count = sum(1 for b in seq_upper[:window] if b == "G" or b == "C")
        last_start = seq_len - window

        for i in range(last_start + 1):
            if i > 0:
                left = seq_upper[i - 1]
                right = seq_upper[i + window - 1]
                if left == "G" or left == "C":
                    gc_count -= 1
                if right == "G" or right == "C":
                    gc_count += 1

            gc = (gc_count / window) * 100.0

            if gc < min_gc or gc > max_gc:
                severity = "high" if gc < 20 or gc > 80 else "medium"
                window_seq = seq[i : i + window]

                violations.append(
                    {
                        "type": "gc_extreme",
                        "position": i,
                        "window_size": window,
                        "gc": round(gc, 1),
                        "context": window_seq,
                        "severity": severity,
                    }
                )

        return violations

    def scan_splice_sites(self, seq: str) -> list[dict[str, Any]]:
        """
        Detect potential splice-site-like patterns

        Scan GT...AG pattern (20-200 bp spacing)
        Plant consensus: GTRAG...YAG (Y=C/T)

        Args:
            seq: DNA sequence

        Returns:
            List of violations

        Raises:
            None.

        Examples:
            >>> engine = RuleEngine()
            >>> engine.scan_splice_sites("GTAG" + "A" * 20 + "CAG")
            [{'type': 'potential_splice_site', ...}]
        """
        violations: list[dict[str, Any]] = []

        # Donor site: GT[AG]AG
        donor_pattern = r"GT[AG]AG"
        # Acceptor site: [CT]AG
        acceptor_pattern = r"[CT]AG"

        donors = [(m.start(), m.group()) for m in re.finditer(donor_pattern, seq)]
        acceptors = [(m.start(), m.group()) for m in re.finditer(acceptor_pattern, seq)]

        # Check 20-200 bp spacing
        for d_pos, d_seq in donors:
            for a_pos, a_seq in acceptors:
                distance = a_pos - d_pos

                if 20 <= distance <= 200:
                    violations.append(
                        {
                            "type": "potential_splice_site",
                            "donor": {"pos": d_pos, "seq": d_seq},
                            "acceptor": {"pos": a_pos, "seq": a_seq},
                            "distance": distance,
                            "severity": "low",
                            "warning": "Potential cryptic splice site",
                        }
                    )

        return violations

    def scan_polya_positive(
        self,
        seq: str,
        required_patterns: set[str] | None = None,
    ) -> list[dict[str, Any]]:
        """
        Positive PolyA validation: check that a region CONTAINS a PolyA signal.

        Used for terminator/3'UTR regions where a PolyA signal must be present
        for proper mRNA polyadenylation.

        Args:
            seq: DNA sequence of the terminator/3'UTR region.
            required_patterns: Set of acceptable PolyA patterns.
                               Defaults to Tier 1 patterns.

        Returns:
            List of violations (non-empty if PolyA signal is missing).
        """
        if required_patterns is None:
            required_patterns = self.POLYA_TIER1_PATTERNS

        for pattern in required_patterns:
            if pattern in seq:
                return []  # At least one PolyA signal found

        return [
            {
                "type": "missing_polya_signal",
                "severity": "high",
                "message": "No PolyA signal found in terminator/3'UTR region.",
                "checked_patterns": sorted(required_patterns),
            }
        ]

    def fix_polya_iterative(
        self,
        seq: str,
        max_rounds: int = 10,
    ) -> dict[str, Any]:
        """
        Iteratively remove all PolyA signals from a CDS via synonymous substitutions.

        Fixing one PolyA violation can create another at a different codon boundary,
        so this method loops until no violations remain or max rounds are reached.

        Args:
            seq: DNA coding sequence (must be divisible by 3).
            max_rounds: Maximum number of fix-scan cycles.

        Returns:
            Dict with success, modified_seq, rounds, and fixes_applied.
        """
        current_seq = seq
        all_fixes: list[dict[str, Any]] = []

        for round_num in range(1, max_rounds + 1):
            violations = self.scan_polya(current_seq)
            # Filter to only polya_signal type (not multiple_polya warnings)
            signal_violations = [v for v in violations if v["type"] == "polya_signal"]

            if not signal_violations:
                return {
                    "success": True,
                    "modified_seq": current_seq,
                    "rounds": round_num - 1,
                    "fixes_applied": all_fixes,
                }

            # Try to fix the first violation
            fix_result = self.fix_violation(current_seq, signal_violations[0])
            if fix_result["success"]:
                current_seq = fix_result["modified_seq"]
                all_fixes.extend(fix_result.get("changes", []))
            else:
                logger.debug(
                    f"Could not fix PolyA at position {signal_violations[0]['position']} "
                    f"in round {round_num}"
                )
                return {
                    "success": False,
                    "modified_seq": current_seq,
                    "rounds": round_num,
                    "fixes_applied": all_fixes,
                    "remaining_violations": len(signal_violations),
                }

        # Max rounds exhausted
        remaining = [v for v in self.scan_polya(current_seq) if v["type"] == "polya_signal"]
        return {
            "success": len(remaining) == 0,
            "modified_seq": current_seq,
            "rounds": max_rounds,
            "fixes_applied": all_fixes,
            "remaining_violations": len(remaining),
        }

    def scan_dinucleotides(
        self,
        seq: str,
        window: int = 50,
        cpg_threshold: float = 0.05,
        tpa_threshold: float = 0.05,
    ) -> list[dict[str, Any]]:
        """
        Detect CpG and TpA dinucleotide-dense regions in CDS.

        CpG dinucleotides trigger methylation-based gene silencing in plants.
        TpA (UpA in RNA) dinucleotides are associated with mRNA instability.

        Args:
            seq: DNA sequence.
            window: Sliding window size (bp).
            cpg_threshold: CpG density (count/window) above which a violation
                is reported. Default 0.05 = >1 CpG per 20 bp.
            tpa_threshold: TpA density threshold (same units).

        Returns:
            List of violation dicts with type, dinucleotide, position,
            density, and severity.

        Examples:
            >>> engine = RuleEngine()
            >>> engine.scan_dinucleotides("ACGACGACG" * 10)  # doctest: +SKIP
            [{'type': 'dinucleotide_hotspot', ...}]
        """
        violations: list[dict[str, Any]] = []
        seq_upper = seq.upper()
        seq_len = len(seq_upper)

        if window < 2 or seq_len < 2:
            return violations

        if seq_len < window:
            # Scan the full sequence as a single window
            cpg_count = count_dinucleotides(seq_upper, "CG")
            tpa_count = count_dinucleotides(seq_upper, "TA")
            cpg_density = cpg_count / seq_len
            if cpg_density > cpg_threshold:
                violations.append(
                    {
                        "type": "dinucleotide_hotspot",
                        "dinucleotide": "CpG",
                        "position": 0,
                        "window_size": seq_len,
                        "count": cpg_count,
                        "density": round(cpg_density, 4),
                        "severity": "high" if cpg_density > cpg_threshold * 2 else "medium",
                    }
                )
            tpa_density = tpa_count / seq_len
            if tpa_density > tpa_threshold:
                violations.append(
                    {
                        "type": "dinucleotide_hotspot",
                        "dinucleotide": "TpA",
                        "position": 0,
                        "window_size": seq_len,
                        "count": tpa_count,
                        "density": round(tpa_density, 4),
                        "severity": "high" if tpa_density > tpa_threshold * 2 else "medium",
                    }
                )
            return violations

        # Rolling dinucleotide counts:
        # a window of `window` bases contains `window - 1` adjacent pairs.
        pair_len = seq_len - 1
        cpg_flags = [1 if seq_upper[i : i + 2] == "CG" else 0 for i in range(pair_len)]
        tpa_flags = [1 if seq_upper[i : i + 2] == "TA" else 0 for i in range(pair_len)]
        pairs_in_window = window - 1
        cpg_count = sum(cpg_flags[:pairs_in_window])
        tpa_count = sum(tpa_flags[:pairs_in_window])
        last_start = seq_len - window

        for i in range(last_start + 1):
            if i > 0:
                # Shift by one base: remove pair at i-1, add pair at i+window-2
                add_idx = i + window - 2
                cpg_count += cpg_flags[add_idx] - cpg_flags[i - 1]
                tpa_count += tpa_flags[add_idx] - tpa_flags[i - 1]

            cpg_density = cpg_count / window
            if cpg_density > cpg_threshold:
                violations.append(
                    {
                        "type": "dinucleotide_hotspot",
                        "dinucleotide": "CpG",
                        "position": i,
                        "window_size": window,
                        "count": cpg_count,
                        "density": round(cpg_density, 4),
                        "severity": "high" if cpg_density > cpg_threshold * 2 else "medium",
                    }
                )

            tpa_density = tpa_count / window
            if tpa_density > tpa_threshold:
                violations.append(
                    {
                        "type": "dinucleotide_hotspot",
                        "dinucleotide": "TpA",
                        "position": i,
                        "window_size": window,
                        "count": tpa_count,
                        "density": round(tpa_density, 4),
                        "severity": "high" if tpa_density > tpa_threshold * 2 else "medium",
                    }
                )

        return violations

    def _calc_cai(self, seq: str) -> float:
        """Calculate CAI using golden set reference weights (Sharp & Li 1987)."""
        _stop = {"TAA", "TAG", "TGA"}
        ws = [
            self._golden_cai_weights.get(seq[i : i + 3], 0.001)
            for i in range(0, len(seq) - 2, 3)
            if seq[i : i + 3] not in _stop
        ]
        if not ws:
            return 0.0
        return math.exp(sum(math.log(w) for w in ws) / len(ws))

    def fix_dinucleotides(
        self,
        seq: str,
        max_rounds: int = 5,
        target_dinucleotides: tuple[str, ...] = ("CG", "TA"),
        mode: str = "balanced",
        cai_floor: float = 0.75,
        max_cai_drop: float | None = None,
    ) -> dict[str, Any]:
        """
        Reduce CpG and TpA dinucleotide density via greedy synonymous substitution.

        Modes:
            aggressive: dinucleotide reduction only; no CAI check (Job 044 behaviour).
            balanced: dinucleotide reduction first; rollback each pass if final CAI
                      drops below cai_floor.
            cai_preserving: rollback each pass if CAI drops more than max_cai_drop
                            below initial CAI.

        Args:
            seq: DNA coding sequence (must be divisible by 3).
            max_rounds: Maximum number of full-sequence passes.
            target_dinucleotides: Dinucleotides to reduce.
            mode: "aggressive" | "balanced" | "cai_preserving".
            cai_floor: Minimum allowed final CAI (balanced mode).
            max_cai_drop: Maximum allowed CAI decrease from initial
                (cai_preserving mode).

        Returns:
            Dict with modified_seq, success, rounds, initial/final counts,
            reduction_pct, mode, cai_before, cai_after.
        """
        seq_upper = seq.upper()
        initial_count = sum(count_dinucleotides(seq_upper, di) for di in target_dinucleotides)
        initial_cai = self._calc_cai(seq_upper)

        if len(seq) % 3 != 0 or len(seq) == 0:
            return {
                "modified_seq": seq,
                "success": False,
                "rounds": 0,
                "initial_count": initial_count,
                "final_count": initial_count,
                "reduction_pct": 0.0,
                "mode": mode,
                "cai_before": round(initial_cai, 4),
                "cai_after": round(initial_cai, 4),
            }

        if mode == "aggressive":
            effective_floor = 0.0
        elif mode == "cai_preserving":
            drop = max_cai_drop if max_cai_drop is not None else 0.003
            effective_floor = initial_cai - drop
        else:
            effective_floor = cai_floor

        targets_set = set(target_dinucleotides)
        current_seq = seq_upper
        round_num = 0

        def _local_count(s: str, codon_start: int) -> int:
            n = len(s) - 1
            total = 0
            for di in range(max(0, codon_start - 1), min(n, codon_start + 3)):
                if s[di : di + 2] in targets_set:
                    total += 1
            return total

        for round_num in range(1, max_rounds + 1):
            improved = False
            pass_start_seq = current_seq

            for codon_idx in range(len(current_seq) // 3):
                codon_start = codon_idx * 3
                original_codon = current_seq[codon_start : codon_start + 3]

                if original_codon not in self.codon_table["codons"]:
                    continue

                aa = self.codon_table["codons"][original_codon]["aa"]
                synonyms = [c for c in self.aa_to_codons.get(aa, []) if c != original_codon]

                if not synonyms:
                    continue

                current_local = _local_count(current_seq, codon_start)
                if current_local == 0:
                    continue

                for alt_codon in synonyms:
                    candidate = (
                        current_seq[:codon_start] + alt_codon + current_seq[codon_start + 3 :]
                    )
                    if _local_count(candidate, codon_start) < current_local:
                        current_seq = candidate
                        improved = True
                        break

            if mode != "aggressive" and improved:
                if self._calc_cai(current_seq) < effective_floor:
                    current_seq = pass_start_seq
                    break

            if not improved:
                break

        final_count = sum(count_dinucleotides(current_seq, di) for di in target_dinucleotides)
        initial_nonzero = initial_count if initial_count > 0 else 1
        reduction_pct = round((initial_count - final_count) / initial_nonzero * 100.0, 1)

        return {
            "modified_seq": current_seq,
            "success": final_count < initial_count,
            "rounds": round_num,
            "initial_count": initial_count,
            "final_count": final_count,
            "reduction_pct": reduction_pct,
            "mode": mode,
            "cai_before": round(initial_cai, 4),
            "cai_after": round(self._calc_cai(current_seq), 4),
        }

    def scan_rare_codon_runs(
        self,
        seq: str,
        min_run: int = 3,
        rarity_threshold: float = 0.3,
    ) -> list[dict[str, Any]]:
        """
        Detect runs of consecutive rare codons (ribosome stalling risk).

        Rare codon: relative adaptiveness w < rarity_threshold.
        A run of >= min_run consecutive rare codons is flagged.
        Reference: Tuller et al. (2010) PMC2565806; Dana & Tuller (2014) PMC4363877.

        Args:
            seq: DNA coding sequence (must be divisible by 3).
            min_run: Minimum consecutive rare codons to flag.
            rarity_threshold: CAI relative adaptiveness w below which a codon is rare.

        Returns:
            List of violation dicts with type, position, codon_index, run_length,
            codons, and severity.
        """
        if len(seq) % 3 != 0 or min_run < 1:
            return []

        violations: list[dict[str, Any]] = []
        seq_upper = seq.upper()
        stop_codons = {"TAA", "TAG", "TGA"}
        current_run_start: int | None = None
        current_run_codons: list[str] = []

        def flush_run() -> None:
            if current_run_start is None or len(current_run_codons) < min_run:
                return
            violations.append(
                {
                    "type": "rare_codon_run",
                    "position": current_run_start * 3,
                    "codon_index": current_run_start,
                    "run_length": len(current_run_codons),
                    "codons": current_run_codons.copy(),
                    "severity": "high" if len(current_run_codons) >= 5 else "medium",
                }
            )

        for codon_index in range(0, len(seq_upper), 3):
            codon = seq_upper[codon_index : codon_index + 3]
            is_rare = codon not in stop_codons and (
                self._rare_codon_weights.get(codon, 0.0) < rarity_threshold
            )

            if is_rare:
                if current_run_start is None:
                    current_run_start = codon_index // 3
                current_run_codons.append(codon)
                continue

            flush_run()
            current_run_start = None
            current_run_codons = []

        flush_run()
        return violations

    def scan_all(
        self,
        seq: str,
        mode: str = "full",
        include: list[str] | None = None,
        exclude: list[str] | None = None,
    ) -> dict[str, list[dict[str, Any]]]:
        """
        Scan all rules

        Args:
            seq: DNA sequence
            mode: Scan mode. "full" runs all scanners; "fast" skips heavier scanners.
            include: Explicit scanner names to run. Overrides mode.
            exclude: Scanner names to exclude from the selected set.

        Returns:
            {
                "polya": [...],
                "are": [...],
                "at_runs": [...],
                "homopolymers": [...],
                "repeats": [...],
                "gc_extremes": [...],
                "splice_sites": [...],
                "dinucleotides": [...],
                "rare_codon_runs": [...]
            }

        Raises:
            None.

        Examples:
            >>> engine = RuleEngine()
            >>> result = engine.scan_all("ATG" * 10)
            >>> "polya" in result
            True
        """
        scanner_map = {
            "polya": self.scan_polya,
            "are": self.scan_are,
            "at_runs": self.scan_at_runs,
            "homopolymers": self.scan_homopolymers,
            "repeats": self.scan_repeats,
            "gc_extremes": self.scan_gc_extremes,
            "splice_sites": self.scan_splice_sites,
            "dinucleotides": self.scan_dinucleotides,
            "rare_codon_runs": self.scan_rare_codon_runs,
        }

        mode_name = mode.lower().strip()
        if include is not None:
            selected_names = [name.strip() for name in include if name.strip()]
        elif mode_name == "full":
            selected_names = list(scanner_map.keys())
        elif mode_name == "fast":
            selected_names = [
                "polya",
                "are",
                "at_runs",
                "homopolymers",
                "gc_extremes",
                "splice_sites",
            ]
        else:
            raise ValueError(f"Unknown scan mode: {mode}. Supported: full, fast")

        if exclude:
            excluded = {name.strip() for name in exclude if name.strip()}
            selected_names = [name for name in selected_names if name not in excluded]

        unknown = sorted({name for name in selected_names if name not in scanner_map})
        if unknown:
            known = ", ".join(scanner_map.keys())
            raise ValueError(f"Unknown scanners: {', '.join(unknown)}. Known scanners: {known}")

        return {name: scanner_map[name](seq) for name in selected_names}  # type: ignore[operator]

    def fix_violation(self, seq: str, violation: dict[str, Any]) -> dict[str, Any]:
        """
        Fix violations via synonymous substitutions

        Args:
            seq: DNA sequence
            violation: Violation entry

        Returns:
            {
                "success": True/False,
                "modified_seq": "...",
                "changes": [{...}],
                "aa_preserved": True/False
            }

        Raises:
            None.

        Examples:
            >>> engine = RuleEngine()
            >>> v = {"type": "polya_signal", "pattern": "AATAAA", "position": 0}
            >>> engine.fix_violation("AATAAA", v)["success"]
            False
        """
        if len(seq) % 3 != 0:
            return {
                "success": False,
                "error": "Sequence length not divisible by 3",
                "aa_preserved": False,
            }

        pos = violation["position"]
        pattern_type = violation["type"]

        # Compute codon range overlapping violation pattern
        if pattern_type == "polya_signal":
            pattern_len = len(violation["pattern"])
        elif pattern_type == "are_element":
            pattern_len = 5  # ATTTA
        elif pattern_type == "at_run":
            pattern_len = violation["length"]
        elif pattern_type == "homopolymer":
            pattern_len = violation["length"]
        else:
            # Other types not supported yet
            return {
                "success": False,
                "error": f"Unsupported violation type: {pattern_type}",
                "aa_preserved": False,
            }

        first_codon_idx = (pos // 3) * 3
        last_codon_idx = ((pos + pattern_len - 1) // 3) * 3

        # Try synonymous substitutions per codon
        modified_seq = list(seq)
        changes: list[dict[str, Any]] = []

        for codon_start in range(first_codon_idx, last_codon_idx + 1, 3):
            if codon_start + 3 > len(seq):
                continue

            original_codon = seq[codon_start : codon_start + 3]

            # Validate amino acid
            if original_codon not in self.codon_table["codons"]:
                continue

            aa = self.codon_table["codons"][original_codon]["aa"]

            # Find synonymous codons
            synonymous_codons = [c for c in self.aa_to_codons.get(aa, []) if c != original_codon]

            if not synonymous_codons:
                continue

            # Try each synonymous codon
            for alt_codon in synonymous_codons:
                # Temporary substitution
                test_seq = modified_seq[:]
                test_seq[codon_start : codon_start + 3] = list(alt_codon)
                test_seq_str = "".join(test_seq)

                # Check if violation pattern is removed
                test_region = test_seq_str[
                    max(0, pos - 10) : min(len(test_seq_str), pos + pattern_len + 10)
                ]

                pattern_removed = False
                if pattern_type == "polya_signal":
                    pattern_removed = violation["pattern"] not in test_region
                elif pattern_type == "are_element":
                    pattern_removed = "ATTTA" not in test_region
                elif pattern_type == "at_run":
                    pattern_removed = not re.search(r"[AT]{6,}", test_region)
                elif pattern_type == "homopolymer":
                    base = violation["base"]
                    pattern_removed = (base * 8) not in test_region

                if pattern_removed:
                    # Success
                    modified_seq = test_seq
                    changes.append(
                        {
                            "pos": codon_start,
                            "original": original_codon,
                            "fixed": alt_codon,
                            "aa": aa,
                        }
                    )

                    return {
                        "success": True,
                        "modified_seq": "".join(modified_seq),
                        "changes": changes,
                        "aa_preserved": True,
                    }

        # Failed to fix
        return {
            "success": False,
            "modified_seq": seq,
            "changes": [],
            "aa_preserved": True,
            "reason": "No synonymous codon available to remove violation",
        }


# --- Usage example ---
if __name__ == "__main__":
    engine = RuleEngine()

    # Test sequence (partial GFP)
    test_seq = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"

    print("=== Scanning for violations ===")
    results = engine.scan_all(test_seq)

    for rule_type, violations in results.items():
        if violations:
            print(f"\n{rule_type.upper()}: {len(violations)} violations")
            for v in violations[:3]:  # Show only the first 3
                print(f"  - Position {v.get('position', 'N/A')}: {v.get('type', 'N/A')}")
