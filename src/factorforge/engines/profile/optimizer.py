"""Profile-based rule optimizer implementation."""

from __future__ import annotations

from typing import Any

from factorforge.core.interfaces import OptimizationResult, OptimizerEngine

from .exporter import SequenceExporter
from .rules.reverse_translator import OptimizationProfile, ReverseTranslator
from .rules.rule_engine import RuleEngine
from .scoring import calculate_composite_score, compute_mfe_evidence
from .validator import InputValidator


class RuleBasedOptimizer(OptimizerEngine):
    """Profile-based rule optimization engine."""

    name = "Profile-based"
    version = "3.1.9"

    def __init__(self) -> None:
        self.validator = InputValidator()
        self.translator = ReverseTranslator()  # Data files use default path
        self.rule_engine = RuleEngine()
        self.exporter = SequenceExporter()

    def optimize(
        self,
        sequence: str,
        profile: str | None = "balanced",
        host: str = "nbenthamiana",
        **kwargs: Any,
    ) -> OptimizationResult:
        """
        Rule-based optimization

        Args:
            sequence: Protein sequence or DNA sequence
            profile: Optimization profile
            host: Host codon table name.
            **kwargs: Additional settings

        Returns:
            OptimizationResult object

        Raises:
            ValueError: If the input sequence is invalid.

        Examples:
            >>> optimizer = RuleBasedOptimizer()
            >>> result = optimizer.optimize("MA", profile="balanced")
            >>> len(result.sequence) == 6
            True
        """
        # 1. Validate input
        val_result = self.validator.validate(sequence)
        if not val_result["valid"]:
            raise ValueError(f"Invalid input sequence: {val_result['errors']}")

        processed_seq = val_result["processed_sequence"]
        seq_type = val_result["type"]
        if seq_type == "fasta":
            seq_type = self.validator.detect_sequence_type(processed_seq).value

        # 2. Normalize profile
        profile_value = (profile or "balanced").lower()
        try:
            opt_profile = OptimizationProfile(profile_value)
        except ValueError as exc:
            supported = ", ".join(p.value for p in OptimizationProfile)
            raise ValueError(
                f"Unknown profile: {profile_value}. Supported profiles: {supported}"
            ) from exc

        if host == "nbenthamiana":
            translator = self.translator
            rule_engine = self.rule_engine
        else:
            translator = ReverseTranslator(host=host)
            rule_engine = RuleEngine(host=host)

        # 3. Reverse-translate (pick the best candidate)
        if seq_type == "dna":
            optimized_dna = processed_seq
            cai = translator.calculate_cai(optimized_dna)
            gc = translator.calculate_gc_content(optimized_dna)
            score = calculate_composite_score(
                cai=cai, gc=gc, sequence=optimized_dna, profile=profile_value
            )
            candidates = [{"sequence": optimized_dna, "cai": cai, "gc": gc, "score": score}]
        else:
            candidates = translator.generate_candidates(
                processed_seq, profile=opt_profile, n=1
            )
            if not candidates:
                raise ValueError("No candidates generated for input sequence.")
            optimized_dna = candidates[0]["sequence"]

        # 4. Rule checks (PolyA, etc.)
        scan_mode = str(kwargs.get("scan_mode", "full"))
        scan_include = kwargs.get("scan_include")
        scan_exclude = kwargs.get("scan_exclude")
        scan_results = rule_engine.scan_all(
            optimized_dna,
            mode=scan_mode,
            include=scan_include,
            exclude=scan_exclude,
        )

        # 5. Build result
        metrics = {
            "cai": candidates[0]["cai"],
            # Keep both names for compatibility across existing tests/callers.
            "gc_content": candidates[0]["gc"],
            "gc_percent": candidates[0]["gc"],
            "score": candidates[0]["score"],
            "violations": sum(len(v) for v in scan_results.values()),
        }
        # MFE provenance: expose whether MFE was actually computed so downstream
        # artifacts (API response, Design Package) never report an uncomputed
        # MFE as a misleading 0.0 (016 audit). Score value is unchanged.
        metrics.update(compute_mfe_evidence(optimized_dna, profile=profile_value))

        return OptimizationResult(
            sequence=optimized_dna,
            metrics=metrics,
            metadata={
                "engine": "profile",
                "profile": profile_value,
                "host": host,
                "scan_mode": scan_mode,
                "scan_results": scan_results,
            },
        )

    def optimize_batch(
        self,
        sequences: list[dict[str, str]] | list[str],
        profile: str | None = "balanced",
        **kwargs: Any,
    ) -> list[OptimizationResult]:
        """Optimize a batch of sequences.

        Args:
            sequences: Either list[str] or list[{"id": str, "sequence": str}].
            profile: Optimization profile.
            **kwargs: Additional optimize options.

        Returns:
            List of OptimizationResult entries in input order.
        """
        results: list[OptimizationResult] = []

        for idx, entry in enumerate(sequences, start=1):
            if isinstance(entry, dict):
                seq = entry.get("sequence", "")
                seq_id = entry.get("id", f"seq{idx}")
            else:
                seq = entry
                seq_id = f"seq{idx}"

            result = self.optimize(seq, profile=profile, **kwargs)
            result.metadata["input_id"] = seq_id
            results.append(result)

        return results

    def validate(self, sequence: str) -> bool:
        """
        Validate input

        Args:
            sequence: Input sequence

        Returns:
            Validity flag

        Raises:
            None.

        Examples:
            >>> optimizer = RuleBasedOptimizer()
            >>> optimizer.validate("MA")
            True
        """
        return bool(self.validator.validate(sequence)["valid"])

    def get_supported_profiles(self) -> list[str]:
        """
        Return list of supported profiles

        Returns:
            List of profile strings

        Raises:
            None.

        Examples:
            >>> optimizer = RuleBasedOptimizer()
            >>> "balanced" in optimizer.get_supported_profiles()
            True
        """
        return [p.value for p in OptimizationProfile]
