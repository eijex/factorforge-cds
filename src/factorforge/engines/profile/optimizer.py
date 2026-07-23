"""Profile-based rule optimizer implementation."""

from __future__ import annotations

import logging
from typing import Any

from factorforge.core.interfaces import OptimizationResult, OptimizerEngine
from factorforge.analysis.metrics import translate_dna

from .exporter import SequenceExporter
from .rules.reverse_translator import OptimizationProfile, ReverseTranslator
from .rules.rule_engine import RuleEngine
from .scoring import calculate_composite_score, compute_mfe_evidence, resolve_host_gc_range
from .validator import InputValidator

logger = logging.getLogger(__name__)


class RuleBasedOptimizer(OptimizerEngine):
    """Profile-based rule optimization engine."""

    name = "Profile-based"
    version = "3.4.0"

    def __init__(
        self,
        codon_table_path: str | None = None,
        generation_reference_id: str | None = None,
    ) -> None:
        """
        Args:
            codon_table_path: Optional path to a codon usage table JSON used for
                design-time reverse translation. When omitted (the default), the
                bundled host reference table is used — this is the only path the
                product `optimize` CLI and the registered engine ever take, so
                their behavior is unchanged. Injection is used by the benchmark
                harness to drive design with alternative source-profile tables.
            generation_reference_id: Explicit manifest reference id when the
                injected design table is also the CAI evaluation reference.
        """
        self.validator = InputValidator()
        self._codon_table_path = codon_table_path
        # When a source-profile table is injected, drive every CAI-dependent
        # profile from it: it is both the design weight table and the golden set
        # that high_cai maximizes against. Otherwise high_cai would keep chasing
        # the bundled golden set and stay invariant to the source profile,
        # making the injection incomplete. gc_target is GC-driven and remains
        # legitimately invariant either way.
        self.translator = ReverseTranslator(
            codon_table_path=codon_table_path,
            golden_set_path=codon_table_path,
            generation_reference_id=generation_reference_id,
            _allow_golden_set_metadata_fallback=codon_table_path is not None,
        )
        self.rule_engine = RuleEngine()
        self.exporter = SequenceExporter()

    def optimize(
        self,
        sequence: str,
        profile: str | None = "balanced",
        host: str = "nbenthamiana",
        seed: int | None = None,
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

        if (
            self._codon_table_path is None
            and host != "nbenthamiana"
            and profile_value == "high_cai"
        ):
            logger.warning(
                "high_cai is anchored to the N. benthamiana golden set and "
                "ignores the requested host=%s; output is host-invariant by design.",
                host,
            )

        if self._codon_table_path is not None or host == "nbenthamiana":
            # An injected design table (benchmark source-profile runs) is
            # authoritative: it already encodes the host it was derived from.
            translator = self.translator
            rule_engine = self.rule_engine
        else:
            translator = ReverseTranslator(host=host)
            rule_engine = RuleEngine(host=host)

        # 3. Optimize DNA/CDS through the translated protein, then restore the
        # original terminal-stop policy.  DNA must never be treated as a
        # protein string (300 bp must not become 300 aa / 900 bp).
        if seq_type == "dna":
            if len(processed_seq) % 3:
                raise ValueError("DNA/CDS input length must be divisible by 3.")
            translated = translate_dna(processed_seq)
            if "*" in translated[:-1]:
                raise ValueError("DNA/CDS input contains an internal stop codon.")
            terminal_stop = processed_seq[-3:] if translated.endswith("*") else ""
            protein = translated.rstrip("*")
            generated = translator.generate_candidates(
                protein, profile=opt_profile, n=1, seed=seed, **{
                    k: v for k, v in kwargs.items()
                    if k not in ("scan_mode", "scan_include", "scan_exclude")
                }
            )
            if not generated:
                raise ValueError("No candidates generated for DNA/CDS input.")
            optimized_dna = generated[0]["sequence"] + terminal_stop
            if len(optimized_dna) != len(processed_seq):
                raise ValueError("CDS invariant failed: optimized nucleotide length differs.")
            if translate_dna(optimized_dna).rstrip("*") != protein:
                raise ValueError("CDS invariant failed: translated protein was not preserved.")
            cai = translator.calculate_cai(optimized_dna)
            gc = translator.calculate_gc_content(optimized_dna)
            score = calculate_composite_score(
                cai=cai, gc=gc, sequence=optimized_dna, profile=profile_value, host=host
            )
            candidates = [{"sequence": optimized_dna, "cai": cai, "gc": gc, "score": score}]
        else:
            translate_kwargs = {
                k: v
                for k, v in kwargs.items()
                if k not in ("scan_mode", "scan_include", "scan_exclude")
            }
            candidates = translator.generate_candidates(
                processed_seq, profile=opt_profile, n=1, seed=seed, **translate_kwargs
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
            "cai_authority": dict(translator.cai_authority),
            # Keep both names for compatibility across existing tests/callers.
            "gc_content": candidates[0]["gc"],
            "gc_percent": candidates[0]["gc"],
            "score": candidates[0]["score"],
            "violations": sum(len(v) for v in scan_results.values()),
        }
        if profile_value == "balanced":
            host_gc_min, host_gc_max = resolve_host_gc_range(host)
            requested_gc_min_percent = float(kwargs.get("target_gc_min", host_gc_min))
            requested_gc_max_percent = float(kwargs.get("target_gc_max", host_gc_max))
            achieved_gc_percent = metrics["gc_percent"]
            metrics.update(
                {
                    "gc_target_reached": (
                        requested_gc_min_percent <= achieved_gc_percent <= requested_gc_max_percent
                    ),
                    "requested_gc_min_percent": requested_gc_min_percent,
                    "requested_gc_max_percent": requested_gc_max_percent,
                }
            )
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
