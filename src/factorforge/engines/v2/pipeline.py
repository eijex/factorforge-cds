"""
Optimization pipeline for FactorForge v2.
Integrates validation, translation, rule scanning, domestication, and construct building.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from factorforge.engines.v2.construct_builder import ConstructBuilder
from factorforge.engines.v2.rules.domesticator import Domesticator
from factorforge.engines.v2.rules.reverse_translator import OptimizationProfile, ReverseTranslator
from factorforge.engines.v2.rules.rule_engine import RuleEngine
from factorforge.engines.v2.scoring import calculate_composite_score
from factorforge.engines.v2.validator import InputValidator
from factorforge.utils.construct_id import generate_construct_id

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord


@dataclass
class PipelineResult:
    """Pipeline output container."""

    sequence: str
    construct: "SeqRecord | None" = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def export_features(self) -> dict[str, Any]:
        """schema.md 호환 피처 dict 반환 (purity_pct 제외 — 실험 후 수동 입력)."""
        metrics = self.metadata.get("metrics", {})
        scan = self.metadata.get("scan_results", {})
        dom = self.metadata.get("domestication", {})

        return {
            "construct_id": self.metadata.get("construct_id", ""),
            "protein_name": "",
            "optimization_profile": self.metadata.get("profile", ""),
            "cai_score": round(metrics.get("cai", 0.0), 4),
            "gc_content_pct": round(metrics.get("gc", 0.0), 2),
            "mfe_kcal_mol": round(metrics.get("mfe", 0.0), 2),
            "polya_signal_count": len(scan.get("polya", [])),
            "domestication_edits": len(dom.get("removed_sites", [])),
            "sequence_length_aa": len(self.sequence) // 3,
            "agro_od600": None,
            "dpi": None,
            "purity_pct": None,
            "yield_mg_per_kg": None,
        }

    def save(self, filepath: Path, format: str = "fasta") -> None:
        """
        Save the result to a file.

        Args:
            filepath: Output file path.
            format: "fasta" or "genbank".

        Raises:
            ImportError: If Biopython is required but not installed.
        """
        format_lower = format.lower()

        if self.construct is not None and format_lower == "genbank":
            try:
                from Bio import SeqIO
            except ImportError as exc:
                raise ImportError("Biopython is required: pip install biopython") from exc

            SeqIO.write(self.construct, str(filepath), "genbank")
            return

        with open(filepath, "w", encoding="utf-8") as handle:
            handle.write(f">optimized\n{self.sequence}\n")


class OptimizationPipeline:
    """End-to-end v2 optimization pipeline."""

    def __init__(
        self,
        profile: str = "balanced",
        construct_template: str | None = None,
        template_dir: Path | None = None,
    ) -> None:
        """
        Args:
            profile: Optimization profile name.
            construct_template: Optional construct template name.
            template_dir: Optional template directory.
        """
        self.profile = profile
        self.construct_template = construct_template
        self.template_dir = template_dir

        self.validator = InputValidator()
        self.translator = ReverseTranslator()
        self.rule_engine = RuleEngine()
        self.domesticator = Domesticator()

        if construct_template:
            if template_dir is None:
                template_dir = Path(__file__).resolve().parents[4] / "data" / "templates"
            self.construct_builder: ConstructBuilder | None = ConstructBuilder(template_dir)
        else:
            self.construct_builder = None

    def run(
        self,
        sequence: str,
        profile: str | None = None,
        construct_template: str | None = None,
        **kwargs: Any,
    ) -> PipelineResult:
        """
        Run the optimization pipeline.

        Args:
            sequence: Input protein or DNA sequence.
            profile: Optional profile override.
            construct_template: Optional template override.
            **kwargs: Additional settings.

        Returns:
            PipelineResult with sequence and metadata.

        Raises:
            ValueError: If input sequence is invalid.
        """
        logger.info(f"Starting optimization pipeline with profile: {profile or self.profile}")

        val_result = self.validator.validate(sequence)
        if not val_result["valid"]:
            logger.error(f"Input validation failed: {val_result['errors']}")
            raise ValueError(f"Invalid input sequence: {val_result['errors']}")

        processed = val_result["processed_sequence"]
        seq_type = val_result["type"]
        logger.debug(f"Detected sequence type: {seq_type}")
        if seq_type == "fasta":
            seq_type = self.validator.detect_sequence_type(processed).value

        effective_profile = (profile or self.profile or "balanced").lower()
        try:
            opt_profile = OptimizationProfile(effective_profile)
        except ValueError as exc:
            supported = ", ".join(p.value for p in OptimizationProfile)
            raise ValueError(
                f"Unknown profile: {effective_profile}. Supported profiles: {supported}"
            ) from exc

        if seq_type == "dna":
            optimized_dna = processed
            cai = self.translator.calculate_cai(optimized_dna)
            gc = self.translator.calculate_gc_content(optimized_dna)
            score = calculate_composite_score(
                cai=cai, gc=gc, sequence=optimized_dna, profile=effective_profile
            )
            candidate_metrics = {"cai": cai, "gc": gc, "score": score}
        else:
            logger.debug(f"Generating candidates with profile: {opt_profile.value}")
            candidates = self.translator.generate_candidates(processed, profile=opt_profile, n=1)
            if not candidates:
                logger.error("No candidates generated for input sequence")
                raise ValueError("No candidates generated for input sequence.")
            optimized_dna = candidates[0]["sequence"]
            candidate_metrics = {
                "cai": candidates[0]["cai"],
                "gc": candidates[0]["gc"],
                "score": candidates[0]["score"],
            }
            logger.info(
                f"Generated optimized sequence: CAI={candidate_metrics['cai']:.3f}, "
                f"GC={candidate_metrics['gc']:.1f}%"
            )

        # Fast pre-check avoids an expensive full rule scan before PolyA fixing.
        has_polya_signal = any(
            pattern in optimized_dna for pattern in self.rule_engine.POLYA_PATTERNS
        )
        if has_polya_signal:
            logger.debug("Potential PolyA signal detected; attempting iterative fix")
            polya_fix = self.rule_engine.fix_polya_iterative(optimized_dna)
            if polya_fix["success"]:
                optimized_dna = polya_fix["modified_seq"]
                logger.info(
                    f"Fixed {len(polya_fix['fixes_applied'])} PolyA violation(s) "
                    f"in {polya_fix['rounds']} round(s)"
                )
            else:
                logger.warning(
                    f"Could not fix all PolyA violations. "
                    f"Remaining: {polya_fix.get('remaining_violations', '?')}"
                )

        # Dinucleotide reduction pass (CpG/TpA greedy synonymous fix)
        if self.rule_engine.scan_dinucleotides(optimized_dna):
            dinu_fix = self.rule_engine.fix_dinucleotides(optimized_dna, mode="balanced")
            if dinu_fix["success"]:
                optimized_dna = dinu_fix["modified_seq"]
                candidate_metrics["cai"] = round(self.translator.calculate_cai(optimized_dna), 4)
                candidate_metrics["gc"] = self.translator.calculate_gc_content(optimized_dna)
                candidate_metrics["score"] = calculate_composite_score(
                    cai=candidate_metrics["cai"],
                    gc=candidate_metrics["gc"],
                    sequence=optimized_dna,
                    profile=effective_profile,
                )
                logger.info(
                    f"Dinucleotide reduction [{dinu_fix['mode']}]: "
                    f"{dinu_fix['initial_count']} -> "
                    f"{dinu_fix['final_count']} ({dinu_fix['reduction_pct']}% reduction, "
                    f"CAI {dinu_fix['cai_before']} -> {dinu_fix['cai_after']}, "
                    f"{dinu_fix['rounds']} round(s))"
                )

        logger.debug("Scanning for final rule violations")
        scan_mode = str(kwargs.get("scan_mode", "full"))
        scan_include = kwargs.get("scan_include")
        scan_exclude = kwargs.get("scan_exclude")
        scan_results = self.rule_engine.scan_all(
            optimized_dna,
            mode=scan_mode,
            include=scan_include,
            exclude=scan_exclude,
        )

        assembly_standard = kwargs.get("assembly_standard", "golden_gate")
        domestication = self.domesticator.domesticate(optimized_dna, standard=assembly_standard)
        domesticated_sequence = domestication.get("domesticated_seq", optimized_dna)

        template_name = construct_template or self.construct_template
        if template_name:
            if self.construct_builder is None:
                template_dir = (
                    self.template_dir or Path(__file__).resolve().parents[4] / "data" / "templates"
                )
                self.construct_builder = ConstructBuilder(template_dir)
            construct_record = self.construct_builder.generate_construct(
                gene_sequence=domesticated_sequence,
                template_name=template_name,
            )
            final_sequence = str(construct_record.seq)
        else:
            construct_record = None
            final_sequence = domesticated_sequence

        metadata: dict[str, Any] = {
            "construct_id": generate_construct_id(),
            "profile": effective_profile,
            "construct_template": template_name,
            "construct_features": len(construct_record.features) if construct_record else 0,
            "validation": val_result,
            "scan_results": scan_results,
            "domestication": domestication,
            "metrics": candidate_metrics,
            "scan_mode": scan_mode,
        }

        return PipelineResult(
            sequence=final_sequence,
            construct=construct_record,
            metadata=metadata,
        )

    def run_batch(
        self,
        sequences: list[dict[str, str]] | list[str],
        profile: str | None = None,
        construct_template: str | None = None,
        **kwargs: Any,
    ) -> list[PipelineResult]:
        """Run the optimization pipeline for a batch of sequences."""
        results: list[PipelineResult] = []
        for idx, entry in enumerate(sequences, start=1):
            if isinstance(entry, dict):
                seq = entry.get("sequence", "")
                seq_id = entry.get("id", f"seq{idx}")
            else:
                seq = entry
                seq_id = f"seq{idx}"
            result = self.run(
                seq,
                profile=profile,
                construct_template=construct_template,
                **kwargs,
            )
            result.metadata["input_id"] = seq_id
            results.append(result)
        return results
