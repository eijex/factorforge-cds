"""Generate structured wet-lab feedback packages for FactorForge designs."""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


@dataclass
class WetLabResult:
    construct_id: str
    factorforge_version: str
    host_profile: str
    profile: Optional[str]
    sequence_hash: str

    protein_name: str
    host_organism: str
    promoter: Optional[str]
    subcellular_targeting: Optional[str]
    expression_system: str
    harvest_timepoint: Optional[str]
    native_control: Optional[str]

    comparison: str
    expression_level: Optional[str]
    notes: Optional[str]

    institution: Optional[str]
    public_listing: bool = True

    submitted_at: Optional[str] = None


class ValidationPackageGenerator:
    """Create validation package files from structured wet-lab result metadata."""

    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)

    def generate(self, result: WetLabResult) -> Path:
        """Write validation package files and return the output directory."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        submitted_at = result.submitted_at or self._utc_now()
        metadata = self._build_metadata(result, submitted_at)

        (self.output_dir / "validation_metadata.json").write_text(
            json.dumps(metadata, indent=2),
            encoding="utf-8",
        )
        (self.output_dir / "validation_summary.txt").write_text(
            self._render_summary(result, submitted_at),
            encoding="utf-8",
        )
        (self.output_dir / "issue_body.md").write_text(
            self._render_issue_body(result, submitted_at),
            encoding="utf-8",
        )
        return self.output_dir

    def _build_metadata(self, result: WetLabResult, submitted_at: str) -> dict:
        return {
            "construct_id": result.construct_id,
            "factorforge_version": result.factorforge_version,
            "host_profile": result.host_profile,
            "profile": result.profile,
            "submitted_at": submitted_at,
            "submitter": {
                "institution": result.institution,
            },
            "experiment": {
                "protein_name": result.protein_name,
                "host_organism": result.host_organism,
                "promoter": result.promoter,
                "subcellular_targeting": result.subcellular_targeting,
                "expression_system": result.expression_system,
                "harvest_timepoint": result.harvest_timepoint,
                "native_control": result.native_control,
            },
            "result": {
                "comparison": result.comparison,
                "expression_level": result.expression_level,
                "notes": result.notes or "",
            },
            "consent": {
                "public_listing": result.public_listing,
            },
            "sequence_hash": result.sequence_hash,
        }

    def _render_summary(self, result: WetLabResult, submitted_at: str) -> str:
        public_listing = "Yes" if result.public_listing else "No"
        return "\n".join(
            [
                "FactorForge Validation Report",
                "==============================",
                f"Construct ID   : {result.construct_id}",
                f"FF Version     : {result.factorforge_version}",
                f"Host Profile   : {result.host_profile}",
                f"Profile        : {self._display(result.profile)}",
                f"Submitted      : {submitted_at}",
                "",
                "--- Experiment ---",
                f"Protein        : {result.protein_name}",
                f"Host Organism  : {result.host_organism}",
                f"Promoter       : {self._display(result.promoter)}",
                f"Targeting      : {self._display(result.subcellular_targeting)}",
                f"Assay          : {result.expression_system}",
                f"Harvest        : {self._display(result.harvest_timepoint)}",
                f"Native Control : {self._display(result.native_control)}",
                "",
                "--- Result ---",
                f"Comparison     : {result.comparison}",
                f"Expression     : {self._display(result.expression_level)}",
                f"Notes          : {self._display(result.notes)}",
                "",
                "--- Sequence ---",
                f"Sequence Hash  : {result.sequence_hash}",
                "",
                "--- Consent ---",
                f"Public Listing : {public_listing}",
                f"Institution    : {self._display(result.institution)}",
                "",
            ]
        )

    def _render_issue_body(self, result: WetLabResult, submitted_at: str) -> str:
        checked = "x" if result.public_listing else " "
        return "\n".join(
            [
                "## Wet-lab Validation Result",
                "",
                "| Field | Value |",
                "|-------|-------|",
                f"| **Construct ID** | {result.construct_id} |",
                f"| **FactorForge Version** | {result.factorforge_version} |",
                f"| **Host Profile** | {result.host_profile} |",
                f"| **Profile** | {self._display(result.profile)} |",
                f"| **Submitted** | {submitted_at} |",
                "",
                "## Experiment Details",
                "",
                "| Field | Value |",
                "|-------|-------|",
                f"| **Protein** | {result.protein_name} |",
                f"| **Host Organism** | {result.host_organism} |",
                f"| **Promoter** | {self._display(result.promoter)} |",
                f"| **Targeting** | {self._display(result.subcellular_targeting)} |",
                f"| **Assay** | {result.expression_system} |",
                f"| **Harvest** | {self._display(result.harvest_timepoint)} |",
                f"| **Native Control** | {self._display(result.native_control)} |",
                "",
                "## Result",
                "",
                "| Field | Value |",
                "|-------|-------|",
                f"| **Comparison** | {result.comparison} |",
                f"| **Expression** | {self._display(result.expression_level)} |",
                f"| **Notes** | {self._display(result.notes)} |",
                "",
                "## Sequence Integrity",
                "",
                f"Sequence Hash: {result.sequence_hash}",
                "",
                "_Raw sequence not included per privacy policy._",
                "",
                "## Consent",
                "",
                f"- [{checked}] Public listing approved",
                f"- Institution: {self._display(result.institution)}",
                "",
            ]
        )

    @staticmethod
    def _display(value: Optional[str]) -> str:
        return value if value else "-"

    @staticmethod
    def _utc_now() -> str:
        return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
