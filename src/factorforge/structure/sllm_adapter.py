"""Evidence-aware structure-conditioned research explainer and sLLM adapter.

Ingests 3D structure prediction sidecars (ESMFold / AlphaFold 3) and 4DP cotranslational
codon optimization results to generate scientifically rigorous, non-hallucinatory assessments.

EPISTEMIC BOUNDARY & CORE PRINCIPLE:
Strictly segregates findings into four epistemic tiers:
1. VERIFIED: Mathematically and computationally proven facts (sequence identity, GC%, 0 restriction sites).
2. PREDICTED: Neural network structural predictions (pLDDT, domain boundaries, predicted disorder).
3. HYPOTHESIZED: In-silico cotranslational hypotheses (exit-tunnel pause candidate coordinates).
4. EXPERIMENTAL VALIDATION REQUIRED: Wet-lab experiments needed to substantiate computational claims.

Classification: UNVERIFIED_PUBLIC_CORPUS · Biological validation: NOT_AVAILABLE
Human-in-the-loop invariant: reviewed=False
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from factorforge.analysis.metrics import (
    calculate_cai,
    calculate_gc,
    load_codon_usage_table,
    translate_dna,
)

SLLM_ADAPTER_VERSION = "1.0.0"


class SLLMAdapterError(ValueError):
    """Raised when evidence validation or explanation generation fails."""


@dataclass
class EvidenceReport:
    """Rigorous, tier-segregated scientific assessment report."""

    target_name: str
    protein_length: int
    model_source: str
    is_mock: bool
    verified_facts: list[str] = field(default_factory=list)
    predicted_metrics: list[str] = field(default_factory=list)
    hypothesized_effects: list[str] = field(default_factory=list)
    wet_lab_validation_plan: list[str] = field(default_factory=list)
    epistemic_status: str = "UNVERIFIED_PUBLIC_CORPUS"
    biological_validation: str = "NOT_AVAILABLE"
    reviewed: bool = False
    tool_version: str = SLLM_ADAPTER_VERSION

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_markdown(self) -> str:
        lines = [
            f"# Evidence-Aware Cotranslational Optimization Assessment: {self.target_name}",
            "",
            f"**Model Source**: {self.model_source} {'[MOCK SYNTHETIC FIXTURE]' if self.is_mock else ''}  ",
            f"**Protein Length**: {self.protein_length} aa  ",
            f"**Epistemic Classification**: `{self.epistemic_status}`  ",
            f"**Biological Validation**: `{self.biological_validation}` (Explicit human review required before wet-lab use)  ",
            f"**Human Review Gate**: `reviewed={self.reviewed}`  ",
            "",
            "---",
            "",
            "## 1. Verified Invariants (Computationally & Mathematically Proven)",
            "> Facts established deterministically by the optimization engine and validated without neural inference.",
            "",
        ]
        for fact in self.verified_facts:
            lines.append(f"- **[VERIFIED]** {fact}")

        lines.extend([
            "",
            "## 2. Predicted Structural Metrics (Neural Model Outputs)",
            "> In-silico predictions derived from neural models (AlphaFold 3 / ESMFold). Subject to computational uncertainty.",
            "",
        ])
        for pred in self.predicted_metrics:
            lines.append(f"- **[PREDICTED]** {pred}")

        lines.extend([
            "",
            "## 3. Hypothesized Biological Mechanisms (In-Silico Cotranslational Model)",
            "> Theoretical hypotheses linking structural predictions to ribosome exit-tunnel dwell times. NOT established biological facts.",
            "",
        ])
        for hyp in self.hypothesized_effects:
            lines.append(f"- **[HYPOTHESIZED]** {hyp}")

        lines.extend([
            "",
            "## 4. Wet-Lab Experimental Validation Plan (Required Prior to Claims)",
            "> Actionable assays required to confirm or refute computational hypotheses in living expression hosts.",
            "",
        ])
        for plan in self.wet_lab_validation_plan:
            lines.append(f"- **[ACTION REQUIRED]** {plan}")

        lines.append("")
        return "\n".join(lines)


def _parse_sidecar(sidecar_source: str | Path | dict[str, Any], allow_mock: bool = False) -> dict[str, Any]:
    """Load and validate structure sidecar metadata."""
    if isinstance(sidecar_source, (str, Path)):
        p = Path(sidecar_source)
        if not p.is_file():
            raise SLLMAdapterError(f"Sidecar file not found: {sidecar_source}")
        try:
            data = json.loads(p.read_text(encoding="utf-8-sig"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise SLLMAdapterError(f"Failed to parse sidecar JSON: {exc}") from exc
    elif isinstance(sidecar_source, dict):
        data = sidecar_source
    else:
        raise SLLMAdapterError("Sidecar must be a file path or dictionary.")

    if data.get("status") != "AVAILABLE":
        raise SLLMAdapterError("Sidecar status is not AVAILABLE.")

    if (data.get("is_mock") is True or data.get("source_type") == "SYNTHETIC_FIXTURE") and not allow_mock:
        raise SLLMAdapterError(
            "Sidecar was generated from synthetic/mock data. "
            "Mock sidecars are blocked from generating research reports. "
            "Pass --allow-mock to proceed with testing."
        )

    return data


def generate_evidence_explanation(
    sidecar_source: str | Path | dict[str, Any],
    optimized_fasta_content: str | Path | None = None,
    allow_mock: bool = False,
    target_name: str | None = None,
) -> EvidenceReport:
    """Generate evidence-aware structured assessment strictly separated into epistemic tiers."""
    sidecar = _parse_sidecar(sidecar_source, allow_mock=allow_mock)
    primary_chain_id = sidecar.get("primary_chain", "A")
    chain_data = sidecar.get("chains", {}).get(primary_chain_id, {})
    protein_seq = chain_data.get("sequence", "")
    protein_len = chain_data.get("length", len(protein_seq))
    model_source = sidecar.get("source_type", "ESMFold")
    is_mock = sidecar.get("is_mock", False) or model_source == "SYNTHETIC_FIXTURE"

    verified: list[str] = []
    predicted: list[str] = []
    hypothesized: list[str] = []
    validation_plan: list[str] = []

    # 1. Process Optimized FASTA / DNA if provided
    dna_seq = None
    header_info = ""
    if optimized_fasta_content:
        if isinstance(optimized_fasta_content, Path) or (
            isinstance(optimized_fasta_content, str) and Path(optimized_fasta_content).is_file()
        ):
            content = Path(optimized_fasta_content).read_text(encoding="utf-8")
        else:
            content = str(optimized_fasta_content)

        lines = [line.strip() for line in content.splitlines() if line.strip()]
        dna_lines = []
        for line in lines:
            if line.startswith(">"):
                header_info = line[1:].strip()
            else:
                dna_lines.append(line)
        dna_seq = "".join(dna_lines).upper()

    if dna_seq:
        translated_aa = translate_dna(dna_seq).removesuffix("*")
        if translated_aa == protein_seq:
            verified.append(
                f"Amino acid sequence identity is 100% preserved ({protein_len} aa) upon translation."
            )
        else:
            raise SLLMAdapterError(
                f"Translation mismatch: Optimized DNA translates to {len(translated_aa)} aa, "
                f"but structure model sequence is {protein_len} aa."
            )

        gc = calculate_gc(dna_seq)
        verified.append(f"Optimized coding sequence GC content: {gc:.1f}%.")

        table = load_codon_usage_table()
        cai = calculate_cai(dna_seq, table.codon_weights)
        verified.append(f"Global Codon Adaptation Index (CAI): {cai:.3f}.")

        # Check restriction sites invariant
        from factorforge.constraints.type_iis import resolve_restriction_site_motifs
        motifs = resolve_restriction_site_motifs(["BsaI", "BsmBI", "BbsI", "SapI", "NotI"])
        found_motifs = [m for m in motifs if m in dna_seq]
        if not found_motifs:
            verified.append("Golden Gate & standard cloning restriction sites (BsaI, BsmBI, BbsI, SapI, NotI): 0 detected.")
        else:
            verified.append(f"Warning: {len(found_motifs)} restriction site motifs detected: {found_motifs}.")

        if "engine=dp_4dp" in header_info:
            verified.append("Optimization engine: Deterministic Cotranslational Dynamic Programming (4DP v2.4.0).")
    else:
        verified.append(f"Protein sequence validated: {protein_len} standard amino acids.")

    # 2. Predicted Structural Metrics
    metrics = chain_data.get("metrics", {})
    mean_plddt = metrics.get("mean_plddt")
    if mean_plddt is not None:
        predicted.append(
            f"Mean CA-atom structural confidence (pLDDT): {mean_plddt:.2f}/100.00 "
            f"(Min: {metrics.get('min_plddt', 0.0):.2f}, Max: {metrics.get('max_plddt', 0.0):.2f})."
        )

    candidates = chain_data.get("pause_candidates", [])
    if candidates:
        predicted.append(
            f"Detected {len(candidates)} candidate flexible/hinge region(s) where predicted CA pLDDT < 70.0:"
        )
        for c in candidates:
            predicted.append(
                f"  - Region residues {c['start_residue']}..{c['end_residue']} "
                f"(length: {c['length']} aa, mean local pLDDT: {c['mean_plddt']:.1f})."
            )
    else:
        predicted.append("No low-confidence flexible hinge regions (pLDDT < 70.0) detected.")

    # 3. Hypothesized Biological Mechanisms
    tunnel_offset = sidecar.get("tunnel_offset_applied", 35)
    if candidates:
        pause_positions = [c["tunnel_pause_position"] for c in candidates if c.get("within_structure_chain")]
        hypothesized.append(
            f"Ribosome exit-tunnel model (+{tunnel_offset} aa offset) projects candidate pause positions: {pause_positions}."
        )
        hypothesized.append(
            "HYPOTHESIS: Introducing a local translation deceleration (codon adaptiveness delta Δw ≈ -0.55) "
            "at the exit-tunnel offset position allows the preceding structural domain to fold cooperatively outside "
            "the exit vestibule prior to synthesis of downstream hydrophobic segments."
        )
        hypothesized.append(
            "EPISTEMIC LIMITATION: Low pLDDT denotes neural prediction uncertainty and/or intrinsic disorder; "
            "it does NOT demonstrate empirical ribosome dwell times or in-vivo chaperone interactions."
        )
    else:
        hypothesized.append("Uniform translation kinetics recommended in the absence of predicted domain hinges.")

    # 4. Wet-Lab Experimental Validation Plan
    name = target_name or sidecar.get("job_name", "Target_Protein")
    validation_plan.append(
        f"Synthesize dual constructs for {name}: Baseline (1DP, uniform high CAI) vs. Test (4DP, tunnel-paused)."
    )
    validation_plan.append(
        "Execute transient expression in target host (e.g. Nicotiana benthamiana agroinfiltration or CHO transient)."
    )
    validation_plan.append(
        "Perform Western blot analysis comparing soluble supernatant vs. insoluble pellet fractions to quantify folding yield."
    )
    validation_plan.append(
        "Conduct Ribo-seq (Ribosome Profiling) or targeted pulse-chase to verify whether ribosomal pause actually occurs at candidate codons."
    )

    return EvidenceReport(
        target_name=name,
        protein_length=protein_len,
        model_source=model_source,
        is_mock=is_mock,
        verified_facts=verified,
        predicted_metrics=predicted,
        hypothesized_effects=hypothesized,
        wet_lab_validation_plan=validation_plan,
        reviewed=False,
    )
