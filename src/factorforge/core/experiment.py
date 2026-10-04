"""
FactorForge Core Experiment and Validation Models (Job 355).

Defines the generalized DesignExperiment, Variant, and EvidenceRecord models
for multi-candidate lineage tracking and structured evidence evaluation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
import hashlib
import json
from typing import Any, Optional


class ValidationLevel(str, Enum):
    """Categorization of evidence confidence and source."""
    DETERMINISTIC_CHECK = "deterministic_check"          # Mathematical invariants (AA identity, Type IIS)
    COMPUTATIONAL_PREDICTION = "computational_prediction" # Advisory models (Splice motifs, RNA MFE)
    EXPERT_DISPOSITION = "expert_disposition"            # Human scientist / partner policy judgement
    WET_LAB_MEASURED = "wet_lab_measured"                # Empirical in-vivo / in-vitro measurement
    NOT_VALIDATED = "not_validated"                      # Unmeasured / pending empirical verification


class Severity(str, Enum):
    """Enforcement severity of findings."""
    HARD_FAIL = "HARD_FAIL"
    WARNING = "WARNING"
    ADVISORY = "ADVISORY"
    INFO = "INFO"


@dataclass(frozen=True)
class EvidenceRecord:
    """A single structured evidence observation for a candidate variant."""
    check_id: str                      # e.g., "assembly.type_iis.bsai"
    category: str                      # "cloning_hygiene", "genetics", "rna_stability", "wet_lab"
    result: str                        # "PASS", "FLAGGED", "WARNING", "NOT_TESTED", "NOT_COMPUTED"
    severity: Severity
    method: str                        # e.g., "ExactPatternScanner v1.0", "HeuristicDonorScanner"
    tool_version: str                  # e.g., "v3.5.4", "NetGene2-bridge"
    confidence: Optional[float]        # 0.0 ~ 1.0 (deterministic = 1.0, untested = None)
    evidence_detail: str               # e.g., "0 sites detected across 708 nt"
    validation_level: ValidationLevel
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def to_dict(self) -> dict[str, Any]:
        return {
            "check_id": self.check_id,
            "category": self.category,
            "result": self.result,
            "severity": self.severity.value if isinstance(self.severity, Severity) else self.severity,
            "method": self.method,
            "tool_version": self.tool_version,
            "confidence": self.confidence,
            "evidence_detail": self.evidence_detail,
            "validation_level": self.validation_level.value if isinstance(self.validation_level, ValidationLevel) else self.validation_level,
            "timestamp": self.timestamp,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> EvidenceRecord:
        return cls(
            check_id=data["check_id"],
            category=data["category"],
            result=data["result"],
            severity=Severity(data["severity"]),
            method=data["method"],
            tool_version=data["tool_version"],
            confidence=data.get("confidence"),
            evidence_detail=data["evidence_detail"],
            validation_level=ValidationLevel(data["validation_level"]),
            timestamp=data.get("timestamp", datetime.now(timezone.utc).isoformat()),
        )


@dataclass(frozen=True)
class Variant:
    """A specific DNA CDS sequence design with lineage, intent, and evidence."""
    variant_id: str
    parent_variant_id: Optional[str]   # None if root baseline
    design_intent: str                 # e.g., "Clean biological baseline", "Watermark provenance"
    interventions: list[str]           # e.g., ["none"], ["watermark"], ["splice_remediation"]
    sequence: str                      # DNA CDS string
    protein_sequence: str              # Translated protein string
    gc_percent: float
    cai_score: Optional[float]         # None if not computed
    policy_version: str                # e.g., "plantform_balanced_codon_v1"
    valid_for_release: bool = True     # Gated: False if any HARD_FAIL check is FAIL
    validation_records: list[EvidenceRecord] = field(default_factory=list)
    provenance_hash: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Check hard failure gating
        has_hard_fail = any(
            r.severity == Severity.HARD_FAIL and r.result == "FAIL"
            for r in self.validation_records
        )
        if has_hard_fail and self.valid_for_release:
            object.__setattr__(self, "valid_for_release", False)

        if not self.provenance_hash:
            # Deterministic SHA-256 digest of sequence and core parameters
            h = hashlib.sha256()
            h.update(self.sequence.encode("utf-8"))
            h.update(self.protein_sequence.encode("utf-8"))
            h.update(self.policy_version.encode("utf-8"))
            h.update(",".join(sorted(self.interventions)).encode("utf-8"))
            h.update(str(self.valid_for_release).encode("utf-8"))
            object.__setattr__(self, "provenance_hash", h.hexdigest())

    def to_dict(self) -> dict[str, Any]:
        return {
            "variant_id": self.variant_id,
            "parent_variant_id": self.parent_variant_id,
            "design_intent": self.design_intent,
            "interventions": list(self.interventions),
            "sequence": self.sequence,
            "protein_sequence": self.protein_sequence,
            "gc_percent": self.gc_percent,
            "cai_score": self.cai_score,
            "policy_version": self.policy_version,
            "valid_for_release": self.valid_for_release,
            "validation_records": [r.to_dict() for r in self.validation_records],
            "provenance_hash": self.provenance_hash,
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Variant:
        val_records = [
            EvidenceRecord.from_dict(r) for r in data.get("validation_records", [])
        ]
        cai_raw = data.get("cai_score")
        cai_score = float(cai_raw) if cai_raw is not None else None
        return cls(
            variant_id=data["variant_id"],
            parent_variant_id=data.get("parent_variant_id"),
            design_intent=data["design_intent"],
            interventions=list(data.get("interventions", [])),
            sequence=data["sequence"],
            protein_sequence=data["protein_sequence"],
            gc_percent=float(data["gc_percent"]),
            cai_score=cai_score,
            policy_version=data["policy_version"],
            valid_for_release=bool(data.get("valid_for_release", True)),
            validation_records=val_records,
            provenance_hash=data.get("provenance_hash", ""),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass(frozen=True)
class ComparisonMatrix:
    """Pairwise analytical comparison between two variants in the experiment."""
    comparison_name: str               # e.g., "Watermark Provenance Confounder Check"
    variant_a_id: str                  # Baseline variant
    variant_b_id: str                  # Comparison variant
    variable_evaluated: str            # e.g., "digital watermark insertion"
    gc_shift_percent: float            # variant_b.gc - variant_a.gc
    sequence_divergence_nt: int        # Absolute mismatch count (Hamming distance)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if abs(self.gc_shift_percent) > 5.0 and "confounder_warning" not in self.metadata:
            meta = dict(self.metadata)
            meta["confounder_warning"] = (
                f"GC shift of {self.gc_shift_percent:+.2f}%p exceeds 5.0%p threshold; "
                f"expression differential risks being confounded by mRNA secondary structure or GC content."
            )
            object.__setattr__(self, "metadata", meta)

    def to_dict(self) -> dict[str, Any]:
        return {
            "comparison_name": self.comparison_name,
            "variant_a_id": self.variant_a_id,
            "variant_b_id": self.variant_b_id,
            "variable_evaluated": self.variable_evaluated,
            "gc_shift_percent": self.gc_shift_percent,
            "sequence_divergence_nt": self.sequence_divergence_nt,
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ComparisonMatrix:
        return cls(
            comparison_name=data["comparison_name"],
            variant_a_id=data["variant_a_id"],
            variant_b_id=data["variant_b_id"],
            variable_evaluated=data["variable_evaluated"],
            gc_shift_percent=float(data["gc_shift_percent"]),
            sequence_divergence_nt=int(data["sequence_divergence_nt"]),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass
class DesignExperiment:
    """A generalized design experiment container holding arbitrary variants and evidence."""
    experiment_id: str
    target_name: str
    host_organism: str
    hypothesis: str
    frozen_regions: list[dict[str, Any]] = field(default_factory=list)
    variants: list[Variant] = field(default_factory=list)
    comparisons: list[ComparisonMatrix] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def add_variant(self, variant: Variant) -> None:
        """Add an arbitrary candidate variant to the experiment."""
        if any(v.variant_id == variant.variant_id for v in self.variants):
            raise ValueError(f"Variant ID '{variant.variant_id}' already exists in experiment.")
        self.variants.append(variant)

    def get_variant(self, variant_id: str) -> Variant:
        """Retrieve a variant by ID."""
        for v in self.variants:
            if v.variant_id == variant_id:
                return v
        raise KeyError(f"Variant '{variant_id}' not found in experiment.")

    def get_lineage(self, variant_id: str) -> list[Variant]:
        """Trace the lineage of a variant up to its root ancestor."""
        chain: list[Variant] = []
        curr_id: Optional[str] = variant_id
        visited: set[str] = set()

        while curr_id:
            if curr_id in visited:
                raise ValueError(f"Circular lineage detected at variant '{curr_id}'.")
            visited.add(curr_id)
            variant = self.get_variant(curr_id)
            chain.append(variant)
            curr_id = variant.parent_variant_id

        return list(reversed(chain))

    def add_comparison(self, comparison: ComparisonMatrix) -> None:
        """Add a pairwise comparison between two variants."""
        self.comparisons.append(comparison)

    def to_dict(self) -> dict[str, Any]:
        """Serialize complete experiment to a dictionary."""
        return {
            "experiment_id": self.experiment_id,
            "target_name": self.target_name,
            "host_organism": self.host_organism,
            "hypothesis": self.hypothesis,
            "frozen_regions": [dict(r) for r in self.frozen_regions],
            "variants": [v.to_dict() for v in self.variants],
            "comparisons": [c.to_dict() for c in self.comparisons],
            "metadata": dict(self.metadata),
        }

    def to_json(self, indent: int = 2) -> str:
        """Serialize complete experiment to JSON."""
        return json.dumps(self.to_dict(), indent=indent, ensure_ascii=False)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> DesignExperiment:
        """Reconstruct experiment from dictionary."""
        exp = cls(
            experiment_id=data["experiment_id"],
            target_name=data["target_name"],
            host_organism=data["host_organism"],
            hypothesis=data["hypothesis"],
            frozen_regions=[dict(r) for r in data.get("frozen_regions", [])],
            metadata=dict(data.get("metadata", {})),
        )
        for v_dict in data.get("variants", []):
            exp.add_variant(Variant.from_dict(v_dict))
        for c_dict in data.get("comparisons", []):
            exp.add_comparison(ComparisonMatrix.from_dict(c_dict))
        return exp

    @classmethod
    def from_json(cls, json_str: str) -> DesignExperiment:
        """Reconstruct experiment from JSON string."""
        return cls.from_dict(json.loads(json_str))
