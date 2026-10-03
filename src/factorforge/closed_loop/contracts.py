"""Data contracts and schemas for Closed-Loop DBTL Lineage."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from typing import Dict, List, Optional
from pydantic import BaseModel, Field


class PolicyRecord(BaseModel):
    """Metadata record of an approved codon distribution policy."""

    policy_id: str
    canonical_sha256: str
    target_distribution: Dict[str, Dict[str, float]]
    region_definition: str = "mature_chain"
    feedback_source_id: str = "doug_2026_10_02_humira_feedback"
    created_at: str = Field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat()
    )


class ConstructRecord(BaseModel):
    """Immutable sequence identity record for a designed construct."""

    construct_id: str
    cds_sha256: str
    protein_sha256: str
    chain: str  # "LC" or "HC"
    watermark_pair_id: str
    watermark_applied: bool
    achieved_distribution: Dict[str, Dict[str, float]]
    policy_id: str
    engine_version: str = "FactorForge-v3.3.0"


class ConstructSetRecord(BaseModel):
    """Co-expressed pair of Light Chain and Heavy Chain constructs."""

    construct_set_id: str
    light_construct_id: str
    heavy_construct_id: str
    watermark_applied: bool
    description: str = ""


class ExperimentRecord(BaseModel):
    """Wet-lab or benchmark experiment metadata."""

    experiment_id: str
    construct_set_id: str
    host: str = "Nicotiana benthamiana"
    expression_system: str = "Agroinfiltration (Transient Leaf Expression)"
    batch_id: str
    experiment_date: str
    protocol_version: str = "PlantForm-Humira-Standard-v1.0"


class MeasurementRecord(BaseModel):
    """Single empirical or synthetic readout record."""

    measurement_id: str
    experiment_id: str
    construct_set_id: str
    biological_replicate_id: str
    technical_replicate_id: str
    measurement_type: str  # "protein_yield", "mrna_abundance", "functional_binding"
    value: Optional[float]
    unit: str  # "mg/L", "relative_fold", "%"
    qc_status: str = "PASS"  # "PASS", "FAIL", "FLAGGED"
    is_synthetic: bool = False  # CRITICAL: If True, flagged as test_only and excluded from real training
    tag: str = "test_only" if True else "empirical"
    notes: str = ""
    recorded_at: str = Field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat()
    )

    @property
    def record_hash(self) -> str:
        payload = {
            "measurement_id": self.measurement_id,
            "experiment_id": self.experiment_id,
            "construct_set_id": self.construct_set_id,
            "biological_replicate_id": self.biological_replicate_id,
            "technical_replicate_id": self.technical_replicate_id,
            "measurement_type": self.measurement_type,
            "value": self.value,
            "unit": self.unit,
            "is_synthetic": self.is_synthetic,
        }
        encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()
