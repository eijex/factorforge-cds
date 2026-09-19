"""Data models and schemas for Prospective Experimental Acquisition & Paired DBTL (Job 283C).

Decouples in silico candidate design, physical experiments, randomized/blinded samples,
raw empirical measurements, and derived outcomes into a tamper-evident paired schema.
"""

from __future__ import annotations
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import hashlib
import json
from typing import Any, Dict, List, Optional

from factorforge.discovery.schemas import TraitVector


@dataclass
class ConstructDesignRecord:
    """In silico construct design snapshot.

    Contains computational traits, hypothesis rationale, and cryptographic sequence digest.
    Sequence DNA is optional in sequence-free memory but preserved in secure reproducibility archive.
    """

    construct_id: str
    target_name: str
    mature_protein_aa_length: int
    construct_aa_length: int
    hypothesis_id: str  # e.g., 'H0', 'H1', 'H2', 'CTRL_POS', 'CTRL_NEG'
    hypothesis_name: str  # e.g., 'Primary Deterministic Optimum'
    generation_contract: str
    sequence_digest: str  # SHA-256 of DNA sequence
    trait_vector: Optional[TraitVector] = None
    sequence_dna: Optional[str] = None
    is_control: bool = False
    control_type: Optional[str] = None  # 'positive_expression', 'negative_empty_vector', etc.
    rationale: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self, include_sequence: bool = True) -> Dict[str, Any]:
        res = asdict(self)
        if self.trait_vector is not None:
            res["trait_vector"] = self.trait_vector.to_dict()
        if not include_sequence:
            res.pop("sequence_dna", None)
        return res


@dataclass
class ExperimentRunRecord:
    """Metadata describing the wet-lab experimental run and environmental context."""

    experiment_id: str
    batch_id: str
    host_organism: str = "Nicotiana benthamiana"
    growth_conditions: str = "24C, 16h light / 8h dark, 60% RH"
    infiltration_od: float = 0.5
    harvest_dpi: int = 4
    protocol_version: str = "AGRO-INFIL-v2.1"
    operator_id: str = "OP-LAB-01"
    run_timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class SampleRecord:
    """Physical biological/technical sample linking blinded plate position to construct.

    Decouples biological sample lineage (plant, leaf, biological replicate)
    from physical assay plate randomization (plate ID, well coordinate).
    """

    sample_id: str
    blinded_plate_position: str  # e.g., 'A01', 'B04' (alias for assay_well_position)
    experiment_id: str
    construct_id: str
    biological_replicate: int  # 1, 2, 3...
    technical_replicate: int = 1
    plant_id: Optional[str] = None  # e.g., 'PLANT-01'
    leaf_id: Optional[str] = None  # e.g., 'LEAF-3'
    biological_replicate_id: Optional[str] = None  # e.g., 'BIO-REP-01'
    assay_plate_id: Optional[str] = "PLATE-01"
    assay_well_position: Optional[str] = None  # e.g., 'A01'
    is_control: bool = False
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.assay_well_position is None:
            self.assay_well_position = self.blinded_plate_position
        if self.blinded_plate_position is None:
            self.blinded_plate_position = self.assay_well_position

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class MeasurementRecord:
    """Raw analytical instrument measurement for a physical sample."""

    measurement_id: str
    sample_id: str
    assay_type: str  # e.g., 'GFP_FLUORESCENCE', 'QUANTITATIVE_ELISA', 'WESTERN_BLOT'
    raw_yield_ug_g_fw: Optional[float] = None
    western_blot_intensity: Optional[float] = None
    solubility_ratio: Optional[float] = None
    chlorosis_score: Optional[int] = None  # 0 (none) to 3 (severe tissue necrosis)
    qc_status: str = "PASS"  # 'PASS', 'WARN', 'FAIL'
    flags: List[str] = field(default_factory=list)
    raw_data_ref: Optional[str] = None
    measurement_source: str = "EMPIRICAL"  # 'EMPIRICAL', 'MOCK_SIMULATED', 'SYNTHETIC_PILOT'
    evidence_tier: str = "TIER_1_VALIDATED"  # 'TIER_1_VALIDATED', 'TEST_ONLY', 'EXPLORATORY'
    training_eligible: bool = True  # strictly False for mock data
    measured_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class DerivedOutcomeRecord:
    """Normalized and comparative outcome derived from raw measurements."""

    outcome_id: str
    sample_id: str
    construct_id: str
    experiment_id: str
    normalized_yield_ug_g_fw: float
    relative_to_h0: Optional[float] = None  # Fold-change relative to H0 within target/batch
    relative_to_pos_ctrl: Optional[float] = None  # Relative to positive expression control
    lod_loq_status: str = "ABOVE_LOQ"  # 'BELOW_LOD', 'BETWEEN_LOD_LOQ', 'ABOVE_LOQ'
    measurement_source: str = "EMPIRICAL"  # 'EMPIRICAL', 'MOCK_SIMULATED', 'SYNTHETIC_PILOT'
    evidence_tier: str = "TIER_1_VALIDATED"  # 'TIER_1_VALIDATED', 'TEST_ONLY', 'EXPLORATORY'
    training_eligible: bool = True  # strictly False for mock data
    notes: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class PairedDBTLDataset:
    """Composite tamper-evident dataset joining in silico designs and wet-lab outcomes."""

    dataset_id: str
    schema_version: str = "0.1.0"
    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    constructs: Dict[str, ConstructDesignRecord] = field(default_factory=dict)
    experiments: Dict[str, ExperimentRunRecord] = field(default_factory=dict)
    samples: Dict[str, SampleRecord] = field(default_factory=dict)
    measurements: List[MeasurementRecord] = field(default_factory=list)
    derived_outcomes: List[DerivedOutcomeRecord] = field(default_factory=list)
    provenance: Dict[str, Any] = field(default_factory=dict)
    archive_sha256: Optional[str] = None

    def compute_sha256(self) -> str:
        """Computes SHA-256 digest over normalized JSON payload."""
        data_str = json.dumps(self.to_dict(include_sequence=False), sort_keys=True)
        return hashlib.sha256(data_str.encode("utf-8")).hexdigest()

    def to_dict(self, include_sequence: bool = True) -> Dict[str, Any]:
        return {
            "$schema": "https://eijex.com/schemas/factorforge/paired-dbtl-v0.1.json",
            "dataset_id": self.dataset_id,
            "schema_version": self.schema_version,
            "created_at": self.created_at,
            "archive_sha256": self.archive_sha256,
            "provenance": self.provenance,
            "constructs": {
                k: v.to_dict(include_sequence=include_sequence) for k, v in self.constructs.items()
            },
            "experiments": {k: v.to_dict() for k, v in self.experiments.items()},
            "samples": {k: v.to_dict() for k, v in self.samples.items()},
            "measurements": [m.to_dict() for m in self.measurements],
            "derived_outcomes": [o.to_dict() for o in self.derived_outcomes],
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> PairedDBTLDataset:
        """Reconstitute dataset from dictionary."""
        constructs = {}
        for cid, cdata in data.get("constructs", {}).items():
            tv = None
            if "trait_vector" in cdata and cdata["trait_vector"]:
                tvd = cdata["trait_vector"]
                tv = TraitVector(
                    cai_golden_set=tvd.get("cai_golden_set", 0.0),
                    global_gc_percent=tvd.get("global_gc_percent", 0.0),
                    initiation_mfe_kcal_mol=tvd.get("initiation_mfe_kcal_mol"),
                    local_50bp_gc_min=tvd.get("local_50bp_gc_min", 0.0),
                    local_50bp_gc_max=tvd.get("local_50bp_gc_max", 0.0),
                    homopolymer_max_run=tvd.get("homopolymer_max_run", 0),
                    synthesis_penalty=tvd.get("synthesis_penalty", 0.0),
                    codon_context_prior=tvd.get("codon_context_prior", 0.0),
                )
            constructs[cid] = ConstructDesignRecord(
                construct_id=cdata["construct_id"],
                target_name=cdata["target_name"],
                mature_protein_aa_length=cdata["mature_protein_aa_length"],
                construct_aa_length=cdata["construct_aa_length"],
                hypothesis_id=cdata["hypothesis_id"],
                hypothesis_name=cdata["hypothesis_name"],
                generation_contract=cdata["generation_contract"],
                sequence_digest=cdata["sequence_digest"],
                trait_vector=tv,
                sequence_dna=cdata.get("sequence_dna"),
                is_control=cdata.get("is_control", False),
                control_type=cdata.get("control_type"),
                rationale=cdata.get("rationale", ""),
                metadata=cdata.get("metadata", {}),
            )

        experiments = {
            eid: ExperimentRunRecord(**edata) for eid, edata in data.get("experiments", {}).items()
        }
        samples = {sid: SampleRecord(**sdata) for sid, sdata in data.get("samples", {}).items()}
        measurements = [MeasurementRecord(**mdata) for mdata in data.get("measurements", [])]
        derived_outcomes = [
            DerivedOutcomeRecord(**odata) for odata in data.get("derived_outcomes", [])
        ]

        return cls(
            dataset_id=data["dataset_id"],
            schema_version=data.get("schema_version", "0.1.0"),
            created_at=data.get("created_at", datetime.now(timezone.utc).isoformat()),
            constructs=constructs,
            experiments=experiments,
            samples=samples,
            measurements=measurements,
            derived_outcomes=derived_outcomes,
            provenance=data.get("provenance", {}),
            archive_sha256=data.get("archive_sha256"),
        )
