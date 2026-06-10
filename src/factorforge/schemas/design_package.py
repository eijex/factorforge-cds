"""Internal API response format for FactorForge CDS artifacts.

This module defines the *internal* Pydantic model used by the optimize API
handler (``api.optimize.handler``). It is NOT the public Open Bio Design
Package contract.

Public contract:
    ``src/factorforge/schemas/design_package.schema.json`` (JSON Schema Draft
    2020-12) is the public output specification. It defines the canonical field
    set (``design_id``, ``claim_boundary``, ``evidence``, etc.) and is tested
    by ``tests/test_design_package_schema.py`` and related files.

Separation rationale:
    The internal model uses API-convenient field names (``construct_id``,
    ``cds_design``) with ``extra="allow"`` for handler flexibility. The public
    schema enforces claim boundary fields and must be validated with jsonschema
    before any output is published or shared externally.
"""

from typing import Any, Optional

from pydantic import BaseModel, ConfigDict, Field


class TargetInfo(BaseModel):
    model_config = ConfigDict(extra="allow")

    name: Optional[str] = None
    uniprot_id: Optional[str] = None
    species: Optional[str] = None
    protein_class: Optional[str] = None


class ConstructPlan(BaseModel):
    model_config = ConfigDict(extra="allow")

    construct_type: Optional[str] = None
    tag: Optional[str] = None
    signal_peptide: Optional[str] = None
    kozak: bool = True


class CdsDesign(BaseModel):
    model_config = ConfigDict(extra="allow")

    engine: str
    product_version: str
    host_profile: str
    objective: str
    profile: Optional[str] = None
    input_length_aa: Optional[int] = None
    output_length_nt: Optional[int] = None
    cai: float
    gc_percent: float


class ConstraintReport(BaseModel):
    model_config = ConfigDict(extra="allow")

    restriction_sites_removed: list[str] = Field(default_factory=list)
    restriction_sites_unresolved: list[str] = Field(default_factory=list)
    polya_warnings: int = 0
    internal_stop_count: int = 0
    codon_rarity_clusters: int = 0
    aa_identity: float = 1.0


class ValidationStatus(BaseModel):
    model_config = ConfigDict(extra="allow")

    in_silico: str = "unchecked"
    aa_identity_check: str = "unchecked"
    gc_check: str = "unchecked"
    polya_check: str = "unchecked"
    moclo_check: str = "unchecked"


class Provenance(BaseModel):
    model_config = ConfigDict(extra="allow")

    input_sequence_hash: str
    output_cds_hash: str
    parameter_hash: str


class WetLabFeedback(BaseModel):
    model_config = ConfigDict(extra="allow")

    status: str = "pending"
    submissions: list[Any] = Field(default_factory=list)


class DesignPackage(BaseModel):
    model_config = ConfigDict(extra="allow")

    design_package_version: str = "1.0"
    construct_id: str
    created_at: str
    target: Optional[TargetInfo] = None
    construct_plan: Optional[ConstructPlan] = None
    cds_design: CdsDesign
    constraint_report: Optional[ConstraintReport] = None
    validation_status: Optional[ValidationStatus] = None
    provenance: Provenance
    wet_lab_feedback: WetLabFeedback = Field(default_factory=WetLabFeedback)
