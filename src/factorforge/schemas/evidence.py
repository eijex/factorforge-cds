from typing import List, Optional, Dict, Any, Literal
from pydantic import BaseModel, Field, ConfigDict
from datetime import datetime

class FindingRecord(BaseModel):
    model_config = ConfigDict(extra='forbid')
    finding_id: str
    rule_id: str
    predictor_name: str
    predictor_version: str
    finding_type: str
    start_pos: int
    end_pos: int
    strand: str
    score: Optional[float] = None
    metadata: Dict[str, Any] = Field(default_factory=dict)

class PredictorRunRecord(BaseModel):
    model_config = ConfigDict(extra='forbid')
    sequence_sha256: str
    predictor_name: str
    predictor_version: str
    predictor_model_context: str
    target_host: str
    run_timestamp: str
    raw_output_sha256: str
    coordinate_system: str
    normalization_version: str
    findings: List[FindingRecord] = Field(default_factory=list)

class ComparisonRecord(BaseModel):
    model_config = ConfigDict(extra='forbid')
    retain_sequence_sha256: str
    resolved_sequence_sha256: str
    reference_finding_id: str
    optimized_finding_id: Optional[str]
    status: Literal['RESOLVED', 'PERSISTENT', 'INTRODUCED']
    match_rule_version: str
    match_distance_nt: int
    match_basis: str
    finding_type: str
    start_pos: int
    end_pos: int
    strand: str
    predictor_name: str

class ConcordanceRecord(BaseModel):
    model_config = ConfigDict(extra='forbid')
    sequence_sha256: str
    finding_type: str
    start_pos: int
    end_pos: int
    strand: str
    concordant_predictors: List[str]
    discordant_predictors: List[str]
    consensus_score: Optional[float] = None
    metadata: Dict[str, Any] = Field(default_factory=dict)

class PolicyDecision(BaseModel):
    model_config = ConfigDict(extra='forbid')
    finding_id: str
    enforcement_level: Literal['HARD_FAIL', 'WARNING', 'IGNORE']
    authorized_action: Literal['BLOCK', 'REPORT_ONLY', 'REGENERATE', 'NONE']
    rule_id: str
    policy_profile_id: str
    policy_version: str
    policy_digest: str
    authorization_source: str