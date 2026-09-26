from typing import List, Optional, Dict, Any, Union
from pydantic import BaseModel, Field
from datetime import datetime

class FindingRecord(BaseModel):
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

class ConcordanceRecord(BaseModel):
    sequence_sha256: str
    finding_type: str
    start_pos: int
    end_pos: int
    strand: str
    concordant_predictors: List[str]
    discordant_predictors: List[str]
    consensus_score: Optional[float] = None
    metadata: Dict[str, Any] = Field(default_factory=dict)

class ComparisonRecord(BaseModel):
    retain_sequence_sha256: str
    resolved_sequence_sha256: str
    status: str  # RESOLVED, PERSISTENT, INTRODUCED
    finding_type: str
    start_pos: int
    end_pos: int
    strand: str
    predictor_name: str

class PolicyDecision(BaseModel):
    sequence_sha256: str
    policy_hash: str
    rule_id: str
    enforcement: str
    authorized_action: str
    decision: str  # BLOCK, WARN, PASS
    message: str