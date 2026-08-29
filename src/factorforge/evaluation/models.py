"""Data models for Shared Evaluation Contract.

This module defines the canonical evaluation structure separating Metrics (raw values)
from Checks (policy application).
"""
from enum import Enum
from typing import List, Optional
from pydantic import BaseModel, Field


class CheckDomain(str, Enum):
    SEQUENCE_INTEGRITY = "sequence_integrity"
    CODON_COMPOSITION = "codon_composition"
    ASSEMBLY = "assembly"
    RNA_RISK = "rna_risk"


class CheckEnforcement(str, Enum):
    HARD_FAIL = "hard_fail"
    GATE = "gate"
    WARNING = "warning"
    INFORMATIONAL = "informational"


class CheckResultValue(str, Enum):
    PASS = "pass"
    FAIL = "fail"
    WARNING = "warning"
    INDETERMINATE = "indeterminate"


class CheckResult(BaseModel):
    check_name: str
    domain: CheckDomain
    enforcement: CheckEnforcement
    result: CheckResultValue
    message: Optional[str] = None


class SequenceIntegrity(BaseModel):
    aa_identity: float
    frame_valid: bool
    internal_stop_count: int


class Metrics(BaseModel):
    cai: Optional[float] = None
    gc_percent: Optional[float] = None
    mfe: Optional[float] = None
    mfe_status: Optional[str] = None
    mfe_reason: Optional[str] = None
    mfe_warning: Optional[str] = None


class EvaluationResult(BaseModel):
    candidate_id: str
    evaluator_version: str = "1.0.0"
    
    sequence_integrity: SequenceIntegrity
    metrics: Metrics
    checks: List[CheckResult] = Field(default_factory=list)
    
    @property
    def passed(self) -> bool:
        """Returns True if no HARD_FAIL or GATE checks returned FAIL."""
        for check in self.checks:
            if check.enforcement in (CheckEnforcement.HARD_FAIL, CheckEnforcement.GATE):
                if check.result == CheckResultValue.FAIL:
                    return False
        return True
