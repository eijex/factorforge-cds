"""Data models for declarative rules, enforcement levels, scope, and authority attribution."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Set


class EnforcementLevel(str, Enum):
    """Rule enforcement tiers in FactorForge & Eijex bio-compiler."""
    HARD_FAIL = "hard_fail"
    WARNING = "warning"
    IGNORE = "ignore"
    ERROR = "error"
    INFORMATIONAL = "informational"

class AuthorizedAction(str, Enum):
    BLOCK = "block"
    REGENERATE = "regenerate"
    REPORT_ONLY = "report_only"

class RuleCategory(str, Enum):
    ASSEMBLY = "assembly"
    RNA_RISK = "rna_risk"
    INTERNAL_POLICY = "internal_policy"
    REGULATORY = "regulatory_policy"

class AuthorityType(str, Enum):
    EIJEX_INTERNAL_POLICY = "EIJEX_INTERNAL_POLICY"
    EXTERNAL_REGULATORY_SOURCE = "EXTERNAL_REGULATORY_SOURCE"
    COMMUNITY_STANDARD = "COMMUNITY_STANDARD"

class EvaluationStage(str, Enum):
    INCREMENTAL = "incremental"
    FINAL_SEQUENCE = "final_sequence"

class EvidenceClass(str, Enum):
    SEQUENCE_SCAN = "sequence_scan"
    MODEL_PREDICTION = "model_prediction"
    GLOBAL_PROPERTY = "global_property"

@dataclass(frozen=True)
class RuleAuthority:
    authority_type: AuthorityType
    source_name: str
    agency: Optional[str] = None
    document_id: Optional[str] = None
    section: Optional[str] = None

@dataclass(frozen=True)
class RuleScope:
    assembly_methods: Set[str] = field(default_factory=lambda: {"*"})
    target_hosts: Set[str] = field(default_factory=lambda: {"*"})
    molecule_types: Set[str] = field(default_factory=lambda: {"dna", "mrna"})

    def matches(self, assembly_method: Optional[str] = None, host: Optional[str] = None) -> bool:
        if assembly_method and "*" not in self.assembly_methods and assembly_method not in self.assembly_methods:
            return False
        if host and "*" not in self.target_hosts and host not in self.target_hosts:
            return False
        return True

@dataclass
class RuleFinding:
    rule_id: str
    start: int
    end: int
    matched_sequence: str
    category: str
    detector_version: str
    evidence_class: str
    message: Optional[str] = None

@dataclass
class RuleDefinition:
    rule_id: str
    name: str
    description: str
    category: RuleCategory
    authority: RuleAuthority
    scope: RuleScope
    # Compatibility fallback for callers that do not provide an SOP profile.
    # Runtime policy may override this value; detectors remain policy-neutral.
    enforcement: EnforcementLevel = EnforcementLevel.IGNORE
    version: str = "1.0.0"
    evaluation_stage: EvaluationStage = EvaluationStage.INCREMENTAL
    evaluator_fn: Optional[Callable[[str, Dict[str, Any]], List[RuleFinding]]] = None

    def evaluate(self, sequence: str, context: Optional[Dict[str, Any]] = None) -> List[RuleFinding]:
        if self.evaluator_fn is None:
            return []
        return self.evaluator_fn(sequence, context or {})

