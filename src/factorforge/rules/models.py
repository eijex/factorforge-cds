# factorforge/src/factorforge/rules/models.py
"""Data models for declarative rules, enforcement levels, scope, and authority attribution."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Set


class EnforcementLevel(str, Enum):
    """Rule enforcement tiers in FactorForge & Eijex bio-compiler."""
    HARD_FAIL = "hard_fail"       # Invariant that must never be violated (e.g. Type IIS restriction site)
    WARNING = "warning"           # High-risk motif (e.g. cryptic splice site, strong hairpin)
    INFORMATIONAL = "informational" # Policy recommendation or heuristic indicator


class RuleCategory(str, Enum):
    """Categorical classification of sequence constraints."""
    ASSEMBLY = "assembly"                 # Physical cloning and assembly constraints
    RNA_RISK = "rna_risk"                 # Structural/translational stability risks
    INTERNAL_POLICY = "internal_policy"   # Eijex engineering policies
    REGULATORY = "regulatory_policy"      # External regulatory agency references


class AuthorityType(str, Enum):
    """Origin and authority attribution of the rule."""
    EIJEX_INTERNAL_POLICY = "EIJEX_INTERNAL_POLICY"
    EXTERNAL_REGULATORY_SOURCE = "EXTERNAL_REGULATORY_SOURCE"
    COMMUNITY_STANDARD = "COMMUNITY_STANDARD"


@dataclass(frozen=True)
class RuleAuthority:
    """Explicit attribution of rule authority."""
    authority_type: AuthorityType
    source_name: str
    agency: Optional[str] = None
    document_id: Optional[str] = None
    section: Optional[str] = None


@dataclass(frozen=True)
class RuleScope:
    """Scope of applicability for the rule."""
    assembly_methods: Set[str] = field(default_factory=lambda: {"*"})
    target_hosts: Set[str] = field(default_factory=lambda: {"*"})
    molecule_types: Set[str] = field(default_factory=lambda: {"dna", "mrna"})

    def matches(self, assembly_method: Optional[str] = None, host: Optional[str] = None) -> bool:
        """Check if scope matches a given execution context."""
        if assembly_method and "*" not in self.assembly_methods and assembly_method not in self.assembly_methods:
            return False
        if host and "*" not in self.target_hosts and host not in self.target_hosts:
            return False
        return True


@dataclass
class RuleDefinition:
    """Declarative definition of a sequence validation and optimization rule."""
    rule_id: str
    name: str
    description: str
    category: RuleCategory
    enforcement: EnforcementLevel
    authority: RuleAuthority
    scope: RuleScope
    version: str = "1.0.0"
    evaluator_fn: Optional[Callable[[str, Dict[str, Any]], Dict[str, Any]]] = None

    def evaluate(self, sequence: str, context: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Evaluate rule on candidate sequence."""
        if self.evaluator_fn is None:
            return {"rule_id": self.rule_id, "passed": True, "details": "No evaluator attached"}
        return self.evaluator_fn(sequence, context or {})
