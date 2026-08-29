from enum import Enum
from typing import Dict, Any, Optional
from pydantic import BaseModel, Field


class LifecycleState(str, Enum):
    RESEARCH = "research"
    CANDIDATE = "candidate"
    BENCHMARKING = "benchmarking"
    REVIEWED = "reviewed"
    APPROVED = "approved"
    PRODUCTION = "production"
    DEPRECATED = "deprecated"


class MethodReference(BaseModel):
    """External reference for a scientific method."""
    origin: str = Field(..., description="e.g., 'external_paper', 'internal_research'")
    reference: str = Field(..., description="DOI / PMID / URL / citation")
    implementation_status: str = Field(default="adapted", description="Status of code adaptation")
    implementation_commit: Optional[str] = Field(default=None, description="Git commit hash of implementation")


class MethodIdentity(BaseModel):
    """Scientific method identity tracked by AgentOS."""
    method_name: str
    method_version: str
    engine_name: str = Field(..., description="The executable engine from EngineRegistry")
    lifecycle_state: LifecycleState = LifecycleState.RESEARCH
    reference: Optional[MethodReference] = None
    description: str = ""


class MethodRegistry:
    """Scientific Method Registry managed by AgentOS.
    
    Separates scientific method identity and lifecycle from physical code execution (EngineRegistry).
    """
    _methods: Dict[str, MethodIdentity] = {}

    @classmethod
    def register(cls, method_id: str, identity: MethodIdentity) -> None:
        cls._methods[method_id] = identity

    @classmethod
    def get_method(cls, method_id: str) -> Optional[MethodIdentity]:
        return cls._methods.get(method_id)

    @classmethod
    def promote_method(cls, method_id: str, new_state: LifecycleState, approved_by: str) -> bool:
        """Promote a method to a new lifecycle state requiring approval."""
        method = cls.get_method(method_id)
        if not method:
            return False
        # In a real DB-backed version, this would log the approval audit trail
        method.lifecycle_state = new_state
        return True

    @classmethod
    def list_methods(cls) -> Dict[str, MethodIdentity]:
        return cls._methods.copy()
