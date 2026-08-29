from datetime import datetime
from typing import Dict, Any, List, Optional
from pydantic import BaseModel, Field

from factorforge.evaluation.models import EvaluationResult


class BenchmarkProvenance(BaseModel):
    """Immutable provenance tracking for benchmark runs."""
    input_sequence_hash: str
    input_type: str = Field(description="protein or CDS")
    terminal_stop_policy: str
    codon_reference_id: str
    constraint_profile_id: str
    evaluator_version: str
    method_version: str
    git_commit: str
    seed: Optional[int] = None
    inference_mode: str
    runtime_environment: str
    timestamp: str = Field(default_factory=lambda: datetime.utcnow().isoformat())


class BenchmarkTarget(BaseModel):
    """A target sequence for benchmarking."""
    target_id: str
    protein_sequence: str
    metadata: Dict[str, Any] = Field(default_factory=dict)


class BenchmarkRunConfig(BaseModel):
    """Configuration for a benchmark execution."""
    suite_name: str = Field(..., description="e.g., Engineering Benchmark Suite")
    host: str = "nbenthamiana"
    target_gc_min_percent: float
    target_gc_max_percent: float
    forbidden_type_iis: List[str] = Field(default_factory=list)
    terminal_stop_policy: str = "append"
    seed: int = 42


class EngineBenchmarkResult(BaseModel):
    """Benchmark result for a single engine on a single target."""
    engine_name: str
    target_id: str
    provenance: BenchmarkProvenance
    optimized_sequence: str
    evaluation: EvaluationResult
    runtime_seconds: float


class BenchmarkSuiteReport(BaseModel):
    """Aggregated benchmark report across all targets and engines."""
    suite_name: str
    run_timestamp: str = Field(default_factory=lambda: datetime.utcnow().isoformat())
    config: BenchmarkRunConfig
    results: List[EngineBenchmarkResult]
    
    # Summary statistics mapped by engine
    pass_counts: Dict[str, str] = Field(description="Raw pass counts e.g., '3/4'")
    average_cai: Dict[str, Optional[float]]
    average_gc: Dict[str, float]
