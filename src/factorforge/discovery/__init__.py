"""FactorForge Discovery Module.

Implements the multi-contract candidate generation, hard constraint filtering,
multi-dimensional trait extraction, Pareto filtering, and diversity-aware
Top-K slate selection (Job 283A).
"""

from factorforge.discovery.acquisition import (
    ConstructDesignRecord,
    DerivedOutcomeRecord,
    ExperimentRunRecord,
    MeasurementRecord,
    PairedDBTLDataset,
    SampleRecord,
)
from factorforge.discovery.acquisition_logger import AcquisitionLogger
from factorforge.discovery.panel_builder import (
    OrthogonalityGate,
    ProspectivePanelBuilder,
    STANDARD_TARGET_PANEL,
)
from factorforge.discovery.schemas import (
    CandidateSlate,
    SlateCandidate,
    SlateSummary,
    TargetMetadata,
    TraitVector,
)
from factorforge.discovery.slate import DiscoverySlateEngine

__all__ = [
    "AcquisitionLogger",
    "CandidateSlate",
    "ConstructDesignRecord",
    "DerivedOutcomeRecord",
    "DiscoverySlateEngine",
    "ExperimentRunRecord",
    "MeasurementRecord",
    "OrthogonalityGate",
    "PairedDBTLDataset",
    "ProspectivePanelBuilder",
    "SampleRecord",
    "SlateCandidate",
    "SlateSummary",
    "STANDARD_TARGET_PANEL",
    "TargetMetadata",
    "TraitVector",
]
