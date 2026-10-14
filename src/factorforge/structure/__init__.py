"""Structure intake, annotation, and recommendation adapters for FactorForge."""

from factorforge.structure.auto_annotator import (
    StructureRecommendationReport,
    recommend_structure_pause_regions,
)

__all__ = [
    "StructureRecommendationReport",
    "recommend_structure_pause_regions",
]
