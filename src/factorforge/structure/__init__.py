"""Structure intake, annotation, and recommendation adapters for FactorForge."""

from factorforge.structure.auto_annotator import (
    CandidateRegion,
    StructureRecommendationReport,
    recommend_structure_pause_regions,
)
from factorforge.structure.esm_annotator import (
    ESMRecommendationReport,
    ESMStructureAnnotator,
    recommend_esm_pause_regions,
)
from factorforge.structure.interpro_sidecar_fetcher import (
    InterProSidecarFetcher,
    fetch_interpro_raw_record,
    fetch_interpro_sidecar,
)
from factorforge.structure.uniprot_sidecar_fetcher import (
    UniProtSidecarFetcher,
    fetch_uniprot_raw_record,
    fetch_uniprot_sidecar,
)
from factorforge.structure.viewer_bridge import (
    COLOR_PALETTE,
    build_structure_position_mapping,
    prepare_viewer_annotations,
)

__all__ = [
    "CandidateRegion",
    "StructureRecommendationReport",
    "recommend_structure_pause_regions",
    "ESMRecommendationReport",
    "ESMStructureAnnotator",
    "recommend_esm_pause_regions",
    "UniProtSidecarFetcher",
    "fetch_uniprot_raw_record",
    "fetch_uniprot_sidecar",
    "InterProSidecarFetcher",
    "fetch_interpro_raw_record",
    "fetch_interpro_sidecar",
    "COLOR_PALETTE",
    "build_structure_position_mapping",
    "prepare_viewer_annotations",
]
