"""Structure intake, annotation, and recommendation adapters for FactorForge."""

from factorforge.structure.alphafold_server_export import (
    AlphaFoldExportError,
    build_alphafold_server_job_payload,
    export_alphafold_server_job,
    parse_protein_fasta,
)
from factorforge.structure.alphafold_server_ingest import (
    AlphaFoldIngestError,
    extract_pause_positions_from_sidecar,
    format_ingest_report,
    ingest_alphafold_server_result,
    parse_af3_cif,
)
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
    "COLOR_PALETTE",
    "AlphaFoldExportError",
    "AlphaFoldIngestError",
    "CandidateRegion",
    "ESMRecommendationReport",
    "ESMStructureAnnotator",
    "InterProSidecarFetcher",
    "StructureRecommendationReport",
    "UniProtSidecarFetcher",
    "build_alphafold_server_job_payload",
    "build_structure_position_mapping",
    "export_alphafold_server_job",
    "extract_pause_positions_from_sidecar",
    "fetch_interpro_raw_record",
    "fetch_interpro_sidecar",
    "fetch_uniprot_raw_record",
    "fetch_uniprot_sidecar",
    "format_ingest_report",
    "ingest_alphafold_server_result",
    "parse_af3_cif",
    "parse_protein_fasta",
    "prepare_viewer_annotations",
    "recommend_esm_pause_regions",
    "recommend_structure_pause_regions",
]
