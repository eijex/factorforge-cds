"""Typed contracts, resource bounds, error codes, and backend capability gates for NB-Fold.

Implements PRD-JOB-398-NB-FOLD-v1.1.0 (Phase A/B local reviewed pipeline):
- Strict size bounds checked before reading full file contents.
- Machine-readable error codes and deterministic nonzero exit codes.
- Canonical SHA-256 digest helpers separating input/reference hashes from output artifact hashes.
- Fail-closed future backend capability stubs (ESMFold, ESM-2) without heavy ML dependencies.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

NBFOLD_PIPELINE_VERSION = "1.2.1"
ANNOTATION_SCHEMA_VERSION = "nbfold-annotation-v1"
REVIEW_SCHEMA_VERSION = "nbfold-review-v1"
OBSERVATIONS_SCHEMA_VERSION = "nbfold-observations-v1"
RUN_MANIFEST_SCHEMA_VERSION = "nbfold-run-manifest-v1"
ERROR_MANIFEST_SCHEMA_VERSION = "nbfold-error-v1"

DOMAIN_CALLER_VERSION = "1.0.1"
TM_SIGNAL_CALLER_VERSION = "1.0.1"
BACKBONE_SS_CALLER_VERSION = "1.0.1"

REVIEW_ACKNOWLEDGEMENT_TYPE = "LOCAL_USER_ACKNOWLEDGEMENT"

# Bounded resource limits (checked before reading full file contents)
MIN_PROTEIN_AA = 3
MAX_PROTEIN_AA = 300
MAX_FASTA_BYTES = 16 * 1024  # 16 KB
MAX_STRUCTURE_BYTES = 250_000  # 250,000 bytes (matches structure_web / auto_annotator)
MAX_SIDECAR_BYTES = 2 * 1024 * 1024  # 2 MB per sidecar
MAX_MANIFEST_JSON_BYTES = 256 * 1024  # 256 KB for annotation/review/sop JSON
MAX_TOTAL_INPUT_BYTES = 5 * 1024 * 1024  # 5 MB total across all input files
DEFAULT_TIMEOUT_SECONDS = 30.0

# Fixed output artifact names (never derived from FASTA headers)
ANNOTATION_FILENAME = "annotation.json"
PAIRED_FASTA_FILENAME = "paired_sequences.fasta"
OBSERVATIONS_FILENAME = "observations.json"
RUN_MANIFEST_FILENAME = "run_manifest.json"
REPORT_MARKDOWN_FILENAME = "NB_FOLD_REPORT.md"
ERROR_MANIFEST_FILENAME = "error_manifest.json"

ANNOTATE_ARTIFACT_NAMES = (
    ANNOTATION_FILENAME,
    ERROR_MANIFEST_FILENAME,
)

OPTIMIZE_ARTIFACT_NAMES = (
    PAIRED_FASTA_FILENAME,
    OBSERVATIONS_FILENAME,
    RUN_MANIFEST_FILENAME,
    REPORT_MARKDOWN_FILENAME,
    ERROR_MANIFEST_FILENAME,
)

# Machine-readable error_code -> numeric CLI exit code mapping
ERROR_EXIT_CODES: dict[str, int] = {
    "INVALID_FASTA": 2,
    "MULTI_RECORD_FASTA": 2,
    "INVALID_SEQUENCE": 2,
    "UNSUPPORTED_INPUT_SIZE": 2,
    "FORMAT_EXTENSION_MISMATCH": 2,
    "INVALID_STRUCTURE_FORMAT": 2,
    "CHAIN_NOT_FOUND": 2,
    "MODEL_NOT_FOUND": 2,
    "DISORDERED_RESIDUE": 2,
    "MODIFIED_OR_AMBIGUOUS_RESIDUE": 2,
    "SEQUENCE_MISMATCH": 2,
    "MISSING_CA_ATOM": 2,
    "INVALID_CONFIDENCE_VALUE": 2,
    "UNKNOWN_CONFIDENCE_SOURCE": 2,
    "SIDECAR_MISMATCH": 2,
    "UNSUPPORTED_HOST": 2,
    "REVIEW_MANIFEST_MISSING": 3,
    "REVIEW_NOT_CONFIRMED": 3,
    "REVIEW_MANIFEST_STALE": 3,
    "INVALID_REVIEW_MANIFEST": 3,
    "REVIEW_CONFLICTS_UNACKNOWLEDGED": 3,
    "INFEASIBLE": 4,
    "REQUIRED_CHECK_UNAVAILABLE": 5,
    "POLICY_BLOCKED": 5,
    "OBSERVATION_ENGINE_ERROR": 5,
    "UNSUPPORTED_BACKEND": 6,
    "OUTPUT_CONFLICT": 7,
    "EXECUTION_TIMEOUT": 8,
    "INTERNAL_ERROR": 1,
}


class NBFoldError(Exception):
    """Structured exception carrying a machine-readable error_code and sanitized message."""

    def __init__(
        self,
        error_code: str,
        message: str,
        *,
        stage: str = "nbfold",
        details: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(message)
        self.error_code = error_code if error_code in ERROR_EXIT_CODES else "INTERNAL_ERROR"
        self.exit_code = ERROR_EXIT_CODES[self.error_code]
        self.message = message
        self.stage = stage
        self.details = details or {}

    def to_manifest_dict(
        self,
        *,
        timestamp_utc: str,
        input_digests: dict[str, str | None] | None = None,
    ) -> dict[str, Any]:
        return {
            "schema_version": ERROR_MANIFEST_SCHEMA_VERSION,
            "pipeline_version": NBFOLD_PIPELINE_VERSION,
            "status": "ERROR",
            "stage": self.stage,
            "error_code": self.error_code,
            "exit_code": self.exit_code,
            "message": self.message,
            "details": self.details,
            "input_digests": input_digests or {},
            "timestamp_utc": timestamp_utc,
            "raw_input_logged": False,
        }


def canonical_json_sha256(payload: Any) -> str:
    """Compute deterministic SHA-256 hex digest over canonical JSON serialization."""
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def sha256_text(text: str) -> str:
    """Compute SHA-256 hex digest of UTF-8 text."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_bytes(data: bytes) -> str:
    """Compute SHA-256 hex digest of raw bytes."""
    return hashlib.sha256(data).hexdigest()


def check_file_size_limit(
    path: Path,
    *,
    max_bytes: int,
    label: str,
    empty_error_code: str = "INVALID_STRUCTURE_FORMAT",
    stage: str = "nbfold",
) -> int:
    """Validate file existence and byte size BEFORE reading file contents into memory."""
    if not isinstance(path, Path):
        path = Path(path)
    if not path.exists() or not path.is_file():
        raise NBFoldError(
            empty_error_code,
            f"{label} file does not exist: {path.name}",
            stage=stage,
        )
    size = path.stat().st_size
    if size == 0:
        raise NBFoldError(
            empty_error_code,
            f"{label} file is empty (0 bytes): {path.name}",
            stage=stage,
        )
    if size > max_bytes:
        raise NBFoldError(
            "UNSUPPORTED_INPUT_SIZE",
            f"{label} file size ({size} bytes) exceeds maximum allowed limit ({max_bytes} bytes).",
            stage=stage,
            details={"file_label": label, "observed_bytes": size, "max_bytes": max_bytes},
        )
    return size


def compute_config_digest(
    *,
    host: str,
    reference_id: str,
    codon_table_sha256: str,
    gc_min_percent: float,
    gc_max_percent: float,
    forbidden_motifs: list[str],
    target_adaptiveness: float,
    penalty_weight: float,
    sop_profile_hash: str | None = None,
) -> str:
    """Compute deterministic SHA-256 digest of resolved reference, policy, and optimizer settings."""
    payload = {
        "host": host,
        "reference_id": reference_id,
        "codon_table_sha256": codon_table_sha256,
        "gc_min_percent": round(float(gc_min_percent), 4),
        "gc_max_percent": round(float(gc_max_percent), 4),
        "forbidden_motifs": sorted(str(m).upper() for m in forbidden_motifs),
        "target_adaptiveness": round(float(target_adaptiveness), 6),
        "penalty_weight": round(float(penalty_weight), 6),
        "sop_profile_hash": sop_profile_hash,
    }
    return canonical_json_sha256(payload)


def compute_annotation_content_sha256(annotation_data: Mapping[str, Any]) -> str:
    """Compute canonical content hash of annotation.json excluding self-hash and helper template."""
    filtered = {
        k: v
        for k, v in annotation_data.items()
        if k not in {"annotation_content_sha256", "review_manifest_template"}
    }
    return canonical_json_sha256(filtered)


@dataclass(frozen=True)
class BackendCapability:
    """Capability descriptor separating file intake from future model inference backends."""

    backend_id: str
    role: str
    available: bool
    phase: str
    requires_external_gpu_or_weights: bool
    status_reason: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


BACKEND_CAPABILITIES: dict[str, BackendCapability] = {
    "alphafold_file": BackendCapability(
        backend_id="alphafold_file",
        role="local_structure_file_intake",
        available=True,
        phase="Phase_A_B",
        requires_external_gpu_or_weights=False,
        status_reason="Active local PDB/mmCIF structure file intake; does not execute structure folding.",
    ),
    "esmfold": BackendCapability(
        backend_id="esmfold",
        role="structure_prediction_inference",
        available=False,
        phase="Phase_D_deferred",
        requires_external_gpu_or_weights=True,
        status_reason="ESMFold structure inference is deferred to Phase D and not available in Phase A/B.",
    ),
    "esm2": BackendCapability(
        backend_id="esm2",
        role="protein_language_model_embedding",
        available=False,
        phase="Phase_D_deferred",
        requires_external_gpu_or_weights=True,
        status_reason="ESM-2 embedding inference is deferred to Phase D and not available in Phase A/B.",
    ),
}


def get_backend_capabilities() -> dict[str, dict[str, Any]]:
    """Return dictionary of backend capability statuses."""
    return {k: v.to_dict() for k, v in BACKEND_CAPABILITIES.items()}


def enforce_backend_gate(
    *,
    structure_backend: str = "alphafold_file",
    enable_esmfold: bool = False,
    enable_esm2: bool = False,
    stage: str = "nbfold",
) -> None:
    """Fail closed with UNSUPPORTED_BACKEND if any unimplemented backend or flag is requested."""
    if enable_esmfold:
        raise NBFoldError(
            "UNSUPPORTED_BACKEND",
            "ESMFold inference backend (enable_esmfold=True) is not implemented in Phase A/B "
            "(available=False). No synthetic or empty structure prediction is returned.",
            stage=stage,
            details={"requested_backend": "esmfold", "available": False},
        )
    if enable_esm2:
        raise NBFoldError(
            "UNSUPPORTED_BACKEND",
            "ESM-2 embedding backend (enable_esm2=True) is not implemented in Phase A/B "
            "(available=False). No synthetic or zero embedding is returned.",
            stage=stage,
            details={"requested_backend": "esm2", "available": False},
        )
    norm_backend = str(structure_backend or "").strip().lower()
    if norm_backend != "alphafold_file":
        raise NBFoldError(
            "UNSUPPORTED_BACKEND",
            f"Requested structure backend '{structure_backend}' is not supported in Phase A/B "
            "(only 'alphafold_file' local file intake is available).",
            stage=stage,
            details={"requested_backend": structure_backend, "available": False},
        )
