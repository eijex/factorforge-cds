"""NB-Fold Phase A/B Orchestration Pipeline, Review Verification, and Atomic Artifact Export.

Implements PRD-JOB-398-NB-FOLD-v1.1.0 Sections 3, 5, 6, 7, 8, 9:
- Phase A (`run_nbfold_annotate`): StructureBridge intake -> `annotation.json` (`reviewed=False`). Never optimizes.
- Phase B (`run_nbfold_optimize`): Explicit `review.json` verification -> paired 1DP (`DPV2Optimizer`) and
  3DP (`StructuralPauseDPOptimizer`) under identical resolved host reference/GC/motif constraints ->
  independent hard QC, sense-CDS plant RNA scan, and 48-nt initiation MFE -> atomic export of
  `paired_sequences.fasta`, `observations.json`, `run_manifest.json`, and `NB_FOLD_REPORT.md`.
"""

from __future__ import annotations

import json
import math
import os
import shutil
import tempfile
import time
import tracemalloc
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from factorforge.analysis.metrics import load_codon_usage_table, translate_dna
from factorforge.constraints.type_iis import get_canonical_forbidden_motifs
from factorforge.engines.profile.scoring import resolve_host_gc_range
from factorforge.guards.initiation_mfe import (
    DEFAULT_INITIATION_WINDOW_NT,
    DEFAULT_TEMPERATURE_C,
    evaluate_initiation_mfe,
)
from factorforge.guards.plant_rna_guards import (
    DETECTOR_VERSION as RNA_DETECTOR_VERSION,
)
from factorforge.guards.plant_rna_guards import (
    TRUTHFULNESS_NOTICE as RNA_TRUTHFULNESS_NOTICE,
)
from factorforge.guards.plant_rna_guards import (
    compare_candidate_plant_rna_risks,
    scan_plant_rna_risks,
    validate_candidate_hard_invariants,
)
from factorforge.nbfold.annotations import (
    evaluate_positions_against_structural_annotations,
)
from factorforge.nbfold.bridge import (
    build_structure_annotation,
    parse_protein_fasta_file,
)
from factorforge.nbfold.contracts import (
    ANNOTATE_ARTIFACT_NAMES,
    ANNOTATION_FILENAME,
    ANNOTATION_SCHEMA_VERSION,
    DEFAULT_TIMEOUT_SECONDS,
    ERROR_MANIFEST_FILENAME,
    MAX_MANIFEST_JSON_BYTES,
    MAX_SIDECAR_BYTES,
    MAX_TOTAL_INPUT_BYTES,
    NBFOLD_PIPELINE_VERSION,
    OBSERVATIONS_FILENAME,
    OBSERVATIONS_SCHEMA_VERSION,
    OPTIMIZE_ARTIFACT_NAMES,
    PAIRED_FASTA_FILENAME,
    REPORT_MARKDOWN_FILENAME,
    REVIEW_ACKNOWLEDGEMENT_TYPE,
    REVIEW_SCHEMA_VERSION,
    RUN_MANIFEST_FILENAME,
    RUN_MANIFEST_SCHEMA_VERSION,
    NBFoldError,
    canonical_json_sha256,
    check_file_size_limit,
    compute_annotation_content_sha256,
    compute_config_digest,
    enforce_backend_gate,
    get_backend_capabilities,
    sha256_bytes,
    sha256_text,
)
from factorforge.registry.versioning import engine_version, product_version
from factorforge.rules.registry import RuleRegistry
from factorforge.sops.config import SopProfile
from factorforge.structure.auto_annotator import RECOMMENDER_VERSION
from factorforge.validation.tradeoff_suite import compute_tradeoff_comparison

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[3]


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _reference_policy_manifest_path() -> Path:
    repo_manifest = REPO_ROOT / "data" / "reference" / "reference_policy_manifest.json"
    if repo_manifest.exists():
        return repo_manifest
    return PACKAGE_ROOT / "data" / "reference" / "reference_policy_manifest.json"


def _resolve_packaged_or_repo_path(relative_path: str) -> Path:
    repo_path = REPO_ROOT / relative_path
    if repo_path.exists():
        return repo_path
    prefix = "src/factorforge/"
    if relative_path.startswith(prefix):
        return PACKAGE_ROOT / relative_path.removeprefix(prefix)
    return PACKAGE_ROOT / relative_path


def resolve_nbfold_configuration(
    *,
    host: str = "nbenthamiana",
    gc_min: float | None = None,
    gc_max: float | None = None,
    target_adaptiveness: float = 0.40,
    penalty_weight: float = 2.0,
    custom_forbidden_motifs: list[str] | None = None,
    sop_profile_path: Path | None = None,
    sop_profile_dict: Mapping[str, Any] | None = None,
    stage: str = "nbfold",
) -> dict[str, Any]:
    """Resolve and checksum-verify host codon reference, GC bounds, motif registry, and SOP policy."""
    if not isinstance(host, str) or host.strip().lower() != "nbenthamiana":
        raise NBFoldError(
            "UNSUPPORTED_HOST",
            f"Unsupported host '{host}'. NB-Fold Phase A/B paired 1DP/3DP supports only 'nbenthamiana'.",
            stage=stage,
            details={"requested_host": host, "supported_hosts": ["nbenthamiana"]},
        )
    internal_host = "nbenthamiana"

    if (
        isinstance(target_adaptiveness, bool)
        or not isinstance(target_adaptiveness, (int, float))
        or not math.isfinite(float(target_adaptiveness))
        or not (0.0 < float(target_adaptiveness) <= 1.0)
    ):
        raise NBFoldError(
            "INVALID_REVIEW_MANIFEST",
            "target_adaptiveness must be a finite number in (0.0, 1.0].",
            stage=stage,
        )
    if (
        isinstance(penalty_weight, bool)
        or not isinstance(penalty_weight, (int, float))
        or not math.isfinite(float(penalty_weight))
        or not (0.0 <= float(penalty_weight) <= 100.0)
    ):
        raise NBFoldError(
            "INVALID_REVIEW_MANIFEST",
            "penalty_weight must be a finite number in [0.0, 100.0].",
            stage=stage,
        )

    default_gc_min, default_gc_max = resolve_host_gc_range(internal_host)
    raw_gc_min = default_gc_min if gc_min is None else gc_min
    raw_gc_max = default_gc_max if gc_max is None else gc_max

    if (
        isinstance(raw_gc_min, bool)
        or isinstance(raw_gc_max, bool)
        or not isinstance(raw_gc_min, (int, float))
        or not isinstance(raw_gc_max, (int, float))
        or not (math.isfinite(float(raw_gc_min)) and math.isfinite(float(raw_gc_max)))
    ):
        raise NBFoldError(
            "INFEASIBLE",
            "GC bounds (gc_min, gc_max) must be finite numbers.",
            stage=stage,
        )

    gc_min_val = float(raw_gc_min)
    gc_max_val = float(raw_gc_max)
    # Normalize fractional (0, 1] input to percentage if both bounds are <= 1.0 and > 0
    if 0.0 <= gc_min_val <= 1.0 and 0.0 < gc_max_val <= 1.0:
        gc_min_pct = round(gc_min_val * 100.0, 4)
        gc_max_pct = round(gc_max_val * 100.0, 4)
    else:
        gc_min_pct = round(gc_min_val, 4)
        gc_max_pct = round(gc_max_val, 4)

    if not (0.0 <= gc_min_pct <= gc_max_pct <= 100.0):
        raise NBFoldError(
            "INFEASIBLE",
            f"Invalid GC percentage range [{gc_min_pct}, {gc_max_pct}]; requires 0 <= gc_min <= gc_max <= 100.",
            stage=stage,
        )

    manifest_path = _reference_policy_manifest_path()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    active_ref_id = str(manifest["active_default_reference_id"])
    ref_entry = next(
        (r for r in manifest.get("references", []) if r.get("reference_id") == active_ref_id),
        None,
    )
    if ref_entry is None:
        raise NBFoldError(
            "INTERNAL_ERROR",
            f"Active default reference '{active_ref_id}' not found in reference policy manifest.",
            stage=stage,
        )

    codon_table_path = _resolve_packaged_or_repo_path(str(ref_entry["codon_table_path"]))
    actual_table_sha256 = sha256_bytes(codon_table_path.read_bytes())
    expected_table_sha256 = str(ref_entry["checksum_sha256"])
    if actual_table_sha256 != expected_table_sha256:
        raise NBFoldError(
            "INTERNAL_ERROR",
            f"Codon reference table checksum mismatch for {active_ref_id}.",
            stage=stage,
        )

    table = load_codon_usage_table(path=codon_table_path)
    codon_weights_sha256 = canonical_json_sha256(table.codon_weights)

    motifs = list(get_canonical_forbidden_motifs())
    if custom_forbidden_motifs:
        for m in custom_forbidden_motifs:
            clean_m = "".join(str(m).upper().split())
            if not clean_m or set(clean_m) - set("ACGT"):
                raise NBFoldError(
                    "INVALID_SEQUENCE",
                    "Custom forbidden motif must contain only unambiguous A/C/G/T bases.",
                    stage=stage,
                )
            if clean_m not in motifs:
                motifs.append(clean_m)
    motifs = sorted(set(motifs))

    sop_profile: SopProfile | None = None
    if sop_profile_path is not None and sop_profile_dict is not None:
        raise NBFoldError(
            "INVALID_REVIEW_MANIFEST",
            "Specify at most one of sop_profile_path or sop_profile_dict.",
            stage=stage,
        )
    if sop_profile_path is not None:
        check_file_size_limit(
            Path(sop_profile_path),
            max_bytes=MAX_MANIFEST_JSON_BYTES,
            label="SOP profile",
            empty_error_code="INVALID_REVIEW_MANIFEST",
            stage=stage,
        )
        try:
            raw_sop = json.loads(Path(sop_profile_path).read_text(encoding="utf-8"))
            known_ids = [r.rule_id for r in RuleRegistry().list_rules()]
            sop_profile = SopProfile.from_dict(raw_sop, known_rule_ids=known_ids)
        except Exception as exc:
            raise NBFoldError(
                "INVALID_REVIEW_MANIFEST",
                f"Invalid SOP profile JSON: {exc}",
                stage=stage,
            ) from exc
    elif sop_profile_dict is not None:
        try:
            known_ids = [r.rule_id for r in RuleRegistry().list_rules()]
            sop_profile = SopProfile.from_dict(sop_profile_dict, known_rule_ids=known_ids)
        except Exception as exc:
            raise NBFoldError(
                "INVALID_REVIEW_MANIFEST",
                f"Invalid SOP profile object: {exc}",
                stage=stage,
            ) from exc

    sop_hash = sop_profile.digest if sop_profile is not None else None
    config_digest = compute_config_digest(
        host=internal_host,
        reference_id=active_ref_id,
        codon_table_sha256=actual_table_sha256,
        gc_min_percent=gc_min_pct,
        gc_max_percent=gc_max_pct,
        forbidden_motifs=motifs,
        target_adaptiveness=float(target_adaptiveness),
        penalty_weight=float(penalty_weight),
        sop_profile_hash=sop_hash,
    )

    return {
        "host": internal_host,
        "reference_id": active_ref_id,
        "reference_tier": str(ref_entry.get("tier", "production_enabled")),
        "reference_policy_version": str(manifest.get("policy_version", "")),
        "reference_claim_boundary": str(ref_entry.get("claim_boundary", "")),
        "codon_table_path": str(codon_table_path),
        "codon_table_sha256": actual_table_sha256,
        "codon_weights": dict(table.codon_weights),
        "codon_weights_sha256": codon_weights_sha256,
        "gc_min_percent": gc_min_pct,
        "gc_max_percent": gc_max_pct,
        "forbidden_motifs": motifs,
        "target_adaptiveness": float(target_adaptiveness),
        "penalty_weight": float(penalty_weight),
        "sop_profile": sop_profile,
        "sop_profile_hash": sop_hash,
        "config_digest": config_digest,
    }


def _check_output_dir_no_conflict(
    output_dir: Path, artifact_names: tuple[str, ...], *, stage: str
) -> None:
    """Refuse to overwrite existing output artifacts by default."""
    out = Path(output_dir)
    if out.exists():
        if not out.is_dir():
            raise NBFoldError(
                "OUTPUT_CONFLICT",
                f"Output path '{out.name}' exists and is not a directory.",
                stage=stage,
            )
        existing = [name for name in artifact_names if (out / name).exists()]
        if existing:
            raise NBFoldError(
                "OUTPUT_CONFLICT",
                f"Refusing to overwrite existing output artifact(s) in '{out.name}': {', '.join(existing)}.",
                stage=stage,
                details={"existing_artifacts": existing},
            )


def write_error_manifest(
    output_dir: Path,
    error: NBFoldError,
    *,
    input_digests: dict[str, str | None] | None = None,
) -> Path | None:
    """Atomically write sanitized error_manifest.json if output_dir has no OUTPUT_CONFLICT."""
    if error.error_code == "OUTPUT_CONFLICT":
        return None
    out = Path(output_dir)
    try:
        out.mkdir(parents=True, exist_ok=True)
        err_target = out / ERROR_MANIFEST_FILENAME
        if err_target.exists():
            return None
        payload = error.to_manifest_dict(
            timestamp_utc=_utc_now_iso(),
            input_digests=input_digests,
        )
        fd, tmp_name = tempfile.mkstemp(prefix=".nbfold_err_", suffix=".tmp", dir=str(out))
        try:
            with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
                json.dump(payload, handle, indent=2, sort_keys=False)
                handle.write("\n")
            os.replace(tmp_name, err_target)
            return err_target
        except (OSError, ValueError, TypeError):
            if os.path.exists(tmp_name):
                os.remove(tmp_name)
            return None
    except (OSError, ValueError, TypeError):
        return None


def verify_review_manifest(
    *,
    review_manifest_path: Path | None,
    annotation_path: Path | None,
    current_annotation: Mapping[str, Any],
    resolved_config: Mapping[str, Any],
    stage: str = "optimize",
) -> dict[str, Any]:
    """Verify explicit user review manifest (review.json) against annotation.json and live inputs."""
    if review_manifest_path is None:
        raise NBFoldError(
            "REVIEW_MANIFEST_MISSING",
            "Optimization requires an explicit --review-manifest (review.json).",
            stage=stage,
        )
    rev_path = Path(review_manifest_path)
    if not rev_path.exists() or not rev_path.is_file():
        raise NBFoldError(
            "REVIEW_MANIFEST_MISSING",
            f"Review manifest file does not exist: {rev_path.name}",
            stage=stage,
        )
    check_file_size_limit(
        rev_path,
        max_bytes=MAX_MANIFEST_JSON_BYTES,
        label="Review manifest",
        empty_error_code="REVIEW_MANIFEST_MISSING",
        stage=stage,
    )

    if annotation_path is None:
        raise NBFoldError(
            "REVIEW_MANIFEST_MISSING",
            "Optimization requires an explicit --annotation (annotation.json).",
            stage=stage,
        )
    ann_path = Path(annotation_path)
    if not ann_path.exists() or not ann_path.is_file():
        raise NBFoldError(
            "REVIEW_MANIFEST_MISSING",
            f"Annotation manifest file does not exist: {ann_path.name}",
            stage=stage,
        )
    check_file_size_limit(
        ann_path,
        max_bytes=MAX_MANIFEST_JSON_BYTES,
        label="Annotation manifest",
        empty_error_code="REVIEW_MANIFEST_MISSING",
        stage=stage,
    )

    ann_raw_bytes = ann_path.read_bytes()
    ann_file_sha256 = sha256_bytes(ann_raw_bytes)
    try:
        loaded_ann = json.loads(ann_raw_bytes.decode("utf-8"))
    except Exception as exc:
        raise NBFoldError(
            "INVALID_REVIEW_MANIFEST",
            "Annotation file (annotation.json) is not valid UTF-8 JSON.",
            stage=stage,
        ) from exc

    if not isinstance(loaded_ann, dict) or loaded_ann.get("schema_version") != ANNOTATION_SCHEMA_VERSION:
        raise NBFoldError(
            "INVALID_REVIEW_MANIFEST",
            f"Annotation file must be a JSON object with schema_version='{ANNOTATION_SCHEMA_VERSION}'.",
            stage=stage,
        )
    if loaded_ann.get("reviewed") is not False:
        raise NBFoldError(
            "INVALID_REVIEW_MANIFEST",
            "annotation.json must retain reviewed=false; approval is recorded only in review.json.",
            stage=stage,
        )

    recomputed_ann_hash = compute_annotation_content_sha256(loaded_ann)
    recorded_ann_hash = loaded_ann.get("annotation_content_sha256")
    if not isinstance(recorded_ann_hash, str) or recorded_ann_hash != recomputed_ann_hash:
        raise NBFoldError(
            "REVIEW_MANIFEST_STALE",
            "annotation.json content was modified after creation (annotation_content_sha256 mismatch).",
            stage=stage,
        )

    # Cross-check loaded annotation.json against current live inputs
    for field_name in (
        "protein_sha256",
        "structure_sha256",
        "sidecar_sha256s",
        "chain_id",
        "model_number",
        "format",
        "provider_family",
        "confidence_metric",
        "source_basis",
        "export_metadata_sha256",
        "export_version",
    ):
        if loaded_ann.get(field_name) != current_annotation.get(field_name):
            raise NBFoldError(
                "REVIEW_MANIFEST_STALE",
                f"annotation.json field '{field_name}' does not match current input files/parameters.",
                stage=stage,
                details={"mismatched_field": field_name},
            )

    rev_raw_bytes = rev_path.read_bytes()
    rev_file_sha256 = sha256_bytes(rev_raw_bytes)
    try:
        review_data = json.loads(rev_raw_bytes.decode("utf-8"))
    except Exception as exc:
        raise NBFoldError(
            "INVALID_REVIEW_MANIFEST",
            "Review manifest (review.json) is not valid UTF-8 JSON.",
            stage=stage,
        ) from exc

    if not isinstance(review_data, dict):
        raise NBFoldError(
            "INVALID_REVIEW_MANIFEST",
            "Review manifest must be a JSON object.",
            stage=stage,
        )

    if review_data.get("schema_version") != REVIEW_SCHEMA_VERSION:
        raise NBFoldError(
            "INVALID_REVIEW_MANIFEST",
            f"Review manifest schema_version must be '{REVIEW_SCHEMA_VERSION}'.",
            stage=stage,
        )

    if review_data.get("review_type") != REVIEW_ACKNOWLEDGEMENT_TYPE:
        raise NBFoldError(
            "INVALID_REVIEW_MANIFEST",
            f"Review manifest review_type must be '{REVIEW_ACKNOWLEDGEMENT_TYPE}'. "
            "Cryptographic signature or server approval claims are not permitted.",
            stage=stage,
        )

    pos_reviewed = review_data.get("positions_reviewed")
    if type(pos_reviewed) is not bool or pos_reviewed is not True:
        raise NBFoldError(
            "REVIEW_NOT_CONFIRMED",
            "Review manifest must explicitly set positions_reviewed=true.",
            stage=stage,
        )

    reviewed_at = review_data.get("reviewed_at")
    if not isinstance(reviewed_at, str) or not reviewed_at.strip():
        raise NBFoldError(
            "INVALID_REVIEW_MANIFEST",
            "Review manifest requires a non-empty 'reviewed_at' timestamp string.",
            stage=stage,
        )

    reviewer_label = review_data.get("reviewer_label")
    if reviewer_label is not None and (
        not isinstance(reviewer_label, str) or not reviewer_label.strip()
    ):
        raise NBFoldError(
            "INVALID_REVIEW_MANIFEST",
            "Optional reviewer_label must be a non-empty string or null.",
            stage=stage,
        )

    # Verify all hashes and parameters against live inputs and annotation.json
    if review_data.get("protein_sha256") != current_annotation["protein_sha256"]:
        raise NBFoldError(
            "REVIEW_MANIFEST_STALE",
            "review.json protein_sha256 does not match the normalized input protein sequence.",
            stage=stage,
        )
    if review_data.get("structure_sha256") != current_annotation["structure_sha256"]:
        raise NBFoldError(
            "REVIEW_MANIFEST_STALE",
            "review.json structure_sha256 does not match the input structure file.",
            stage=stage,
        )
    if review_data.get("sidecar_sha256s") != current_annotation["sidecar_sha256s"]:
        raise NBFoldError(
            "REVIEW_MANIFEST_STALE",
            "review.json sidecar_sha256s do not match the current sidecar inputs.",
            stage=stage,
        )
    if (
        review_data.get("chain_id") != current_annotation["chain_id"]
        or type(review_data.get("model_number")) is not int
        or review_data.get("model_number") != current_annotation["model_number"]
    ):
        raise NBFoldError(
            "REVIEW_MANIFEST_STALE",
            "review.json chain_id or model_number does not match the selected structure chain/model.",
            stage=stage,
        )

    rev_ann_hash = review_data.get("annotation_sha256")
    if not isinstance(rev_ann_hash, str) or rev_ann_hash not in {
        recomputed_ann_hash,
        ann_file_sha256,
    }:
        raise NBFoldError(
            "REVIEW_MANIFEST_STALE",
            "review.json annotation_sha256 does not match annotation.json.",
            stage=stage,
        )

    if review_data.get("config_digest") != resolved_config["config_digest"]:
        raise NBFoldError(
            "REVIEW_MANIFEST_STALE",
            "review.json config_digest does not match the resolved reference/policy/optimizer configuration.",
            stage=stage,
        )

    protein_length = int(current_annotation["protein_length"])
    positions = review_data.get("reviewed_positions")
    if (
        not isinstance(positions, list)
        or not positions
        or any(type(p) is not int or not (1 <= p <= protein_length) for p in positions)
    ):
        raise NBFoldError(
            "INVALID_REVIEW_MANIFEST",
            f"reviewed_positions must be a non-empty list of 1-based integers within [1, {protein_length}].",
            stage=stage,
        )
    if len(positions) != len(set(positions)):
        raise NBFoldError(
            "INVALID_REVIEW_MANIFEST",
            "reviewed_positions must not contain duplicate positions.",
            stage=stage,
        )

    sorted_positions = sorted(positions)
    ann_recommended = sorted(loaded_ann.get("recommended_positions", []))
    ann_curated = loaded_ann.get("user_curated_positions")
    mod_rationale = review_data.get("modification_rationale")

    if sorted_positions != ann_recommended:
        has_valid_curated = (
            isinstance(ann_curated, list) and sorted(ann_curated) == sorted_positions
        )
        has_valid_rationale = isinstance(mod_rationale, str) and bool(mod_rationale.strip())
        if not (has_valid_curated or has_valid_rationale):
            raise NBFoldError(
                "REVIEW_MANIFEST_STALE",
                "reviewed_positions differ from annotation.json recommended_positions without "
                "updated annotation user_curated_positions or explicit modification_rationale.",
                stage=stage,
            )

    tools_block = current_annotation.get("independent_structural_tools", {})
    active_conflicts = evaluate_positions_against_structural_annotations(
        positions=sorted_positions,
        domain_caller_report=tools_block.get("domain_boundary_assignment", {}),
        tm_signal_caller_report=tools_block.get("transmembrane_topology", {}),
        backbone_ss_report=tools_block.get("backbone_secondary_structure", {}),
    )
    conflicts_ack = review_data.get("conflicts_acknowledged")
    conflict_rationale = review_data.get("conflict_resolution_rationale")
    if active_conflicts and (
        type(conflicts_ack) is not bool
        or conflicts_ack is not True
        or not isinstance(conflict_rationale, str)
        or not conflict_rationale.strip()
    ):
        raise NBFoldError(
            "REVIEW_CONFLICTS_UNACKNOWLEDGED",
            "Structural annotation conflicts exist for the selected positions or sidecars; "
            "review.json must set conflicts_acknowledged=true and provide a non-empty "
            "conflict_resolution_rationale.",
            stage=stage,
            details={
                "conflict_count": len(active_conflicts),
                "conflict_codes": [c.get("conflict_code") for c in active_conflicts],
            },
        )

    return {
        "schema_version": REVIEW_SCHEMA_VERSION,
        "review_type": REVIEW_ACKNOWLEDGEMENT_TYPE,
        "positions_reviewed": True,
        "reviewed_at": reviewed_at.strip(),
        "reviewer_label": reviewer_label.strip() if isinstance(reviewer_label, str) else None,
        "reviewed_positions": sorted_positions,
        "annotation_sha256": recomputed_ann_hash,
        "annotation_file_sha256": ann_file_sha256,
        "review_manifest_sha256": rev_file_sha256,
        "config_digest": resolved_config["config_digest"],
        "modification_rationale": mod_rationale.strip() if isinstance(mod_rationale, str) else None,
        "conflicts_present": bool(active_conflicts),
        "conflict_count": len(active_conflicts),
        "conflicts_acknowledged": bool(conflicts_ack) if active_conflicts else True,
        "conflict_resolution_rationale": (
            conflict_rationale.strip() if isinstance(conflict_rationale, str) and conflict_rationale.strip() else None
        ),
        "active_conflicts": active_conflicts,
    }


def evaluate_candidate_observations(
    *,
    candidate_id: str,
    engine_name: str,
    engine_ver: str,
    cds: str,
    expected_protein: str,
    cai: float,
    gc_fraction: float,
    structure_fit: float,
    gc_min_pct: float,
    gc_max_pct: float,
    forbidden_motifs: list[str],
    left_flank: str = "",
    right_flank: str = "",
    upstream_context: str | None = None,
    initiation_window_nt: int = DEFAULT_INITIATION_WINDOW_NT,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    sop_profile: SopProfile | None = None,
    unmapped_construct: bool = False,
) -> dict[str, Any]:
    """Run independent hard checks, sense-CDS plant RNA scan, initiation MFE, and rule policy."""
    hard_qc = validate_candidate_hard_invariants(
        cds=cds,
        protein=expected_protein,
        gc_min=gc_min_pct / 100.0,
        gc_max=gc_max_pct / 100.0,
        forbidden_motifs=forbidden_motifs,
        left_flank=left_flank,
        right_flank=right_flank,
    )

    clean_cds = "".join(cds.upper().split())
    clean_prot = "".join(expected_protein.upper().split()).rstrip("*")
    translated = translate_dna(clean_cds) if len(clean_cds) % 3 == 0 else ""
    if translated.endswith("*") and not clean_prot.endswith("*"):
        translated = translated[:-1]

    sense_cds_established = (not unmapped_construct) and (translated == clean_prot)

    if not sense_cds_established:
        rna_dict: dict[str, Any] = {
            "sequence_length": len(clean_cds),
            "detector_version": RNA_DETECTOR_VERSION,
            "evidence_class": "SEQUENCE_SCAN",
            "scope": "not_established",
            "summary_status": "NOT_AVAILABLE",
            "reason": "Expanded or unverified construct lacks confirmed sense-CDS coordinate mapping.",
            "donor_count": None,
            "polya_count": None,
            "are_count": None,
            "are_clustered": None,
            "cluster_details": [],
            "total_findings": None,
            "findings": [],
            "truthfulness_notice": RNA_TRUTHFULNESS_NOTICE,
        }
        mfe_dict: dict[str, Any] = {
            "status": "NOT_AVAILABLE",
            "evidence_class": "MODEL_PREDICTION",
            "engine": "ViennaRNA_RNAfold",
            "input_scope": "not_established",
            "window_start": 0,
            "window_end": 0,
            "folded_sequence_length": 0,
            "temperature_c": float(temperature_c),
            "mfe_kcal_mol": None,
            "secondary_structure": None,
            "runtime_ms": 0.0,
            "reason": "Sense-CDS start boundary is not established; refusing to fold unmapped construct flanks.",
            "advisory_flag": "UNMEASURED",
        }
        mfe_context_val = None
    else:
        rna_res = scan_plant_rna_risks(clean_cds)
        rna_dict = rna_res.to_dict(include_sequences=False)
        mfe_res = evaluate_initiation_mfe(
            cds=clean_cds,
            upstream_context=upstream_context,
            window_nt=initiation_window_nt,
            temperature_c=temperature_c,
        )
        mfe_dict = mfe_res.to_dict()
        mfe_context_val = mfe_res.mfe_kcal_mol if mfe_res.status == "AVAILABLE" else None

    registry = RuleRegistry(sop_profile=sop_profile)
    rule_eval = registry.evaluate_sequence(
        clean_cds,
        host="Nicotiana benthamiana",
        context={"initiation_mfe": mfe_context_val},
    )

    if mfe_dict["status"] == "ERROR":
        obs_status = "ERROR"
    elif not sense_cds_established or rule_eval["evaluation_status"] == "INCOMPLETE":
        obs_status = "INCOMPLETE"
    elif rule_eval["evaluation_status"] == "FAILED" or not hard_qc["passed"]:
        obs_status = "FAILED"
    elif mfe_dict["status"] == "NOT_AVAILABLE":
        obs_status = "COMPLETE_WITH_OPTIONAL_UNAVAILABLE"
    else:
        obs_status = "COMPLETE"

    return {
        "candidate_id": candidate_id,
        "engine": engine_name,
        "engine_version": engine_ver,
        "cds_sha256": sha256_text(clean_cds),
        "length_nt": len(clean_cds),
        "cai": cai,
        "gc_fraction": gc_fraction,
        "gc_percent": round(gc_fraction * 100.0, 2),
        "structure_fit": structure_fit,
        "codon_weight_target_conformity": structure_fit,
        "scope": "sense_cds_verified" if sense_cds_established else "not_established",
        "observation_status": obs_status,
        "hard_constraints": hard_qc,
        "plant_rna_scan": rna_dict,
        "initiation_mfe": mfe_dict,
        "rule_evaluation": rule_eval,
        "yield_prediction": {
            "status": "NOT_AVAILABLE",
            "value": None,
            "units": None,
            "model_version": None,
            "reason": (
                "Phase A/B does not train or execute a 14D YieldPredict model; "
                "13D checkpoints are not repurposed."
            ),
        },
    }


def _format_paired_fasta(
    *,
    cds_1dp: str,
    cds_3dp: str,
    cai_1dp: float,
    cai_3dp: float,
    gc_1dp_pct: float,
    gc_3dp_pct: float,
    fit_1dp: float,
    fit_3dp: float,
) -> str:
    """Format paired 1DP and 3DP CDS sequences with fixed candidate headers."""

    def wrap80(seq: str) -> str:
        return "\n".join(seq[i : i + 80] for i in range(0, len(seq), 80))

    h1 = (
        f">baseline_1dp|engine=dp_v2|role=comparator|length_nt={len(cds_1dp)}|"
        f"cai={cai_1dp:.4f}|gc_percent={gc_1dp_pct:.2f}|"
        f"codon_weight_target_conformity={fit_1dp:.4f}"
    )
    h3 = (
        f">structure_3dp|engine=structural_pause|role=experimental_candidate|length_nt={len(cds_3dp)}|"
        f"cai={cai_3dp:.4f}|gc_percent={gc_3dp_pct:.2f}|"
        f"codon_weight_target_conformity={fit_3dp:.4f}"
    )
    return f"{h1}\n{wrap80(cds_1dp)}\n{h3}\n{wrap80(cds_3dp)}\n"


def _build_markdown_report(
    *,
    annotation_data: Mapping[str, Any],
    review_provenance: Mapping[str, Any],
    resolved_config: Mapping[str, Any],
    observations: Mapping[str, Any],
    runtime_ms: float,
) -> str:
    """Render truthful, non-certifying NB_FOLD_REPORT.md."""
    comp = observations["paired_comparison"]
    c1 = observations["candidates"]["baseline_1dp"]
    c3 = observations["candidates"]["structure_3dp"]
    pol = observations["policy_summary"]
    mfe1 = c1["initiation_mfe"]
    mfe3 = c3["initiation_mfe"]
    rna1 = c1["plant_rna_scan"]
    rna3 = c3["plant_rna_scan"]

    mfe1_str = (
        f"{mfe1['mfe_kcal_mol']:.2f} kcal/mol ({mfe1['status']})"
        if mfe1["mfe_kcal_mol"] is not None
        else f"null ({mfe1['status']}: {mfe1.get('reason')})"
    )
    mfe3_str = (
        f"{mfe3['mfe_kcal_mol']:.2f} kcal/mol ({mfe3['status']})"
        if mfe3["mfe_kcal_mol"] is not None
        else f"null ({mfe3['status']}: {mfe3.get('reason')})"
    )

    lines = [
        "# NB-Fold Phase A/B Paired CDS Evaluation Report",
        "",
        f"- **Pipeline Version:** `{NBFOLD_PIPELINE_VERSION}` (`product_version={product_version()}`)",
        "- **Execution Mode:** `local_cli_phase_ab` (stateless local file intake; no external network calls)",
        f"- **Runtime:** `{runtime_ms:.2f} ms`",
        "",
        "## 1. Scientific & Causal Claim Boundary",
        "",
        "- AlphaFold predictions depend on amino acid sequence. Because `baseline_1dp` and `structure_3dp` encode the identical amino acid sequence via synonymous codons, AlphaFold structure or pLDDT cannot evaluate or compare post-design CDS folding or expression.",
        "- Structure confidence (`CA_atom_confidence` for AF3, `residue_plddt` for AF2, or `experimental_bfactor`) is used strictly as an upstream structural hypothesis for human review.",
        "- `structure_fit` (`codon_weight_target_conformity = 1.0 - mean(|w_i - target_adaptiveness|)`) is an in-silico codon-weight proxy, not measured ribosome elongation velocity or physical co-translational folding accuracy.",
        "- Plant RNA scan findings (`SEQUENCE_SCAN`) report sense-transcript motif matches; zero hits do not guarantee absence of in-vivo splicing or polyadenylation.",
        "- 5' Initiation MFE (`MODEL_PREDICTION`) is a local thermodynamic proxy over the 48-nt CDS start window.",
        "- **Yield Prediction:** `NOT_AVAILABLE` (`value: null`). Phase A/B does not train or run a 14D `YieldPredict` model.",
        "",
        "## 2. Structure Intake & Human Review Provenance",
        "",
        f"- **Protein Length:** `{annotation_data['protein_length']} aa` (`protein_sha256: {annotation_data['protein_sha256']}`)",
        f"- **Structure Format / Chain / Model:** `{annotation_data['format']}` · Chain `{annotation_data['chain_id']}` · Model `{annotation_data['model_number']}` (`structure_sha256: {annotation_data['structure_sha256']}`)",
        f"- **Provider Family & Source Basis:** `{annotation_data['provider_family']}` (`source_basis: {annotation_data['source_basis']}`, `{annotation_data['provenance_verification']}`)",
        f"- **Confidence Metric & Policy:** `{annotation_data['confidence_metric']}` (`{annotation_data['selected_atom_policy']}`)",
        f"- **Domain-Relative Confidence (PAE):** `{annotation_data['domain_relative_confidence']['status']}` (Domain segmentation: `{annotation_data['domain_relative_confidence']['domain_segmentation_status']}`)",
        f"- **Independent Structural Callers:** Domain `{annotation_data['independent_structural_tools']['domain_boundary_assignment']['status']}` · TM/Signal `{annotation_data['independent_structural_tools']['transmembrane_topology']['status']}` · Backbone SS `{annotation_data['independent_structural_tools']['backbone_secondary_structure']['status']}` · External DSSP `{annotation_data['independent_structural_tools']['secondary_structure_dssp']['status']}`",
        f"- **Structural Review Conflicts:** `count={review_provenance.get('conflict_count', 0)}` · `acknowledged={review_provenance.get('conflicts_acknowledged', True)}` · `rationale={review_provenance.get('conflict_resolution_rationale')}`",
        f"- **User Review Contract:** `{review_provenance['review_type']}` (`positions_reviewed={review_provenance['positions_reviewed']}`, `reviewed_at={review_provenance['reviewed_at']}`)",
        f"- **Reviewed 1-Based Positions ({len(review_provenance['reviewed_positions'])}):** `{review_provenance['reviewed_positions']}`",
        "",
        "## 3. Shared Reference & Constraint Envelope",
        "",
        f"- **Host:** `{resolved_config['host']}`",
        f"- **Codon Reference ID:** `{resolved_config['reference_id']}` (`tier={resolved_config['reference_tier']}`, `sha256={resolved_config['codon_table_sha256']}`)",
        f"- **Resolved GC Band:** `{resolved_config['gc_min_percent']:.2f}% - {resolved_config['gc_max_percent']:.2f}%`",
        f"- **Canonical Forbidden Motifs:** `{', '.join(resolved_config['forbidden_motifs'])}`",
        f"- **3DP Target Adaptiveness / Penalty Weight:** `target_adaptiveness={resolved_config['target_adaptiveness']}` · `penalty_weight={resolved_config['penalty_weight']}`",
        f"- **Configuration Digest:** `{resolved_config['config_digest']}`",
        "",
        "## 4. Paired 1DP (`dp_v2`) vs 3DP (`structural_pause`) Comparison",
        "",
        "| Metric | `baseline_1dp` (`dp_v2`) | `structure_3dp` (`structural_pause`) | Delta (3DP - 1DP) |",
        "|---|---:|---:|---:|",
        f"| **CAI** | `{c1['cai']:.4f}` | `{c3['cai']:.4f}` | `{comp['delta_1d_cai']:+.4f}` |",
        f"| **GC %** | `{c1['gc_percent']:.2f}%` | `{c3['gc_percent']:.2f}%` | `{comp['delta_gc_percent']:+.2f}%` |",
        f"| **`structure_fit` (`codon_weight_target_conformity`)** | `{c1['structure_fit']:.4f}` | `{c3['structure_fit']:.4f}` | `{comp['delta_3d_gain']:+.4f}` |",
        f"| **Hard Constraints (`AA / GC / Type IIS / Stop`)** | `{'PASS' if c1['hard_constraints']['passed'] else 'FAIL'}` | `{'PASS' if c3['hard_constraints']['passed'] else 'FAIL'}` | `--` |",
        f"| **Total Codon Edits** | `--` | `{comp['total_edits']}` | `--` |",
        f"| **Target-Zone Edits (±1 codon flank)** | `--` | `{comp['target_zone_edits']}` | `--` |",
        f"| **Outside-Zone Edits** | `--` | `{comp['outside_edits']}` | `--` |",
        f"| **Edit Localization Ratio** | `--` | `{comp['localization_ratio']*100:.1f}%` | `--` |",
        "",
        "## 5. Sense-CDS Observations & Policy Status",
        "",
        f"- **Observation Scope:** `{observations['scope']}`",
        f"- **`baseline_1dp` Plant RNA Scan (`SEQUENCE_SCAN`):** `{rna1['summary_status']}` (Donors: `{rna1['donor_count']}`, PolyA: `{rna1['polya_count']}`, ARE: `{rna1['are_count']}`, Clustered: `{rna1['are_clustered']}`)",
        f"- **`structure_3dp` Plant RNA Scan (`SEQUENCE_SCAN`):** `{rna3['summary_status']}` (Donors: `{rna3['donor_count']}`, PolyA: `{rna3['polya_count']}`, ARE: `{rna3['are_count']}`, Clustered: `{rna3['are_clustered']}`)",
        f"- **`baseline_1dp` 5' Initiation MFE (`MODEL_PREDICTION`):** `{mfe1_str}`",
        f"- **`structure_3dp` 5' Initiation MFE (`MODEL_PREDICTION`):** `{mfe3_str}`",
        f"- **Policy Status / Decision:** `policy_status={pol['policy_status']}` · `policy_decision={pol['policy_decision']}` · `automated_decision={pol['automated_decision']}`",
        "- **Yield Prediction:** `status=NOT_AVAILABLE` (`value=null`)",
        "",
    ]
    return "\n".join(lines)


def run_nbfold_annotate(
    *,
    protein_fasta_path: Path,
    structure_file_path: Path,
    output_dir: Path,
    declared_format: str | None = None,
    chain_id: str = "A",
    model_number: int = 1,
    provider_family: str | None = None,
    confidence_sidecar_path: Path | None = None,
    pae_sidecar_path: Path | None = None,
    export_metadata_path: Path | None = None,
    interpro_sidecar_path: Path | None = None,
    uniprot_sidecar_path: Path | None = None,
    host: str = "nbenthamiana",
    gc_min: float | None = None,
    gc_max: float | None = None,
    target_adaptiveness: float = 0.40,
    penalty_weight: float = 2.0,
    custom_forbidden_motifs: list[str] | None = None,
    sop_profile_path: Path | None = None,
    sop_profile_dict: Mapping[str, Any] | None = None,
    min_region_len: int = 3,
    max_region_len: int = 15,
    structure_backend: str = "alphafold_file",
    enable_esmfold: bool = False,
    enable_esm2: bool = False,
) -> dict[str, Any]:
    """Execute Phase A/C (`factorforge nbfold annotate`) and atomically write `annotation.json`.

    Never executes 1DP or 3DP optimization and always sets `reviewed=False`.
    """
    out = Path(output_dir)
    input_digests: dict[str, str | None] = {}
    staging_dir: Path | None = None

    try:
        _check_output_dir_no_conflict(out, ANNOTATE_ARTIFACT_NAMES, stage="annotate")
        enforce_backend_gate(
            structure_backend=structure_backend,
            enable_esmfold=enable_esmfold,
            enable_esm2=enable_esm2,
            stage="annotate",
        )
        resolved_cfg = resolve_nbfold_configuration(
            host=host,
            gc_min=gc_min,
            gc_max=gc_max,
            target_adaptiveness=target_adaptiveness,
            penalty_weight=penalty_weight,
            custom_forbidden_motifs=custom_forbidden_motifs,
            sop_profile_path=sop_profile_path,
            sop_profile_dict=sop_profile_dict,
            stage="annotate",
        )
        config_preview = {
            "host": resolved_cfg["host"],
            "reference_id": resolved_cfg["reference_id"],
            "codon_table_sha256": resolved_cfg["codon_table_sha256"],
            "gc_min_percent": resolved_cfg["gc_min_percent"],
            "gc_max_percent": resolved_cfg["gc_max_percent"],
            "forbidden_motifs": resolved_cfg["forbidden_motifs"],
            "target_adaptiveness": resolved_cfg["target_adaptiveness"],
            "penalty_weight": resolved_cfg["penalty_weight"],
            "sop_profile_hash": resolved_cfg["sop_profile_hash"],
        }

        annotation_payload = build_structure_annotation(
            protein_fasta_path=Path(protein_fasta_path),
            structure_file_path=Path(structure_file_path),
            declared_format=declared_format,
            chain_id=chain_id,
            model_number=model_number,
            provider_family=provider_family,
            confidence_sidecar_path=confidence_sidecar_path,
            pae_sidecar_path=pae_sidecar_path,
            export_metadata_path=export_metadata_path,
            interpro_sidecar_path=interpro_sidecar_path,
            uniprot_sidecar_path=uniprot_sidecar_path,
            min_region_len=min_region_len,
            max_region_len=max_region_len,
            default_config_digest=resolved_cfg["config_digest"],
            resolved_config_preview=config_preview,
            stage="annotate",
        )
        input_digests["protein_sha256"] = annotation_payload.get("protein_sha256")
        input_digests["structure_sha256"] = annotation_payload.get("structure_sha256")

        out.mkdir(parents=True, exist_ok=True)
        staging_dir = Path(tempfile.mkdtemp(prefix=".nbfold_staging_", dir=str(out)))
        staged_ann = staging_dir / ANNOTATION_FILENAME
        staged_ann.write_text(
            json.dumps(annotation_payload, indent=2, sort_keys=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        final_ann = out / ANNOTATION_FILENAME
        os.replace(staged_ann, final_ann)
        shutil.rmtree(staging_dir, ignore_errors=True)
        staging_dir = None

        return {
            "annotation_path": str(final_ann),
            "annotation": annotation_payload,
        }
    except NBFoldError as err:
        if staging_dir is not None:
            shutil.rmtree(staging_dir, ignore_errors=True)
        write_error_manifest(out, err, input_digests=input_digests)
        raise
    except Exception as exc:
        if staging_dir is not None:
            shutil.rmtree(staging_dir, ignore_errors=True)
        wrapped = NBFoldError(
            "INTERNAL_ERROR",
            f"Unexpected internal error during annotate: {type(exc).__name__}",
            stage="annotate",
        )
        write_error_manifest(out, wrapped, input_digests=input_digests)
        raise wrapped from exc


def run_nbfold_optimize(
    *,
    protein_fasta_path: Path,
    structure_file_path: Path,
    annotation_path: Path,
    review_manifest_path: Path,
    output_dir: Path,
    declared_format: str | None = None,
    chain_id: str | None = None,
    model_number: int | None = None,
    provider_family: str | None = None,
    confidence_sidecar_path: Path | None = None,
    pae_sidecar_path: Path | None = None,
    export_metadata_path: Path | None = None,
    interpro_sidecar_path: Path | None = None,
    uniprot_sidecar_path: Path | None = None,
    host: str = "nbenthamiana",
    gc_min: float | None = None,
    gc_max: float | None = None,
    target_adaptiveness: float = 0.40,
    penalty_weight: float = 2.0,
    custom_forbidden_motifs: list[str] | None = None,
    left_flank: str = "",
    right_flank: str = "",
    upstream_context: str | None = None,
    initiation_window_nt: int = DEFAULT_INITIATION_WINDOW_NT,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    sop_profile_path: Path | None = None,
    sop_profile_dict: Mapping[str, Any] | None = None,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    structure_backend: str = "alphafold_file",
    enable_esmfold: bool = False,
    enable_esm2: bool = False,
) -> dict[str, Any]:
    """Execute Phase B/C (`factorforge nbfold optimize`) after verifying explicit `review.json`."""
    out = Path(output_dir)
    input_digests: dict[str, str | None] = {}
    staging_dir: Path | None = None
    started_clock = time.perf_counter()
    tracemalloc.start()

    try:
        _check_output_dir_no_conflict(out, OPTIMIZE_ARTIFACT_NAMES, stage="optimize")
        enforce_backend_gate(
            structure_backend=structure_backend,
            enable_esmfold=enable_esmfold,
            enable_esm2=enable_esm2,
            stage="optimize",
        )
        if (
            isinstance(timeout_seconds, bool)
            or not isinstance(timeout_seconds, (int, float))
            or not math.isfinite(float(timeout_seconds))
            or float(timeout_seconds) <= 0.0
        ):
            raise NBFoldError(
                "EXECUTION_TIMEOUT",
                "timeout_seconds must be a positive finite number.",
                stage="optimize",
            )

        # Pre-check total input byte size across FASTA + structure + annotation + review + sidecars
        total_input_bytes = 0
        for p_item, max_b, lbl, code in (
            (protein_fasta_path, MAX_MANIFEST_JSON_BYTES, "Protein FASTA", "INVALID_FASTA"),
            (structure_file_path, MAX_SIDECAR_BYTES, "Structure", "INVALID_STRUCTURE_FORMAT"),
            (annotation_path, MAX_MANIFEST_JSON_BYTES, "Annotation manifest", "REVIEW_MANIFEST_MISSING"),
            (review_manifest_path, MAX_MANIFEST_JSON_BYTES, "Review manifest", "REVIEW_MANIFEST_MISSING"),
        ):
            if p_item is None:
                raise NBFoldError(code, f"{lbl} path is required.", stage="optimize")
            sz = check_file_size_limit(
                Path(p_item), max_bytes=max_b, label=lbl, empty_error_code=code, stage="optimize"
            )
            total_input_bytes += sz

        for sp in (
            confidence_sidecar_path,
            pae_sidecar_path,
            export_metadata_path,
            interpro_sidecar_path,
            uniprot_sidecar_path,
            sop_profile_path,
        ):
            if sp is not None:
                total_input_bytes += check_file_size_limit(
                    Path(sp),
                    max_bytes=MAX_SIDECAR_BYTES,
                    label="Sidecar/SOP input",
                    empty_error_code="SIDECAR_MISMATCH",
                    stage="optimize",
                )
        if total_input_bytes > MAX_TOTAL_INPUT_BYTES:
            raise NBFoldError(
                "UNSUPPORTED_INPUT_SIZE",
                f"Total input size ({total_input_bytes} bytes) exceeds 5 MB limit ({MAX_TOTAL_INPUT_BYTES} bytes).",
                stage="optimize",
            )

        # Read chain/model/provider defaults from annotation.json if not explicitly overridden on CLI
        try:
            ann_preview = json.loads(Path(annotation_path).read_text(encoding="utf-8"))
        except Exception as exc:
            raise NBFoldError(
                "INVALID_REVIEW_MANIFEST",
                "Annotation manifest (annotation.json) is not valid UTF-8 JSON.",
                stage="optimize",
            ) from exc
        if not isinstance(ann_preview, dict):
            raise NBFoldError(
                "INVALID_REVIEW_MANIFEST",
                "Annotation manifest must be a JSON object.",
                stage="optimize",
            )

        eff_chain = chain_id if chain_id is not None else ann_preview.get("chain_id", "A")
        eff_model = model_number if model_number is not None else ann_preview.get("model_number", 1)
        eff_provider = (
            provider_family
            if provider_family is not None
            else ann_preview.get("provider_family")
        )
        eff_format = declared_format if declared_format is not None else ann_preview.get("format")

        resolved_cfg = resolve_nbfold_configuration(
            host=host,
            gc_min=gc_min,
            gc_max=gc_max,
            target_adaptiveness=target_adaptiveness,
            penalty_weight=penalty_weight,
            custom_forbidden_motifs=custom_forbidden_motifs,
            sop_profile_path=sop_profile_path,
            sop_profile_dict=sop_profile_dict,
            stage="optimize",
        )

        # Re-run StructureBridge verification on current inputs
        current_annotation = build_structure_annotation(
            protein_fasta_path=Path(protein_fasta_path),
            structure_file_path=Path(structure_file_path),
            declared_format=eff_format,
            chain_id=eff_chain,
            model_number=eff_model,
            provider_family=eff_provider,
            confidence_sidecar_path=confidence_sidecar_path,
            pae_sidecar_path=pae_sidecar_path,
            export_metadata_path=export_metadata_path,
            interpro_sidecar_path=interpro_sidecar_path,
            uniprot_sidecar_path=uniprot_sidecar_path,
            default_config_digest=resolved_cfg["config_digest"],
            stage="optimize",
        )
        input_digests["protein_sha256"] = current_annotation.get("protein_sha256")
        input_digests["structure_sha256"] = current_annotation.get("structure_sha256")

        # Verify explicit user review contract
        review_provenance = verify_review_manifest(
            review_manifest_path=Path(review_manifest_path),
            annotation_path=Path(annotation_path),
            current_annotation=current_annotation,
            resolved_config=resolved_cfg,
            stage="optimize",
        )
        input_digests["annotation_sha256"] = review_provenance["annotation_sha256"]
        input_digests["review_manifest_sha256"] = review_provenance["review_manifest_sha256"]

        parsed_protein = parse_protein_fasta_file(Path(protein_fasta_path), stage="optimize")
        reviewed_positions = review_provenance["reviewed_positions"]

        # Execute paired 1DP and 3DP under identical resolved reference and constraints
        try:
            comparison = compute_tradeoff_comparison(
                protein=parsed_protein.sequence,
                pause_positions=reviewed_positions,
                codon_weights=resolved_cfg["codon_weights"],
                target_gc_min=resolved_cfg["gc_min_percent"] / 100.0,
                target_gc_max=resolved_cfg["gc_max_percent"] / 100.0,
                forbidden_motifs=resolved_cfg["forbidden_motifs"],
                lambda_weight=resolved_cfg["penalty_weight"],
                target_adaptiveness=resolved_cfg["target_adaptiveness"],
                left_flank=left_flank,
                right_flank=right_flank,
            )
        except (RuntimeError, ValueError) as exc:
            raise NBFoldError(
                "INFEASIBLE",
                f"Paired 1DP/3DP optimization is infeasible under the resolved constraints: {exc}",
                stage="optimize",
                details={
                    "gc_min_percent": resolved_cfg["gc_min_percent"],
                    "gc_max_percent": resolved_cfg["gc_max_percent"],
                    "protein_length_aa": len(parsed_protein.sequence),
                },
            ) from exc

        if (time.perf_counter() - started_clock) > float(timeout_seconds):
            raise NBFoldError(
                "EXECUTION_TIMEOUT",
                f"Optimization exceeded timeout budget ({timeout_seconds}s).",
                stage="optimize",
            )

        # Evaluate independent QC, RNA scan, MFE, and SOP policy on both candidates
        obs_1dp = evaluate_candidate_observations(
            candidate_id="baseline_1dp",
            engine_name="dp_v2",
            engine_ver=engine_version("dp"),
            cds=comparison.cds_1dp,
            expected_protein=parsed_protein.sequence,
            cai=comparison.cai_1dp,
            gc_fraction=comparison.gc_1dp,
            structure_fit=comparison.structure_fit_1dp,
            gc_min_pct=resolved_cfg["gc_min_percent"],
            gc_max_pct=resolved_cfg["gc_max_percent"],
            forbidden_motifs=resolved_cfg["forbidden_motifs"],
            left_flank=left_flank,
            right_flank=right_flank,
            upstream_context=upstream_context,
            initiation_window_nt=initiation_window_nt,
            temperature_c=temperature_c,
            sop_profile=resolved_cfg["sop_profile"],
        )
        obs_3dp = evaluate_candidate_observations(
            candidate_id="structure_3dp",
            engine_name="structural_pause",
            engine_ver="not_independently_versioned",
            cds=comparison.cds_3dp,
            expected_protein=parsed_protein.sequence,
            cai=comparison.cai_3dp,
            gc_fraction=comparison.gc_3dp,
            structure_fit=comparison.structure_fit_3dp,
            gc_min_pct=resolved_cfg["gc_min_percent"],
            gc_max_pct=resolved_cfg["gc_max_percent"],
            forbidden_motifs=resolved_cfg["forbidden_motifs"],
            left_flank=left_flank,
            right_flank=right_flank,
            upstream_context=upstream_context,
            initiation_window_nt=initiation_window_nt,
            temperature_c=temperature_c,
            sop_profile=resolved_cfg["sop_profile"],
        )

        if (
            not comparison.hard_constraints_pass
            or not obs_1dp["hard_constraints"]["passed"]
            or not obs_3dp["hard_constraints"]["passed"]
        ):
            raise NBFoldError(
                "INFEASIBLE",
                "Generated candidate(s) failed independent full-CDS hard constraint verification.",
                stage="optimize",
            )

        if (
            obs_1dp["initiation_mfe"]["status"] == "ERROR"
            or obs_3dp["initiation_mfe"]["status"] == "ERROR"
        ):
            raise NBFoldError(
                "OBSERVATION_ENGINE_ERROR",
                "Initiation MFE engine failed or returned invalid folding observations.",
                stage="optimize",
                details={
                    "baseline_1dp_mfe_status": obs_1dp["initiation_mfe"]["status"],
                    "structure_3dp_mfe_status": obs_3dp["initiation_mfe"]["status"],
                },
            )

        req_unavail = (
            obs_1dp["rule_evaluation"].get("required_unavailable_count", 0)
            + obs_3dp["rule_evaluation"].get("required_unavailable_count", 0)
        )
        if req_unavail > 0:
            raise NBFoldError(
                "REQUIRED_CHECK_UNAVAILABLE",
                "Required policy observation(s) are unavailable; execution is INCOMPLETE/BLOCKED.",
                stage="optimize",
                details={"required_unavailable_count": req_unavail},
            )

        any_hard_fail = (
            not obs_1dp["rule_evaluation"]["all_passed"]
            or not obs_3dp["rule_evaluation"]["all_passed"]
        )
        if any_hard_fail:
            raise NBFoldError(
                "POLICY_BLOCKED",
                "One or more candidates failed a required policy rule (policy_decision=BLOCK).",
                stage="optimize",
                details={
                    "baseline_1dp_hard_fails": obs_1dp["rule_evaluation"]["hard_fail_count"],
                    "structure_3dp_hard_fails": obs_3dp["rule_evaluation"]["hard_fail_count"],
                },
            )

        elapsed_s = time.perf_counter() - started_clock
        if elapsed_s > float(timeout_seconds):
            raise NBFoldError(
                "EXECUTION_TIMEOUT",
                f"Optimization exceeded timeout budget ({timeout_seconds}s).",
                stage="optimize",
            )

        rna_comparison = compare_candidate_plant_rna_risks(
            comparison.cds_1dp, comparison.cds_3dp, protein_name="paired_nbfold"
        )
        rna_comparison.pop("summary_table", None)

        has_optional_unavail = (
            obs_1dp["initiation_mfe"]["status"] == "NOT_AVAILABLE"
            or obs_3dp["initiation_mfe"]["status"] == "NOT_AVAILABLE"
        )
        has_warnings = (
            obs_1dp["rule_evaluation"]["warning_count"] > 0
            or obs_3dp["rule_evaluation"]["warning_count"] > 0
        )
        policy_status = (
            "COMPLETE_WITH_OPTIONAL_UNAVAILABLE" if has_optional_unavail else "COMPLETE"
        )
        policy_decision = "PASS"
        automated_decision = "CONDITIONAL_PASS" if has_warnings else "PASS"

        observations_payload: dict[str, Any] = {
            "schema_version": OBSERVATIONS_SCHEMA_VERSION,
            "pipeline_version": NBFOLD_PIPELINE_VERSION,
            "scope": "sense_cds_verified",
            "candidates": {
                "baseline_1dp": obs_1dp,
                "structure_3dp": obs_3dp,
            },
            "paired_comparison": {
                "cai_1dp": comparison.cai_1dp,
                "cai_3dp": comparison.cai_3dp,
                "delta_1d_cai": comparison.delta_1d_cai,
                "gc_1dp_fraction": comparison.gc_1dp,
                "gc_3dp_fraction": comparison.gc_3dp,
                "gc_1dp_percent": round(comparison.gc_1dp * 100.0, 2),
                "gc_3dp_percent": round(comparison.gc_3dp * 100.0, 2),
                "delta_gc_fraction": comparison.delta_gc,
                "delta_gc_percent": round(comparison.delta_gc * 100.0, 2),
                "structure_fit_1dp": comparison.structure_fit_1dp,
                "structure_fit_3dp": comparison.structure_fit_3dp,
                "codon_weight_target_conformity_1dp": comparison.structure_fit_1dp,
                "codon_weight_target_conformity_3dp": comparison.structure_fit_3dp,
                "delta_3d_gain": comparison.delta_3d_gain,
                "delta_codon_weight_target_conformity": comparison.delta_3d_gain,
                "conformity_formula": "1.0 - mean(|w_i - target_adaptiveness|) over reviewed_positions",
                "target_adaptiveness": resolved_cfg["target_adaptiveness"],
                "penalty_weight": resolved_cfg["penalty_weight"],
                "reviewed_positions": reviewed_positions,
                "total_edits": comparison.total_edits,
                "target_zone_edits": comparison.target_zone_edits,
                "outside_edits": comparison.outside_edits,
                "localization_ratio": comparison.localization_ratio,
                "hard_constraints_pass": comparison.hard_constraints_pass,
                "plant_rna_comparison": rna_comparison,
            },
            "structural_annotation_summary": {
                "domain_relative_confidence_status": current_annotation["domain_relative_confidence"]["status"],
                "domain_segmentation_status": current_annotation["domain_relative_confidence"]["domain_segmentation_status"],
                "domain_boundary_assignment_status": current_annotation["independent_structural_tools"]["domain_boundary_assignment"]["status"],
                "transmembrane_topology_status": current_annotation["independent_structural_tools"]["transmembrane_topology"]["status"],
                "backbone_secondary_structure_status": current_annotation["independent_structural_tools"]["backbone_secondary_structure"]["status"],
                "secondary_structure_dssp_status": current_annotation["independent_structural_tools"]["secondary_structure_dssp"]["status"],
                "conflict_count": review_provenance["conflict_count"],
                "conflicts_acknowledged": review_provenance["conflicts_acknowledged"],
            },
            "policy_summary": {
                "policy_status": policy_status,
                "policy_decision": policy_decision,
                "automated_decision": automated_decision,
                "ruleset_digest": obs_3dp["rule_evaluation"]["ruleset_digest"],
                "sop_profile": (
                    resolved_cfg["sop_profile"].provenance()
                    if resolved_cfg["sop_profile"] is not None
                    else "default_rule_registry"
                ),
            },
            "yield_prediction": {
                "status": "NOT_AVAILABLE",
                "value": None,
                "units": None,
                "model_version": None,
                "reason": (
                    "Phase A/B does not train or execute a 14D YieldPredict model; "
                    "13D checkpoints are not repurposed."
                ),
            },
        }

        runtime_ms = round((time.perf_counter() - started_clock) * 1000.0, 2)
        _, peak_memory_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        fasta_content = _format_paired_fasta(
            cds_1dp=comparison.cds_1dp,
            cds_3dp=comparison.cds_3dp,
            cai_1dp=comparison.cai_1dp,
            cai_3dp=comparison.cai_3dp,
            gc_1dp_pct=round(comparison.gc_1dp * 100.0, 2),
            gc_3dp_pct=round(comparison.gc_3dp * 100.0, 2),
            fit_1dp=comparison.structure_fit_1dp,
            fit_3dp=comparison.structure_fit_3dp,
        )
        observations_text = json.dumps(observations_payload, indent=2, sort_keys=False) + "\n"
        report_md_text = _build_markdown_report(
            annotation_data=current_annotation,
            review_provenance=review_provenance,
            resolved_config=resolved_cfg,
            observations=observations_payload,
            runtime_ms=runtime_ms,
        )

        artifact_hashes = {
            PAIRED_FASTA_FILENAME: sha256_text(fasta_content),
            OBSERVATIONS_FILENAME: sha256_text(observations_text),
            REPORT_MARKDOWN_FILENAME: sha256_text(report_md_text),
        }

        run_manifest_payload: dict[str, Any] = {
            "schema_version": RUN_MANIFEST_SCHEMA_VERSION,
            "pipeline_version": NBFOLD_PIPELINE_VERSION,
            "product_version": product_version(),
            "execution_mode": "local_cli_phase_ab",
            "status": "COMPLETED",
            "timestamp_utc": _utc_now_iso(),
            "runtime_ms": runtime_ms,
            "timeout_budget_seconds": float(timeout_seconds),
            "peak_memory_bytes": int(peak_memory_bytes),
            "input_provenance": {
                "protein_length_aa": len(parsed_protein.sequence),
                "protein_sha256": parsed_protein.protein_sha256,
                "protein_fasta_file_sha256": parsed_protein.file_sha256,
                "structure_file_sha256": current_annotation["structure_sha256"],
                "structure_format": current_annotation["format"],
                "chain_id": current_annotation["chain_id"],
                "model_number": current_annotation["model_number"],
                "provider_family": current_annotation["provider_family"],
                "source_basis": current_annotation["source_basis"],
                "export_version": current_annotation["export_version"],
                "provenance_verification": "unverified_local_claim",
                "sidecar_sha256s": current_annotation["sidecar_sha256s"],
            },
            "review_provenance": {
                **review_provenance,
                "authority_disclaimer": (
                    "LOCAL_USER_ACKNOWLEDGEMENT only; not a cryptographic signature or "
                    "ResearchOS/AgentOS/ValidationHub server approval."
                ),
            },
            "reference_and_constraints": {
                "host": resolved_cfg["host"],
                "reference_id": resolved_cfg["reference_id"],
                "reference_tier": resolved_cfg["reference_tier"],
                "reference_policy_version": resolved_cfg["reference_policy_version"],
                "codon_table_sha256": resolved_cfg["codon_table_sha256"],
                "codon_weights_sha256": resolved_cfg["codon_weights_sha256"],
                "reference_claim_boundary": resolved_cfg["reference_claim_boundary"],
                "resolved_gc_min_percent": resolved_cfg["gc_min_percent"],
                "resolved_gc_max_percent": resolved_cfg["gc_max_percent"],
                "forbidden_motifs": resolved_cfg["forbidden_motifs"],
                "target_adaptiveness": resolved_cfg["target_adaptiveness"],
                "penalty_weight": resolved_cfg["penalty_weight"],
                "left_flank_length_nt": len(left_flank),
                "right_flank_length_nt": len(right_flank),
                "config_digest": resolved_cfg["config_digest"],
            },
            "engines_and_backends": {
                "comparator_1dp": {
                    "engine": "dp_v2",
                    "version": engine_version("dp"),
                    "objective": "feasibility_best_log_adaptiveness",
                },
                "candidate_3dp": {
                    "engine": "structural_pause",
                    "version": "not_independently_versioned",
                    "objective": "log_adaptiveness_minus_position_penalty",
                },
                "structure_recommender": {
                    "tool": "auto_annotator",
                    "version": RECOMMENDER_VERSION,
                },
                "rna_scanner": {
                    "tool": "plant_rna_guards",
                    "version": RNA_DETECTOR_VERSION,
                },
                "initiation_mfe": {
                    "engine": obs_3dp["initiation_mfe"]["engine"],
                    "status": obs_3dp["initiation_mfe"]["status"],
                    "temperature_c": float(temperature_c),
                    "window_nt": int(initiation_window_nt),
                },
                "future_backends": get_backend_capabilities(),
                "watermark_enabled": False,
                "external_network_calls": False,
                "server_persistence": "not_connected",
            },
            "output_artifact_sha256s": artifact_hashes,
        }
        run_manifest_text = json.dumps(run_manifest_payload, indent=2, sort_keys=False) + "\n"

        # Atomic staging and finalization inside output_dir
        out.mkdir(parents=True, exist_ok=True)
        staging_dir = Path(tempfile.mkdtemp(prefix=".nbfold_staging_", dir=str(out)))

        staged_files = {
            PAIRED_FASTA_FILENAME: fasta_content,
            OBSERVATIONS_FILENAME: observations_text,
            REPORT_MARKDOWN_FILENAME: report_md_text,
            # Publish completion evidence only after every referenced artifact.
            RUN_MANIFEST_FILENAME: run_manifest_text,
        }
        for fname, content in staged_files.items():
            (staging_dir / fname).write_text(content, encoding="utf-8", newline="\n")

        for fname in staged_files:
            os.replace(staging_dir / fname, out / fname)

        shutil.rmtree(staging_dir, ignore_errors=True)
        staging_dir = None

        return {
            "output_dir": str(out),
            "artifacts": {
                name: str(out / name) for name in staged_files
            },
            "observations": observations_payload,
            "run_manifest": run_manifest_payload,
        }

    except NBFoldError as err:
        if tracemalloc.is_tracing():
            tracemalloc.stop()
        if staging_dir is not None:
            shutil.rmtree(staging_dir, ignore_errors=True)
        write_error_manifest(out, err, input_digests=input_digests)
        raise
    except Exception as exc:
        if tracemalloc.is_tracing():
            tracemalloc.stop()
        if staging_dir is not None:
            shutil.rmtree(staging_dir, ignore_errors=True)
        wrapped = NBFoldError(
            "INTERNAL_ERROR",
            f"Unexpected internal error during optimize: {type(exc).__name__}",
            stage="optimize",
        )
        write_error_manifest(out, wrapped, input_digests=input_digests)
        raise wrapped from exc
