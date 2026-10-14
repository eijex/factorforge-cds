"""Click CLI subcommands for NB-Fold (`factorforge nbfold annotate` and `factorforge nbfold optimize`)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from factorforge.nbfold.contracts import DEFAULT_TIMEOUT_SECONDS, NBFoldError
from factorforge.nbfold.pipeline import run_nbfold_annotate, run_nbfold_optimize


@click.group("nbfold")
def nbfold_cli() -> None:
    """NB-Fold local reviewed AlphaFold-file to FactorForge 1DP/3DP pipeline (Phase A/B)."""


@nbfold_cli.command("annotate")
@click.option("--protein-fasta", required=True, type=click.Path(path_type=Path), help="Single-record protein FASTA file (3-300 aa).")
@click.option("--structure-file", required=True, type=click.Path(path_type=Path), help="Local PDB or mmCIF structure file (<= 250 KB).")
@click.option("--format", "declared_format", default=None, type=str, help="Structure format ('pdb', 'cif', or 'mmcif').")
@click.option("--chain", "chain_id", default="A", type=str, help="Explicit structure chain ID.")
@click.option("--model", "model_number", default=1, type=int, help="1-based file model order.")
@click.option("--provider-family", default=None, type=str, help="Confidence provider family ('af2', 'af3', 'experimental', 'unknown').")
@click.option("--confidence-sidecar", default=None, type=click.Path(path_type=Path), help="Optional verified confidence JSON sidecar.")
@click.option("--pae-sidecar", default=None, type=click.Path(path_type=Path), help="Optional verified PAE JSON sidecar.")
@click.option("--export-metadata", default=None, type=click.Path(path_type=Path), help="Optional export metadata JSON sidecar.")
@click.option("--interpro-sidecar", default=None, type=click.Path(path_type=Path), help="Optional verified InterPro domain JSON sidecar.")
@click.option("--uniprot-sidecar", default=None, type=click.Path(path_type=Path), help="Optional verified UniProt feature JSON sidecar.")
@click.option("--host", default="nbenthamiana", type=str, help="Target expression host (default: nbenthamiana).")
@click.option("--gc-min", default=None, type=float, help="Optional minimum GC percentage override.")
@click.option("--gc-max", default=None, type=float, help="Optional maximum GC percentage override.")
@click.option("--target-adaptiveness", default=0.40, type=float, help="3DP target codon adaptiveness in (0, 1].")
@click.option("--penalty-weight", default=2.0, type=float, help="3DP position penalty weight in [0, 100].")
@click.option("--min-region-len", default=3, type=int, help="Minimum contiguous candidate region length.")
@click.option("--max-region-len", default=15, type=int, help="Maximum contiguous candidate region length.")
@click.option("--structure-backend", default="alphafold_file", type=str, help="Structure intake backend ('alphafold_file').")
@click.option("--enable-esmfold", is_flag=True, default=False, help="Future ESMFold flag (fails closed in Phase A/B/C).")
@click.option("--enable-esm2", is_flag=True, default=False, help="Future ESM-2 flag (fails closed in Phase A/B/C).")
@click.option("--output-dir", required=True, type=click.Path(path_type=Path), help="Target directory for annotation.json.")
def annotate_cmd(
    protein_fasta: Path,
    structure_file: Path,
    declared_format: str | None,
    chain_id: str,
    model_number: int,
    provider_family: str | None,
    confidence_sidecar: Path | None,
    pae_sidecar: Path | None,
    export_metadata: Path | None,
    interpro_sidecar: Path | None,
    uniprot_sidecar: Path | None,
    host: str,
    gc_min: float | None,
    gc_max: float | None,
    target_adaptiveness: float,
    penalty_weight: float,
    min_region_len: int,
    max_region_len: int,
    structure_backend: str,
    enable_esmfold: bool,
    enable_esm2: bool,
    output_dir: Path,
) -> None:
    """Phase A/C: Validate structure mapping, observe confidence & structural annotations, and write annotation.json."""
    try:
        res = run_nbfold_annotate(
            protein_fasta_path=protein_fasta,
            structure_file_path=structure_file,
            output_dir=output_dir,
            declared_format=declared_format,
            chain_id=chain_id,
            model_number=model_number,
            provider_family=provider_family,
            confidence_sidecar_path=confidence_sidecar,
            pae_sidecar_path=pae_sidecar,
            export_metadata_path=export_metadata,
            interpro_sidecar_path=interpro_sidecar,
            uniprot_sidecar_path=uniprot_sidecar,
            host=host,
            gc_min=gc_min,
            gc_max=gc_max,
            target_adaptiveness=target_adaptiveness,
            penalty_weight=penalty_weight,
            min_region_len=min_region_len,
            max_region_len=max_region_len,
            structure_backend=structure_backend,
            enable_esmfold=enable_esmfold,
            enable_esm2=enable_esm2,
        )
        ann = res["annotation"]
        rev_sum = ann.get("structural_review_summary", {})
        summary = {
            "status": ann["status"],
            "workflow_state": ann["workflow_state"],
            "reviewed": ann["reviewed"],
            "provider_family": ann["provider_family"],
            "confidence_metric": ann["confidence_metric"],
            "recommended_positions": ann["recommended_positions"],
            "conflict_count": rev_sum.get("conflict_count", 0),
            "requires_explicit_conflict_acknowledgement": rev_sum.get(
                "requires_explicit_conflict_acknowledgement", False
            ),
            "annotation_path": res["annotation_path"],
            "annotation_content_sha256": ann["annotation_content_sha256"],
        }
        click.echo(json.dumps(summary))
    except NBFoldError as err:
        err_payload = {
            "status": "ERROR",
            "stage": err.stage,
            "error_code": err.error_code,
            "exit_code": err.exit_code,
            "message": err.message,
        }
        click.echo(json.dumps(err_payload), err=True)
        sys.exit(err.exit_code)


@nbfold_cli.command("optimize")
@click.option("--protein-fasta", required=True, type=click.Path(path_type=Path), help="Single-record protein FASTA file (3-300 aa).")
@click.option("--structure-file", required=True, type=click.Path(path_type=Path), help="Local PDB or mmCIF structure file (<= 250 KB).")
@click.option("--annotation", "annotation_path", required=True, type=click.Path(path_type=Path), help="Path to annotation.json from Phase A.")
@click.option("--review-manifest", "review_manifest_path", required=True, type=click.Path(path_type=Path), help="Path to explicit user review.json.")
@click.option("--format", "declared_format", default=None, type=str, help="Optional structure format ('pdb', 'cif', or 'mmcif').")
@click.option("--chain", "chain_id", default=None, type=str, help="Optional explicit structure chain ID.")
@click.option("--model", "model_number", default=None, type=int, help="Optional 1-based file model order.")
@click.option("--provider-family", default=None, type=str, help="Optional provider family ('af2', 'af3', 'experimental').")
@click.option("--confidence-sidecar", default=None, type=click.Path(path_type=Path), help="Optional verified confidence JSON sidecar.")
@click.option("--pae-sidecar", default=None, type=click.Path(path_type=Path), help="Optional verified PAE JSON sidecar.")
@click.option("--export-metadata", default=None, type=click.Path(path_type=Path), help="Optional export metadata JSON sidecar.")
@click.option("--interpro-sidecar", default=None, type=click.Path(path_type=Path), help="Optional verified InterPro domain JSON sidecar.")
@click.option("--uniprot-sidecar", default=None, type=click.Path(path_type=Path), help="Optional verified UniProt feature JSON sidecar.")
@click.option("--host", default="nbenthamiana", type=str, help="Target expression host (default: nbenthamiana).")
@click.option("--gc-min", default=None, type=float, help="Optional minimum GC percentage override.")
@click.option("--gc-max", default=None, type=float, help="Optional maximum GC percentage override.")
@click.option("--target-adaptiveness", default=0.40, type=float, help="3DP target codon adaptiveness in (0, 1].")
@click.option("--penalty-weight", default=2.0, type=float, help="3DP position penalty weight in [0, 100].")
@click.option("--sop-profile", "sop_profile_path", default=None, type=click.Path(path_type=Path), help="Optional SOP profile JSON path.")
@click.option("--upstream-context", default=None, type=str, help="Optional 5' upstream nucleotide flank for initiation MFE.")
@click.option("--temperature-c", default=37.0, type=float, help="Initiation MFE folding temperature in Celsius.")
@click.option("--timeout-seconds", default=DEFAULT_TIMEOUT_SECONDS, type=float, help="Execution timeout budget in seconds.")
@click.option("--structure-backend", default="alphafold_file", type=str, help="Structure intake backend ('alphafold_file').")
@click.option("--enable-esmfold", is_flag=True, default=False, help="Future ESMFold flag (fails closed in Phase A/B/C).")
@click.option("--enable-esm2", is_flag=True, default=False, help="Future ESM-2 flag (fails closed in Phase A/B/C).")
@click.option("--output-dir", required=True, type=click.Path(path_type=Path), help="Target directory for paired design artifacts.")
def optimize_cmd(
    protein_fasta: Path,
    structure_file: Path,
    annotation_path: Path,
    review_manifest_path: Path,
    declared_format: str | None,
    chain_id: str | None,
    model_number: int | None,
    provider_family: str | None,
    confidence_sidecar: Path | None,
    pae_sidecar: Path | None,
    export_metadata: Path | None,
    interpro_sidecar: Path | None,
    uniprot_sidecar: Path | None,
    host: str,
    gc_min: float | None,
    gc_max: float | None,
    target_adaptiveness: float,
    penalty_weight: float,
    sop_profile_path: Path | None,
    upstream_context: str | None,
    temperature_c: float,
    timeout_seconds: float,
    structure_backend: str,
    enable_esmfold: bool,
    enable_esm2: bool,
    output_dir: Path,
) -> None:
    """Phase B/C: Verify review.json, run paired 1DP/3DP optimization, and export artifacts atomically."""
    try:
        res = run_nbfold_optimize(
            protein_fasta_path=protein_fasta,
            structure_file_path=structure_file,
            annotation_path=annotation_path,
            review_manifest_path=review_manifest_path,
            output_dir=output_dir,
            declared_format=declared_format,
            chain_id=chain_id,
            model_number=model_number,
            provider_family=provider_family,
            confidence_sidecar_path=confidence_sidecar,
            pae_sidecar_path=pae_sidecar,
            export_metadata_path=export_metadata,
            interpro_sidecar_path=interpro_sidecar,
            uniprot_sidecar_path=uniprot_sidecar,
            host=host,
            gc_min=gc_min,
            gc_max=gc_max,
            target_adaptiveness=target_adaptiveness,
            penalty_weight=penalty_weight,
            sop_profile_path=sop_profile_path,
            upstream_context=upstream_context,
            temperature_c=temperature_c,
            timeout_seconds=timeout_seconds,
            structure_backend=structure_backend,
            enable_esmfold=enable_esmfold,
            enable_esm2=enable_esm2,
        )
        obs = res["observations"]
        summary = {
            "status": "COMPLETED",
            "policy_status": obs["policy_summary"]["policy_status"],
            "policy_decision": obs["policy_summary"]["policy_decision"],
            "automated_decision": obs["policy_summary"]["automated_decision"],
            "yield_prediction_status": obs["yield_prediction"]["status"],
            "conflict_count": res["run_manifest"]["review_provenance"]["conflict_count"],
            "artifacts": res["artifacts"],
            "runtime_ms": res["run_manifest"]["runtime_ms"],
        }
        click.echo(json.dumps(summary))
    except NBFoldError as err:
        err_payload = {
            "status": "ERROR",
            "stage": err.stage,
            "error_code": err.error_code,
            "exit_code": err.exit_code,
            "message": err.message,
        }
        click.echo(json.dumps(err_payload), err=True)
        sys.exit(err.exit_code)
