"""Click CLI subcommands for FactorForge structure intake, sidecars, and Mol* viewer integration."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import click

from factorforge.structure.interpro_sidecar_fetcher import fetch_interpro_sidecar
from factorforge.structure.uniprot_sidecar_fetcher import fetch_uniprot_sidecar
from factorforge.structure.viewer_bridge import (
    build_structure_position_mapping,
    prepare_viewer_annotations,
)


@click.group("structure")
def structure_cli() -> None:
    """Structure intake, UniProt/InterPro automated sidecars, and Mol* viewer integration."""


@structure_cli.command("fetch-uniprot")
@click.option("--accession", "-a", required=True, help="UniProt accession (e.g. P21926).")
@click.option("--sequence", "-s", default=None, help="Optional protein sequence to verify identity.")
@click.option("--output", "-o", type=click.Path(path_type=Path), help="Output path for JSON sidecar.")
@click.option("--cache-dir", type=click.Path(path_type=Path), default=None, help="Local cache directory.")
def fetch_uniprot_cmd(
    accession: str,
    sequence: Optional[str],
    output: Optional[Path],
    cache_dir: Optional[Path],
) -> None:
    """Fetch and normalize UniProt experimental/curated features into FactorForge sidecar."""
    result = fetch_uniprot_sidecar(
        accession=accession,
        target_protein_sequence=sequence,
        cache_dir=cache_dir,
    )
    payload_str = json.dumps(result, indent=2)
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(payload_str, encoding="utf-8")
        click.echo(f"UniProt sidecar saved to {output} (status: {result.get('status')})")
    else:
        click.echo(payload_str)


@structure_cli.command("fetch-interpro")
@click.option("--accession", "-a", required=True, help="UniProt accession for InterPro entry search.")
@click.option("--sequence", "-s", required=True, help="Target protein sequence for coordinate bounds.")
@click.option("--output", "-o", type=click.Path(path_type=Path), help="Output path for JSON sidecar.")
@click.option("--cache-dir", type=click.Path(path_type=Path), default=None, help="Local cache directory.")
def fetch_interpro_cmd(
    accession: str,
    sequence: str,
    output: Optional[Path],
    cache_dir: Optional[Path],
) -> None:
    """Fetch and normalize InterPro domain families into FactorForge sidecar."""
    result = fetch_interpro_sidecar(
        accession=accession,
        target_protein_sequence=sequence,
        cache_dir=cache_dir,
    )
    payload_str = json.dumps(result, indent=2)
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(payload_str, encoding="utf-8")
        click.echo(f"InterPro sidecar saved to {output} (status: {result.get('status')})")
    else:
        click.echo(payload_str)


@structure_cli.command("annotate")
@click.option("--accession", "-a", required=True, help="UniProt accession (e.g. P21926).")
@click.option("--structure-file", "--cif", type=click.Path(path_type=Path, exists=True), default=None, help="3D structure file (.cif, .mmcif, .pdb).")
@click.option("--chain", default="A", help="Structure chain ID (default: A).")
@click.option("--sequence", "-s", default=None, help="Optional protein sequence.")
@click.option("--interpro/--no-interpro", default=True, help="Whether to fetch InterPro domains.")
@click.option("--output", "-o", type=click.Path(path_type=Path), help="Output JSON path for consolidated annotations and viewer payload.")
@click.option("--cache-dir", type=click.Path(path_type=Path), default=None, help="Cache directory.")
def annotate_cmd(
    accession: str,
    structure_file: Optional[Path],
    chain: str,
    sequence: Optional[str],
    interpro: bool,
    output: Optional[Path],
    cache_dir: Optional[Path],
) -> None:
    """Consolidate UniProt, InterPro, and 3D structure coordinates into verified Mol* viewer payload."""
    click.echo(f"Fetching UniProt features for '{accession}'...")
    uniprot_res = fetch_uniprot_sidecar(accession=accession, target_protein_sequence=sequence, cache_dir=cache_dir)

    target_seq = sequence
    if not target_seq and uniprot_res.get("status") == "AVAILABLE":
        target_seq = uniprot_res.get("sequence", {}).get("value")

    interpro_res = None
    if interpro and target_seq:
        click.echo(f"Fetching InterPro domains for '{accession}'...")
        interpro_res = fetch_interpro_sidecar(accession=accession, target_protein_sequence=target_seq, cache_dir=cache_dir)

    # Merge candidate regions
    all_candidates = []
    if uniprot_res.get("status") == "AVAILABLE":
        all_candidates.extend(uniprot_res.get("candidate_regions", []))
    if interpro_res and interpro_res.get("status") == "AVAILABLE":
        all_candidates.extend(interpro_res.get("candidate_regions", []))

    # Resolve 3D coordinates against structure if provided
    position_mapping = None
    if structure_file and target_seq:
        click.echo(f"Mapping coordinates against structure '{structure_file.name}' chain '{chain}'...")
        struct_text = structure_file.read_text(encoding="utf-8")
        fmt = structure_file.suffix.lstrip(".").lower()
        try:
            position_mapping = build_structure_position_mapping(
                protein_sequence=target_seq,
                structure_content=struct_text,
                format=fmt,
                chain_id=chain,
            )
            click.echo(f"Successfully mapped {len(position_mapping)} residues to structure auth_seq_id.")
        except Exception as exc:
            click.echo(f"Warning: Structure coordinate mapping failed ({exc}). Falling back to unverified sequence indices.")

    viewer_annotations = prepare_viewer_annotations(
        {"candidate_regions": all_candidates},
        auth_asym_id=chain,
        position_mapping=position_mapping,
    )

    consolidated = {
        "accession": accession,
        "protein_length": len(target_seq) if target_seq else None,
        "structure_file": str(structure_file) if structure_file else None,
        "chain_id": chain,
        "coordinate_mapping_verified": position_mapping is not None,
        "uniprot": uniprot_res,
        "interpro": interpro_res,
        "candidate_regions_total": len(all_candidates),
        "viewer_annotations": viewer_annotations,
        "reviewed": False,  # Mandatory human review invariant
        "epistemic_status": "UNVERIFIED_PUBLIC_CORPUS",
        "biological_validation": "NOT_AVAILABLE",
    }

    out_json = json.dumps(consolidated, indent=2)
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(out_json, encoding="utf-8")
        click.echo(f"Consolidated annotation report & viewer payload saved to {output}")
    else:
        click.echo(out_json)
