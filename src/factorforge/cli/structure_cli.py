"""Click CLI subcommands for FactorForge structure intake, sidecars, and Mol* viewer integration."""

from __future__ import annotations

import json
from pathlib import Path

import click

from factorforge.structure.alphafold_server_export import (
    AlphaFoldExportError,
    export_alphafold_server_job,
    parse_protein_fasta,
)
from factorforge.structure.alphafold_server_ingest import (
    AlphaFoldIngestError,
    format_ingest_report,
    ingest_alphafold_server_result,
)
from factorforge.structure.esmfold_client import (
    ESMFoldClientError,
    format_esm_report,
    predict_and_ingest_esm,
)
from factorforge.structure.interpro_sidecar_fetcher import fetch_interpro_sidecar
from factorforge.structure.sllm_adapter import (
    SLLMAdapterError,
    generate_evidence_explanation,
)
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
    sequence: str | None,
    output: Path | None,
    cache_dir: Path | None,
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
    output: Path | None,
    cache_dir: Path | None,
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
    structure_file: Path | None,
    chain: str,
    sequence: str | None,
    interpro: bool,
    output: Path | None,
    cache_dir: Path | None,
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
        except (ValueError, KeyError, OSError, RuntimeError) as exc:
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


@structure_cli.command("export-af3")
@click.option("--fasta", "-f", "fasta_file", type=click.Path(path_type=Path, exists=True, dir_okay=False), default=None, help="Input FASTA file with protein sequence(s).")
@click.option("--sequence", "-s", default=None, help="Direct protein sequence (for single chain).")
@click.option("--name", "-n", default=None, help="Optional non-confidential job name included in the manifest.")
@click.option("--counts", "-c", default=None, help="Comma-separated copy counts for chains (e.g. '1,1' or '2').")
@click.option("--seed", default=1, type=int, help="Model seed integer (default: 1).")
@click.option("--output", "-o", type=click.Path(path_type=Path), default=None, help="Output JSON path (e.g. job_request.json).")
@click.option("--overwrite", is_flag=True, default=False, help="Allow overwriting existing output file.")
def export_af3_cmd(
    fasta_file: Path | None,
    sequence: str | None,
    name: str | None,
    counts: str | None,
    seed: int,
    output: Path | None,
    overwrite: bool,
) -> None:
    """Export offline AlphaFold Server (alphafoldserver.com) job JSON payload.

    Writes protein sequences to JSON locally. Review the manifest and ensure
    permission to share its contents before uploading it to an external service.
    """
    if (fasta_file is None) == (sequence is None):
        raise click.UsageError("Provide either --fasta / -f or --sequence / -s, exactly one.")

    try:
        if fasta_file is not None and output is not None and (
            output.resolve() == fasta_file.resolve() or (output.exists() and output.samefile(fasta_file))
        ):
            raise click.UsageError("The output path must differ from the input FASTA.")
        raw_sequences = parse_protein_fasta(fasta_file) if fasta_file is not None else [sequence]
        parsed_counts = [1] * len(raw_sequences)
        if counts is not None:
            try:
                parsed_counts = [int(part) for part in counts.split(",")]
            except ValueError as exc:
                raise click.BadParameter(
                    "Invalid chain count. Provide positive integers separated by commas.",
                    param_hint="--counts",
                ) from exc
            if len(parsed_counts) != len(raw_sequences) or any(count < 1 for count in parsed_counts):
                raise click.BadParameter(
                    "Invalid chain count. Provide exactly one positive integer per input chain.",
                    param_hint="--counts",
                )
        chains = list(zip(raw_sequences, parsed_counts, strict=True))

        _, json_str = export_alphafold_server_job(
            chains=chains,
            output_path=output,
            job_name=name or "factorforge_af3_job",
            model_seeds=seed,
            overwrite=overwrite,
        )

        if output:
            click.echo(
                f"AlphaFold Server job JSON saved to {output} "
                f"({len(chains)} chain definition(s)).\n"
                "Epistemic: UNVERIFIED_PUBLIC_CORPUS · Biological validation: NOT_AVAILABLE"
            )
        else:
            click.echo(json_str, nl=False)

    except FileExistsError as fee:
        raise click.ClickException(str(fee))
    except AlphaFoldExportError as afe:
        raise click.ClickException(f"AlphaFold Server validation error: {afe}")
    except (OSError, UnicodeError) as exc:
        raise click.ClickException("Could not read the FASTA or write the output JSON.") from exc




@structure_cli.command("ingest-af3")
@click.option("--cif", "-c", "cif_file", type=click.Path(path_type=Path, exists=True), required=True, help="AlphaFold 3 mmCIF structure file (*_model_0.cif).")
@click.option("--chain", default="A", help="Chain ID to analyze (e.g. 'A', 'B', or 'ALL', default: A).")
@click.option("--summary-json", "-s", type=click.Path(path_type=Path, exists=True), default=None, help="Optional summary_confidences_0.json file.")
@click.option("--pause-window", default=35, type=click.IntRange(min=0), help="Hypothetical model residue offset (default: 35).")
@click.option("--output", "-o", type=click.Path(path_type=Path), default=None, help="Output JSON path for structured sidecar.")
@click.option("--json", "json_output", is_flag=True, default=False, help="Output structured result as JSON to stdout.")
def ingest_af3_cmd(
    cif_file: Path,
    chain: str,
    summary_json: Path | None,
    pause_window: int,
    output: Path | None,
    json_output: bool,
) -> None:
    """Ingest AlphaFold Server 3D structure (.cif) and confidence metrics.

    Reports CA confidence and hypothetical offsets in model coordinates.
    The report contains sequence and does not establish precursor or CDS mapping.
    """
    try:
        result = ingest_alphafold_server_result(
            cif_path=cif_file,
            chain=chain,
            summary_json_path=summary_json,
            output_path=output,
            tunnel_offset=pause_window,
        )

        if json_output:
            click.echo(json.dumps(result, indent=2))
        else:
            click.echo(format_ingest_report(result))

        if output and not json_output:
            click.echo(f"\nIngest sidecar saved to: {output}")

    except AlphaFoldIngestError as aie:
        raise click.ClickException(f"AlphaFold Server ingest error: {aie}")
    except FileExistsError as exc:
        raise click.ClickException("Output already exists. Choose a different path.") from exc
    except (OSError, UnicodeError) as exc:
        raise click.ClickException("Could not read the AF3 input or write the sidecar.") from exc


@structure_cli.command("predict-esm")
@click.argument("target", required=False, default=None)
@click.option("--fasta", "-f", "fasta_file", type=click.Path(path_type=Path, exists=True), default=None, help="Input protein FASTA file.")
@click.option("--sequence", "-s", default=None, help="Input protein amino acid sequence string.")
@click.option("--output-sidecar", "-o", type=click.Path(path_type=Path), default=None, help="Output JSON path for structured pause sidecar.")
@click.option("--output-pdb", type=click.Path(path_type=Path), default=None, help="Optional output path to save predicted 3D structure (.pdb).")
@click.option("--provider", type=click.Choice(["esmatlas", "nvidia", "mock"], case_sensitive=False), default="esmatlas", help="ESMFold inference backend provider (default: esmatlas).")
@click.option("--api-key", default=None, help="Credential for the configured NVIDIA endpoint.")
@click.option("--allow-remote", is_flag=True, help="Allow sending the input sequence to the selected provider.")
@click.option("--pause-window", default=35, type=click.IntRange(min=0), help="Ribosome exit-tunnel residue offset (default: 35).")
@click.option("--timeout", default=30.0, type=float, help="API request timeout in seconds (default: 30.0).")
@click.option("--json", "json_output", is_flag=True, default=False, help="Output structured result as JSON to stdout.")
def predict_esm_cmd(
    target: str | None,
    fasta_file: Path | None,
    sequence: str | None,
    output_sidecar: Path | None,
    output_pdb: Path | None,
    provider: str,
    api_key: str | None,
    allow_remote: bool,
    pause_window: int,
    timeout: float,
    json_output: bool,
) -> None:
    """Generate an experimental sidecar from ESMFold or an offline fixture."""
    if sum(value is not None for value in (target, fasta_file, sequence)) != 1:
        raise click.UsageError("Provide exactly one of TARGET, --fasta, or --sequence.")
    if provider != "mock" and not allow_remote:
        raise click.UsageError("Use --allow-remote to send the sequence to the selected provider.")

    try:
        target_path = Path(target) if target and len(target) <= 240 else None
        input_path = fasta_file or (target_path if target_path and target_path.is_file() else None)
        if input_path is not None:
            seqs = parse_protein_fasta(input_path)
            if len(seqs) != 1:
                raise click.UsageError("Provide a FASTA file with exactly one protein.")
            raw_seq = seqs[0]
        else:
            raw_seq = target if target is not None else sequence
        result, _pdb_text = predict_and_ingest_esm(
            sequence=raw_seq,
            provider=provider.lower(),
            api_key=api_key,
            tunnel_offset=pause_window,
            output_sidecar_path=output_sidecar,
            output_pdb_path=output_pdb,
            timeout=timeout,
        )

        if json_output:
            click.echo(json.dumps(result, indent=2))
        else:
            click.echo(format_esm_report(result))

        if output_sidecar and not json_output:
            click.echo(f"\nSidecar saved to: {output_sidecar}")
        if output_pdb and not json_output:
            click.echo(f"3D Structure saved to: {output_pdb}")

    except (ESMFoldClientError, AlphaFoldExportError) as exc:
        raise click.ClickException(f"ESMFold error: {exc}") from exc
    except (OSError, UnicodeError) as exc:
        raise click.ClickException(f"File I/O error: {exc}") from exc


@structure_cli.command("explain")
@click.option("--sidecar", "-s", "sidecar_path", type=click.Path(path_type=Path, exists=True), required=True, help="Structure sidecar JSON (from ingest-af3 or predict-esm).")
@click.option("--fasta", "-f", "fasta_path", type=click.Path(path_type=Path, exists=True), default=None, help="Optional optimized DNA FASTA to verify invariants against.")
@click.option("--allow-mock", is_flag=True, default=False, help="Allow synthetic/mock sidecars for testing.")
@click.option("--output", "-o", "output_path", type=click.Path(path_type=Path), default=None, help="Output Markdown report path.")
@click.option("--json", "json_output", is_flag=True, default=False, help="Output structured report as JSON to stdout.")
def explain_cmd(
    sidecar_path: Path,
    fasta_path: Path | None,
    allow_mock: bool,
    output_path: Path | None,
    json_output: bool,
) -> None:
    """Generate evidence-aware assessment strictly separated into Verified/Predicted/Hypothesized tiers."""
    try:
        report = generate_evidence_explanation(
            sidecar_source=sidecar_path,
            optimized_fasta_content=fasta_path,
            allow_mock=allow_mock,
        )

        if json_output:
            click.echo(json.dumps(report.to_dict(), indent=2))
        else:
            click.echo(report.to_markdown())

        if output_path:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(report.to_markdown(), encoding="utf-8")
            if not json_output:
                click.echo(f"\nReport saved to: {output_path}")

    except SLLMAdapterError as exc:
        raise click.ClickException(f"Evidence adapter error: {exc}") from exc
    except (OSError, UnicodeError) as exc:
        raise click.ClickException(f"File I/O error: {exc}") from exc

