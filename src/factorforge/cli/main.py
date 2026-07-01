"""
FactorForge CLI

Usage:
  factorforge optimize input.fasta -e profile -p balanced -o output.fasta
  factorforge optimize input.fasta -e profile -p balanced --template standard_expression -o output.gb --format genbank
  factorforge list-engines
"""

from pathlib import Path
import hashlib
import json
import sys

import click

from factorforge import __version__
from factorforge.engines.registry import EngineRegistry
from factorforge.engines.profile.utils import parse_fasta_records

HOST_MAP = {"nbenthamiana": "nbenthamiana", "by2": "ntabacum"}
HOST_TAXIDS = {"nbenthamiana": 4100, "ntabacum": 4097}
PACKAGE_ROOT = Path(__file__).resolve().parents[1]
FACTORFORGE_REPO_ROOT = Path(__file__).resolve().parents[3]
REFERENCE_POLICY_MANIFEST_PATH = (
    FACTORFORGE_REPO_ROOT / "data" / "reference" / "reference_policy_manifest.json"
)
BUNDLED_REFERENCE_POLICY_MANIFEST_PATH = PACKAGE_ROOT / "data" / "reference" / "reference_policy_manifest.json"


def _reference_policy_manifest_path() -> Path:
    """Return the repo manifest during development, or the packaged manifest in wheels."""
    if REFERENCE_POLICY_MANIFEST_PATH.exists():
        return REFERENCE_POLICY_MANIFEST_PATH
    return BUNDLED_REFERENCE_POLICY_MANIFEST_PATH


def _resolve_packaged_or_repo_path(relative_path: str) -> Path:
    """Resolve manifest paths in a source checkout or an installed wheel."""
    repo_path = FACTORFORGE_REPO_ROOT / relative_path
    if repo_path.exists():
        return repo_path

    package_prefix = "src/factorforge/"
    if relative_path.startswith(package_prefix):
        return PACKAGE_ROOT / relative_path.removeprefix(package_prefix)
    return PACKAGE_ROOT / relative_path


def _load_reference_policy_manifest() -> dict:
    """Load the checksum/tier policy manifest for expert CLI reference selection."""
    return json.loads(_reference_policy_manifest_path().read_text(encoding="utf-8"))


def _reference_entries_by_id() -> dict[str, dict]:
    manifest = _load_reference_policy_manifest()
    return {entry["reference_id"]: entry for entry in manifest["references"]}


def _reference_id_choices() -> tuple[str, ...]:
    return tuple(_reference_entries_by_id())


REFERENCE_ID_CHOICES = _reference_id_choices()


def _reference_entry_by_id(reference_id: str) -> dict:
    entries = _reference_entries_by_id()
    try:
        return entries[reference_id]
    except KeyError as exc:
        choices = ", ".join(sorted(entries))
        raise click.UsageError(
            f"Unknown reference_id {reference_id!r}. Supported values: {choices}"
        ) from exc


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resolve_reference_by_id(reference_id: str) -> Path:
    """Resolve a manifest reference_id to a checksum-verified codon table path."""
    entry = _reference_entry_by_id(reference_id)
    codon_table_path = _resolve_packaged_or_repo_path(str(entry["codon_table_path"]))
    expected = str(entry["checksum_sha256"])
    if not codon_table_path.exists():
        raise click.UsageError(
            f"Codon table file for {reference_id} does not exist: {codon_table_path}"
        )
    actual = _sha256_file(codon_table_path)
    if actual != expected:
        raise click.UsageError(
            "Checksum mismatch for "
            f"{reference_id} at {codon_table_path}: expected {expected}, actual {actual}"
        )

    if entry["tier"] != "production_enabled":
        limitations = "; ".join(entry.get("known_limitations", []))
        warning = (
            f"Warning: reference_id={reference_id} has tier={entry['tier']}; "
            f"{entry['claim_boundary']}"
        )
        if limitations:
            warning = f"{warning} Known limitations: {limitations}"
        click.echo(warning, err=True)

    return codon_table_path


def _validate_reference_host(reference_id: str, internal_host: str) -> None:
    entry = _reference_entry_by_id(reference_id)
    expected_taxid = HOST_TAXIDS[internal_host]
    actual_taxid = int(entry["ncbi_taxid"])
    if actual_taxid != expected_taxid:
        raise click.UsageError(
            f"reference_id={reference_id} targets {entry['organism']} "
            f"(NCBI taxid {actual_taxid}) and is incompatible with "
            f"--host {internal_host} (expected NCBI taxid {expected_taxid})."
        )


def _configure_stdio() -> None:
    """Best-effort UTF-8 for Windows consoles."""
    for stream in (sys.stdout, sys.stderr):
        try:
            reconfigure = getattr(stream, "reconfigure", None)
            if callable(reconfigure):
                reconfigure(encoding="utf-8")
        except Exception:
            pass


def _parse_csv_option(value):
    """Parse comma-separated option values."""
    if not value:
        return None
    parsed = [item.strip() for item in value.split(",") if item.strip()]
    return parsed or None


def _wrap_sequence(sequence, width=80):
    """Wrap sequence to fixed-width lines."""
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def _build_dp_result(
    sequence: str,
    objective: str,
    gc_min: float,
    gc_max: float,
    codon_table_path: Path | None = None,
):
    """Run the constraint-based DP feasibility engine for a single protein sequence."""
    if objective != "feasibility_best":
        raise ValueError("DP engine currently supports --objective feasibility_best.")
    if gc_min > gc_max:
        raise ValueError("--gc-min must be <= --gc-max.")

    from factorforge.analysis.metrics import load_codon_usage_table
    from factorforge.analysis.feasibility import analyze_feasibility

    table = load_codon_usage_table(path=codon_table_path)
    result = analyze_feasibility(
        sequence,
        table.codon_weights,
        target_gc_low=gc_min,
        target_gc_high=gc_max,
    )
    best = result["target"]["best_candidate"]
    feasible = best is not None
    if best is None:
        best = result["best_candidate_without_gc"]
    if best is None:
        raise ValueError("No DP candidate generated.")

    reason = (
        f"Maximum CAI under GC {gc_min:g}-{gc_max:g}%"
        if feasible
        else "Maximum CAI without GC constraint; requested GC range was infeasible"
    )
    return best, result, reason


def _format_dp_fasta(sequence_id: str, dna_sequence: str, cai: float, gc: float) -> str:
    """Format a DP result as FASTA."""
    header = f">{sequence_id}|engine=dp|objective=feasibility_best|cai={cai:.3f}|gc={gc:.2f}"
    return f"{header}\n{_wrap_sequence(dna_sequence)}\n"


def _option_was_explicitly_set(option_name: str) -> bool:
    """Return whether an option was provided on the command line."""
    ctx = click.get_current_context(silent=True)
    if ctx is None:
        return False
    get_parameter_source = getattr(ctx, "get_parameter_source", None)
    if not callable(get_parameter_source):
        return False
    return get_parameter_source(option_name) == click.core.ParameterSource.COMMANDLINE


def _engine_option_was_explicitly_set() -> bool:
    """Return whether --engine/-e was provided on the command line."""
    return _option_was_explicitly_set("engine")


def _format_profile_fasta(sequence_id: str, profile: str, result) -> str:
    """Format a profile optimization result as FASTA."""
    cai = float(result.metrics.get("cai", 0.0))
    gc = float(result.metrics.get("gc_percent", result.metrics.get("gc_content", 0.0)))
    score = float(result.metrics.get("score", 0.0))
    header = f">{sequence_id}|profile={profile}|cai={cai:.3f}|gc={gc:.2f}|score={score:.3f}"
    return f"{header}\n{_wrap_sequence(result.sequence)}\n"


def _format_profile_comparison_table(profile_results) -> str:
    """Format profile optimization metrics as a comparison table."""
    divider = "─" * 45
    lines = [
        "Profile comparison results:",
        divider,
        f"{'Profile':<18}{'CAI':>7} {'GC%':>7} {'Score':>8}",
        divider,
    ]
    for profile_name, result in profile_results:
        cai = float(result.metrics.get("cai", 0.0))
        gc = float(result.metrics.get("gc_percent", result.metrics.get("gc_content", 0.0)))
        score = float(result.metrics.get("score", 0.0))
        lines.append(f"{profile_name:<18}{cai:>7.3f} {gc:>7.2f} {score:>8.3f}")
    lines.append(divider)
    return "\n".join(lines)


@click.group()
@click.version_option(version=__version__)
def cli():
    """FactorForge - Codon Optimization Platform"""
    _configure_stdio()


@cli.command()
def list_engines():
    """List available optimization engines"""
    engines = EngineRegistry.list_engines()

    click.echo("\nAvailable Engines:\n")
    for name, info in engines.items():
        click.echo(f"  - {name}: {info['name']} v{info['version']}")
    click.echo()


@cli.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "--engine",
    "-e",
    default="dp",
    type=click.Choice(["dp", "profile"], case_sensitive=False),
    help="Engine (dp, profile)",
)
@click.option(
    "--host",
    default="nbenthamiana",
    type=click.Choice(["nbenthamiana", "by2"], case_sensitive=False),
    help="Expression host: nbenthamiana (default) or by2 (Tobacco BY-2 / N. tabacum)",
)
@click.option("--profile", "-p", default="balanced", help="Optimization profile")
@click.option(
    "--objective",
    default="feasibility_best",
    type=click.Choice(["feasibility_best"], case_sensitive=False),
    help="DP objective",
)
@click.option("--gc-min", type=float, default=40.0, help="Minimum target GC percentage")
@click.option("--gc-max", type=float, default=47.0, help="Maximum target GC percentage")
@click.option("--template", "construct_template", help="Construct template name")
@click.option("--output", "-o", help="Output file")
@click.option("--format", "output_format", default="fasta", help="Output format (fasta, genbank)")
@click.option(
    "--reference-id",
    type=click.Choice(REFERENCE_ID_CHOICES, case_sensitive=False),
    default=None,
    help="Expert/research codon-reference ID; checksum-validated.",
)
@click.option(
    "--compare-profiles",
    help=(
        "Comma-separated profiles to compare "
        "(e.g. balanced,high_cai,gc_target). Implies --engine profile."
    ),
)
@click.option(
    "--scan-mode",
    default="full",
    type=click.Choice(["full", "fast"], case_sensitive=False),
    help="Rule scan mode",
)
@click.option("--scan-include", help="Comma-separated scanner names to include")
@click.option("--scan-exclude", help="Comma-separated scanner names to exclude")
def optimize(
    input_file,
    engine,
    host,
    profile,
    objective,
    gc_min,
    gc_max,
    construct_template,
    output,
    output_format,
    reference_id,
    compare_profiles,
    scan_mode,
    scan_include,
    scan_exclude,
):
    """Optimize protein sequence"""
    compare_profile_list = _parse_csv_option(compare_profiles)
    engine = engine.lower()
    host_value = host.lower()
    internal_host = HOST_MAP[host_value]
    host_was_explicit = _option_was_explicitly_set("host")
    reference_id = reference_id.lower() if reference_id else None
    reference_table_path = None

    if reference_id is not None:
        _validate_reference_host(reference_id, internal_host)
        reference_table_path = resolve_reference_by_id(reference_id)

    if host_was_explicit and engine == "dp" and _engine_option_was_explicitly_set():
        raise click.UsageError("--host is only supported with --engine profile")
    if host_was_explicit and engine == "dp" and internal_host != "nbenthamiana":
        engine = "profile"

    if host_was_explicit and internal_host != "nbenthamiana":
        if profile == "high_cai" or "high_cai" in (compare_profile_list or []):
            raise click.UsageError(
                "high_cai requires the N. benthamiana golden-set reference "
                "and is not available for non-default hosts."
            )

    if compare_profile_list:
        if engine == "dp" and _engine_option_was_explicitly_set():
            raise click.UsageError("--compare-profiles cannot be used with --engine dp.")
        engine = "profile"

    if reference_table_path is not None and construct_template:
        raise click.UsageError("--reference-id is not supported with --template mode.")

    try:
        # Read file
        with open(input_file, encoding="utf-8") as f:
            raw_input = f.read()

        scan_include_list = _parse_csv_option(scan_include)
        scan_exclude_list = _parse_csv_option(scan_exclude)

        fasta_records = None
        sequence = raw_input.strip()
        if raw_input.lstrip().startswith(">"):
            fasta_records = parse_fasta_records(raw_input)
            if len(fasta_records) == 1:
                sequence = fasta_records[0][1]

        if compare_profile_list:
            if fasta_records is not None and len(fasta_records) > 1:
                raise ValueError("Profile comparison requires a single input sequence.")
            if construct_template:
                raise ValueError("Profile comparison does not support --template mode.")
            if output_format.lower() != "fasta":
                raise ValueError("Profile comparison only supports FASTA output.")

            if reference_table_path is not None:
                from factorforge.engines.profile.optimizer import RuleBasedOptimizer

                optimizer = RuleBasedOptimizer(codon_table_path=str(reference_table_path))
            else:
                optimizer = EngineRegistry.get("profile")
            profile_results = []
            for profile_name in compare_profile_list:
                result = optimizer.optimize(
                    sequence,
                    profile=profile_name,
                    host=internal_host,
                    scan_mode=scan_mode,
                    scan_include=scan_include_list,
                    scan_exclude=scan_exclude_list,
                )
                profile_results.append((profile_name, result))

            if output:
                first_profile, first_result = profile_results[0]
                sequence_id = Path(input_file).stem or "factorforge_profile"
                fasta = _format_profile_fasta(sequence_id, first_profile, first_result)
                with open(output, "w", encoding="utf-8") as f:
                    f.write(fasta)
                click.echo(f"Saved to: {output}")

            click.echo(_format_profile_comparison_table(profile_results))
            return

        if fasta_records is not None and len(fasta_records) > 1:
            if engine == "dp":
                raise ValueError("Multi-FASTA input requires --engine profile.")
            if construct_template:
                raise ValueError("Multi-FASTA input does not support --template mode.")
            if output_format.lower() != "fasta":
                raise ValueError("Multi-FASTA input only supports FASTA output.")

            if reference_table_path is not None and engine == "profile":
                from factorforge.engines.profile.optimizer import RuleBasedOptimizer

                optimizer = RuleBasedOptimizer(codon_table_path=str(reference_table_path))
            else:
                optimizer = EngineRegistry.get(engine)
            payload = [{"id": seq_id, "sequence": seq} for seq_id, seq in fasta_records]
            if hasattr(optimizer, "optimize_batch"):
                results = optimizer.optimize_batch(
                    payload,
                    profile=profile,
                    host=internal_host,
                    scan_mode=scan_mode,
                    scan_include=scan_include_list,
                    scan_exclude=scan_exclude_list,
                )
            else:
                results = [
                    optimizer.optimize(
                        seq,
                        profile=profile,
                        host=internal_host,
                        scan_mode=scan_mode,
                        scan_include=scan_include_list,
                        scan_exclude=scan_exclude_list,
                    )
                    for _id, seq in fasta_records
                ]

            combined_fasta = []
            for idx, result in enumerate(results):
                seq_id = payload[idx]["id"]
                cai = result.metrics.get("cai", 0.0)
                gc = result.metrics.get("gc_percent", result.metrics.get("gc_content", 0.0))
                score = result.metrics.get("score", 0.0)
                header = (
                    f">{seq_id}|profile={profile}|cai={float(cai):.3f}|"
                    f"gc={float(gc):.2f}|score={float(score):.3f}"
                )
                combined_fasta.append(f"{header}\n{_wrap_sequence(result.sequence)}")
            out_content = "\n".join(combined_fasta) + "\n"

            if output:
                with open(output, "w", encoding="utf-8") as f:
                    f.write(out_content)
                click.echo(f"Saved batch FASTA to: {output}")
            else:
                click.echo(f"\n{out_content}")
            click.echo(f"Batch optimized: {len(results)} sequences")
            return

        if engine == "dp":
            if construct_template:
                raise ValueError("DP engine does not support --template mode.")
            if output_format.lower() != "fasta":
                raise ValueError("DP engine only supports FASTA output.")

            best, feasibility, recommendation_reason = _build_dp_result(
                sequence,
                objective=objective,
                gc_min=gc_min,
                gc_max=gc_max,
                codon_table_path=reference_table_path,
            )
            dna_sequence = best["dna_sequence"]
            cai = float(best["cai"])
            gc = float(best["gc"])
            sequence_id = Path(input_file).stem or "factorforge_dp"
            fasta = _format_dp_fasta(sequence_id, dna_sequence, cai, gc)

            click.echo("Optimizing with DP feasibility engine...")
            if output:
                with open(output, "w", encoding="utf-8") as f:
                    f.write(fasta)
                click.echo(f"Saved to: {output}")
            else:
                click.echo(f"\n{fasta}")

            click.echo("Metrics:")
            click.echo(f"  - cai: {cai:.3f}")
            click.echo(f"  - gc_percent: {gc:.2f}")
            click.echo(f"  - target_gc_min: {float(feasibility['target']['gc_low']):.2f}")
            click.echo(f"  - target_gc_max: {float(feasibility['target']['gc_high']):.2f}")
            click.echo(f"  - target_feasible: {bool(feasibility['target']['best_candidate'])}")
            click.echo(f"  - recommendation_reason: {recommendation_reason}")
            return

        if engine == "profile" and construct_template:
            from factorforge.engines.profile.pipeline import OptimizationPipeline

            pipeline = OptimizationPipeline(
                profile=profile,
                construct_template=construct_template,
                host=internal_host,
            )
            result = pipeline.run(
                sequence,
                host=internal_host,
                scan_mode=scan_mode,
                scan_include=scan_include_list,
                scan_exclude=scan_exclude_list,
            )

            if output_format.lower() == "genbank" and not output:
                raise ValueError("GenBank output requires --output file path.")

            if output:
                result.save(Path(output), format=output_format)
                click.echo(f"Saved to: {output}")
            else:
                click.echo(f"\n{result.sequence}\n")

            click.echo("Metrics:")
            for key, value in result.metadata.get("metrics", {}).items():
                click.echo(f"  - {key}: {value}")
        else:
            if output_format.lower() != "fasta":
                raise ValueError("Non-FASTA output requires --template with profile pipeline.")

            # Get engine
            if reference_table_path is not None and engine == "profile":
                from factorforge.engines.profile.optimizer import RuleBasedOptimizer

                optimizer = RuleBasedOptimizer(codon_table_path=str(reference_table_path))
            else:
                optimizer = EngineRegistry.get(engine)

            # Optimize
            click.echo(f"Optimizing with {optimizer.name} v{optimizer.version}...")
            result = optimizer.optimize(
                sequence,
                profile=profile,
                host=internal_host,
                scan_mode=scan_mode,
                scan_include=scan_include_list,
                scan_exclude=scan_exclude_list,
            )

            # Output results
            if output:
                with open(output, "w", encoding="utf-8") as f:
                    f.write(result.sequence)
                click.echo(f"Saved to: {output}")
            else:
                click.echo(f"\n{result.sequence}\n")

            # Output metrics
            click.echo("Metrics:")
            for key, value in result.metrics.items():
                click.echo(f"  - {key}: {value}")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()


if __name__ == "__main__":
    cli()
