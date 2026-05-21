"""
FactorForge CLI

Usage:
  factorforge optimize input.fasta -e v2 -p balanced -o output.fasta
  factorforge optimize input.fasta -e v2 -p balanced --template standard_expression -o output.gb --format genbank
  factorforge list-engines
"""

from pathlib import Path
import sys

import click

from factorforge import __version__
from factorforge.engines.registry import EngineRegistry
from factorforge.engines.v2.utils import parse_fasta_records


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
    default="v2",
    type=click.Choice(["v2"], case_sensitive=False),
    help="Engine (v2)",
)
@click.option("--profile", "-p", default="balanced", help="Optimization profile")
@click.option("--template", "construct_template", help="Construct template name")
@click.option("--output", "-o", help="Output file")
@click.option("--format", "output_format", default="fasta", help="Output format (fasta, genbank)")
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
    profile,
    construct_template,
    output,
    output_format,
    scan_mode,
    scan_include,
    scan_exclude,
):
    """Optimize protein sequence"""
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

        if fasta_records is not None and len(fasta_records) > 1:
            if construct_template:
                raise ValueError("Multi-FASTA input does not support --template mode.")
            if output_format.lower() != "fasta":
                raise ValueError("Multi-FASTA input only supports FASTA output.")

            optimizer = EngineRegistry.get(engine)
            payload = [{"id": seq_id, "sequence": seq} for seq_id, seq in fasta_records]
            if hasattr(optimizer, "optimize_batch"):
                results = optimizer.optimize_batch(
                    payload,
                    profile=profile,
                    scan_mode=scan_mode,
                    scan_include=scan_include_list,
                    scan_exclude=scan_exclude_list,
                )
            else:
                results = [
                    optimizer.optimize(
                        seq,
                        profile=profile,
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

        if engine == "v2" and construct_template:
            from factorforge.engines.v2.pipeline import OptimizationPipeline

            pipeline = OptimizationPipeline(profile=profile, construct_template=construct_template)
            result = pipeline.run(
                sequence,
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
                raise ValueError("Non-FASTA output requires --template with v2 pipeline.")

            # Get engine
            optimizer = EngineRegistry.get(engine)

            # Optimize
            click.echo(f"Optimizing with {optimizer.name} v{optimizer.version}...")
            result = optimizer.optimize(
                sequence,
                profile=profile,
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
