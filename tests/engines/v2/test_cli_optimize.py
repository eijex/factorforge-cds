"""CLI tests for v2 optimize command."""

from __future__ import annotations

import sys
from pathlib import Path

from click.testing import CliRunner

# Add project src to path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.cli.main import cli


def test_cli_optimize_multifasta_batch(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = tmp_path / "input.fasta"
    output_file = tmp_path / "output.fasta"
    input_file.write_text(">a\nMAKL\n>b\nMAKL\n", encoding="utf-8")

    result = runner.invoke(
        cli,
        [
            "optimize",
            str(input_file),
            "--engine",
            "v2",
            "--profile",
            "balanced",
            "--scan-mode",
            "fast",
            "--output",
            str(output_file),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Batch optimized: 2 sequences" in result.output
    content = output_file.read_text(encoding="utf-8")
    assert content.count(">") == 2
    assert ">a|profile=balanced|" in content
    assert ">b|profile=balanced|" in content


def test_cli_multifasta_template_rejected(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = tmp_path / "input.fasta"
    input_file.write_text(">a\nMAKL\n>b\nMAKL\n", encoding="utf-8")

    result = runner.invoke(
        cli,
        [
            "optimize",
            str(input_file),
            "--engine",
            "v2",
            "--profile",
            "balanced",
            "--template",
            "standard_expression",
        ],
    )

    assert result.exit_code != 0
    assert "Multi-FASTA input does not support --template mode." in result.output


def test_cli_optimize_defaults_to_dp_engine(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = tmp_path / "input.fasta"
    input_file.write_text(">test\nMSKGEELFTGVVPILVELD\n", encoding="utf-8")

    result = runner.invoke(cli, ["optimize", str(input_file)])

    assert result.exit_code == 0, result.output
    assert "Optimizing with DP feasibility engine" in result.output
    assert ">input|engine=dp|objective=feasibility_best|" in result.output
    assert "recommendation_reason" in result.output


def test_cli_optimize_v2_profile_still_works(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = tmp_path / "input.fasta"
    input_file.write_text(">test\nMSKGEELFTGVVPILVELD\n", encoding="utf-8")

    result = runner.invoke(
        cli,
        [
            "optimize",
            str(input_file),
            "--engine",
            "v2",
            "--profile",
            "gc_target",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Optimizing with Rule-based v3.1.3" in result.output
    assert "Metrics:" in result.output
