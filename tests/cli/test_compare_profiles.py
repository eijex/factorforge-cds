"""CLI tests for profile comparison mode."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from factorforge.cli.main import cli


SAMPLE_PROTEIN = "MSKGEELFTGVVPILVELD"


def test_cli_compare_profiles_outputs_table(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = tmp_path / "input.fasta"
    input_file.write_text(f">test\n{SAMPLE_PROTEIN}\n", encoding="utf-8")

    result = runner.invoke(
        cli,
        [
            "optimize",
            str(input_file),
            "--compare-profiles",
            "balanced,high_cai",
            "--scan-mode",
            "fast",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Profile comparison results:" in result.output
    assert "Profile" in result.output
    assert "CAI" in result.output
    assert "GC%" in result.output
    assert "Score" in result.output
    assert "balanced" in result.output
    assert "high_cai" in result.output


def test_cli_compare_profiles_rejects_explicit_dp_engine(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = tmp_path / "input.fasta"
    input_file.write_text(SAMPLE_PROTEIN, encoding="utf-8")

    result = runner.invoke(
        cli,
        [
            "optimize",
            str(input_file),
            "--engine",
            "dp",
            "--compare-profiles",
            "balanced,high_cai",
        ],
    )

    assert result.exit_code != 0
    assert "--compare-profiles cannot be used with --engine dp." in result.output


def test_cli_compare_profiles_rejects_invalid_profile(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = tmp_path / "input.fasta"
    input_file.write_text(SAMPLE_PROTEIN, encoding="utf-8")

    result = runner.invoke(
        cli,
        [
            "optimize",
            str(input_file),
            "--compare-profiles",
            "balanced,bad_profile",
            "--scan-mode",
            "fast",
        ],
    )

    assert result.exit_code != 0
    assert "Unknown profile: bad_profile" in result.output


def test_cli_compare_profiles_writes_first_profile_fasta(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = tmp_path / "input.fasta"
    output_file = tmp_path / "output.fasta"
    input_file.write_text(SAMPLE_PROTEIN, encoding="utf-8")

    result = runner.invoke(
        cli,
        [
            "optimize",
            str(input_file),
            "--compare-profiles",
            "balanced,high_cai",
            "--scan-mode",
            "fast",
            "--output",
            str(output_file),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Profile comparison results:" in result.output
    assert output_file.exists()
    assert output_file.read_text(encoding="utf-8").startswith(
        ">input|engine=profile|profile=balanced|"
    )
