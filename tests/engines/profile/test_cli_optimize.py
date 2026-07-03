"""CLI tests for profile optimize command."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

from click.testing import CliRunner

# Add project src to path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.cli.main import cli
from factorforge.engines.profile.utils import parse_fasta_records


def _patch_profile_optimizer(monkeypatch, sequence: str) -> None:
    class FakeProfileOptimizer:
        name = "Profile-based"
        version = "test"

        def optimize(self, *_args, **_kwargs):
            return SimpleNamespace(
                sequence=sequence,
                metrics={"cai": 0.9, "gc_percent": 50.0, "score": 0.8},
            )

    monkeypatch.setattr(
        "factorforge.cli.main.EngineRegistry.get",
        lambda engine: FakeProfileOptimizer(),
    )


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
            "profile",
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
            "profile",
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
    assert "target_cai=" not in result.output
    assert "  - target_cai:" not in result.output
    assert "recommendation_reason" in result.output


def test_cli_dp_explicit_cai_target_is_echoed(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = tmp_path / "input.fasta"
    input_file.write_text(">test\nMSKGEELFTGVVPILVELD\n", encoding="utf-8")

    result = runner.invoke(
        cli,
        [
            "optimize",
            str(input_file),
            "--engine",
            "dp",
            "--cai-target",
            "0.99",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "target_cai=0.990" in result.output
    assert "  - target_cai: 0.990" in result.output


def test_cli_optimize_profile_engine(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = tmp_path / "input.fasta"
    input_file.write_text(">test\nMSKGEELFTGVVPILVELD\n", encoding="utf-8")

    result = runner.invoke(
        cli,
        [
            "optimize",
            str(input_file),
            "--engine",
            "profile",
            "--profile",
            "gc_target",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Optimizing with Profile-based v3.3.0" in result.output
    assert "Metrics:" in result.output


def test_cli_profile_single_output_writes_parseable_fasta_with_input_id(
    monkeypatch,
    tmp_path: Path,
) -> None:
    expected_sequence = "ATGAAACCCGGGTTT"
    _patch_profile_optimizer(monkeypatch, expected_sequence)
    runner = CliRunner()
    input_file = tmp_path / "input.fasta"
    output_file = tmp_path / "output.fasta"
    input_file.write_text(">source_record description\nMSKGEELFTGVVPILVELD\n", encoding="utf-8")

    result = runner.invoke(
        cli,
        [
            "optimize",
            str(input_file),
            "--engine",
            "profile",
            "--profile",
            "balanced",
            "--scan-mode",
            "fast",
            "--output",
            str(output_file),
        ],
    )

    assert result.exit_code == 0, result.output
    content = output_file.read_text(encoding="utf-8")
    assert content.startswith(">source_record|profile=balanced|")
    assert "|cai=" in content.splitlines()[0]
    assert "|gc=" in content.splitlines()[0]
    assert "|score=" in content.splitlines()[0]
    records = parse_fasta_records(content)
    assert len(records) == 1
    assert records[0][1] == expected_sequence


def test_cli_profile_single_raw_input_output_uses_file_stem_id(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = tmp_path / "raw_protein.txt"
    output_file = tmp_path / "output.fasta"
    input_file.write_text("MSKGEELFTGVVPILVELD", encoding="utf-8")

    result = runner.invoke(
        cli,
        [
            "optimize",
            str(input_file),
            "--engine",
            "profile",
            "--profile",
            "balanced",
            "--scan-mode",
            "fast",
            "--output",
            str(output_file),
        ],
    )

    assert result.exit_code == 0, result.output
    assert output_file.read_text(encoding="utf-8").startswith(">raw_protein|profile=balanced|")


def test_cli_profile_agentops_output_contract_is_fasta(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = tmp_path / "input.faa"
    output_file = tmp_path / "output.fasta"
    input_file.write_text(">agentops_input\nMSKGEELFTGVVPILVELD\n", encoding="utf-8")

    result = runner.invoke(
        cli,
        [
            "optimize",
            str(input_file),
            "--engine",
            "profile",
            "--profile",
            "assembly_friendly",
            "--scan-mode",
            "fast",
            "--output",
            str(output_file),
        ],
    )

    assert result.exit_code == 0, result.output
    content = output_file.read_text(encoding="utf-8")
    assert content.startswith(">agentops_input|profile=assembly_friendly|")
    records = parse_fasta_records(content)
    assert len(records) == 1
    assert records[0][1]
