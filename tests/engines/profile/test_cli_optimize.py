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


def _header_fields(header_line: str) -> dict[str, str]:
    assert header_line.startswith(">")
    parts = header_line[1:].split("|")
    fields = {"id": parts[0]}
    fields.update(part.split("=", 1) for part in parts[1:])
    return fields


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
    assert ">a|engine=profile|profile=balanced|" in content
    assert ">b|engine=profile|profile=balanced|" in content


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
    assert "Optimizing with Profile-based v3.4.1" in result.output
    assert "Metrics:" in result.output


def test_profile_cli_explicit_reference_reports_same_reference_not_fallback(
    tmp_path: Path,
) -> None:
    runner = CliRunner()
    input_file = tmp_path / "input.fasta"
    input_file.write_text(">test\nMSKGEELF\n", encoding="utf-8")

    result = runner.invoke(
        cli,
        [
            "optimize",
            str(input_file),
            "--engine",
            "profile",
            "--reference-id",
            "nbenthamiana_nbev11_hc_v2",
            "--scan-mode",
            "fast",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "'reference_id': 'nbenthamiana_nbev11_hc_v2'" in result.output
    assert "'reference_relationship': 'same_as_generation_reference'" in result.output
    assert "'fallback_used': False" in result.output
    assert "fallback_to_generation_reference" not in result.output


def test_cli_template_profile_metrics_include_cai_authority(tmp_path: Path) -> None:
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
            "balanced",
            "--template",
            "standard_expression",
            "--scan-mode",
            "fast",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Metrics:" in result.output
    assert "  - cai_authority:" in result.output
    assert "distinct_from_generation_reference" in result.output


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
    first_line = content.splitlines()[0]
    assert first_line == (
        ">source_record|engine=profile|profile=balanced|"
        "cai=0.900|gc=50.00|score=0.800"
    )
    fields = _header_fields(first_line)
    assert fields == {
        "id": "source_record",
        "engine": "profile",
        "profile": "balanced",
        "cai": "0.900",
        "gc": "50.00",
        "score": "0.800",
    }
    assert "objective" not in fields
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
    assert output_file.read_text(encoding="utf-8").startswith(
        ">raw_protein|engine=profile|profile=balanced|"
    )


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
    assert content.startswith(">agentops_input|engine=profile|profile=assembly_friendly|")
    fields = _header_fields(content.splitlines()[0])
    assert fields["engine"] == "profile"
    assert fields["profile"] == "assembly_friendly"
    assert {"cai", "gc", "score"} <= fields.keys()
    assert "objective" not in fields
    records = parse_fasta_records(content)
    assert len(records) == 1
    assert records[0][1]
