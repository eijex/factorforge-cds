"""CLI tests for factorforge slate command (Job 283A)."""

import json
from click.testing import CliRunner
from factorforge.cli.main import cli


def test_cli_slate_stdout_summary():
    runner = CliRunner()
    result = runner.invoke(cli, ["slate", "-s", "DIQMTQSPSSLSASVGDRVT", "--top-k", "3"])
    assert result.exit_code == 0
    assert "FactorForge Discovery Slate" in result.output
    assert "Top-K Candidates" in result.output
    assert "Diversity Index" in result.output


def test_cli_slate_json_output():
    runner = CliRunner()
    result = runner.invoke(
        cli, ["slate", "-s", "EVQLVESGGGLVQPGRSLRL", "--top-k", "2", "--json-output"]
    )
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data["$schema"] == "https://eijex.com/schemas/factorforge/candidate-slate-v0.1.json"
    assert data["slate_summary"]["top_k_count"] <= 2
    assert len(data["candidates"]) >= 1
