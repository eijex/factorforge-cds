"""End-to-end and CLI tests for `factorforge structure predict-esm`."""

import json
from pathlib import Path

from Bio.Seq import Seq
from click.testing import CliRunner

from factorforge.cli.main import cli

SAMPLE_PROTEIN = "M" + "KIF" * 12 + "LMM"


def test_cli_predict_esm_with_sequence(tmp_path: Path):
    runner = CliRunner()
    sidecar_out = tmp_path / "sidecar.json"
    pdb_out = tmp_path / "model.pdb"

    res = runner.invoke(
        cli,
        [
            "structure",
            "predict-esm",
            SAMPLE_PROTEIN,
            "--provider",
            "mock",
            "-o",
            str(sidecar_out),
            "--output-pdb",
            str(pdb_out),
            "--pause-window",
            "35",
        ],
    )

    assert res.exit_code == 0, res.output
    assert "ESMFold Structure Report" in res.output
    assert f"Protein Length: {len(SAMPLE_PROTEIN)} aa" in res.output
    assert sidecar_out.is_file()
    assert pdb_out.is_file()

    data = json.loads(sidecar_out.read_text(encoding="utf-8"))
    assert data["status"] == "AVAILABLE"
    assert data["source_type"] == "SYNTHETIC_FIXTURE"
    assert data["coordinate_system"] == "PDB_res_seq_1based"
    assert data["chains"]["A"]["sequence"] == SAMPLE_PROTEIN


def test_cli_predict_esm_json_output(tmp_path: Path):
    runner = CliRunner()
    res = runner.invoke(
        cli,
        [
            "structure",
            "predict-esm",
            SAMPLE_PROTEIN,
            "--provider",
            "mock",
            "--json",
        ],
    )

    assert res.exit_code == 0, res.output
    data = json.loads(res.output)
    assert data["status"] == "AVAILABLE"
    assert data["source_type"] == "SYNTHETIC_FIXTURE"


def test_cli_predict_esm_fasta_file(tmp_path: Path):
    runner = CliRunner()
    fasta_file = tmp_path / "target.fasta"
    fasta_file.write_text(f">target_protein\n{SAMPLE_PROTEIN}\n", encoding="utf-8")
    sidecar_out = tmp_path / "sidecar.json"

    res = runner.invoke(
        cli,
        [
            "structure",
            "predict-esm",
            "-f",
            str(fasta_file),
            "--provider",
            "mock",
            "-o",
            str(sidecar_out),
        ],
    )

    assert res.exit_code == 0, res.output
    assert sidecar_out.is_file()


def test_cli_predict_esm_e2e_piping_to_optimize_4dp(tmp_path: Path):
    """Verify that predict-esm sidecar seamlessly drives `factorforge optimize --engine dp --pause-sidecar`."""
    runner = CliRunner()

    fasta_file = tmp_path / "input.fasta"
    fasta_file.write_text(f">input\n{SAMPLE_PROTEIN}\n", encoding="utf-8")

    sidecar_file = tmp_path / "esm_sidecar.json"
    optimized_fasta = tmp_path / "optimized.fasta"

    # Step 1: Predict structure and emit a synthetic sidecar
    predict_res = runner.invoke(
        cli,
        [
            "structure",
            "predict-esm",
            "-f",
            str(fasta_file),
            "--provider",
            "mock",
            "-o",
            str(sidecar_file),
            "--pause-window",
            "35",
        ],
    )
    assert predict_res.exit_code == 0, predict_res.output
    assert sidecar_file.is_file()

    # Step 2a: Verify mock sidecar is BLOCKED by default without --allow-mock
    blocked_res = runner.invoke(
        cli,
        [
            "optimize",
            str(fasta_file),
            "--engine",
            "dp",
            "--pause-sidecar",
            str(sidecar_file),
            "-o",
            str(optimized_fasta),
        ],
    )
    assert blocked_res.exit_code != 0
    assert "synthetic/mock data" in blocked_res.output
    assert not optimized_fasta.exists()

    # Step 2b: Pass --allow-mock to authorize fixture for testing
    opt_res = runner.invoke(
        cli,
        [
            "optimize",
            str(fasta_file),
            "--engine",
            "dp",
            "--pause-sidecar",
            str(sidecar_file),
            "--allow-mock",
            "--gc-min",
            "0",
            "--gc-max",
            "100",
            "-o",
            str(optimized_fasta),
        ],
    )
    assert opt_res.exit_code == 0, opt_res.output
    assert "sequence-matched model offset positions: [38]" in opt_res.output
    assert optimized_fasta.is_file()

    header, *dna_lines = optimized_fasta.read_text(encoding="utf-8").splitlines()
    assert "engine=dp_4dp" in header
    assert "pauses=38" in header

    translated = str(Seq("".join(dna_lines)).translate()).removesuffix("*")
    assert translated == SAMPLE_PROTEIN


def test_remote_request_requires_explicit_opt_in(monkeypatch, tmp_path):
    def unexpected_call(**kwargs):
        raise AssertionError("Remote path must not be called without consent")

    monkeypatch.setattr("factorforge.cli.structure_cli.predict_and_ingest_esm", unexpected_call)
    output = tmp_path / "sidecar.json"
    result = CliRunner().invoke(cli, ["structure", "predict-esm", SAMPLE_PROTEIN, "-o", str(output)])
    assert result.exit_code == 2
    assert "--allow-remote" in result.output
    assert not output.exists()


def test_cli_rejects_multiple_fasta_records(tmp_path):
    path = tmp_path / "input.fasta"
    path.write_text(">one\nMKIF\n>two\nMLLF\n", encoding="utf-8")
    result = CliRunner().invoke(cli, ["structure", "predict-esm", str(path), "--provider", "mock"])
    assert result.exit_code == 2
    assert "exactly one protein" in result.output
