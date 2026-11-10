"""Unit and CLI tests for evidence-aware structure-conditioned sLLM adapter."""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from factorforge.cli.main import cli
from factorforge.structure.esmfold_client import predict_and_ingest_esm
from factorforge.structure.sllm_adapter import (
    SLLMAdapterError,
    generate_evidence_explanation,
)

SAMPLE_PROTEIN = "M" + "KIF" * 12 + "LMM"


@pytest.fixture
def synthetic_assets(tmp_path: Path):
    fasta_file = tmp_path / "protein.fasta"
    fasta_file.write_text(f">sample_protein\n{SAMPLE_PROTEIN}\n", encoding="utf-8")
    sidecar_file = tmp_path / "sidecar.json"
    predict_and_ingest_esm(
        sequence=SAMPLE_PROTEIN,
        provider="mock",
        output_sidecar_path=sidecar_file,
        tunnel_offset=35,
    )
    return fasta_file, sidecar_file


def test_mock_sidecar_blocked_without_allow_mock(synthetic_assets):
    _, sidecar = synthetic_assets
    with pytest.raises(SLLMAdapterError, match="synthetic/mock data"):
        generate_evidence_explanation(sidecar_source=sidecar, allow_mock=False)


def test_evidence_explanation_tiers(synthetic_assets, tmp_path: Path):
    fasta, sidecar = synthetic_assets

    # Generate an optimized DNA file using 4DP
    optimized_dna_file = tmp_path / "opt.fasta"
    runner = CliRunner()
    opt_res = runner.invoke(
        cli,
        [
            "optimize",
            str(fasta),
            "--engine",
            "dp",
            "--pause-sidecar",
            str(sidecar),
            "--allow-mock",
            "--gc-min",
            "0",
            "--gc-max",
            "100",
            "-o",
            str(optimized_dna_file),
        ],
    )
    assert opt_res.exit_code == 0, opt_res.output

    # Generate evidence-aware report
    report = generate_evidence_explanation(
        sidecar_source=sidecar,
        optimized_fasta_content=optimized_dna_file,
        allow_mock=True,
        target_name="Test_Protein",
    )

    # 1. Verified invariants
    assert any("100% preserved" in f for f in report.verified_facts)
    assert any("GC content" in f for f in report.verified_facts)
    assert any("CAI" in f for f in report.verified_facts)
    assert any("0 detected" in f for f in report.verified_facts)

    # 2. Predicted metrics
    assert any("pLDDT" in p for p in report.predicted_metrics)

    # 3. Hypothesized effects (explicitly marked as hypothesis)
    assert any("HYPOTHESIS" in h for h in report.hypothesized_effects)
    assert any("EPISTEMIC LIMITATION" in h for h in report.hypothesized_effects)

    # 4. Actionable wet-lab plan
    assert any("Western blot" in w for w in report.wet_lab_validation_plan)
    assert any("Ribo-seq" in w for w in report.wet_lab_validation_plan)

    # Invariants
    assert report.reviewed is False
    assert report.biological_validation == "NOT_AVAILABLE"

    # Markdown rendering
    md = report.to_markdown()
    assert "## 1. Verified Invariants" in md
    assert "## 2. Predicted Structural Metrics" in md
    assert "## 3. Hypothesized Biological Mechanisms" in md
    assert "## 4. Wet-Lab Experimental Validation Plan" in md


def test_mismatched_translation_rejected(synthetic_assets, tmp_path: Path):
    _, sidecar = synthetic_assets
    wrong_dna_file = tmp_path / "wrong.fasta"
    # DNA translating to 'AAAA'
    wrong_dna_file.write_text(">wrong\nGCCGCCGCCGCC\n", encoding="utf-8")

    with pytest.raises(SLLMAdapterError, match="Translation mismatch"):
        generate_evidence_explanation(
            sidecar_source=sidecar,
            optimized_fasta_content=wrong_dna_file,
            allow_mock=True,
        )


def test_cli_explain_command(synthetic_assets, tmp_path: Path):
    _, sidecar = synthetic_assets
    report_out = tmp_path / "assessment.md"
    runner = CliRunner()

    res = runner.invoke(
        cli,
        [
            "structure",
            "explain",
            "-s",
            str(sidecar),
            "--allow-mock",
            "-o",
            str(report_out),
        ],
    )
    assert res.exit_code == 0, res.output
    assert report_out.is_file()
    content = report_out.read_text(encoding="utf-8")
    assert "Evidence-Aware Cotranslational Optimization Assessment" in content
    assert "## 1. Verified Invariants" in content
    assert "## 4. Wet-Lab Experimental Validation Plan" in content


def test_cli_explain_json_output(synthetic_assets):
    _, sidecar = synthetic_assets
    runner = CliRunner()
    res = runner.invoke(
        cli,
        [
            "structure",
            "explain",
            "-s",
            str(sidecar),
            "--allow-mock",
            "--json",
        ],
    )
    assert res.exit_code == 0, res.output
    data = json.loads(res.output)
    assert "verified_facts" in data
    assert "hypothesized_effects" in data
    assert data["reviewed"] is False
