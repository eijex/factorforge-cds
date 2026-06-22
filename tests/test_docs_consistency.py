"""Docs-as-Code consistency guardrails.

Checks that key documentation files stay in sync with their source-of-truth
artifacts. Fails CI if claim wording, ablation layer definitions, or
reproducibility anchors drift.
"""
from __future__ import annotations
import hashlib
import json
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
BENCH_DIR = ROOT / "reproducibility" / "benchmark_v0.5.1"
ABLATION_SPEC = ROOT / "benchmarks" / "ablation" / "ablation_spec.yaml"
BENCHMARK_README = BENCH_DIR / "README.md"
MANIFEST_PATH = BENCH_DIR / "MANIFEST.json"
FROZEN_SUMMARY = BENCH_DIR / "data" / "benchmark_summary.frozen.json"


# ---------------------------------------------------------------------------
# Reproducibility anchor files
# ---------------------------------------------------------------------------

def test_manifest_json_exists_and_valid():
    assert MANIFEST_PATH.exists(), f"Missing: {MANIFEST_PATH}"
    data = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    for key in ("factorforge_version", "git_commit", "scoring_contract", "inputs", "outputs"):
        assert key in data, f"MANIFEST.json missing key: {key}"


def test_benchmark_summary_frozen_exists():
    assert FROZEN_SUMMARY.exists(), f"Missing: {FROZEN_SUMMARY}"


def test_manifest_scoring_contract_v1_1():
    data = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert data["scoring_contract"] == "v1.1", (
        f"Expected scoring_contract='v1.1', got {data['scoring_contract']!r}"
    )


def test_manifest_inputs_sha256_non_empty():
    data = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    for name, entry in data["inputs"].items():
        sha = entry.get("sha256", "")
        assert sha and len(sha) == 64, (
            f"MANIFEST.json inputs[{name!r}].sha256 is missing or invalid"
        )


def test_manifest_input_sha256_matches_files():
    # Hash the committed git blob, not local working-tree bytes: on Windows,
    # core.autocrlf=true plus this repo's .gitattributes (eol=lf for
    # json/yaml/yml) means the working tree can hold CRLF while git stores
    # LF, so path.read_bytes() can produce a value that mismatches what is
    # actually committed (the root cause of two prior MANIFEST drift incidents).
    data = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    for name, entry in data["inputs"].items():
        blob = subprocess.run(
            ["git", "show", f"HEAD:{entry['path']}"],
            cwd=ROOT,
            capture_output=True,
            check=True,
        ).stdout
        actual = hashlib.sha256(blob).hexdigest()
        assert actual == entry["sha256"], f"MANIFEST.json hash drift for {name}: {entry['path']}"


def test_manifest_commands_use_reproducible_entrypoints():
    data = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    commands = "\n".join(data["commands"])
    assert "python scripts/benchmark.py" not in commands
    assert "python -m benchmarks.run_benchmark" in commands
    assert "--dataset nbenthamiana_full --mode formal --seed 320" in commands
    assert "python -m benchmarks.ablation.run_ablation" in commands


def test_manifest_separates_software_and_benchmark_dois():
    data = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert data["archives"]["software_release"]["doi"] == "10.5281/zenodo.20640931"
    assert data["archives"]["corrected_benchmark_dataset"]["doi"] == "10.5281/zenodo.20676276"


# ---------------------------------------------------------------------------
# Ablation spec ↔ benchmark README consistency
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def ablation_spec():
    return yaml.safe_load(ABLATION_SPEC.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def benchmark_readme():
    return BENCHMARK_README.read_text(encoding="utf-8")


def test_ablation_spec_layers_in_benchmark_readme(ablation_spec, benchmark_readme):
    """Every layer key (L0–L5) from ablation_spec.yaml must appear in the benchmark README."""
    layers = ablation_spec.get("layers", {})
    assert layers, "ablation_spec.yaml has no 'layers' section"
    for level in layers:
        assert level in benchmark_readme, (
            f"ablation_spec layer '{level}' not documented in benchmark README"
        )


def test_ablation_spec_layer_names_in_benchmark_readme(ablation_spec, benchmark_readme):
    """Every layer 'name' value from ablation_spec.yaml must appear in the benchmark README."""
    layers = ablation_spec.get("layers", {})
    for level, cfg in layers.items():
        name = cfg.get("name", "")
        assert name and name in benchmark_readme, (
            f"Layer {level} name '{name}' not found in benchmark README"
        )


# ---------------------------------------------------------------------------
# Evidence boundary / claim wording
# ---------------------------------------------------------------------------

def test_main_readme_evidence_boundary():
    """Main README must contain in-silico or computational disclaimer."""
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "in-silico" in readme or "computational" in readme, (
        "Main README missing evidence boundary disclaimer (in-silico / computational)"
    )


def test_benchmark_readme_evidence_boundary(benchmark_readme):
    """Benchmark README must contain in-silico or computational disclaimer."""
    assert "in-silico" in benchmark_readme or "computational" in benchmark_readme, (
        "Benchmark README missing evidence boundary disclaimer"
    )


def test_benchmark_readme_no_wet_lab_claim(benchmark_readme):
    """Benchmark README must not contain unsupported wet-lab outcome claims."""
    forbidden = [
        "predicts expression",
        "predicts yield",
        "wet-lab success",
        "guarantees cloning",
        "guarantees synthesis",
        "predicts synthesis success",
    ]
    for phrase in forbidden:
        assert phrase.lower() not in benchmark_readme.lower(), (
            f"Benchmark README contains forbidden claim phrase: {phrase!r}"
        )
