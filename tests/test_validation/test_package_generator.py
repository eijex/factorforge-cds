"""ValidationPackageGenerator tests."""

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SRC_PATH = str(_REPO_ROOT / "src")

# Derive factorforge source path from the currently-imported module so that
# subprocess calls work regardless of which Python executable is used.
import factorforge as _ff
_FF_SRC_PATH = str(Path(_ff.__file__).resolve().parent.parent)

from factorforge.validation import ValidationPackageGenerator
from factorforge.validation.package_generator import WetLabResult


@pytest.fixture
def sample_result():
    seq = "ATGAAACCC"
    return WetLabResult(
        construct_id="CF-20260525-000000",
        factorforge_version="3.2.4",
        host_profile="nbenthamiana",
        profile="balanced",
        sequence_hash="sha256:" + hashlib.sha256(seq.encode()).hexdigest(),
        protein_name="TestProtein",
        host_organism="N. benthamiana",
        promoter="35S",
        subcellular_targeting="Cytosol",
        expression_system="Agroinfiltration",
        harvest_timepoint="5 dpi",
        native_control="Yes",
        comparison="Improved",
        expression_level="Detected",
        notes="Test note",
        institution="Test Lab",
        public_listing=True,
    )


def test_generate_creates_three_files(sample_result, tmp_path):
    gen = ValidationPackageGenerator(tmp_path)
    output = gen.generate(sample_result)

    assert (output / "validation_metadata.json").exists()
    assert (output / "validation_summary.txt").exists()
    assert (output / "issue_body.md").exists()


def test_metadata_no_raw_sequence(sample_result, tmp_path):
    gen = ValidationPackageGenerator(tmp_path)
    output = gen.generate(sample_result)
    meta = json.loads((output / "validation_metadata.json").read_text(encoding="utf-8"))

    assert "sequence" not in meta
    assert meta["sequence_hash"].startswith("sha256:")


def test_metadata_construct_id(sample_result, tmp_path):
    gen = ValidationPackageGenerator(tmp_path)
    output = gen.generate(sample_result)
    meta = json.loads((output / "validation_metadata.json").read_text(encoding="utf-8"))

    assert meta["construct_id"] == "CF-20260525-000000"


def test_issue_body_contains_construct_id(sample_result, tmp_path):
    gen = ValidationPackageGenerator(tmp_path)
    output = gen.generate(sample_result)
    body = (output / "issue_body.md").read_text(encoding="utf-8")

    assert "CF-20260525-000000" in body


def test_summary_readable(sample_result, tmp_path):
    gen = ValidationPackageGenerator(tmp_path)
    output = gen.generate(sample_result)
    summary = (output / "validation_summary.txt").read_text(encoding="utf-8")

    assert "FactorForge Validation Report" in summary
    assert "CF-20260525-000000" in summary


def test_cli_generates_package_without_raw_sequence(tmp_path):
    output_dir = tmp_path / "validation_package"

    env = os.environ.copy()
    existing_pp = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = _FF_SRC_PATH + (os.pathsep + existing_pp if existing_pp else "")

    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "factorforge.validation.cli",
            "--construct-id",
            "CF-20260525-000001",
            "--version",
            "3.2.4",
            "--profile",
            "balanced",
            "--sequence",
            "ATGAAACCC",
            "--protein-name",
            "TestProtein",
            "--expression-system",
            "Agroinfiltration",
            "--comparison",
            "Improved",
            "--expression-level",
            "Detected",
            "--output-dir",
            str(output_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )

    combined_output = "\n".join(
        path.read_text(encoding="utf-8") for path in output_dir.iterdir() if path.is_file()
    )
    meta = json.loads((output_dir / "validation_metadata.json").read_text(encoding="utf-8"))

    assert "Validation package generated" in completed.stdout
    assert "ATGAAACCC" not in combined_output
    assert meta["sequence_hash"].startswith("sha256:")
