"""Benchmark summary codon table metadata guardrails."""
from __future__ import annotations
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FORMAL_SUMMARY_JSON = ROOT / "benchmarks" / "results" / "v3.2.0" / "benchmark_summary.json"
FORMAL_SUMMARY_MD = ROOT / "benchmarks" / "results" / "v3.2.0" / "benchmark_summary.md"
DATASETS_DIR = ROOT / "benchmarks" / "datasets"
GENOME_EXTENSIONS = {".fasta", ".fa", ".fna", ".gz", ".fastq"}


# --- formal benchmark_summary.json ---

def test_formal_benchmark_summary_json_exists():
    assert FORMAL_SUMMARY_JSON.exists(), f"Missing: {FORMAL_SUMMARY_JSON}"


@pytest.fixture(scope="module")
def formal_summary():
    return json.loads(FORMAL_SUMMARY_JSON.read_text(encoding="utf-8"))


def test_formal_summary_codon_table_fields_present(formal_summary):
    required = [
        "codon_table_id",
        "codon_table_sha256",
        "codon_table_source_status",
        "codon_table_build_path_status",
        "codon_table_reference_note",
    ]
    for field in required:
        assert field in formal_summary, f"Missing codon table field in benchmark_summary.json: {field}"


def test_formal_summary_codon_table_id(formal_summary):
    assert formal_summary["codon_table_id"] == "nbenthamiana_legacy_kazusa_sgn_v101"


def test_formal_summary_codon_table_sha256_present(formal_summary):
    sha = formal_summary.get("codon_table_sha256")
    assert sha and len(sha) == 64, (
        f"codon_table_sha256 must be a 64-char hex string, got: {sha!r}"
    )


def test_formal_summary_codon_table_source_status(formal_summary):
    assert formal_summary["codon_table_source_status"] == "legacy_metadata_only"


def test_formal_summary_codon_table_build_path_status(formal_summary):
    assert formal_summary["codon_table_build_path_status"] == "incomplete"


def test_formal_summary_is_formal_dataset(formal_summary):
    """Ensure committed summary is the formal run, not a synthetic/smoke run."""
    assert formal_summary.get("dataset_id") == "formal", (
        f"benchmark_summary.json must be from the formal dataset, got: {formal_summary.get('dataset_id')!r}"
    )


# --- formal benchmark_summary.md ---

def test_formal_benchmark_summary_md_exists():
    assert FORMAL_SUMMARY_MD.exists(), f"Missing: {FORMAL_SUMMARY_MD}"


def test_formal_summary_md_contains_codon_table_note():
    content = FORMAL_SUMMARY_MD.read_text(encoding="utf-8")
    assert "legacy_metadata_only" in content or "codon table" in content.lower(), (
        "benchmark_summary.md must contain codon table provenance note"
    )


def test_formal_summary_md_contains_required_wording():
    content = FORMAL_SUMMARY_MD.read_text(encoding="utf-8")
    assert "legacy FactorForge reference" in content, (
        "benchmark_summary.md must include required provenance wording"
    )


# --- smoke run: codon table fields in generated JSON ---

def test_smoke_summary_contains_codon_table_fields(tmp_path):
    """Smoke run must produce benchmark_summary.json with codon_table_* fields."""
    from benchmarks.run_benchmark import run

    out_csv = tmp_path / "results.csv"
    out_md = tmp_path / "summary.md"
    run(
        dataset="synthetic",
        mode="smoke",
        out_csv=out_csv,
        out_md=out_md,
        proteins_fasta=ROOT / "tests" / "fixtures" / "small_proteins.fasta",
        native_fasta=ROOT / "tests" / "fixtures" / "small_native_cds.fasta",
    )
    summary_json = tmp_path / "benchmark_summary.json"
    assert summary_json.exists(), "smoke run must produce benchmark_summary.json"

    data = json.loads(summary_json.read_text(encoding="utf-8"))
    for field in [
        "codon_table_id",
        "codon_table_sha256",
        "codon_table_source_status",
        "codon_table_build_path_status",
        "codon_table_reference_note",
    ]:
        assert field in data, f"Missing codon table field in smoke summary JSON: {field}"

    assert data["codon_table_id"] == "nbenthamiana_legacy_kazusa_sgn_v101"


# --- raw FASTA / genome files must not be committed ---

def test_no_raw_genome_fasta_committed():
    """Raw genome/CDS/protein FASTA files must not be committed to the repo."""
    forbidden_patterns = ["*.fasta", "*.fa", "*.fna", "*.fastq"]
    allowed_fixtures = ROOT / "tests" / "fixtures"

    import subprocess
    result = subprocess.run(
        ["git", "ls-files", "--", "*.fasta", "*.fa", "*.fna", "*.fastq"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    allowed_prefixes = (
        "tests/fixtures/",   # unit test fixtures
        "archive/",          # historical archived files
        "examples/",         # public example files
    )
    committed = [
        p for p in result.stdout.strip().splitlines()
        if p and not any(p.startswith(pfx) for pfx in allowed_prefixes)
    ]
    assert not committed, (
        f"Raw FASTA/genome files must not be committed outside allowed locations:\n"
        + "\n".join(committed)
    )
