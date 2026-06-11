"""Focused regression locks for the v3.2.0 benchmark contract."""

import csv
import json
from pathlib import Path

import pytest
import yaml

from benchmarks.config import load_benchmark_config
from benchmarks.run_benchmark import run
from factorforge.engines.profile.rules.domesticator import Domesticator
from factorforge.registry.registry_loader import load_registry, resolve_ref

ROOT = Path(__file__).resolve().parents[1]
BASELINE_METHODS = {
    "random_synonymous",
    "greedy_cai",
    "native_reference",
}
RESULT_REQUIRED_FIELDS = {
    "method",
    "method_type",
    "sequence_id",
    "aa_identity",
    "internal_stop_count",
    "invalid_codon_count",
    "length_multiple_of_three",
    "cai",
    "gc_percent",
    "gc_in_target_range",
    "forbidden_type_iis_site_count",
    "biological_pass",
    "assembly_pass",
    "multi_constraint_pass",
    "runtime_seconds",
    "replicate",
    "seed",
}
SUMMARY_REQUIRED_FIELDS = {
    "run_id",
    "timestamp",
    "factorforge_version",
    "registry_version",
    "registry_hash",
    "spec_hash",
    "runtime_seconds",
    "dataset_id",
    "dataset_n",
    "methods",  # per-method breakdown; pass_rate_* are now nested under methods.<name>
}
RAW_SEQUENCE_FIELDS = {"sequence", "raw_sequence", "protein", "cds", "output_cds", "input_sequence"}


@pytest.fixture(scope="module")
def benchmark_outputs(tmp_path_factory: pytest.TempPathFactory) -> tuple[list[dict], dict, str]:
    output_dir = tmp_path_factory.mktemp("benchmark-regression")
    out_csv = output_dir / "benchmark_results.csv"
    out_md = output_dir / "benchmark_summary.md"
    run(
        dataset="synthetic",
        mode="regression",
        out_csv=out_csv,
        out_md=out_md,
        proteins_fasta=ROOT / "tests/fixtures/small_proteins.fasta",
        native_fasta=ROOT / "tests/fixtures/small_native_cds.fasta",
    )
    rows = list(csv.DictReader(out_csv.open(encoding="utf-8")))
    summary = json.loads((output_dir / "benchmark_summary.json").read_text(encoding="utf-8"))
    return rows, summary, out_md.read_text(encoding="utf-8")


def test_active_enzyme_set_matches_registry_production_and_canonical_set() -> None:
    registry_type_iis = resolve_ref(
        load_registry(), "parameters.constraints.assembly.type_iis_enzymes.value"
    )
    config_type_iis = load_benchmark_config().forbidden_type_iis
    canonical = {"BsaI", "BpiI", "BsmBI"}
    assert set(registry_type_iis) == set(config_type_iis) == set(Domesticator.GOLDEN_GATE_ENZYMES)
    assert set(registry_type_iis) == canonical
    assert "BbsI" not in registry_type_iis


def test_baseline_names_and_random_seeds_are_stable(benchmark_outputs) -> None:
    rows, _, _ = benchmark_outputs
    assert BASELINE_METHODS.issubset({row["method"] for row in rows})
    random_seeds = {int(row["seed"]) for row in rows if row["method"] == "random_synonymous"}
    assert random_seeds == {320, 321, 322}


def test_result_schema_has_no_raw_sequence_columns(benchmark_outputs) -> None:
    rows, _, _ = benchmark_outputs
    assert RESULT_REQUIRED_FIELDS.issubset(rows[0])
    assert RAW_SEQUENCE_FIELDS.isdisjoint(rows[0])


def test_summary_schema_contains_registry_provenance(benchmark_outputs) -> None:
    _, summary, _ = benchmark_outputs
    assert SUMMARY_REQUIRED_FIELDS.issubset(summary)
    assert summary["registry_version"]
    assert summary["registry_hash"].startswith("sha256:")


def test_benchmark_contract_does_not_make_performance_claims(benchmark_outputs) -> None:
    _, _, summary_md = benchmark_outputs
    prohibited_claims = (
        "guarantees expression",
        "increases yield",
        "wet-lab proven",
        "clinical performance improved",
    )
    lowered = summary_md.lower()
    assert all(claim not in lowered for claim in prohibited_claims)


def test_benchmark_spec_is_in_silico_only_and_output_contract_is_stable() -> None:
    spec = yaml.safe_load((ROOT / "benchmarks/benchmark_spec.yaml").read_text(encoding="utf-8"))
    assert spec["benchmark"]["status"] == "in_silico_only"
    assert spec["random"] == {"seed": 320, "replicates": 3}
    assert set(spec["methods"]["baselines"]) == BASELINE_METHODS
