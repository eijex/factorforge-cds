import csv
from pathlib import Path
from benchmarks.run_benchmark import run

ROOT = Path(__file__).resolve().parents[1]


def test_smoke_runs_and_writes_outputs(tmp_path):
    out_csv = tmp_path / "results.csv"
    out_md = tmp_path / "summary.md"
    run(dataset="synthetic", mode="smoke", out_csv=out_csv, out_md=out_md,
        proteins_fasta=ROOT / "tests" / "fixtures" / "small_proteins.fasta",
        native_fasta=ROOT / "tests" / "fixtures" / "small_native_cds.fasta")
    assert out_csv.exists() and out_md.exists()
    rows = list(csv.DictReader(out_csv.open(encoding="utf-8")))
    methods = {r["method"] for r in rows}
    # every method produced at least one row
    assert "random_synonymous" in methods
    assert "factorforge_balanced" in methods
    assert "multi_constraint_pass" in rows[0]
