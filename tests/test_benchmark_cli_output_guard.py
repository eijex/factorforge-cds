"""CLI output-directory safety tests for the benchmark runner."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest

import benchmarks.run_benchmark as run_benchmark


def _invoke_main(monkeypatch: pytest.MonkeyPatch, args: list[str], fake_run) -> None:
    monkeypatch.setattr(sys, "argv", ["run_benchmark.py", *args])
    monkeypatch.setattr(run_benchmark, "run", fake_run)
    run_benchmark.main()


def test_cli_default_output_dir_stays_v320(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict[str, Path] = {}

    def fake_run(**kwargs) -> None:
        captured["out_csv"] = kwargs["out_csv"]
        captured["out_md"] = kwargs["out_md"]

    _invoke_main(monkeypatch, ["--force"], fake_run)

    assert captured["out_csv"] == run_benchmark.DEFAULT_RESULTS_DIR.resolve() / "benchmark_results.csv"
    assert captured["out_md"] == run_benchmark.DEFAULT_RESULTS_DIR.resolve() / "benchmark_summary.md"


def test_cli_refuses_existing_results_without_force(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "existing"
    output_dir.mkdir()
    existing_csv = output_dir / "benchmark_results.csv"
    existing_json = output_dir / "benchmark_summary.json"
    existing_csv.write_text("do-not-overwrite\n", encoding="utf-8")
    existing_json.write_text('{"status": "do-not-overwrite"}\n', encoding="utf-8")
    before_csv_hash = hashlib.sha256(existing_csv.read_bytes()).hexdigest()
    before_json_hash = hashlib.sha256(existing_json.read_bytes()).hexdigest()
    before_csv_mtime = existing_csv.stat().st_mtime_ns
    before_json_mtime = existing_json.stat().st_mtime_ns

    def fail_run(**_kwargs) -> None:
        raise AssertionError("run() must not be called when overwrite guard fails")

    with pytest.raises(SystemExit) as exc_info:
        _invoke_main(monkeypatch, ["--out-dir", str(output_dir)], fail_run)

    message = str(exc_info.value)
    assert "Pass --force" in message
    assert "benchmark_results.csv" in message
    assert "benchmark_summary.json" in message
    assert hashlib.sha256(existing_csv.read_bytes()).hexdigest() == before_csv_hash
    assert hashlib.sha256(existing_json.read_bytes()).hexdigest() == before_json_hash
    assert existing_csv.stat().st_mtime_ns == before_csv_mtime
    assert existing_json.stat().st_mtime_ns == before_json_mtime


def test_cli_force_allows_existing_results_overwrite(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "existing"
    output_dir.mkdir()
    (output_dir / "benchmark_results.csv").write_text("old\n", encoding="utf-8")
    (output_dir / "benchmark_summary.json").write_text('{"status": "old"}\n', encoding="utf-8")

    def fake_run(**kwargs) -> None:
        kwargs["out_csv"].write_text("new\n", encoding="utf-8")
        kwargs["out_md"].write_text("# new\n", encoding="utf-8")
        (kwargs["out_md"].parent / "benchmark_summary.json").write_text(
            '{"status": "new"}\n',
            encoding="utf-8",
        )

    _invoke_main(monkeypatch, ["--out-dir", str(output_dir), "--force"], fake_run)

    assert (output_dir / "benchmark_results.csv").read_text(encoding="utf-8") == "new\n"
    assert (output_dir / "benchmark_summary.md").read_text(encoding="utf-8") == "# new\n"
    assert (output_dir / "benchmark_summary.json").read_text(encoding="utf-8") == (
        '{"status": "new"}\n'
    )
