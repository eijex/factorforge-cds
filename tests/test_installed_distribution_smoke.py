"""Protect distribution smoke against accidentally testing the checkout."""

import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "smoke_installed_distribution.py"
SPEC = importlib.util.spec_from_file_location("distribution_smoke", SCRIPT)
smoke = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(smoke)


@pytest.fixture
def environment(tmp_path):
    module = tmp_path / "lib" / "site-packages" / "factorforge" / "__init__.py"
    module.parent.mkdir(parents=True)
    module.touch()
    metadata = tmp_path / "conda-meta" / "factorforge-cds-3.5.11-test.json"
    metadata.parent.mkdir()
    record = {
        "name": "factorforge-cds",
        "version": "3.5.11",
        "channel": "https://conda.anaconda.org/bioconda/noarch",
        "files": [module.relative_to(tmp_path).as_posix()],
    }
    metadata.write_text(json.dumps(record), encoding="utf-8")
    return tmp_path, module, metadata, record


def test_reject_checkout_outside_environment(environment):
    prefix, _, _, _ = environment
    with pytest.raises(ValueError, match="outside"):
        smoke.verify_origin(SCRIPT, prefix, "3.5.11")


@pytest.mark.parametrize(
    "change",
    [
        {"version": "0.0.0"},
        {"channel": "https://evil.example/bioconda"},
        {"channel": "https://conda.anaconda.org/bioconda-impostor"},
        {"files": []},
    ],
)
def test_reject_wrong_artifact(environment, change):
    prefix, module, metadata, record = environment
    record.update(change)
    metadata.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(ValueError):
        smoke.verify_origin(module, prefix, "3.5.11")


def test_reject_missing_conda_record(environment):
    prefix, module, metadata, _ = environment
    metadata.unlink()
    with pytest.raises(ValueError, match="exactly one"):
        smoke.verify_origin(module, prefix, "3.5.11")


def test_failure_is_retained_for_ci(monkeypatch, tmp_path):
    def fail(report):
        report["calls"].append({"exit_code": 1})
        raise ValueError("synthetic launcher failure")

    monkeypatch.setattr(smoke, "exercise", fail)
    assert smoke.main(["--output", str(tmp_path)]) == 1
    report = json.loads((tmp_path / "smoke-result.json").read_text())
    assert report["status"] == "FAIL"
    assert "synthetic launcher failure" in report["error"]
    assert report["calls"] == [{"exit_code": 1}]
