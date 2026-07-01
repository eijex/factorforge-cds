"""CLI tests for checksum/tier-gated codon-reference expert mode."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import click
import pytest
from click.testing import CliRunner

import factorforge.cli.main as cli_main


SAMPLE_PROTEIN = "MSKGEELFTGVVPILVELD"


@pytest.fixture(scope="module")
def reference_entries() -> dict[str, dict]:
    return cli_main._reference_entries_by_id()


def _write_input(tmp_path: Path) -> Path:
    input_file = tmp_path / "input.fasta"
    input_file.write_text(f">test\n{SAMPLE_PROTEIN}\n", encoding="utf-8")
    return input_file


def test_reference_id_choice_values_match_manifest(reference_entries: dict[str, dict]) -> None:
    assert set(cli_main.REFERENCE_ID_CHOICES) == set(reference_entries)
    assert len(cli_main.REFERENCE_ID_CHOICES) == 5


def test_resolve_reference_by_id_loads_all_manifest_references(
    reference_entries: dict[str, dict],
) -> None:
    for reference_id, entry in reference_entries.items():
        path = cli_main.resolve_reference_by_id(reference_id)
        assert path.exists(), reference_id
        assert path == cli_main.FACTORFORGE_REPO_ROOT / entry["codon_table_path"]


def test_non_production_reference_ids_print_manifest_claim_boundary(
    reference_entries: dict[str, dict],
    capsys: pytest.CaptureFixture[str],
) -> None:
    non_production_entries = [
        entry for entry in reference_entries.values() if entry["tier"] != "production_enabled"
    ]
    assert len(non_production_entries) == 4

    observed_claim_boundaries = set()
    for entry in non_production_entries:
        capsys.readouterr()
        cli_main.resolve_reference_by_id(entry["reference_id"])
        err = capsys.readouterr().err
        assert f"reference_id={entry['reference_id']}" in err
        assert f"tier={entry['tier']}" in err
        assert entry["claim_boundary"] in err
        observed_claim_boundaries.add(entry["claim_boundary"])

    assert len(observed_claim_boundaries) == len(non_production_entries)


def test_packaged_install_fallback_resolves_bundled_manifest_and_tables(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    package_root = tmp_path / "factorforge"
    bundled_manifest = package_root / "data" / "reference" / "reference_policy_manifest.json"
    packaged_table = package_root / "data" / "profiles" / "codons.json"
    bundled_manifest.parent.mkdir(parents=True)
    packaged_table.parent.mkdir(parents=True)
    packaged_table.write_text('{"codons": {}}', encoding="utf-8")
    expected_sha = hashlib.sha256(packaged_table.read_bytes()).hexdigest()
    bundled_manifest.write_text(
        json.dumps(
            {
                "references": [
                    {
                        "reference_id": "packaged_reference",
                        "tier": "production_enabled",
                        "codon_table_path": "src/factorforge/data/profiles/codons.json",
                        "checksum_sha256": expected_sha,
                        "organism": "Nicotiana benthamiana",
                        "ncbi_taxid": 4100,
                        "claim_boundary": "Packaged test boundary.",
                        "known_limitations": [],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(cli_main, "FACTORFORGE_REPO_ROOT", tmp_path / "missing-repo-root")
    monkeypatch.setattr(
        cli_main,
        "REFERENCE_POLICY_MANIFEST_PATH",
        tmp_path / "missing-repo-root" / "data" / "reference" / "reference_policy_manifest.json",
    )
    monkeypatch.setattr(cli_main, "PACKAGE_ROOT", package_root)
    monkeypatch.setattr(cli_main, "BUNDLED_REFERENCE_POLICY_MANIFEST_PATH", bundled_manifest)

    assert cli_main._reference_policy_manifest_path() == bundled_manifest
    assert cli_main.resolve_reference_by_id("packaged_reference") == packaged_table


def test_resolve_reference_by_id_checksum_mismatch_hard_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    table_path = tmp_path / "codons.json"
    table_path.write_text('{"codons": {}}', encoding="utf-8")
    actual_sha = hashlib.sha256(table_path.read_bytes()).hexdigest()
    expected_sha = "0" * 64
    manifest_path = tmp_path / "reference_policy_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "references": [
                    {
                        "reference_id": "fake_reference",
                        "tier": "production_enabled",
                        "codon_table_path": "codons.json",
                        "checksum_sha256": expected_sha,
                        "organism": "Nicotiana benthamiana",
                        "ncbi_taxid": 4100,
                        "claim_boundary": "Test boundary.",
                        "known_limitations": [],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(cli_main, "FACTORFORGE_REPO_ROOT", tmp_path)
    monkeypatch.setattr(cli_main, "REFERENCE_POLICY_MANIFEST_PATH", manifest_path)

    with pytest.raises(click.UsageError) as exc_info:
        cli_main.resolve_reference_by_id("fake_reference")

    message = str(exc_info.value)
    assert "Checksum mismatch" in message
    assert expected_sha in message
    assert actual_sha in message


def test_cli_rejects_host_reference_mismatch(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = _write_input(tmp_path)

    result = runner.invoke(
        cli_main.cli,
        [
            "optimize",
            str(input_file),
            "--host",
            "by2",
            "--reference-id",
            "nbenthamiana_nbev11_hc_v2",
        ],
    )

    assert result.exit_code != 0
    assert "incompatible with --host ntabacum" in result.output
    assert "NCBI taxid 4100" in result.output
    assert "expected NCBI taxid 4097" in result.output


def test_cli_reference_id_omitted_keeps_default_path(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_if_called(reference_id: str) -> Path:
        raise AssertionError(f"unexpected reference resolution: {reference_id}")

    monkeypatch.setattr(cli_main, "resolve_reference_by_id", fail_if_called)
    runner = CliRunner()
    input_file = _write_input(tmp_path)

    result = runner.invoke(cli_main.cli, ["optimize", str(input_file)])

    assert result.exit_code == 0, result.output
    assert "Optimizing with DP feasibility engine" in result.output


def test_cli_accepts_production_reference_id_for_dp(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = _write_input(tmp_path)

    result = runner.invoke(
        cli_main.cli,
        [
            "optimize",
            str(input_file),
            "--reference-id",
            "nbenthamiana_nbev11_hc_v2",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Optimizing with DP feasibility engine" in result.output
    assert "Warning: reference_id=" not in result.output


def test_cli_non_production_reference_id_warns_and_runs(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = _write_input(tmp_path)

    result = runner.invoke(
        cli_main.cli,
        [
            "optimize",
            str(input_file),
            "--reference-id",
            "nbenthamiana_qld183_v103",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Warning: reference_id=nbenthamiana_qld183_v103" in result.output
    assert "Research comparator only" in result.output
    assert "Optimizing with DP feasibility engine" in result.output


def test_cli_reference_id_with_template_is_rejected(tmp_path: Path) -> None:
    runner = CliRunner()
    input_file = _write_input(tmp_path)

    result = runner.invoke(
        cli_main.cli,
        [
            "optimize",
            str(input_file),
            "--engine",
            "profile",
            "--reference-id",
            "nbenthamiana_nbev11_hc_v2",
            "--template",
            "standard_expression",
        ],
    )

    assert result.exit_code != 0
    assert "--reference-id is not supported with --template mode." in result.output
