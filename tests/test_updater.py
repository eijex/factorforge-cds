"""Unit tests for FactorForge explicit read-only updater check."""

from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from factorforge.cli.main import cli
from factorforge.cli.updater import (
    check_factorforge_update,
    format_update_report,
    get_installed_distribution_info,
    query_pypi_latest_version,
)


def test_get_installed_distribution_info():
    """Verify package manager provenance detection for active environment."""
    info = get_installed_distribution_info()
    assert "name" in info
    assert "installed_version" in info
    assert "installer" in info
    assert isinstance(info["is_editable"], bool)
    assert isinstance(info["is_externally_managed"], bool)


def test_query_pypi_latest_version_filtering():
    """Verify prerelease, yanked, and incompatible python version filtering."""
    mock_payload = {
        "info": {"version": "3.6.0"},
        "releases": {
            "3.4.0": [{"yanked": False, "requires_python": ">=3.9"}],
            "3.5.0": [{"yanked": False, "requires_python": ">=3.9"}],
            "3.6.0b1": [{"yanked": False, "requires_python": ">=3.9"}],  # Prerelease -> skip
            "3.6.0rc1": [{"yanked": False, "requires_python": ">=3.9"}],  # Prerelease -> skip
            "3.6.0": [{"yanked": False, "requires_python": ">=3.9"}],    # Highest stable compatible
            "3.7.0": [{"yanked": True, "requires_python": ">=3.9"}],     # Yanked -> skip
            "3.8.0": [{"yanked": False, "requires_python": ">=3.99"}],   # Incompatible Python -> skip
        },
    }

    mock_resp = MagicMock()
    for artifacts in mock_payload["releases"].values():
        for artifact in artifacts:
            artifact["packagetype"] = "sdist"
    mock_payload["releases"]["9.0.0"] = [{
        "packagetype": "bdist_wheel", "filename": "factorforge_cds-9.0.0-cp311-cp311-madeup_arch.whl",
    }]
    mock_resp.status = 200
    mock_resp.read.return_value = json.dumps(mock_payload).encode("utf-8")
    mock_resp.__enter__.return_value = mock_resp
    mock_resp.__exit__.return_value = None

    with patch("urllib.request.urlopen", return_value=mock_resp):
        res = query_pypi_latest_version("factorforge-cds", target_python="3.11.0")

    assert res["available"] is True
    assert res["latest_version"] == "3.6.0"
    assert res["error"] is None


def test_query_pypi_network_error_fail_safe():
    """Verify fail-open contract on timeout or connection error."""
    with patch("urllib.request.urlopen", side_effect=TimeoutError("Connection timed out")):
        res = query_pypi_latest_version("factorforge-cds", timeout=0.1)

    assert res["available"] is False
    assert res["latest_version"] is None
    assert "lookup failed" in str(res["error"])


def test_update_unverified_install_never_queries_public_index():
    with (
        patch("factorforge.cli.updater.get_installed_distribution_info", return_value={
            "name": "factorforge-cds", "installed_version": "3.5.12",
            "installer": "pip", "direct_url": {}, "is_editable": False,
            "is_externally_managed": False,
        }),
        patch("factorforge.cli.updater.query_pypi_latest_version") as lookup,
    ):
        assert check_factorforge_update()["status"] == "CHECK_UNAVAILABLE"
        lookup.assert_not_called()


def test_check_factorforge_update_editable_handling():
    """Verify that development/editable installs are routed without PyPI query."""
    with (
        patch(
            "factorforge.cli.updater.get_installed_distribution_info",
            return_value={
                "name": "factorforge-cds",
                "installed_version": "3.5.12",
                "installer": "editable",
                "is_editable": True,
                "direct_url": {"dir_info": {"editable": True}},
                "is_externally_managed": False,
                "sys_prefix": "C:/dummy",
                "python_version": "3.11.0",
            },
        ),
        patch("factorforge.cli.updater.query_pypi_latest_version") as mock_pypi,
    ):
        res = check_factorforge_update()
        # PyPI lookup must NOT be called for editable install
        mock_pypi.assert_not_called()

    assert res["status"] == "DEVELOPMENT_MODE"
    assert res["upgrade_command"] is None


def test_check_factorforge_update_conda_handling():
    """Verify that Conda installations are identified and isolated from PyPI."""
    with (
        patch(
            "factorforge.cli.updater.get_installed_distribution_info",
            return_value={
                "name": "factorforge-cds",
                "installed_version": "3.5.0",
                "installer": "conda",
                "is_editable": False,
                "direct_url": None,
                "is_externally_managed": False,
                "sys_prefix": "C:/dummy/conda",
                "python_version": "3.11.0",
            },
        ),
        patch("factorforge.cli.updater.query_pypi_latest_version") as mock_pypi,
    ):
        res = check_factorforge_update()
        mock_pypi.assert_not_called()

    assert res["status"] == "CHECK_UNAVAILABLE"
    assert res["upgrade_command"] is None


def test_check_factorforge_update_pip_states():
    """Verify UPDATE_AVAILABLE and UP_TO_DATE comparisons under pip."""
    base_info = {
        "name": "factorforge-cds",
        "installed_version": "3.5.0",
        "installer": "pip",
        "is_editable": False,
        "direct_url": None,
        "is_externally_managed": False,
        "sys_prefix": "C:/dummy/venv",
        "python_version": "3.11.0",
    }

    # Case 1: Newer version on PyPI -> UPDATE_AVAILABLE
    with patch("factorforge.cli.updater.get_installed_distribution_info", return_value=base_info):
        mock_pypi = {"latest_version": "3.6.0", "available": True, "error": None}
        res_update = check_factorforge_update(mock_pypi_result=mock_pypi)

    assert res_update["status"] == "UPDATE_AVAILABLE"
    assert res_update["upgrade_command"] is None
    assert res_update["latest_version"] == "3.6.0"

    # Case 2: Matching/newer installed version -> UP_TO_DATE
    base_info_current = dict(base_info, installed_version="3.6.0")
    with patch("factorforge.cli.updater.get_installed_distribution_info", return_value=base_info_current):
        mock_pypi = {"latest_version": "3.6.0", "available": True, "error": None}
        res_current = check_factorforge_update(mock_pypi_result=mock_pypi)

    assert res_current["status"] == "UP_TO_DATE"


def test_cli_self_update_text_and_json():
    """Verify CLI command output in human-readable and JSON format."""
    runner = CliRunner()

    # Plain text run
    res_text = runner.invoke(cli, ["self-update", "--check"])
    assert res_text.exit_code == 0
    assert "FactorForge Distribution:" in res_text.output
    assert "Epistemic: UNVERIFIED_PUBLIC_CORPUS" in res_text.output

    # JSON run
    assert runner.invoke(cli, ["self-update", "--json"]).exit_code == 2
    res_json = runner.invoke(cli, ["self-update", "--check", "--json"])
    assert res_json.exit_code == 0
    parsed = json.loads(res_json.output)
    assert "package_name" in parsed
    assert "installed_version" in parsed
    assert "status" in parsed
    assert "epistemic" in parsed


def test_zero_mutation_guarantee():
    """Verify that self-update never invokes subprocess execution or installer calls."""
    runner = CliRunner()
    with patch("subprocess.run") as mock_sub_run, \
         patch("subprocess.Popen") as mock_sub_popen, \
         patch("os.system") as mock_os_system:
        res = runner.invoke(cli, ["self-update", "--check"])
        assert res.exit_code == 0

        mock_sub_run.assert_not_called()
        mock_sub_popen.assert_not_called()
        mock_os_system.assert_not_called()


def test_format_update_report():
    """Verify human-readable format generation."""
    sample = {
        "package_name": "factorforge-cds",
        "installed_version": "3.5.0",
        "installer": "pip",
        "status": "UPDATE_AVAILABLE",
        "message": "A new version of FactorForge is available",
        "latest_version": "3.6.0",
        "upgrade_command": "python -m pip install --upgrade factorforge-cds",
        "epistemic": "UNVERIFIED_PUBLIC_CORPUS · Biological validation: NOT_AVAILABLE",
    }
    text = format_update_report(sample)
    assert "Installed Version:        3.5.0" in text
    assert "Latest Remote Version:    3.6.0" in text
    assert "To upgrade manually, run:" in text
    assert "python -m pip install --upgrade factorforge-cds" in text
