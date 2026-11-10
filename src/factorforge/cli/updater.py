"""FactorForge explicit, read-only package updater checker.

Inspects the running interpreter's package manager metadata (PEP 610 direct_url,
INSTALLER, PEP 668 EXTERNALLY-MANAGED) and checks remote distribution channels
strictly in read-only mode without executing any installation subprocesses.

EPISTEMIC NOTICE:
UNVERIFIED_PUBLIC_CORPUS · Biological validation: NOT_AVAILABLE
Explicit read-only update checker · Zero mutation guarantee
"""

from __future__ import annotations

import importlib.metadata
import json
import math
import sys
import sysconfig
import urllib.request
from pathlib import Path
from typing import Any

from packaging.specifiers import InvalidSpecifier, SpecifierSet
from packaging.tags import sys_tags
from packaging.utils import InvalidWheelFilename, parse_wheel_filename
from packaging.version import InvalidVersion, Version
from packaging.version import parse as parse_version

DEFAULT_PYPI_TIMEOUT = 2.0
PRIMARY_DIST_NAME = "factorforge-cds"
FALLBACK_DIST_NAME = "factorforge"


def get_installed_distribution_info(
    dist_name: str = PRIMARY_DIST_NAME,
) -> dict[str, Any]:
    """Inspect the running environment and distribution metadata for FactorForge.

    Reads from importlib.metadata, checks PEP 610 direct_url.json, PEP 376 INSTALLER,
    and PEP 668 EXTERNALLY-MANAGED.

    Returns:
        Dictionary containing resolved metadata:
            - name: Distribution name
            - installed_version: Version string from metadata or __version__
            - installer: 'pip', 'conda', 'uv', 'editable', or 'unknown'
            - is_editable: True if installed in editable mode (PEP 610)
            - direct_url: Parsed PEP 610 direct_url payload if present
            - is_externally_managed: True if EXTERNALLY-MANAGED marker is present
            - sys_prefix: Path of current sys.prefix
    """
    resolved_name = dist_name
    dist: importlib.metadata.Distribution | None = None

    try:
        dist = importlib.metadata.Distribution.from_name(dist_name)
    except importlib.metadata.PackageNotFoundError:
        try:
            dist = importlib.metadata.Distribution.from_name(FALLBACK_DIST_NAME)
            resolved_name = FALLBACK_DIST_NAME
        except importlib.metadata.PackageNotFoundError:
            dist = None

    installed_ver: str = "unknown"
    if dist is not None:
        installed_ver = dist.version
    else:
        try:
            import factorforge
            installed_ver = getattr(factorforge, "__version__", "unknown")
        except Exception:
            installed_ver = "unknown"

    installer = "unknown"
    is_editable = False
    direct_url_data: dict[str, Any] | None = None

    if dist is not None:
        # Check PEP 610 direct_url.json
        try:
            direct_url_text = dist.read_text("direct_url.json")
            if direct_url_text is not None:
                direct_url_data = {}
                direct_url_data = json.loads(direct_url_text)
                if not isinstance(direct_url_data, dict):
                    direct_url_data = {}
                else:
                    dir_info = direct_url_data.get("dir_info", {})
                    if dir_info.get("editable") is True:
                        is_editable = True
                        installer = "editable"
        except Exception:
            pass

        # Check PEP 376 INSTALLER if not editable
        if not is_editable:
            try:
                installer_text = dist.read_text("INSTALLER")
                if installer_text:
                    installer = installer_text.strip().lower()
            except Exception:
                pass

    # Check conda environment indicators
    if installer != "editable":
        conda_meta = Path(sys.prefix) / "conda-meta"
        if conda_meta.is_dir():
            matches = list(conda_meta.glob(f"{resolved_name}-*.json"))
            if matches:
                installer = "conda"

    # Check PEP 668 externally managed environment
    is_ext_managed = False
    prefix_path = Path(sys.prefix)
    if sys.prefix == sys.base_prefix and (Path(sysconfig.get_path("stdlib")) / "EXTERNALLY-MANAGED").is_file():
        is_ext_managed = True

    return {
        "name": resolved_name,
        "installed_version": installed_ver,
        "installer": installer,
        "is_editable": is_editable,
        "direct_url": direct_url_data,
        "is_externally_managed": is_ext_managed,
        "sys_prefix": str(prefix_path),
        "python_version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
    }


def query_pypi_latest_version(
    package_name: str = PRIMARY_DIST_NAME,
    timeout: float = DEFAULT_PYPI_TIMEOUT,
    target_python: str | None = None,
) -> dict[str, Any]:
    """Query PyPI JSON API for the latest stable, compatible release.

    Performs bounded network query with strict timeout and fail-open contract.
    Filters out prereleases and releases incompatible with the target Python.

    Returns:
        Dictionary:
            - latest_version: Semantic version string or None
            - available: Boolean indicating query success
            - error: Error description if unavailable
    """
    if not math.isfinite(timeout) or timeout <= 0:
        return {"latest_version": None, "available": False, "error": "Invalid lookup timeout"}
    url = f"https://pypi.org/pypi/{package_name}/json"
    req = urllib.request.Request(
        url,
        headers={"User-Agent": f"factorforge-self-update-checker/1.0 ({sys.platform})"},
    )

    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            if response.status != 200:
                return {
                    "latest_version": None,
                    "available": False,
                    "error": f"PyPI returned HTTP status {response.status}",
                }
            payload = json.loads(response.read().decode("utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return {
            "latest_version": None,
            "available": False,
            "error": "PyPI lookup failed or returned invalid data",
        }

    releases = payload.get("releases") if isinstance(payload, dict) else None
    if not isinstance(releases, dict) or not releases:
        return {"latest_version": None, "available": False, "error": "No verifiable release artifacts"}

    curr_py_str = target_python or f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    curr_py_ver = parse_version(curr_py_str)

    valid_versions: list[Version] = []
    supported_tags = set(sys_tags())

    for ver_str, artifacts in releases.items():
        if not isinstance(artifacts, list):
            continue
        try:
            ver_obj = parse_version(ver_str)
        except InvalidVersion:
            continue

        # Skip prereleases and development releases under standard policy
        if ver_obj.is_prerelease or ver_obj.is_devrelease:
            continue

        # Check Python version compatibility against release artifacts
        compatible = False
        for art in artifacts:
            if not isinstance(art, dict) or art.get("yanked", False):
                continue
            artifact_type = art.get("packagetype")
            if artifact_type == "bdist_wheel":
                filename = art.get("filename")
                if not isinstance(filename, str):
                    continue
                try:
                    _, _, _, wheel_tags = parse_wheel_filename(filename)
                except InvalidWheelFilename:
                    continue
                if not supported_tags.intersection(wheel_tags):
                    continue
            elif artifact_type != "sdist":
                continue
            req_py = art.get("requires_python")
            if req_py:
                try:
                    spec = SpecifierSet(req_py)
                    if curr_py_ver not in spec:
                        continue
                except (InvalidSpecifier, TypeError):
                    continue
            compatible = True
            break

        if compatible:
            valid_versions.append(ver_obj)

    if not valid_versions:
        return {
            "latest_version": None,
            "available": False,
            "error": "No compatible stable release found on PyPI",
        }

    valid_versions.sort()
    latest_obj = valid_versions[-1]

    return {
        "latest_version": str(latest_obj),
        "available": True,
        "error": None,
    }


def check_factorforge_update(
    dist_name: str = PRIMARY_DIST_NAME,
    timeout: float = DEFAULT_PYPI_TIMEOUT,
    mock_pypi_result: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Execute an explicit, read-only check for FactorForge distribution updates.

    Determines ownership and queries the corresponding channel strictly in
    read-only mode. Never invokes installer subprocesses.

    Returns:
        Structured update result dictionary.
    """
    info = get_installed_distribution_info(dist_name)
    current_ver = info["installed_version"]
    installer = info["installer"]
    is_editable = info["is_editable"]

    result: dict[str, Any] = {
        "package_name": info["name"],
        "installed_version": current_ver,
        "installer": installer,
        "is_editable": is_editable,
        "is_externally_managed": info["is_externally_managed"],
        "status": "CHECK_UNAVAILABLE",
        "latest_version": None,
        "upgrade_command": None,
        "message": "",
        "epistemic": "UNVERIFIED_PUBLIC_CORPUS · Biological validation: NOT_AVAILABLE",
    }

    # 1. Editable / Development mode handling
    if is_editable:
        result["status"] = "DEVELOPMENT_MODE"
        result["message"] = (
            "FactorForge is installed in development/editable mode from source. "
            "Index updates are not applicable."
        )
        return result

    # 2. Conda / Bioconda package handling
    if installer == "conda":
        result["status"] = "CHECK_UNAVAILABLE"
        result["message"] = (
            "Installed via Conda/Bioconda channel. Querying PyPI is disallowed "
            "under package ownership separation."
        )
        return result

    if (
        installer not in {"pip", "uv"} or info.get("direct_url") is not None
        or info["is_externally_managed"] or info["name"] != PRIMARY_DIST_NAME
    ):
        result["message"] = "A public PyPI comparison is unavailable for this installation."
        return result

    pypi_data = mock_pypi_result if mock_pypi_result is not None else query_pypi_latest_version(
        package_name=info["name"],
        timeout=timeout,
    )

    if not pypi_data.get("available"):
        result["status"] = "CHECK_UNAVAILABLE"
        result["message"] = f"Public PyPI comparison unavailable: {pypi_data.get('error', 'Unknown error')}"
        return result

    latest_ver_str = pypi_data["latest_version"]
    result["latest_version"] = latest_ver_str

    try:
        c_ver = parse_version(current_ver)
        l_ver = parse_version(latest_ver_str)

        if c_ver < l_ver:
            result["status"] = "UPDATE_AVAILABLE"
            result["message"] = f"Public PyPI reports a newer version: v{current_ver} -> v{latest_ver_str}"
        else:
            result["status"] = "UP_TO_DATE"
            result["message"] = f"No newer matching release on public PyPI (installed v{current_ver})."
    except (InvalidVersion, TypeError):
        result["status"] = "CHECK_UNAVAILABLE"
        result["message"] = "Version comparison failed."

    return result


def format_update_report(result: dict[str, Any]) -> str:
    """Format structured update result into human-readable terminal text."""
    lines = [
        f"FactorForge Distribution: {result['package_name']}",
        f"Installed Version:        {result['installed_version']}",
        f"Detected Installer:       {result['installer']}",
        f"Status:                   {result['status']}",
        f"Details:                  {result['message']}",
    ]

    if result.get("latest_version"):
        lines.append(f"Latest Remote Version:    {result['latest_version']}")

    if result.get("upgrade_command"):
        lines.append("")
        lines.append("To upgrade manually, run:")
        lines.append(f"  {result['upgrade_command']}")

    lines.append("")
    lines.append(f"Epistemic: {result['epistemic']}")
    return "\n".join(lines)
