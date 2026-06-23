#!/usr/bin/env python3
"""
Update recipes/meta.yaml's sha256 after a release's sdist is live on PyPI.

`release.py`'s bump()/auto_release() never touch this sha256 anymore (PyPI
doesn't have the new version yet at that point in the pipeline — the fetch
always failed there). This script is the explicit, separate, post-publish
step: run it only after the version is confirmed live on PyPI (e.g. after
scripts/post_publish_smoke.py / the publish.yml smoke-test job has passed).

It never touches the already-tagged release commit: it creates a dedicated
branch, commits the sha256 update there, and prints the push command —
it does not push or modify any other branch.

Usage:
    python scripts/update_conda_recipe_sha.py 3.2.5
"""

import argparse
import json
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).parent.parent
PACKAGE = "factorforge-cds"


def fetch_pypi_sha256(package: str, version: str) -> str | None:
    """Fetch the sdist sha256 for a given PyPI package version."""
    url = f"https://pypi.org/pypi/{package}/{version}/json"
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = json.loads(resp.read())
    except Exception as exc:
        print(f"ERROR: could not reach PyPI for {package}=={version}: {exc}")
        return None
    for entry in data.get("urls", []):
        if entry.get("packagetype") == "sdist":
            return entry["digests"]["sha256"]
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("version", help="Released version already live on PyPI, e.g. 3.2.5")
    args = parser.parse_args()
    version = args.version

    sha256 = fetch_pypi_sha256(PACKAGE, version)
    if not sha256:
        print(
            f"ERROR: {PACKAGE}=={version} sdist sha256 not found on PyPI yet.\n"
            "Run this only after the post-publish smoke test confirms the version is live."
        )
        sys.exit(1)

    meta_path = ROOT / "recipes" / "meta.yaml"
    content = meta_path.read_text(encoding="utf-8")
    updated = re.sub(r"(sha256:\s*)[0-9a-f]{64}", lambda m: f"{m.group(1)}{sha256}", content)
    version_updated = re.sub(
        r'\{% set version = "[^"]+" %\}', f'{{% set version = "{version}" %}}', updated
    )
    if version_updated == content:
        print(f"Nothing to update — recipes/meta.yaml already matches {version}/{sha256}.")
        sys.exit(0)

    branch = f"update-conda-recipe-sha-v{version}"
    subprocess.run(["git", "checkout", "-b", branch], cwd=str(ROOT), check=True)
    meta_path.write_text(version_updated, encoding="utf-8")
    subprocess.run(["git", "add", "recipes/meta.yaml"], cwd=str(ROOT), check=True)
    subprocess.run(
        ["git", "commit", "-m", f"chore: update conda recipe sha256 for v{version}"],
        cwd=str(ROOT),
        check=True,
    )

    print(f"\nrecipes/meta.yaml updated on branch '{branch}' (sha256: {sha256}).")
    print(f"This commit does NOT touch the v{version} release tag/commit.")
    print(f"Push and open a PR (or push to your bioconda-recipes fork) manually:")
    print(f"  git push origin {branch}")


if __name__ == "__main__":
    main()
