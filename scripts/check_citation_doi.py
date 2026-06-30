"""Check that CITATION.cff's exact-version Zenodo DOI matches the latest
release minted under the concept DOI.

Background: `scripts/release.py` supports `--zenodo-doi` to update
CITATION.cff, but that flag can only be run *after* Zenodo asynchronously
mints a new DOI following a GitHub Release — it cannot happen inside the
same release command. That manual follow-up step has been forgotten on
three consecutive releases (v3.2.3, v3.2.4, v3.2.5). This script detects
the drift in CI instead of relying on someone noticing by chance.

Warn-only by default (`--strict` to fail CI on mismatch) — querying
Zenodo over the network is a new external dependency for the CI pipeline,
and a Zenodo outage should not block unrelated merges while this check is
new and unproven.
"""

from __future__ import annotations

import argparse
import re
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CITATION_PATH = ROOT / "CITATION.cff"
CONCEPT_DOI = "10.5281/zenodo.20407330"
DOI_RESOLVER_URL = f"https://doi.org/{CONCEPT_DOI}"


def read_citation_doi(path: Path) -> str | None:
    text = path.read_text(encoding="utf-8")
    m = re.search(r'^doi:\s*"(10\.5281/zenodo\.\d+)"', text, re.MULTILINE)
    return m.group(1) if m else None


def resolve_latest_doi(timeout: float = 10.0) -> str | None:
    """Follow the concept DOI's redirect chain to the latest version's
    exact record URL (e.g. https://zenodo.org/records/20838848) and
    extract the DOI from the record ID. Returns None on any network
    failure (caller decides whether that's fatal)."""
    req = urllib.request.Request(DOI_RESOLVER_URL, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            final_url = resp.geturl()
    except urllib.error.URLError:
        return None
    m = re.search(r"zenodo\.org/records/(\d+)", final_url)
    if not m:
        return None
    return f"10.5281/zenodo.{m.group(1)}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--strict", action="store_true",
        help="Exit 1 on mismatch (default: warn-only, exit 0)",
    )
    parser.add_argument(
        "--citation-path", type=Path, default=CITATION_PATH,
        help="Path to CITATION.cff (for testing with a fixture)",
    )
    args = parser.parse_args()

    current_doi = read_citation_doi(args.citation_path)
    if current_doi is None:
        print(f"WARN: no doi field found in {args.citation_path}")
        return 0

    latest_doi = resolve_latest_doi()
    if latest_doi is None:
        print("WARN: could not resolve Zenodo concept DOI (network issue or "
              "Zenodo unavailable) — skipping comparison, not failing CI.")
        return 0

    if current_doi == latest_doi:
        print(f"OK: CITATION.cff doi ({current_doi}) matches Zenodo's latest release.")
        return 0

    message = (
        f"CITATION.cff doi is stale.\n"
        f"  CITATION.cff has: {current_doi}\n"
        f"  Zenodo latest is: {latest_doi}\n"
        f"  Fix with: python scripts/release.py <version> --zenodo-doi {latest_doi}"
    )
    if args.strict:
        print(f"FAIL: {message}")
        return 1
    print(f"WARN: {message}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
