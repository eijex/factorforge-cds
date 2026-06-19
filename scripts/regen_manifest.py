"""Recompute reproducibility/benchmark_v0.5.1/MANIFEST.json input SHA-256 values.

Hashes are computed from the committed git blob (`git show HEAD:<path>`), not
from local working-tree bytes. On Windows, `core.autocrlf=true` combined with
this repo's `.gitattributes` (`eol=lf` for json/yaml/yml) means the working
tree can hold CRLF bytes while git stores LF — hashing `path.read_bytes()`
then produces a value that mismatches what is actually committed (the root
cause of the Job 129 / Job 136 MANIFEST drift incidents).

Usage:
    python scripts/regen_manifest.py            # report diff, exit 1 if any
    python scripts/regen_manifest.py --write     # rewrite MANIFEST.json in place
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "reproducibility" / "benchmark_v0.5.1" / "MANIFEST.json"


def git_blob_sha256(rel_path: str) -> str:
    """SHA-256 of a file's committed content at HEAD, not its working-tree bytes."""
    content = subprocess.run(
        ["git", "show", f"HEAD:{rel_path}"],
        cwd=ROOT,
        capture_output=True,
        check=True,
    ).stdout
    return hashlib.sha256(content).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="rewrite MANIFEST.json with recomputed hashes")
    args = parser.parse_args()

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    mismatches = []
    for name, entry in manifest["inputs"].items():
        recomputed = git_blob_sha256(entry["path"])
        if recomputed != entry["sha256"]:
            mismatches.append((name, entry["path"], entry["sha256"], recomputed))
            if args.write:
                entry["sha256"] = recomputed

    if not mismatches:
        print(f"OK - all {len(manifest['inputs'])} input hashes match git HEAD content.")
        return 0

    print(f"DRIFT - {len(mismatches)} input(s) mismatch committed content:")
    for name, path, old, new in mismatches:
        print(f"  {name} ({path})")
        print(f"    old: {old}")
        print(f"    new: {new}")

    if args.write:
        MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(f"\nRewrote {MANIFEST_PATH.relative_to(ROOT)} with recomputed hashes.")
        return 0

    print("\nRun with --write to apply, or investigate why content changed without a manifest update.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
