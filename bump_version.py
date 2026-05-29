#!/usr/bin/env python3
"""
Bump FactorForge version across all version-bearing files.

Usage:
    python bump_version.py 3.1.4
    python bump_version.py 3.2.0 --dry-run

After running, manually:
  1. Add a new [X.Y.Z] entry in CHANGELOG.md (move from [Unreleased])
  2. Add a new entry in web/index.html changelog panel (version-specific HTML)
  3. Add a new entry in docs/changelog.md (summarized)
  4. git commit -m "chore: release vX.Y.Z"
  5. git tag -a vX.Y.Z -m "Release vX.Y.Z" && git push && git push --tags
  6. GitHub Actions automatically creates GitHub Release from the tag

JOSS paper review (paper.md — keep in sync with feature changes):
  - Patch release: no update needed unless claims changed
  - Feature release (minor, e.g. v3.2 → v3.3): review Software Design section
  - Algorithm added (tAI, codon pair bias, 5' UTR MFE): update State of the Field differentiators
  - New host added: update Statement of Need + State of the Field
  - Zenodo DOI updated: update paper.bib @kim2026factorforge doi field
"""

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).parent

# (file_path, [(exact_old_string, exact_new_string), ...])
# Strings are literal — no regex. All occurrences of old_string are replaced.
TARGETS: list[tuple[str, list[tuple[str, str]]]] = []


def build_targets(old: str, new: str) -> list[tuple[str, list[tuple[str, str]]]]:
    return [
        ("pyproject.toml", [
            (f'version = "{old}"', f'version = "{new}"'),
        ]),
        ("CITATION.cff", [
            (f"version: {old}", f"version: {new}"),
            (f'date-released: "{_today()}"', f'date-released: "{_today()}"'),  # no-op, keeps date current
        ]),
        ("src/factorforge/__init__.py", [
            (f'__version__ = "{old}"', f'__version__ = "{new}"'),
        ]),
        ("src/factorforge/engines/__init__.py", [
            (f'"version": "{old}"', f'"version": "{new}"'),
        ]),
        ("src/factorforge/engines/profile/__init__.py", [
            (f'__version__ = "{old}"', f'__version__ = "{new}"'),
        ]),
        ("src/factorforge/engines/profile/optimizer.py", [
            (f'version = "{old}"', f'version = "{new}"'),
        ]),
        ("api/optimize.py", [
            (f"Product Version: {old}", f"Product Version: {new}"),
            (f'"product": "{old}"', f'"product": "{new}"'),
            (f'"rule_engine": "{old}"', f'"rule_engine": "{new}"'),
            (f'"dp_engine": "{old}"', f'"dp_engine": "{new}"'),
        ]),
        ("web/index.html", [
            (f"v{old} Release Notes", f"v{new} Release Notes"),
        ]),
        ("web/js/app.js", [
            (f"FactorForge v{old} Engaged", f"FactorForge v{new} Engaged"),
            (f"?.product || '{old}'", f"?.product || '{new}'"),
        ]),
        ("ROADMAP.md", [
            (f"**v{old}**", f"**v{new}**"),
        ]),
        ("tests/api/test_optimize_contract.py", [
            (f'"product"] == "{old}"', f'"product"] == "{new}"'),
        ]),
        ("README.md", [
            (f"FactorForge v{old} (", f"FactorForge v{new} ("),
        ]),
        ("tests/engines/profile/test_cli_optimize.py", [
            (f"Profile-based v{old}", f"Profile-based v{new}"),
        ]),
        ("docs/tutorials/gfp-nbenthamiana.md", [
            (f"Profile-based v{old}", f"Profile-based v{new}"),
        ]),
        ("tests/test_validation/test_package_generator.py", [
            (f'factorforge_version="{old}"', f'factorforge_version="{new}"'),
            # "--version" arg appears on its own line: '            "3.1.3",'
            (f'            "{old}",\n            "--profile"', f'            "{new}",\n            "--profile"'),
        ]),
        ("tests/test_schemas/test_design_package.py", [
            (f'"product_version": "{old}"', f'"product_version": "{new}"'),
        ]),
    ]


def _today() -> str:
    from datetime import date
    return date.today().isoformat()


def _update_citation_date(path: Path, new: str, dry_run: bool) -> list[str]:
    """Update CITATION.cff date-released to today."""
    today = _today()
    content = path.read_text(encoding="utf-8")
    updated = re.sub(r'date-released: "\d{4}-\d{2}-\d{2}"', f'date-released: "{today}"', content)
    if content == updated:
        return []
    if not dry_run:
        path.write_text(updated, encoding="utf-8")
    return [f"  date-released → {today}"]


def bump(old: str, new: str, dry_run: bool = False) -> int:
    targets = build_targets(old, new)
    errors = 0
    total_changes = 0

    for rel_path, replacements in targets:
        path = ROOT / rel_path
        if not path.exists():
            print(f"  SKIP (not found): {rel_path}")
            continue

        content = path.read_text(encoding="utf-8")
        original = content
        changes = []

        for old_str, new_str in replacements:
            if old_str == new_str:
                continue
            if old_str in content:
                content = content.replace(old_str, new_str)
                changes.append(f"  {old_str!r} → {new_str!r}")
            else:
                # String not found — warn unless it's an optional pattern
                print(f"  WARN: pattern not found in {rel_path}: {old_str!r}")

        # Special: update CITATION.cff date
        if rel_path == "CITATION.cff":
            date_changes = _update_citation_date(path, new, dry_run=True)
            if date_changes:
                today = _today()
                content = re.sub(r'date-released: "\d{4}-\d{2}-\d{2}"', f'date-released: "{today}"', content)
                changes.extend(date_changes)

        if content != original:
            total_changes += 1
            print(f"\n{'[DRY RUN] ' if dry_run else ''}Updated: {rel_path}")
            for c in changes:
                print(c)
            if not dry_run:
                path.write_text(content, encoding="utf-8")

    print(f"\n{'[DRY RUN] ' if dry_run else ''}{'─' * 40}")
    print(f"Files modified: {total_changes}")
    if dry_run:
        print("No files were written (--dry-run).")
    else:
        print(f"Version bumped: {old} → {new}")
        print("\nNext steps (manual):")
        print("  1. Move [Unreleased] entries to a new [X.Y.Z] block in CHANGELOG.md")
        print(f"  2. Add v{new} entry to web/index.html changelog panel")
        print("     → Copy the previous version block, update version/date/bullet points")
        print("     → Set previous 'Current' block: remove emerald classes, remove Current badge")
        print("     → Set new block: border-emerald-500, dot bg-emerald-500, add Current badge")
        print("  3. Add a summary entry in docs/changelog.md")
        print(f"  4. git commit -m 'chore: release v{new}'")
        print(f"  5. git tag -a v{new} -m 'Release v{new}' && git push && git push --tags")
        print("  6. GitHub Release is created automatically by GitHub Actions")
    return errors


def main() -> None:
    parser = argparse.ArgumentParser(description="Bump FactorForge version strings.")
    parser.add_argument("new_version", help="New version, e.g. 3.1.4")
    parser.add_argument("--from", dest="old_version", default=None,
                        help="Old version to replace (auto-detected from pyproject.toml if omitted)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would change without writing files")
    args = parser.parse_args()

    new = args.new_version
    if not re.fullmatch(r"\d+\.\d+\.\d+", new):
        print(f"ERROR: version must be X.Y.Z, got: {new!r}")
        sys.exit(1)

    old = args.old_version
    if old is None:
        toml = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
        m = re.search(r'^version = "(\d+\.\d+\.\d+)"', toml, re.MULTILINE)
        if not m:
            print("ERROR: could not auto-detect current version from pyproject.toml")
            sys.exit(1)
        old = m.group(1)

    if old == new:
        print(f"Nothing to do — already at {new}")
        sys.exit(0)

    print(f"Bumping {old} → {new}{' (dry run)' if args.dry_run else ''}\n")
    errors = bump(old, new, dry_run=args.dry_run)
    sys.exit(errors)


if __name__ == "__main__":
    main()
