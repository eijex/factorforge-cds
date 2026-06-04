#!/usr/bin/env python3
"""
Bump FactorForge version across all version-bearing files.

Usage:
    python release.py 3.1.7
    python release.py 3.2.0 --dry-run
    python release.py 3.2.0 --strict   # fail on any missing pattern

Publication sync (keep documentation in sync with feature changes):
  - Patch release: no update needed unless claims changed
  - Feature release (minor): review Software Design documentation
  - Algorithm added (tAI, codon pair bias, 5' UTR MFE): update State of the Field documentation
  - New host added: update Statement of Need + supported hosts documentation
  - Zenodo DOI updated: update citation metadata
"""

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent

# (file_path, [(exact_old_string, exact_new_string)], required=True)
# required=True → missing pattern increments errors (fails with --strict)
# required=False → missing pattern is a WARN only


def build_workspace_targets(old: str, new: str, workspace: Path) -> list[tuple[Path, list[tuple[str, str]], bool]]:
    """Optional cross-repo targets in the eijex-workspace (planning/tracking repo)."""
    return [
        (workspace / "README.md", [
            (f"FactorForge v{old}", f"FactorForge v{new}"),
        ], False),
        (workspace / "ROADMAP.md", [
            (f"FactorForge v{old}", f"FactorForge v{new}"),
        ], False),
        (workspace / "CLAUDE.md", [
            (f"pyproject.toml: 현재 v{old}", f"pyproject.toml: 현재 v{new}"),
        ], False),
        (workspace / "refs/papers/SUBMISSION_STRATEGY.md", [
            (f"v{old} 공개", f"v{new} 공개"),
        ], False),
        (workspace / "refs/papers/factorforge-joss/README.md", [
            (f"(v{old})", f"(v{new})"),
        ], False),
    ]


def build_mcp_targets(old: str, new: str, mcp: Path) -> list[tuple[Path, list[tuple[str, str]], bool]]:
    """Optional cross-repo targets in the eijex-mcp repo."""
    return [
        (mcp / "src/app/_lib/mcp-tools.ts", [
            (f"FactorForge v{old} stable design path", f"FactorForge v{new} stable design path"),
            (f"DP feasibility design (v{old})", f"DP feasibility design (v{new})"),
        ], False),
    ]


def build_targets(old: str, new: str) -> list[tuple[str, list[tuple[str, str]], bool]]:
    return [
        ("pyproject.toml", [
            (f'version = "{old}"', f'version = "{new}"'),
        ], True),
        ("CITATION.cff", [
            (f"version: {old}", f"version: {new}"),
        ], True),
        ("src/factorforge/__init__.py", [
            (f'__version__ = "{old}"', f'__version__ = "{new}"'),
        ], True),
        ("src/factorforge/engines/__init__.py", [
            (f'"version": "{old}"', f'"version": "{new}"'),
        ], True),
        ("src/factorforge/engines/profile/__init__.py", [
            (f'__version__ = "{old}"', f'__version__ = "{new}"'),
        ], True),
        ("src/factorforge/engines/profile/optimizer.py", [
            (f'version = "{old}"', f'version = "{new}"'),
        ], True),
        ("api/optimize.py", [
            (f"Product Version: {old}", f"Product Version: {new}"),
            (f'"product": "{old}"', f'"product": "{new}"'),
            (f'"rule_engine": "{old}"', f'"rule_engine": "{new}"'),
            (f'"dp_engine": "{old}"', f'"dp_engine": "{new}"'),
        ], True),
        ("web/index.html", [
            (f"v{old} Release Notes", f"v{new} Release Notes"),
        ], True),
        ("web/js/app.js", [
            (f"FactorForge v{old} Engaged", f"FactorForge v{new} Engaged"),
            (f"?.product || '{old}'", f"?.product || '{new}'"),
        ], True),
        ("ROADMAP.md", [
            (f"**v{old}**", f"**v{new}**"),
        ], True),
        ("tests/api/test_optimize_contract.py", [
            (f'"product"] == "{old}"', f'"product"] == "{new}"'),
        ], True),
        ("README.md", [
            (f"FactorForge v{old} (", f"FactorForge v{new} ("),
        ], True),
        ("tests/engines/profile/test_cli_optimize.py", [
            (f"Profile-based v{old}", f"Profile-based v{new}"),
        ], True),
        ("docs/tutorials/gfp-nbenthamiana.md", [
            (f"Profile-based v{old}", f"Profile-based v{new}"),
        ], True),
        ("tests/test_validation/test_package_generator.py", [
            (f'factorforge_version="{old}"', f'factorforge_version="{new}"'),
            (f'            "{old}",\n            "--profile"', f'            "{new}",\n            "--profile"'),
        ], True),
        ("tests/test_schemas/test_design_package.py", [
            (f'"product_version": "{old}"', f'"product_version": "{new}"'),
        ], True),
    ]


def _today() -> str:
    from datetime import date
    return date.today().isoformat()


def _update_citation_date(path: Path, dry_run: bool) -> list[str]:
    """Update CITATION.cff date-released to today."""
    today = _today()
    content = path.read_text(encoding="utf-8")
    updated = re.sub(r'date-released: "\d{4}-\d{2}-\d{2}"', f'date-released: "{today}"', content)
    if content == updated:
        return []
    if not dry_run:
        path.write_text(updated, encoding="utf-8")
    return [f"  date-released → {today}"]


def _check_residual(old: str, dry_run: bool) -> int:
    """Scan for leftover occurrences of old version string after bump."""
    residual_errors = 0
    skip_dirs = {'.git', '__pycache__', '.pytest_cache', 'build', 'dist', 'archive', '.venv', 'venv'}
    skip_files = {'CHANGELOG.md', 'docs/changelog.md'}

    for path in ROOT.rglob('*'):
        if path.is_dir():
            continue
        if any(p in path.parts for p in skip_dirs):
            continue
        if path.name in skip_files:
            continue
        suffix = path.suffix.lower()
        if suffix not in {'.py', '.toml', '.yml', '.yaml', '.md', '.html', '.js', '.cff', '.json', '.txt'}:
            continue
        try:
            text = path.read_text(encoding='utf-8', errors='ignore')
            if old in text:
                rel = path.relative_to(ROOT)
                print(f"  RESIDUAL: {rel} still contains {old!r}")
                residual_errors += 1
        except Exception:
            pass
    return residual_errors


def bump(old: str, new: str, dry_run: bool = False, strict: bool = False, workspace: Path | None = None, mcp: Path | None = None) -> int:
    targets = build_targets(old, new)
    errors = 0
    total_changes = 0

    for rel_path, replacements, required in targets:
        path = ROOT / rel_path
        if not path.exists():
            msg = f"  {'ERROR' if required else 'SKIP'} (not found): {rel_path}"
            print(msg)
            if required:
                errors += 1
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
                if required:
                    print(f"  {'ERROR' if strict else 'WARN'}: pattern not found in {rel_path}: {old_str!r}")
                    if strict:
                        errors += 1
                else:
                    print(f"  WARN: pattern not found in {rel_path}: {old_str!r}")

        # Update CITATION.cff date-released
        if rel_path == "CITATION.cff":
            date_changes = _update_citation_date(path, dry_run=True)
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

    # Workspace targets (PlantFormOrg or other tracking repo)
    if workspace is not None:
        if not workspace.is_dir():
            print(f"  WARN: --workspace path not found: {workspace}")
        else:
            print(f"\n{'[DRY RUN] ' if dry_run else ''}Workspace: {workspace}")
            for abs_path, replacements, required in build_workspace_targets(old, new, workspace):
                if not abs_path.exists():
                    print(f"  {'ERROR' if required else 'SKIP'} (not found): {abs_path.name}")
                    if required:
                        errors += 1
                    continue
                content = abs_path.read_text(encoding="utf-8")
                original = content
                changes = []
                for old_str, new_str in replacements:
                    if old_str in content:
                        content = content.replace(old_str, new_str)
                        changes.append(f"  {old_str!r} → {new_str!r}")
                    else:
                        print(f"  WARN: pattern not found in {abs_path.name}: {old_str!r}")
                if content != original:
                    total_changes += 1
                    print(f"\n{'[DRY RUN] ' if dry_run else ''}Updated: {abs_path.relative_to(workspace)}")
                    for c in changes:
                        print(c)
                    if not dry_run:
                        abs_path.write_text(content, encoding="utf-8")

    # MCP targets (eijex-mcp repo)
    if mcp is not None:
        if not mcp.is_dir():
            print(f"  WARN: --mcp path not found: {mcp}")
        else:
            print(f"\n{'[DRY RUN] ' if dry_run else ''}MCP: {mcp}")
            for abs_path, replacements, required in build_mcp_targets(old, new, mcp):
                if not abs_path.exists():
                    print(f"  {'ERROR' if required else 'SKIP'} (not found): {abs_path.name}")
                    if required:
                        errors += 1
                    continue
                content = abs_path.read_text(encoding="utf-8")
                original = content
                changes = []
                for old_str, new_str in replacements:
                    if old_str in content:
                        content = content.replace(old_str, new_str)
                        changes.append(f"  {old_str!r} → {new_str!r}")
                    else:
                        print(f"  WARN: pattern not found in {abs_path.name}: {old_str!r}")
                if content != original:
                    total_changes += 1
                    print(f"\n{'[DRY RUN] ' if dry_run else ''}Updated: {abs_path.relative_to(mcp)}")
                    for c in changes:
                        print(c)
                    if not dry_run:
                        abs_path.write_text(content, encoding="utf-8")

    if errors:
        print(f"Errors: {errors} (patterns not found in required files)")

    if dry_run:
        print("No files were written (--dry-run).")
        return errors

    print(f"Version bumped: {old} → {new}")

    # Residual check: scan for leftover old version strings
    print(f"\nChecking for residual {old!r} strings...")
    residual = _check_residual(old, dry_run=False)
    if residual:
        print(f"  WARNING: {residual} file(s) still contain old version string — check manually")
    else:
        print(f"  OK: no residual {old!r} strings found")

    print("\nNext steps:")
    print("  --- Pre-release gate (if not done) ---")
    print("  0a. git status --short  (must be clean before bump)")
    print("  0b. python -m ruff check .")
    print("  0c. python -m pytest tests/ -v --tb=short")
    print()
    print("  --- Changelog & docs (manual) ---")
    print("  1. Move [Unreleased] entries to a new [X.Y.Z] block in CHANGELOG.md; update comparison links")
    print(f"  2. Add v{new} entry to web/index.html changelog panel")
    print("     → Copy the previous version block, update version/date/bullet points")
    print("     → Set previous 'Current' block: remove emerald classes, remove Current badge")
    print("     → Set new block: border-emerald-500, dot bg-emerald-500, add Current badge")
    print("  3. Add a summary entry in docs/changelog.md")
    print()
    print("  --- Commit & CI gate (before tagging) ---")
    print(f"  4. git commit -m 'chore: release v{new}'")
    print("  5. git push")
    print("  6. Wait for CI to pass (github.com/eijex/factorforge-cds/actions)")
    print()
    print("  --- Tag & publish ---")
    print(f"  7. git tag -a v{new} -m 'Release v{new}' && git push --tags")
    print("  8. GitHub Actions publishes PyPI + Docker + GitHub Release + Zenodo automatically")
    print()
    print("  --- Post-release verification ---")
    print(f"  9.  pip install factorforge-cds=={new} && factorforge --help  (PyPI smoke test)")
    print(f" 10.  docker run ghcr.io/eijex/factorforge-cds:v{new} factorforge --help  (Docker smoke test)")
    print(" 11.  Confirm Zenodo DOI: https://zenodo.org/doi/10.5281/zenodo.20407331")
    print(" 12.  Update Bioconda recipes/meta.yaml (version + SHA256) → push fork branch")
    print(" 13.  Close completed GitHub Issues; close milestone if all done")

    return errors


def main() -> None:
    parser = argparse.ArgumentParser(description="Bump FactorForge version strings.")
    parser.add_argument("new_version", help="New version, e.g. 3.1.7")
    parser.add_argument("--from", dest="old_version", default=None,
                        help="Old version to replace (auto-detected from pyproject.toml if omitted)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would change without writing files")
    parser.add_argument("--strict", action="store_true",
                        help="Exit with error if any required pattern is not found")
    parser.add_argument("--workspace", default=None,
                        help="Path to eijex-workspace to also bump cross-repo docs (default: C:\\Work\\eijex\\eijex-workspace)")
    parser.add_argument("--mcp", default=None,
                        help="Path to eijex-mcp repo to also bump MCP tool version strings (e.g. C:\\Work\\eijex\\eijex-mcp)")
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

    workspace = Path(args.workspace) if args.workspace else None
    mcp = Path(args.mcp) if args.mcp else None
    print(f"Bumping {old} → {new}{' (dry run)' if args.dry_run else ''}{' [strict]' if args.strict else ''}\n")
    errors = bump(old, new, dry_run=args.dry_run, strict=args.strict, workspace=workspace, mcp=mcp)
    sys.exit(errors)


if __name__ == "__main__":
    main()
