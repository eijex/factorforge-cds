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
import urllib.request
import json
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
    ]


def build_web_targets(old: str, new: str, web: Path) -> list[tuple[Path, list[tuple[str, str]], bool]]:
    """Optional cross-repo targets in the eijex-web repo."""
    return [
        (web / "src/app/components/StatsBar.tsx", [
            (f'"v{old}"', f'"v{new}"'),
        ], False),
    ]


def _update_changelog_current(root: Path, old: str, new: str, dry_run: bool) -> list[str]:
    """Toggle the CURRENT badge in web/index.html changelog panel."""
    path = root / "web/index.html"
    if not path.exists():
        return []
    today = _today()
    content = path.read_text(encoding="utf-8")
    original = content
    h3_old = f'<h3 class="font-bold text-slate-800 dark:text-white">v{old}</h3>'
    h3_new = f'<h3 class="font-bold text-slate-800 dark:text-white">v{new}</h3>'
    content = content.replace(h3_old, h3_new)
    date_pattern = re.compile(
        r'(<h3 class="font-bold text-slate-800 dark:text-white">v' + re.escape(new) + r'</h3>\s*'
        r'<span class="text-slate-400 text-\[10px\]">)(\d{4}-\d{2}-\d{2})(</span>)'
    )
    content = date_pattern.sub(lambda m: f"{m.group(1)}{today}{m.group(3)}", content)
    if content == original:
        return []
    changes = [f"  web/index.html: CURRENT v{old} → v{new}, date → {today}"]
    if not dry_run:
        path.write_text(content, encoding="utf-8")
    return changes


def build_mcp_targets(old: str, new: str, mcp: Path) -> list[tuple[Path, list[tuple[str, str]], bool]]:
    """Optional cross-repo targets in the eijex-mcp repo."""
    return [
        (mcp / "src/app/_lib/mcp-tools.ts", [
            (f"FactorForge v{old} stable design path", f"FactorForge v{new} stable design path"),
            (f"DP feasibility design (v{old})", f"DP feasibility design (v{new})"),
        ], False),
        (mcp / "src/app/api/mcp/route.ts", [
            (f"FactorForge CDS v{old}", f"FactorForge CDS v{new}"),
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
        (".github/ISSUE_TEMPLATE/wet_lab_result.yml", [
            (f'"e.g. {old}"', f'"e.g. {new}"'),
        ], True),
        ("recipes/meta.yaml", [
            (f'{{% set version = "{old}" %}}', f'{{% set version = "{new}" %}}'),
        ], True),
    ]


def _fetch_pypi_sha256(package: str, version: str) -> str | None:
    """Fetch the sdist sha256 for a given PyPI package version."""
    url = f"https://pypi.org/pypi/{package}/{version}/json"
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = json.loads(resp.read())
        for entry in data.get("urls", []):
            if entry.get("packagetype") == "sdist":
                return entry["digests"]["sha256"]
        return None
    except Exception:
        return None


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


def _update_citation_doi(path: Path, new_doi: str, dry_run: bool) -> list[str]:
    """Update CITATION.cff doi field to the new Zenodo version DOI."""
    content = path.read_text(encoding="utf-8")
    updated = re.sub(r'doi: "10\.5281/zenodo\.\d+"', f'doi: "{new_doi}"', content)
    if content == updated:
        return []
    if not dry_run:
        path.write_text(updated, encoding="utf-8")
    return [f"  doi → {new_doi}"]


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


def bump(old: str, new: str, dry_run: bool = False, strict: bool = False, workspace: Path | None = None, mcp: Path | None = None, web: Path | None = None, zenodo_doi: str | None = None) -> int:
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

        # Update recipes/meta.yaml sha256 via PyPI fetch
        if rel_path == "recipes/meta.yaml" and not dry_run:
            print(f"  Fetching sha256 for factorforge-cds {new} from PyPI...")
            sha256 = _fetch_pypi_sha256("factorforge-cds", new)
            if sha256:
                old_sha_pattern = re.compile(r'(sha256:\s*)[0-9a-f]{64}')
                updated = old_sha_pattern.sub(lambda m: f"{m.group(1)}{sha256}", content)
                if updated != content:
                    content = updated
                    changes.append(f"  sha256 → {sha256}")
                    print(f"  sha256 updated: {sha256}")
            else:
                placeholder = "FIXME-fetch-failed"
                old_sha_pattern = re.compile(r'(sha256:\s*)[0-9a-f]{64}')
                updated = old_sha_pattern.sub(lambda m: f"{m.group(1)}{placeholder}", content)
                if updated != content:
                    content = updated
                    changes.append(f"  sha256 → {placeholder} (fetch failed — update manually)")
                print(f"  WARN: could not fetch sha256 from PyPI for factorforge-cds {new} — set placeholder")

        # Update CITATION.cff date-released and optional doi
        if rel_path == "CITATION.cff":
            date_changes = _update_citation_date(path, dry_run=True)
            if date_changes:
                today = _today()
                content = re.sub(r'date-released: "\d{4}-\d{2}-\d{2}"', f'date-released: "{today}"', content)
                changes.extend(date_changes)
            if zenodo_doi:
                doi_changes = _update_citation_doi(path, zenodo_doi, dry_run=True)
                if doi_changes:
                    content = re.sub(r'doi: "10\.5281/zenodo\.\d+"', f'doi: "{zenodo_doi}"', content)
                    changes.extend(doi_changes)

        if content != original:
            total_changes += 1
            print(f"\n{'[DRY RUN] ' if dry_run else ''}Updated: {rel_path}")
            for c in changes:
                print(c)
            if not dry_run:
                path.write_text(content, encoding="utf-8")

    # Update web/index.html changelog CURRENT version/date
    changelog_changes = _update_changelog_current(ROOT, old, new, dry_run=dry_run)
    if changelog_changes:
        total_changes += 1
        print(f"\n{'[DRY RUN] ' if dry_run else ''}Updated: web/index.html (changelog CURRENT)")
        for c in changelog_changes:
            print(c)

    print(f"\n{'[DRY RUN] ' if dry_run else ''}{'─' * 40}")
    print(f"Files modified: {total_changes}")

    # Workspace targets (eijex-workspace or other tracking repo)
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

    # Web targets (eijex-web repo)
    if web is not None:
        if not web.is_dir():
            print(f"  WARN: --web path not found: {web}")
        else:
            print(f"\n{'[DRY RUN] ' if dry_run else ''}Web: {web}")
            for abs_path, replacements, required in build_web_targets(old, new, web):
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
                    print(f"\n{'[DRY RUN] ' if dry_run else ''}Updated: {abs_path.relative_to(web)}")
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
    print(f"  2. Add v{new} entry to web/index.html changelog panel  ← NEW block with bullet points (manual)")
    print("     → Copy the v{old} block above it; write v{new} bullet points")
    print("     → New block: border-emerald-500, dot bg-emerald-500, Current badge")
    print("     → Old block: border-slate-200, dot bg-slate-300, no Current badge")
    print("     (Version number + date in the existing CURRENT block are auto-bumped by this script)")
    print("  3. Add a summary entry in docs/changelog.md")
    print()
    print("  --- Commit & CI gate (before tagging) ---")
    print(f"  4. git commit -m 'chore: release v{new}'")
    print("  5. git push")
    print("  6. Wait for CI to pass (github.com/eijex/factorforge-cds/actions)")
    print()
    print("  --- Public surface audit (before tagging) ---")
    print(f"  6a. python ~/.codex/skills/factorforge-public-surface-audit/scripts/audit_public_surface.py \\")
    print(f"        --workspace C:\\Work\\eijex --live")
    print("        → fix any findings before tagging")
    print()
    print("  --- Tag & publish ---")
    print(f"  7. git tag -a v{new} -m 'Release v{new}' && git push --tags")
    print("  8. GitHub Actions publishes PyPI + Docker + GitHub Release + Zenodo automatically")
    print()
    print("  --- Post-release verification ---")
    print(f"  9.  pip install factorforge-cds=={new} && factorforge --help  (PyPI smoke test)")
    print(f" 10.  docker run ghcr.io/eijex/factorforge-cds:v{new} factorforge --help  (Docker smoke test)")
    print(" 11.  Wait for Zenodo to mint the new DOI (triggered by GitHub Release)")
    print("       Then update CITATION.cff doi:")
    print(f"         python scripts/release.py {new} --zenodo-doi 10.5281/zenodo.<NEW_RECORD_ID>")
    print(" 12.  Push Bioconda fork branch (recipes/meta.yaml already bumped by this script)")
    print("       → git push origin add-factorforge-cds  (or open/update PR on bioconda/bioconda-recipes)")
    print(" 13.  Close completed GitHub Issues; close milestone if all done")
    print()
    print("  --- Post-release external audit ---")
    print(f" 14.  python ~/.codex/skills/factorforge-public-surface-audit/scripts/audit_public_surface.py \\")
    print(f"        --workspace C:\\Work\\eijex --external \\")
    print(f"        --url https://pypi.org/project/factorforge-cds/{new}/")
    print("        → confirms no policy violations leaked into published surfaces")

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
                        help="Path to the internal planning/tracking workspace repo to also bump cross-repo docs")
    parser.add_argument("--mcp", default=None,
                        help="Path to the eijex-mcp repo to also bump MCP tool version strings")
    parser.add_argument("--web", default=None,
                        help="Path to the eijex-web repo to also bump StatsBar version string")
    parser.add_argument("--zenodo-doi", dest="zenodo_doi", default=None,
                        help="Zenodo software version DOI for this release (e.g. 10.5281/zenodo.20640931). "
                             "Run after Zenodo mints the DOI post-tagging.")
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
    web = Path(args.web) if args.web else None
    print(f"Bumping {old} → {new}{' (dry run)' if args.dry_run else ''}{' [strict]' if args.strict else ''}\n")
    errors = bump(old, new, dry_run=args.dry_run, strict=args.strict, workspace=workspace, mcp=mcp, web=web, zenodo_doi=args.zenodo_doi)
    sys.exit(errors)


if __name__ == "__main__":
    main()
