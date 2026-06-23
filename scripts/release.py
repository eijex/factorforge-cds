#!/usr/bin/env python3
"""
Bump FactorForge version across all version-bearing files.

Usage:
    python release.py 3.2.2                    # version bump only (manual release)
    python release.py 3.2.2 --auto             # full automated release pipeline
    python release.py 3.2.2 --auto --dry-run   # preview without writing
    python release.py 3.2.2 --dry-run          # preview bump only
    python release.py 3.2.2 --strict           # fail on any missing pattern
    python release.py 3.2.2 --zenodo-doi 10.5281/zenodo.XYZ  # post-tag: update DOI

--auto pipeline (fully automated):
    [0] Pre-flight  — clean git tree, [Unreleased] has content
    [1] Version bump — all version-bearing files
    [2] Freeze      — examples/worked_example/run_example.py --freeze
    [3] Audit       — scripts/audit_public_surface.py --live (if --audit-script given)
    [4] Tests       — pytest tests/ -q
    [5] CHANGELOG   — [Unreleased] → [X.Y.Z] — date + comparison links
    [6] What's New  — auto-fill web/index.html bullets from CHANGELOG entries
    [7] Commit      — git add -A && git commit
    [8] Tag         — git tag -a vX.Y.Z
    [9] Push        — git push && git push --tags
    → GitHub Actions publishes PyPI + Docker + GitHub Release + Zenodo

Publication sync (keep documentation in sync with feature changes):
  - Patch release: no update needed unless claims changed
  - Feature release (minor): review Software Design documentation
  - Algorithm added (tAI, codon pair bias, 5' UTR MFE): update State of the Field documentation
  - New host added: update Statement of Need + supported hosts documentation
  - Zenodo DOI updated: update citation metadata
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

from packaging.version import Version

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

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
    """Demote old CURRENT block to history and insert new CURRENT placeholder block."""
    path = root / "web/index.html"
    if not path.exists():
        return []
    today = _today()
    content = path.read_text(encoding="utf-8")
    original = content
    changes = []

    # 1. Navbar button
    if f"v{old} Release Notes" in content:
        content = content.replace(f"v{old} Release Notes", f"v{new} Release Notes")
        changes.append(f"  navbar button: v{old} → v{new}")

    # 2. Demote old CURRENT block: emerald → slate (replace first occurrence only)
    emerald_outer = 'border-l-2 border-emerald-500">'
    slate_outer   = 'border-l-2 border-slate-200 dark:border-slate-700">'
    if emerald_outer in content:
        content = content.replace(emerald_outer, slate_outer, 1)
        changes.append("  demoted CURRENT outer border: emerald → slate")

    emerald_dot = 'bg-emerald-500 rounded-full border-4 border-white dark:border-slate-900 shadow-[0_0_10px_rgba(16,185,129,0.5)]">'
    slate_dot   = 'bg-slate-300 dark:bg-slate-600 rounded-full border-4 border-white dark:border-slate-900">'
    if emerald_dot in content:
        content = content.replace(emerald_dot, slate_dot, 1)
        changes.append("  demoted CURRENT dot: emerald → slate")

    badge = re.compile(
        r'\s*<span\s[^>]*bg-emerald-100[^>]*>Current</span>\s*'
    )
    content, n = badge.subn(lambda m: "\n                        ", content, count=1)
    if n:
        changes.append("  removed Current badge")

    # 3. Insert new CURRENT block before <!-- Version {old} -->
    marker = f"<!-- Version {old} -->"
    if marker in content:
        new_block = (
            f"<!-- Version {new} -->\n"
            f'                <div class="relative pl-8 border-l-2 border-emerald-500">\n'
            f'                    <div\n'
            f'                        class="absolute -left-[9px] top-0 w-4 h-4 bg-emerald-500 rounded-full border-4 border-white dark:border-slate-900 shadow-[0_0_10px_rgba(16,185,129,0.5)]">\n'
            f'                    </div>\n'
            f'                    <div class="flex items-center space-x-2 mb-2">\n'
            f'                        <span\n'
            f'                            class="px-2 py-0.5 bg-emerald-100 dark:bg-emerald-900/30 text-emerald-700 dark:text-emerald-400 text-[10px] font-bold rounded-md uppercase">Current</span>\n'
            f'                        <h3 class="font-bold text-slate-800 dark:text-white">v{new}</h3>\n'
            f'                        <span class="text-slate-400 text-[10px]">{today}</span>\n'
            f'                    </div>\n'
            f'                    <ul class="text-xs text-slate-600 dark:text-slate-400 space-y-2 list-disc ml-4">\n'
            f'                        <!-- TODO: add v{new} release notes here -->\n'
            f'                    </ul>\n'
            f'                </div>\n'
            f'                '
        )
        content = content.replace(marker, new_block + marker, 1)
        changes.append(f"  inserted new CURRENT block for v{new} (fill in release notes)")

    if content == original:
        return []
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
        ("docs/benchmark.md", [
            (f"FactorForge v{old} provides", f"FactorForge v{new} provides"),
        ], True),
        ("docs/validation/RELEASE_GATE.md", [
            (f"# FactorForge v{old} Release Gate", f"# FactorForge v{new} Release Gate"),
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
        ("src/factorforge/registry/current_parameter_registry.yaml", [
            (f'version: "{old}"', f'version: "{new}"'),
        ], True),
        ("tests/conftest.py", [
            (f'"registry_version": "{old}"', f'"registry_version": "{new}"'),
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
    skip_dirs = {'.git', '__pycache__', '.pytest_cache', 'build', 'dist', 'archive', '.venv', 'venv', 'node_modules'}
    # skip_names: matched against path.name (filename only)
    skip_names = {'CHANGELOG.md', 'changelog.md'}
    # skip_rel: matched against path relative to ROOT (posix form)
    skip_rel = {'benchmarks/results', 'reproducibility/'}

    for path in ROOT.rglob('*'):
        if path.is_dir():
            continue
        if any(p in path.parts for p in skip_dirs):
            continue
        if path.name in skip_names:
            continue
        rel_posix = path.relative_to(ROOT).as_posix()
        if any(rel_posix.startswith(s) for s in skip_rel):
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

        # NOTE: recipes/meta.yaml's sha256 is intentionally NOT touched here.
        # PyPI does not have the new version yet at this point in the pipeline,
        # so a fetch here always fails (see scripts/update_conda_recipe_sha.py,
        # which updates the sha256 only after PyPI confirms the version is live).

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
    print("  --- Regenerate frozen example (required every release) ---")
    print("  3a. python examples/worked_example/run_example.py --freeze")
    print("      → updates evidence.registry_version and evidence.registry_hash in design_package.json")
    print("      → git add examples/worked_example/output/")
    print()
    print("  --- Commit & CI gate (before tagging) ---")
    print(f"  4. git commit -m 'chore: release v{new}'")
    print("  5. git push")
    print("  6. Wait for CI to pass (github.com/eijex/factorforge-cds/actions)")
    print()
    print("  --- Public surface audit (before tagging) ---")
    print(f"  6a. python scripts/audit_public_surface.py --live")
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
    print(f" 12.  python scripts/update_conda_recipe_sha.py {new}  (only after PyPI confirms the sdist is live)")
    print("       → commits on a dedicated branch, then push it: git push origin update-conda-recipe-sha-v" + new)
    print("       → open/update PR on bioconda/bioconda-recipes")
    print(" 13.  Close completed GitHub Issues; close milestone if all done")
    print()
    print("  --- Post-release external audit ---")
    print(f" 14.  python scripts/audit_public_surface.py --external \\")
    print(f"        --url https://pypi.org/project/factorforge-cds/{new}/")
    print("        → confirms no policy violations leaked into published surfaces")

    return errors


# ─────────────────────────────────────────────────────────────────────────────
# Auto-release helpers
# ─────────────────────────────────────────────────────────────────────────────

def _run_step(label: str, cmd: list[str]) -> None:
    """Print step header, run command, abort on non-zero exit."""
    print(f"  $ {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, cwd=str(ROOT))
    if result.returncode != 0:
        print(f"  ABORT: '{label}' failed (exit {result.returncode})")
        sys.exit(result.returncode)


def _changelog_move_unreleased(root: Path, old: str, new: str, dry_run: bool) -> list[str]:
    """
    In CHANGELOG.md:
    - Rename '## [Unreleased]' → '## [new] — YYYY-MM-DD'
    - Update [Unreleased] comparison link: v{old}...HEAD → v{new}...HEAD
    - Insert new [new]: v{old}...v{new} comparison link
    """
    path = root / "CHANGELOG.md"
    content = path.read_text(encoding="utf-8")
    today = _today()
    changes: list[str] = []

    # 1. Rename [Unreleased] header (bracket format only — leaves '## Unreleased' intact)
    old_header = "## [Unreleased]"
    new_header = f"## [{new}] — {today}"
    if old_header in content:
        content = content.replace(old_header, new_header, 1)
        changes.append(f"  {old_header} -> {new_header}")

    # 2. Update comparison links
    #    [Unreleased]: https://.../compare/vX.Y.Z...HEAD
    link_pat = re.compile(
        r'(\[Unreleased\]: https://[^\s]+/compare/v)[^\s.]+(\.[^\s.]+\.[^\s]+)(\.\.\.HEAD)'
    )
    m = link_pat.search(content)
    if m:
        base = m.group(1)
        suffix = m.group(3)
        old_link = m.group(0)
        new_unreleased = f"[Unreleased]: {base.split('[Unreleased]: ')[1]}{new}{suffix}"
        version_link = f"[{new}]: {base.split('[Unreleased]: ')[1]}{old}...v{new}"
        # reconstruct properly using the full match
        full_base_url = re.search(r'https://[^\s]+/compare/', old_link).group(0)
        new_unreleased = f"[Unreleased]: {full_base_url}v{new}...HEAD"
        version_link = f"[{new}]: {full_base_url}v{old}...v{new}"
        content = content.replace(old_link, f"{new_unreleased}\n{version_link}", 1)
        changes.append(f"  comparison links updated")

    if not changes:
        return []
    if not dry_run:
        path.write_text(content, encoding="utf-8")
    return changes


def _parse_changelog_section(root: Path, version: str) -> list[tuple[str, str]]:
    """
    Parse CHANGELOG.md for [version] section.
    Returns list of (category, item_text) for top-level '- ' bullets only.
    """
    path = root / "CHANGELOG.md"
    content = path.read_text(encoding="utf-8")

    pat = re.compile(
        rf'^## \[{re.escape(version)}\][^\n]*\n(.*?)(?=^## \[|\Z)',
        re.MULTILINE | re.DOTALL,
    )
    m = pat.search(content)
    if not m:
        return []

    entries: list[tuple[str, str]] = []
    category = "Changed"
    for line in m.group(1).splitlines():
        stripped = line.strip()
        if stripped.startswith("### "):
            category = stripped[4:].strip()
        elif stripped.startswith("- ") or stripped.startswith("* "):
            text = stripped[2:].strip()
            # Collapse bold markdown: **Foo** → Foo
            text = re.sub(r'\*\*([^*]+)\*\*', r'\1', text)
            # Trim at first colon if the lead looks like a label "Label — detail"
            if " — " in text:
                text = text.split(" — ")[0].strip()
            entries.append((category, text))
    return entries


def _update_web_whats_new_bullets(
    root: Path, new: str, entries: list[tuple[str, str]], dry_run: bool
) -> list[str]:
    """Fill in the <!-- TODO: add v{new} release notes here --> placeholder."""
    path = root / "web/index.html"
    if not path.exists():
        return []
    content = path.read_text(encoding="utf-8")
    placeholder = f"<!-- TODO: add v{new} release notes here -->"
    if placeholder not in content:
        return []

    indent = "                        "
    if entries:
        li_lines = [f"{indent}<li><b>{cat}:</b> {text}</li>" for cat, text in entries]
        bullets = "\n".join(li_lines)
    else:
        bullets = f"{indent}<li>See CHANGELOG.md for details.</li>"

    updated = content.replace(placeholder, bullets)
    if content == updated:
        return []
    if not dry_run:
        path.write_text(updated, encoding="utf-8")
    return [f"  web/index.html: inserted {len(entries)} bullet(s) for v{new}"]


def auto_release(
    old: str,
    new: str,
    workspace: "Path | None",
    mcp: "Path | None",
    web_repo: "Path | None",
    audit_script: "Path | None" = None,
    skip_tests: bool = False,
    dry_run: bool = False,
) -> int:
    """Full automated release pipeline. Returns 0 on success."""
    bar = "=" * 56
    print(f"\n{bar}")
    print(f"  AUTO RELEASE  {old} -> {new}{' [DRY RUN]' if dry_run else ''}")
    print(f"{bar}")

    # ── [0] PRE-FLIGHT ──────────────────────────────────────
    print("\n[0] Pre-flight checks")

    status = subprocess.run(
        ["git", "status", "--porcelain"], capture_output=True, text=True, cwd=str(ROOT)
    )
    if status.stdout.strip():
        print(f"  ABORT: working tree is not clean -- commit or stash first\n{status.stdout}")
        sys.exit(1)
    print("  git status: clean")

    cl_text = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    unrel = re.search(r"## \[Unreleased\](.*?)(?=\n## \[|\Z)", cl_text, re.DOTALL)
    if not unrel or not unrel.group(1).strip():
        print("  ABORT: CHANGELOG.md [Unreleased] is empty — write release notes first")
        sys.exit(1)
    print("  CHANGELOG [Unreleased]: has content")

    # ── [1]-[6]: BUMP, FREEZE, AUDIT, TEST, CHANGELOG, WHAT'S NEW ───────────
    # Wrapped together: steps [1]/[1b]/[2]/[5]/[6] write files, and any later
    # failure in this pre-commit window must leave the working tree exactly as
    # it was before this function started (verified by `git status --porcelain`
    # in the AC1 dry-run/failing-test scenario).
    try:
        # ── [1] VERSION BUMP ────────────────────────────────
        print(f"\n[1] Version bump: {old} -> {new}")
        if not dry_run:
            errors = bump(
                old, new, dry_run=False, strict=False,
                workspace=workspace, mcp=mcp, web=web_repo,
            )
            if errors:
                print(f"  ABORT: {errors} bump error(s)")
                sys.exit(1)
        else:
            print("  (skipped — dry run)")

        # ── [1b] REGENERATE MANIFEST HASHES ─────────────────
        # The version bump above rewrites `current_parameter_registry.yaml`'s
        # `version:` field, which is one of MANIFEST.json's hashed inputs — every
        # release otherwise drifts this hash (this recurred across Jobs 129, 136,
        # and again during the v3.2.3 release itself).
        print("\n[1b] Regenerate MANIFEST.json provenance hashes")
        if not dry_run:
            _run_step("regen_manifest", ["python", "scripts/regen_manifest.py", "--write"])
        else:
            print("  (skipped — dry run)")

        # ── [2] FREEZE EXAMPLE ──────────────────────────────
        print("\n[2] Regenerate frozen example")
        if not dry_run:
            _run_step("freeze", ["python", "examples/worked_example/run_example.py", "--freeze"])
        else:
            print("  (skipped — dry run)")

        # ── [3] PUBLIC SURFACE AUDIT ─────────────────────────
        print("\n[3] Public surface audit")
        if audit_script and audit_script.exists():
            if not dry_run:
                r = subprocess.run(["python", str(audit_script), "--live"], cwd=str(ROOT))
                if r.returncode != 0:
                    print("  ABORT: audit found critical issues — fix before releasing")
                    sys.exit(1)
                print("  audit: passed")
            else:
                print(f"  (dry run) would run: python {audit_script} --live")
        else:
            print("  SKIP: --audit-script not provided")
            print("  Reminder: run audit manually before tagging")

        # ── [4] TESTS ─────────────────────────────────────────
        if not skip_tests:
            print("\n[4] Test suite")
            if not dry_run:
                _run_step("pytest", ["python", "-m", "pytest", "tests/", "-q", "--tb=short"])
            else:
                print("  (skipped — dry run)")
        else:
            print("\n[4] Tests — skipped (--skip-tests)")

        # ── [5] CHANGELOG ─────────────────────────────────────
        print("\n[5] CHANGELOG.md: [Unreleased] -> [" + new + "]")
        if not dry_run:
            cl_changes = _changelog_move_unreleased(ROOT, old, new, dry_run=False)
            for c in cl_changes:
                print(c)
        else:
            print(f"  (dry run) would rename and update comparison links")

        # ── [6] WHAT'S NEW AUTO-FILL ──────────────────────────
        print("\n[6] What's New bullets (web/index.html)")
        if not dry_run:
            entries = _parse_changelog_section(ROOT, new)
            wn = _update_web_whats_new_bullets(ROOT, new, entries, dry_run=False)
            if wn:
                for c in wn:
                    print(c)
            else:
                print("  NOTE: TODO placeholder not found — bullets may already be written")
        else:
            print("  (dry run) would auto-fill from CHANGELOG entries")
    except SystemExit as exc:
        if not dry_run:
            print("\n  Cleaning up working tree after failure (git checkout -- .)")
            subprocess.run(["git", "checkout", "--", "."], cwd=str(ROOT))
        raise exc

    # ── [7] COMMIT ──────────────────────────────────────────
    print("\n[7] Git commit")
    if not dry_run:
        _run_step("git add", ["git", "add", "-A"])
        _run_step("git commit", ["git", "commit", "-m", f"chore: release v{new}"])
    else:
        print(f"  (dry run) git add -A && git commit -m 'chore: release v{new}'")

    # ── [8] TAG ─────────────────────────────────────────────
    print("\n[8] Git tag")
    if not dry_run:
        _run_step("git tag", ["git", "tag", "-a", f"v{new}", "-m", f"Release v{new}"])
    else:
        print(f"  (dry run) git tag -a v{new} -m 'Release v{new}'")

    # ── [9] PUSH ────────────────────────────────────────────
    print("\n[9] Push")
    if not dry_run:
        _run_step("git push", ["git", "push"])
        _run_step("git push tags", ["git", "push", "--tags"])
    else:
        print("  (dry run) git push && git push --tags")

    # ── DONE ────────────────────────────────────────────────
    print(f"\n{bar}")
    print(f"  Release v{new} {'[DRY RUN — nothing written]' if dry_run else 'COMPLETE'}")
    print(f"{bar}")

    if not dry_run:
        print("\nGitHub Actions (automatic): PyPI + Docker + GitHub Release + Zenodo")
        print("\nPost-release (3 manual items):")
        print(f"  1. docs/changelog.md — add a one-line summary entry for v{new}")
        print(f"  2. After Zenodo mints the new DOI:")
        print(f"       python scripts/release.py {new} --zenodo-doi 10.5281/zenodo.<ID>")
        print(f"  3. Bioconda: python scripts/update_conda_recipe_sha.py {new}  (after PyPI confirms the sdist is live), then push the branch it creates")
        cross = [r for r in [workspace, mcp, web_repo] if r]
        if cross:
            print("\nCross-repo commits (bump already applied, need separate push):")
            for r in cross:
                print(f"    cd {r}")
                print(f"    git add -A && git commit -m 'chore: FactorForge v{new}' && git push")

    return 0


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
    parser.add_argument("--auto", action="store_true",
                        help="Full automated release pipeline: bump + freeze + audit + test + commit + tag + push")
    parser.add_argument("--audit-script", dest="audit_script", default=None,
                        help="Path to audit_public_surface.py (optional; skipped if not provided)")
    parser.add_argument("--skip-tests", action="store_true",
                        help="Skip pytest in --auto mode (use only when CI already confirmed green)")
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

    if Version(new) < Version(old):
        print(f"ERROR: new version {new} is older than current version {old} — refusing to downgrade")
        sys.exit(1)

    workspace = Path(args.workspace) if args.workspace else None
    mcp = Path(args.mcp) if args.mcp else None
    web = Path(args.web) if args.web else None
    audit = Path(args.audit_script) if args.audit_script else None

    if args.auto:
        sys.exit(auto_release(
            old, new,
            workspace=workspace, mcp=mcp, web_repo=web,
            audit_script=audit,
            skip_tests=args.skip_tests,
            dry_run=args.dry_run,
        ))

    print(f"Bumping {old} -> {new}{' (dry run)' if args.dry_run else ''}{' [strict]' if args.strict else ''}\n")
    errors = bump(old, new, dry_run=args.dry_run, strict=args.strict, workspace=workspace, mcp=mcp, web=web, zenodo_doi=args.zenodo_doi)
    sys.exit(errors)


if __name__ == "__main__":
    main()
