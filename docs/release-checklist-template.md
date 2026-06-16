# FactorForge Release Checklist — vX.Y.Z

> Copy this file for each release. Fill in version and date.
> **Template location**: `docs/release-checklist-template.md`

**Release version**: vX.Y.Z
**Date**: YYYY-MM-DD
**Executor**: Maintainer

---

## Before Release (manual — do this during development)

- [ ] `CHANGELOG.md` `## [Unreleased]` section has content describing this release
- [ ] `web/index.html` What's New TODO placeholder is **optional** — `--auto` fills it from CHANGELOG automatically

---

## Release Command (one command, fully automated)

```bash
python scripts/release.py X.Y.Z --auto \
  --audit-script "C:\Users\munky\.codex\skills\factorforge-public-surface-audit\scripts\audit_public_surface.py" \
  --workspace C:\Work\eijex\eijex-workspace \
  --mcp C:\Work\eijex\eijex-mcp \
  --web C:\Work\eijex\eijex-web
```

**`--auto` handles (in order):**

| Step | What it does |
|------|-------------|
| [0] Pre-flight | git clean check + [Unreleased] content check |
| [1] Version bump | All 22 version-bearing files updated |
| [2] Freeze | `examples/worked_example/run_example.py --freeze` |
| [3] Audit | `audit_public_surface.py --live` (abort on critical) |
| [4] Tests | `pytest tests/ -q --tb=short` (abort on failure) |
| [5] CHANGELOG | `[Unreleased]` → `[X.Y.Z] — date` + comparison links |
| [6] What's New | Auto-fills `web/index.html` bullets from CHANGELOG entries |
| [7] Commit | `git add -A && git commit -m "chore: release vX.Y.Z"` |
| [8] Tag | `git tag -a vX.Y.Z` |
| [9] Push | `git push && git push --tags` |

GitHub Actions then auto-publishes: **PyPI + Docker + GitHub Release + Zenodo**

- [ ] `--auto` completed without errors
- [ ] GitHub Actions: PyPI publish ✅
- [ ] GitHub Actions: Docker publish ✅
- [ ] GitHub Actions: GitHub Release created ✅
- [ ] GitHub Actions: Zenodo DOI updated ✅

---

## Post-release (3 items — still manual)

- [ ] `docs/changelog.md` — add one-line summary entry for vX.Y.Z

- [ ] After Zenodo mints the new version DOI:
  ```bash
  python scripts/release.py X.Y.Z --zenodo-doi 10.5281/zenodo.<NEW_RECORD_ID>
  git add CITATION.cff && git commit -m "chore: update Zenodo DOI for vX.Y.Z" && git push
  ```

- [ ] Bioconda (`recipes/meta.yaml` already bumped by `--auto`):
  ```bash
  cd <bioconda-recipes-fork> && git push origin add-factorforge-cds
  ```

---

## Cross-repo commits (eijex-mcp, eijex-web, eijex-workspace)

Files were updated by `--auto`, but these are separate git repos that need their own commits:

```bash
# eijex-mcp
cd C:\Work\eijex\eijex-mcp
git add -A && git commit -m "chore: FactorForge vX.Y.Z" && git push

# eijex-web
cd C:\Work\eijex\eijex-web
git add -A && git commit -m "chore: FactorForge vX.Y.Z" && git push

# eijex-workspace
cd C:\Work\eijex\eijex-workspace
git add -A && git commit -m "chore: FactorForge vX.Y.Z" && git push
```

- [ ] eijex-mcp committed + pushed
- [ ] eijex-web committed + pushed
- [ ] eijex-workspace committed + pushed

---

## Quick verification (spot-check after push)

- [ ] `https://factorforge.eijex.com` — version header shows vX.Y.Z, What's New CURRENT block correct
- [ ] `https://pypi.org/project/factorforge-cds/X.Y.Z/` — version live
- [ ] `https://zenodo.org/doi/10.5281/zenodo.20407330` — concept DOI resolves to new record

---

## Registry/UI Consistency (only if registry parameters changed)

```bash
python -c "
import yaml, pathlib
r = yaml.safe_load(pathlib.Path('src/factorforge/registry/current_parameter_registry.yaml').read_text())
gc = r['parameters']['optimization']['gc_range_nbenthamiana_global']['value']
print(f'Registry GC range: {gc[0]}-{gc[1]}%')
print('Verify matches: docs/index.md / web/index.html / web/js/app.js')
"
```
