# FactorForge Release Checklist — vX.Y.Z

> Copy this file for each release. Fill in version, date, and check off items.
> **Template location**: `docs/release-checklist-template.md`

**Release version**: vX.Y.Z
**Date**: YYYY-MM-DD
**Executor**: Claude / user

---

## Pre-release Gate

- [ ] `git status --short` — working tree clean
- [ ] `python -m ruff check .` — no lint errors
- [ ] `python -m pytest tests/ -v --tb=short` — all pass (Ubuntu + Windows)
- [ ] `[Unreleased]` section in `CHANGELOG.md` has content to release

---

## Step 1 — Version Bump (automated)

```bash
python scripts/release.py X.Y.Z \
  --workspace C:\Work\eijex\eijex-workspace \
  --mcp C:\Work\eijex\eijex-mcp \
  --web C:\Work\eijex\eijex-web
```

Expected files updated automatically:
- [ ] `pyproject.toml`
- [ ] `CITATION.cff` (version + date-released)
- [ ] `src/factorforge/__init__.py`
- [ ] `src/factorforge/engines/__init__.py`
- [ ] `src/factorforge/engines/profile/__init__.py`
- [ ] `src/factorforge/engines/profile/optimizer.py`
- [ ] `api/optimize.py`
- [ ] `web/index.html` (Release Notes heading + changelog CURRENT version/date)
- [ ] `web/js/app.js`
- [ ] `ROADMAP.md`
- [ ] `tests/` (version assertion strings)
- [ ] `.github/ISSUE_TEMPLATE/wet_lab_result.yml`
- [ ] `recipes/meta.yaml` (version + PyPI sha256 auto-fetch)
- [ ] eijex-mcp: `mcp-tools.ts`, `route.ts`
- [ ] eijex-web: `StatsBar.tsx`

Residual check: confirm `"No residual X.Y.Z-old strings found"` in output.

---

## Step 2 — Manual Changelog & Docs

- [ ] `CHANGELOG.md`: move `[Unreleased]` → `[X.Y.Z] — YYYY-MM-DD`; update comparison links
- [ ] `web/index.html`: add new `vX.Y.Z` CURRENT block with bullet points (above old CURRENT)
  - New block: `border-emerald-500`, dot `bg-emerald-500`, **Current** badge
  - Old block: `border-slate-200`, dot `bg-slate-300`, no badge
- [ ] `docs/changelog.md`: add summary entry

---

## Step 3 — Commit & CI

- [ ] `git add -A && git commit -m "chore: release vX.Y.Z"`
- [ ] `git push`
- [ ] CI passes — **all 8 matrix jobs green** (4 Python × Ubuntu + Windows)

---

## Step 4 — Public Surface Audit (pre-tag)

```bash
python ~/.codex/skills/factorforge-public-surface-audit/scripts/audit_public_surface.py \
  --workspace C:\Work\eijex --live
```

- [ ] No new findings (false positives documented and accepted)

---

## Step 5 — Tag & Publish

```bash
git tag -a vX.Y.Z -m "Release vX.Y.Z" && git push --tags
```

- [ ] Tag pushed
- [ ] GitHub Actions: PyPI publish ✅
- [ ] GitHub Actions: Docker publish ✅
- [ ] GitHub Actions: GitHub Release created ✅
- [ ] GitHub Actions: Zenodo DOI updated ✅

---

## Step 6 — Post-release Verification

- [ ] `pip install factorforge-cds==X.Y.Z && factorforge --help` (PyPI smoke)
- [ ] `docker run ghcr.io/eijex/factorforge-cds:vX.Y.Z factorforge --help` (Docker smoke)
- [ ] Zenodo DOI: https://zenodo.org/doi/10.5281/zenodo.20407331
- [ ] Bioconda PR: push fork branch `add-factorforge-cds`
  ```bash
  cd <bioconda-recipes-fork> && git push origin add-factorforge-cds
  ```
  (recipes/meta.yaml already bumped by release.py — update PR title if needed)

---

## Step 7 — Post-release External Audit

```bash
python ~/.codex/skills/factorforge-public-surface-audit/scripts/audit_public_surface.py \
  --workspace C:\Work\eijex --live --external \
  --url https://pypi.org/project/factorforge-cds/X.Y.Z/
```

- [ ] No new findings

---

## Step 8 — Cross-repo Commits

- [ ] eijex-mcp: commit + push (`mcp-tools.ts`, `route.ts`)
- [ ] eijex-web: commit + push (`StatsBar.tsx`)
- [ ] eijex-workspace: Done 루틴 (job Done + ROADMAP + done-job.py + commit + push)

---

## Step 9 — Cleanup

- [ ] Close completed GitHub Issues / milestone
- [ ] Verify factorforge.eijex.com shows new CURRENT version in "What's New" modal
- [ ] Verify mcp.eijex.com tool descriptions show new version
