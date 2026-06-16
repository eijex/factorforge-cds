# FactorForge Release Checklist — vX.Y.Z

> Copy this file for each release. Fill in version, date, and check off items.
> **Template location**: `docs/release-checklist-template.md`

**Release version**: vX.Y.Z
**Date**: YYYY-MM-DD
**Executor**: Maintainer

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
  --workspace <private-workspace-path> \
  --mcp <mcp-repository-path> \
  --web <website-repository-path>
```

Script handles 16 version-bearing files automatically (pyproject.toml, CITATION.cff, __init__.py files, api/optimize.py, web/, ROADMAP.md, tests/, issue template, recipes/meta.yaml, eijex-mcp, eijex-web StatsBar.tsx).

- [ ] Script ran without errors
- [ ] Residual check: `"No residual X.Y.Z-old strings found"` in output
- [ ] `CITATION.cff` version + date-released correct (spot-check)

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

Run the automated audit first, then verify each public surface manually.

```bash
python scripts/audit_public_surface.py --live --external
```

- [ ] No new findings (false positives documented and accepted)

### 4-A. factorforge.eijex.com (web app)
- [ ] Page title / header shows vX.Y.Z
- [ ] "What's New" modal: vX.Y.Z block is CURRENT (emerald border), previous version demoted
- [ ] Disclaimer text: no guarantee/yield/predict-expression claims
- [ ] Footer links: **Share Wet-lab Results (GitHub)** + **Share Wet-lab Results (Form)** both present
- [ ] Results panel button label: **Share Wet-lab Results (GitHub)**

### 4-B. www.eijex.com (eijex-web)
- [ ] StatsBar: mean CAI / AA identity / CDS benchmarked N= / version all reflect new release
- [ ] Hero copy: scope matches claim boundary (no "expression outcome" language)
- [ ] Footer Resources: **Share Wet-lab Results (Form)** link present

### 4-C. eijex.github.io/factorforge-cds/ (MkDocs docs)
- [ ] Homepage benchmark table (CAI, GC%, AA identity, N=) matches new release
- [ ] Homepage "Share Wet-lab Results" section: Form + GitHub + email all present
- [ ] Validation page submission links: Form + GitHub + email all present
- [ ] Tutorials: no outdated claim language ("codon optimization" → "CDS design" where applicable)

### 4-D. github.com/eijex/factorforge-cds (README + repo)
- [ ] README version badge / citation line reflects vX.Y.Z
- [ ] README wet-lab section: **Share Wet-lab Results (Form)** + **(GitHub)** both present
- [ ] CONTRIBUTING.md wet-lab section: Form + GitHub both present
- [ ] VALIDATION.md submission links: Form + GitHub + email all present
- [ ] CITATION.cff: version + date-released updated

### 4-E. mcp.eijex.com / github.com/eijex/eijex-mcp
- [ ] `mcp-tools.ts`, `route.ts` version strings updated (by release.py)
- [ ] `https://mcp.eijex.com` tool descriptions show vX.Y.Z

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
- [ ] `docker run ghcr.io/eijex/factorforge-cds:X.Y.Z factorforge --help` (Docker smoke)
- [ ] Zenodo software — concept DOI auto-updated: https://zenodo.org/doi/10.5281/zenodo.20407330
- [ ] Zenodo software — new record URL visible: https://zenodo.org/records/20407331 (verify vX.Y.Z)
- [ ] Zenodo benchmark data — record URL: https://zenodo.org/records/20676276 (benchmark_results.csv v2, scoring_contract v1.1)
- [ ] Bioconda PR: push fork branch `add-factorforge-cds`
  ```bash
  cd <bioconda-recipes-fork> && git push origin add-factorforge-cds
  ```
  (recipes/meta.yaml already bumped by release.py — update PR title if needed)

---

## Step 7 — Post-release External Audit (deployed verification)

Repeat Step 4 checks against live deployed URLs (not local files).

- [ ] `https://factorforge.eijex.com` — version header, What's New modal, footer links
- [ ] `https://www.eijex.com` — StatsBar numbers, Hero copy, Footer Form link
- [ ] `https://eijex.github.io/factorforge-cds/` — benchmark table, submission links
- [ ] `https://github.com/eijex/factorforge-cds` — README version, wet-lab links
- [ ] `https://github.com/eijex/eijex-mcp` — version strings in mcp-tools.ts, route.ts
- [ ] `https://mcp.eijex.com` — tool descriptions show vX.Y.Z
- [ ] `https://pypi.org/project/factorforge-cds/X.Y.Z/` — version published, description accurate
- [ ] `https://zenodo.org/records/20407331` — new version record visible
- [ ] Automated audit: no new findings
  ```bash
  python scripts/audit_public_surface.py --live --external
  ```

---

## Step 8 — Cross-repo Commits

- [ ] eijex-mcp: commit + push (`mcp-tools.ts`, `route.ts`)
- [ ] eijex-web: commit + push (`StatsBar.tsx`)
- [ ] Private release tracking: mark the release task complete and commit/push the tracking update

---

## Step 9 — Cleanup

- [ ] Close completed GitHub Issues / milestone

---

## Step 10 — Registry/UI Consistency Check

Prevent parameter drift: verify that hardcoded values in docs and UI match the registry source of truth.

```bash
python -c "
import yaml, pathlib
r = yaml.safe_load(pathlib.Path('src/factorforge/registry/current_parameter_registry.yaml').read_text())
gc = r['parameters']['optimization']['gc_range_nbenthamiana_global']['value']
print(f'Registry GC range: {gc[0]}–{gc[1]}%')
print('Verify this matches:')
print('  docs/index.md — GC in target range Target column')
print('  web/index.html — Target Zone label (line ~593)')
print('  web/js/app.js — Target Max/Min chart lines')
"
```

- [ ] `docs/index.md` benchmark table GC target = registry value
- [ ] `web/index.html` "Target Zone X–Y%" = registry value
- [ ] `web/js/app.js` `Array(n).fill(X)` Target Max/Min = registry value
