# FactorForge Agent Operating Rules

> Implementation rules for any coding agent working in this repo.

---

## 1. Repo Role

FactorForge is the **public implementation repo**.

- Source code: `src/factorforge/`
- Scripts: `scripts/`
- Tests: `tests/`
- Docs: `docs/` (public user documentation)
- Configs: `configs/`

## 2. Package Name

The PyPI distribution name is **`factorforge-cds`** — install with `pip install factorforge-cds`.
The Python import name is **`factorforge`**.

```python
from factorforge.engines.profile import RuleBasedOptimizer
from factorforge.engines.profile.pipeline import OptimizationPipeline
```

## 3. Branch & Commit Conventions

Branch naming:
```
feat/short-description
fix/short-description
docs/short-description
chore/short-description
```

Commit message format:
```
feat: add Wolffia globosa codon table
fix: correct GC% calculation for short sequences
docs: update README installation section
chore: bump version to 3.1.1
```

## 4. Source of Truth

Before editing, read in this order:

1. Current repo code
2. public user documentation under `docs/`
3. `VALIDATION.md` and `docs/validation.md` — claims policy

## 5. Scope Guardrails

Do not change without an explicit job:
- production optimizer behavior (`src/factorforge/engines/profile/`)
- CodonTokenizer vocab or special token definitions
- Codon table values

Treat as read-mostly unless the active job targets them:
- `experiments/` — results storage only
- `models/` — checkpoints; do not commit `.pt`/`.ckpt`/`.bin` files (see `.gitignore`)

## 6. Claims Policy

Do not add any of the following to code comments, docstrings, or docs:

- "increases yield"
- "guarantees expression"
- "validated expression optimizer"
- "wet-lab proven"
- "validated expression optimizer"

FactorForge is an **in-silico CDS design assistant**. Public claims must stay limited to
deterministic CDS design, sequence metrics, validation checks, and reproducible output
paths unless wet-lab evidence is explicitly documented.

## 7. Code Style & Testing

Before submitting a PR:

```bash
ruff check .      # linting
ruff format .     # formatting
pytest tests/ -v  # all tests must pass
```

New features must include tests. Bug fixes should add a regression test.

## 8. Security

Do **not** open a public GitHub issue for security vulnerabilities.
Report privately to: eijex.lab@gmail.com

## 9. Do Not Commit

- `models/*.pt`, `models/*.ckpt`, `models/*.bin` — use Google Drive or HuggingFace
- `data/raw/`, `data/embeddings/` — large files, Drive only
- `CLAUDE.md`, `.claude/` — internal workflow files, not for public repo
- `.env`, secrets of any kind

## 10. Completion Report

After completing a job, print to terminal:

```
=== JOB NNN COMPLETE ===
Files created/modified: [list]
Test results: [results or N/A]
Notes: [if any, otherwise none]
========================
```

## 11. Public Document Checklist

After any change, update the relevant public-facing files before pushing:

| Change type | Files to update |
|-------------|----------------|
| New feature / API change | `README.md`, `CHANGELOG.md` |
| Bug fix | `CHANGELOG.md` |
| Version bump | Run `python scripts/release.py X.Y.Z` — updates all 16 version-bearing files automatically. Then manually update `web/index.html` changelog panel and `docs/changelog.md`. **Then run Step 4 public surface audit** (see `docs/release-checklist-template.md`). |
| Wet-lab submission link change | Update all 8 surfaces: `web/index.html`, `web/js/app.js`, `README.md`, `CONTRIBUTING.md`, `VALIDATION.md`, `docs/index.md`, `docs/validation.md`, `docs/feedback-inbox.md`. Label convention: **Share Wet-lab Results (GitHub)** / **Share Wet-lab Results (Form)**. |
| New distribution method | `README.md` (Installation), `CHANGELOG.md` |
| Web UI change | `web/index.html`, `web/README.md` — also run: `npx playwright test` |
| New engine / model | `docs/model_cards/`, `CHANGELOG.md` |
| Wet-lab result added | `VALIDATION.md` |
| New API endpoint | `docs/cli.md` (API section), `CHANGELOG.md`, **eijex-mcp sync** (see `CHANGELOG.md` Step 14) |

## 12. Public Repo Structure Checklist

For every audit, release, or repo-structure change, check the repository as a
user- and paper-facing artifact:

- `README.md`, `docs/`, and `web/` should describe only supported public
  runtime paths first.
- `src/factorforge/` should contain installable, supported package code; historical
  engines and experiments should stay under `archive/`.
- `archive/` is public provenance, not an active import surface.
- Claims must stay aligned with `VALIDATION.md` and `docs/validation.md`.

## 13. Public Surface Coverage (all releases)

Every release must pass Step 4 of `docs/release-checklist-template.md` before tagging.
The 8 public surfaces and their required checks:

| Surface | URL | Key checks |
|---------|-----|-----------|
| Web app | https://factorforge.eijex.com | version, What's New modal, wet-lab links, disclaimer |
| Corporate site | https://www.eijex.com | StatsBar numbers, Hero copy, Footer Form link |
| Docs | https://eijex.github.io/factorforge-cds/ | benchmark table, submission links |
| GitHub repo | https://github.com/eijex/factorforge-cds | README version, wet-lab links |
| MCP GitHub | https://github.com/eijex/eijex-mcp | version strings in mcp-tools.ts, route.ts |
| MCP site | https://mcp.eijex.com | tool descriptions version |
| PyPI | https://pypi.org/project/factorforge-cds/ | version, description |
| Zenodo | https://zenodo.org/records/20407330 | new version record visible |
- Run or report proportional checks for CLI, API, package build metadata, and
  untracked-file hygiene before push.
