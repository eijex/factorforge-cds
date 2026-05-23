# FactorForge Agent Operating Rules

> Implementation rules for any coding agent working in this repo.
> Planning, job tracking, and workflow rules live in the PlantFormOrg private repo.

---

## 1. Repo Role

FactorForge is the **public implementation repo**.

- Source code: `src/factorforge/`
- Scripts: `scripts/`
- Tests: `tests/`
- Docs: `docs/` (design, reports, model cards)
- Configs: `configs/`

Planning docs, job specs, and internal workflow tracking are in a separate private repo. Do not duplicate them here.

## 2. Package Name

The PyPI distribution name is **`factorforge-cds`** — install with `pip install factorforge-cds`.
The Python import name is **`factorforge`** — not `codonforge`.

```python
from factorforge.ml.metrics import calculate_cai
from factorforge.engines.v3.synonym_mask import build_synonym_mask
```

Do not use `codonforge` in any import path, file path, or script.

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
chore: bump version to 3.1.0
```

## 4. Source of Truth

Before editing, read in this order:

1. Current repo code
2. `docs/reports/v3_audit_report.md` — v2/v3 audit findings
3. `docs/design/v3_alpha_product_boundary.md` — claims policy

## 5. Scope Guardrails

Do not change without an explicit job:
- v2 optimizer production behavior (`src/factorforge/engines/v2/`)
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
- "production replacement for v2"

v3-alpha is an **in-silico CDS design assistant and evidence layer**.
v2 remains the production fallback until v3 shows consistent advantage on defined benchmarks.

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
| Version bump | `pyproject.toml`, `CHANGELOG.md`, `src/factorforge/__init__.py` |
| New distribution method | `README.md` (Installation), `CHANGELOG.md` |
| Web UI change | `web/index.html`, `web/README.md` |
| New engine / model | `docs/model_cards/`, `CHANGELOG.md` |
| Wet-lab result added | `VALIDATION.md` |
