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

Planning docs, job specs, and internal workflow tracking are in the PlantFormOrg private repo (`C:\Work\PlantFormOrg`). Do not duplicate them here.

## 2. Package Name

The Python package name is **`factorforge`** — not `codonforge`.

```python
from factorforge.ml.metrics import calculate_cai
from factorforge.engines.v3.synonym_mask import build_synonym_mask
```

Do not use `codonforge` in any import path, file path, or script.

## 3. Source of Truth

Before editing, read in this order:

1. Current repo code
2. `docs/reports/v3_audit_report.md` — v2/v3 audit findings
3. `docs/design/run4_experiment_plan.md` — Run 4 strategy
4. `docs/design/v3_alpha_product_boundary.md` — claims policy
5. `configs/v3_training_config_alpha_run1.yml` — alpha_run1 training config

## 4. Scope Guardrails

Do not change without an explicit job:
- v2 optimizer production behavior (`src/factorforge/engines/v2/`)
- CodonTokenizer vocab or special token definitions
- Codon table values

Treat as read-mostly unless the active job targets them:
- `experiments/` — results storage only
- `models/` — checkpoints; do not commit `.pt`/`.ckpt`/`.bin` files (see `.gitignore`)

## 5. Claims Policy

Do not add any of the following to code comments, docstrings, or docs:

- "increases yield"
- "guarantees expression"
- "validated expression optimizer"
- "wet-lab proven"
- "production replacement for v2"

v3-alpha is an **in-silico CDS design assistant and evidence layer**.
v2 remains the production fallback until v3 shows consistent advantage on defined benchmarks.

## 6. Do Not Commit

- `models/*.pt`, `models/*.ckpt`, `models/*.bin` — use Google Drive or HuggingFace
- `data/raw/`, `data/embeddings/` — large files, Drive only
- `CLAUDE.md`, `.claude/` — internal workflow files, not for public repo
- `.env`, secrets of any kind

## 7. Completion Report

After completing a job, print to terminal:

```
=== JOB NNN 완료 보고 ===
생성/수정 파일: [목록]
테스트 결과: [결과 또는 N/A]
특이사항: [있으면 기재, 없으면 없음]
========================
```
