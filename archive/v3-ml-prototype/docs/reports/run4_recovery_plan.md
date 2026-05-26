# Run 4 Recovery Plan

> Date: 2026-05-20  
> Scope: Run 4B recovery design for v3-alpha smoke testing  
> Boundary: in-silico sequence design only; no expression, yield, wet-lab, or production-replacement claims

## Why Run 4 Failed GC

Run 4 completed training and preserved the hard safety properties:

- validator pass rate was 100% in the quick constrained eval;
- amino acid identity was 100%;
- no internal stop codons or invalid codons were observed;
- synonym mask coverage was 1.0;
- v3 CAI stayed effectively identical to the v2 teacher.

The failure was the GC objective. The v2 `high_cai` teacher distribution had low global GC in the analyzed eval subset, with v2 mean GC around 33.23%. Run 4 learned this teacher closely, so v3 decoded mean GC stayed around 33.24%, below the required 40% lower bound.

The loss log showed the bounded GC term was active but too weak. At step 20,000, expected GC averaged about 0.338 and bounded GC penalty remained positive. With `gc_weight=0.2`, the penalty contribution was too small relative to CE-to-pseudo-label to move the model into the 0.40-0.55 target band.

## v2 Profile Availability

The actual v2 optimizer profiles are:

- `balanced`
- `high_cai`
- `gc_target`
- `assembly_friendly`
- `ramp`
- `viral_delivery`

The requested `local_gc_balanced` and `motif_clean` profiles do not exist in the v2 optimizer. Run 4B therefore uses `gc_target` when available and falls back to exact feasibility analysis for a 40-55% GC candidate when a v2 candidate is unavailable or not GC-feasible.

## Why Mixed Teacher Is Needed

Pure `high_cai` pseudo-labels teach the model to reproduce a low-GC teacher. A stronger GC loss helps, but the training target should also expose valid synonymous choices that satisfy the GC band when feasible.

Run 4B generates multiple candidates per protein:

1. v2 `high_cai` candidate.
2. v2 `gc_target` candidate.
3. Feasibility-derived 40-55% GC candidate if needed.

Each candidate is validated, scored, and logged. The selected teacher uses a simple multi-objective score:

```text
CAI
+ GC bonus if global GC is in 40-55%
- forbidden motif penalty
- repeat/homopolymer penalty
```

This keeps the teacher expression-oriented while preventing the model from inheriting low-GC `high_cai` labels when a better GC-feasible synonymous sequence exists.

## Smoke Grid

Run 4B is smoke-only. Full training remains out of scope until a smoke configuration passes the acceptance gate.

| Config | gc_weight | gc_lambda_low | expected_log_cai_weight | max_steps |
|---|---:|---:|---:|---:|
| A | 1.0 | 1.0 | 0.05 | 2,000 |
| B | 2.0 | 1.0 | 0.05 | 2,000 |
| C | 2.0 | 2.0 | 0.05 | 2,000 |
| D | 2.0 | 2.0 | 0.10 | 2,000 |

The grid isolates two recovery levers:

- stronger global GC objective (`gc_weight`);
- stronger low-GC correction (`gc_lambda_low`);
- modestly stronger CAI preservation (`expected_log_cai_weight=0.10`) for the highest-GC-pressure setting.

## Full Training Gate

Do not run full Run 4B training unless a smoke run passes:

- validator pass rate = 100%;
- amino acid identity = 100%;
- decoded global GC mean >= 40%;
- decoded global GC max <= 60%;
- expected_GC >= 0.40;
- v3/v2 CAI ratio mean >= 0.97.

If no smoke configuration passes, revise the teacher selection strategy before increasing training scale.

## Outputs

Run 4B smoke outputs should be saved per config:

- `experiments/results/run4b_smoke/config_A/`
- `experiments/results/run4b_smoke/config_B/`
- `experiments/results/run4b_smoke/config_C/`
- `experiments/results/run4b_smoke/config_D/`

Each directory should include:

- `run4_loss_log.csv`
- `run4_comparison.csv`
- `run4_summary.json`

Use `scripts/4_evaluation/run4b_smoke_compare.py` to summarize pass/fail criteria across configs.
