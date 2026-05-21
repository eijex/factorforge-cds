# Run 4B Full Training Analysis

> Date: 2026-05-21  
> Scope: FactorForge v3-alpha Run 4B full eval, Config B  
> Config: `configs/v3_training_config_alpha_run2.yml`  
> Boundary: in-silico sequence design only; no expression, yield, wet-lab, or production-replacement claims

## Inputs

Run 4B full evaluation artifacts:

- `experiments/results/alpha_run2/alpha_run2_summary.json`
- `experiments/results/alpha_run2/alpha_run2_comparison.csv`
- `experiments/results/alpha_run2/alpha_run2_loss_log.csv`

Baseline references:

- Run 4 failed GC target at approximately 33.24% decoded global GC.
- Run 4B smoke Config B passed at approximately 42.52% decoded global GC on 200 eval samples.

## Summary

Run 4B Config B passes the full eval gate on all hard safety criteria and the primary optimization criteria.

| Metric | Run 4 (026) | Run 4B (028) | v2 baseline in Run 4B eval |
|---|---:|---:|---:|
| eval sequences | 100-sample quick eval | 3,876 | 3,876 |
| decoded global GC mean | ~33.24% | 42.5440% | 42.5440% |
| decoded global GC min | 23.17% | 41.9948% | 40.3646% |
| decoded global GC max | 38.29% | 53.8077% | 53.8077% |
| v3/v2 CAI ratio mean | ~1.0000 | 1.0000 | N/A |
| v3/v2 CAI ratio min | 0.9971 | 0.9682 | N/A |
| validator pass rate | 100% | 100% | N/A |
| amino acid identity min | 1.0 | 1.0 | N/A |

The low-GC failure from Run 4 is resolved. Run 4B decoded mean GC is 42.5440%, an improvement of approximately +9.30 percentage points over the Run 4 decoded mean.

## Safety Gate

| Criterion | Result | Verdict |
|---|---:|---|
| validator pass rate = 100% | 100% | pass |
| amino acid identity = 100% | min 1.0 | pass |
| no internal stop codons | 0 total | pass |
| no invalid codons | 0 total | pass |
| forbidden motif count | 0 total | pass |

## Optimization Gate

| Criterion | Result | Verdict |
|---|---:|---|
| decoded global GC mean >= 40% | 42.5440% | pass |
| decoded global GC max <= 60% | 53.8077% | pass |
| v3/v2 CAI ratio mean >= 0.97 | 1.0000 | pass |
| expected_GC >= 0.40 at final eval | 0.4274 mean | pass |
| bounded_GC_penalty remains 0 | 0.0000 mean | pass |

The minimum individual v3/v2 CAI ratio is 0.9682, and 2 of 3,876 evaluated sequences fall below 0.97 individually. The acceptance criterion is mean CAI ratio, which passes comfortably.

## Loss Log

The downloaded loss log shows the final step 20,000 eval components:

| Component | Mean |
|---|---:|
| CE_to_pseudo_label | 2.0007 |
| expected_GC | 0.4274 |
| bounded_GC_penalty | 0.0000 |
| total_loss | 2.0120 |
| synonym_mask_coverage | 1.0000 |

At the first eval step 500, expected_GC was already 0.4264 and bounded_GC_penalty was 0.0000. This supports the Run 4B recovery hypothesis: the mixed teacher distribution, not only stronger loss weighting, corrected the GC target behavior.

## v3 vs v2 Comparison

Run 4B primarily learns to reproduce the mixed v2 teacher while preserving hard biological validity.

| Metric | v2 mean | v3 mean | Delta |
|---|---:|---:|---:|
| CAI | 0.802525 | 0.802534 | +0.000009 |
| global GC | 42.543985% | 42.543970% | -0.000015 pp |
| local GC min | 38.9839% | 39.0041% | +0.0202 pp |
| local GC max | 46.1008% | 46.0922% | -0.0086 pp |

Additional v3-only aggregate metrics:

- homopolymer count mean: 0.3088
- repeat count mean: 3.7763
- forbidden motif count: 0 total

The current `run4_comparison.csv` does not include v2 homopolymer, repeat, or forbidden motif columns, so those deltas cannot be computed from the saved artifact.

## Recommendation

Proceed with Run 4B as the successful GC recovery run for v3-alpha evidence/reranking work. The result should not be framed as v3 production replacement for v2; it is best interpreted as a validated ML reproduction layer over a GC-corrected mixed teacher.

Recommended next phase:

1. Keep v2 as the production optimizer.
2. Use Run 4B as a v3-alpha candidate for evidence generation and reranking experiments.
3. Add richer v2-vs-v3 comparison columns in future eval artifacts for repeat, homopolymer, forbidden motif, and synthesis-risk deltas.
