# Run 4 Full Training Analysis

> Date: 2026-05-20  
> Scope: Run 4 Kaggle full training, constrained quick eval on 100 eval sequences  
> Boundary: in-silico sequence metrics only; no expression, yield, wet-lab, or production-replacement claims

## Summary

Run 4 training completed on Kaggle and produced a final model checkpoint plus loss log. Constrained decoding on the first 100 eval sequences preserved amino acid identity and validator safety checks, and v3 CAI stayed effectively identical to the v2 pseudo-label teacher.

The run does not meet the Run 4 success criteria because decoded global GC remains below the required 40% lower bound. The loss log confirms the GC objective was active but too weak relative to CE-to-pseudo-label: final eval expected GC averaged 0.338, below `gc_low=0.40`, with a positive bounded GC penalty still present at step 20,000.

Recommendation: stop this Run 4 candidate and run a recovery experiment before any broader v3-alpha claim. The likely recovery path is to strengthen the low-GC penalty and/or change the teacher target from pure `high_cai` v2 pseudo-labels to `gc_target` or mixed teacher labels.

## Inputs

- Config: `configs/v3_training_config_alpha_run1.yml`
- Training rows: 34,878 v2 pseudo-labels
- Eval rows: 3,876 v2 pseudo-labels
- Kaggle model checkpoint: `/kaggle/working/factorforge_run4_drive_export/model/pytorch_model.pt`
- Local results:
  - `experiments/results/alpha_run1/alpha_run1_loss_log.csv`
  - `experiments/results/alpha_run1/run4_comparison.csv`
  - `experiments/results/alpha_run1/run4_summary.json`

## Training Result

- Training completed: yes
- Steps: 20,000
- Final eval loss: 1.9964 mean total loss in the logged final eval rows
- Checkpoints saved at 5,000 / 10,000 / 15,000 / 20,000
- Final model saved by Kaggle
- Synonym mask coverage at final eval: 1.0 minimum

Final eval loss-component means at step 20,000:

| metric | value |
|---|---:|
| expected_GC | 0.337993 |
| bounded_GC_penalty | 0.062007 |
| total_loss | 1.996351 |
| synonym_mask_coverage_min | 1.000000 |

Interpretation: the synonym mask and bounded GC penalty are wired, but GC remains below the target interval. With `gc_weight=0.2`, a typical low-GC penalty near 0.06 contributes about 0.012 to total loss, which is too small compared with CE values around 1.7-2.0.

## Constrained Eval Result

The downloaded result is a 100-sequence quick eval (`decoded_count=100`), not the full 3,876-sequence eval. The quick eval is sufficient for a stop decision because every decoded sequence remains below the 40% global GC criterion.

| metric | value |
|---|---:|
| decoded_count | 100 |
| validator_pass_rate | 1.0000 |
| min_amino_acid_identity | 1.0000 |
| mean v2 CAI | 0.732558 |
| mean v3 CAI | 0.732585 |
| mean v3 CAI delta | 0.0035% |
| min v3 CAI delta | -0.2893% |
| mean v2 global GC | 33.2339% |
| mean v3 global GC | 33.2370% |
| min v3 global GC | 23.1678% |
| max v3 global GC | 38.2869% |
| mean v2 local GC min | 23.9833% |
| mean v3 local GC min | 23.9833% |
| mean v2 local GC max | 42.4500% |
| mean v3 local GC max | 42.4500% |
| internal stop codons | 0 |
| invalid codons | 0 |
| forbidden motifs | 0 |

## Criteria

### Stage 1 Safety Criteria

| criterion | result |
|---|---|
| amino acid identity = 100% | PASS |
| validator pass rate = 100% | PASS |
| no internal stop codons | PASS |
| no invalid codons | PASS |
| v3 CAI not catastrophically below v2 | PASS |

### Stage 2 Run 4 Criteria

| criterion | result |
|---|---|
| global GC >= 40% and <= 60% | FAIL |
| local GC improved vs v2 | FAIL |
| repeat/homopolymer reduced vs v2 | not demonstrated |
| forbidden motif reduced vs v2 | not demonstrated |
| CAI loss limited | PASS |

## Diagnosis

Run 4 learned the v2 pseudo-label teacher closely. This is visible in the comparison CSV: v2 and v3 CAI/GC/local GC are effectively identical for the quick eval set. That behavior confirms constrained decoding and teacher imitation, but it also means Run 4 inherits the low-GC behavior of the current `high_cai` pseudo-labels.

The final loss log shows expected GC remained in the 0.328-0.352 range. Since `bounded_GC_penalty` remained positive at step 20,000, the model had not moved into the 0.40-0.55 target interval. The current loss weighting does not overcome the CE pressure toward the low-GC teacher distribution.

## Recommendation

Stop this Run 4 candidate. Do not proceed to v3-alpha next-phase claims based on this run.

Recommended recovery job:

1. Run a small recovery smoke with `gc_weight=1.0` and optionally `gc_lambda_low=2.0`.
2. Generate a second teacher set using `gc_target` or a mixed `high_cai`/`gc_target` pseudo-label strategy.
3. Before another full training run, require a 500-2,000 step smoke to show `expected_GC >= 0.40` and decoded global GC >= 40% on a representative eval subset.
4. Keep the current validator and synonym-mask checks unchanged.

## Files

- `alpha_run1_loss_log.csv`
- `alpha_run1_comparison.csv`
- `alpha_run1_summary.json`
- `alpha_run1_analysis_report.md`
