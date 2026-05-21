# FactorForge v3 Loss Debug Report

> Date: 2026-05-20  
> Scope: STEP 5 differentiable GC/CAI loss audit and unit-tested components

## Current Loss Location

The active v3 ESM2+BART training loss is implemented in:

- `scripts/3_training/train_v3_esm2_bart.py`

There is no separate `src/factorforge/engines/v3/training/loss.py` module in the current repository.

## Pre-STEP 5 Findings

- Cross-entropy is computed over codon token labels with `ignore_index=-100`.
- GC penalty was already differentiable because it used `softmax(logits)` and a codon-level GC vector.
- GC penalty was a single-point absolute penalty: `abs(expected_gc - gc_target)`.
- Padding labels were excluded by `labels != -100`.
- Special-token label positions were excluded by `labels >= 5`.
- Special-token probability mass at valid codon positions was not renormalized away; special tokens had GC value `0.0`.
- Expected CAI / expected log CAI was not implemented.
- Synonymous codon masks were not implemented for training or decoding.

## STEP 5 Changes

Implemented in `scripts/3_training/train_v3_esm2_bart.py`:

- `_codon_token_mask()`
- `_codon_log_cai_vector()`
- `_valid_label_positions()`
- `_masked_codon_probs()`
- `_compute_expected_gc()`
- `_compute_expected_log_cai()`
- `_compute_bounded_gc_penalty()`

The new expected-metric functions mask non-codon tokens before softmax, so special-token probability mass is excluded and codon probabilities are renormalized over valid codon tokens. Padding positions remain excluded through the label mask.

## Bounded GC Penalty

The training loss now supports:

```text
loss_gc = lambda_low * relu(gc_low - expected_gc)
        + lambda_high * relu(expected_gc - gc_high)
```

Default bounds:

- `gc_low = 0.40`
- `gc_high = 0.55`

Config keys:

- `gc_low`
- `gc_high`
- `gc_lambda_low`
- `gc_lambda_high`
- `gc_weight`

## Expected Log CAI

Expected log CAI is computed as:

```text
mean_i sum_c p_i(c) * log(w_c)
```

where `p_i(c)` is the masked codon probability at position `i`, and `w_c` is the codon relative adaptiveness weight from the current v3 codon usage table. It is opt-in through `cai_weight` or `expected_log_cai_weight`.

## Remaining Gaps

- Synonymous codon masks are now produced by the dataset/collate path when JSONL records include `protein_sequence`, `protein_seq`, or `sequence`.
- Expected GC and expected log CAI now receive the batch synonym mask when available.
- A constrained greedy decoding helper now restricts codon selection to synonymous codons for each amino-acid position.
- No full training run was started.
- Run 4 should not proceed until STEP 6 review decides whether the current audit/test coverage is sufficient.

## Pre-Run-4 Safety Gate Follow-up

Added after STEP 6 review:

- `src/factorforge/engines/v3/synonym_mask.py`
- `src/factorforge/engines/v3/inference/constrained_decoder.py`
- dataset/collate propagation of `synonym_mask`
- training loss use of `synonym_mask`
- constrained v3 decode with validator and v2 fallback helper

Remaining before full Run 4:

- approved protein FASTA selection and actual train/eval pseudo-label generation;
- decide whether to apply synonym masks to CE logits directly, not only metric losses;
- run a tiny training smoke only after explicit approval.

Run 4 config `configs/v3_training_config_alpha_run1.yml` now sets bounded GC and expected log CAI weights, and `--dry-run` model initialization succeeds without loading training data or starting optimization.

## Tests

Added `tests/engines/v3/test_training_loss.py`.

Covered:

- deterministic expected GC;
- expected GC changes with high-GC vs low-GC logits;
- bounded GC penalty inside/below/above range;
- expected log CAI increases for preferred codons;
- gradients flow back to logits;
- special tokens and padding do not affect expected GC/CAI.
- synonym masks allow only synonymous codons;
- constrained decoding preserves amino acid sequence;
- invalid v3 candidates fall back to v2.
