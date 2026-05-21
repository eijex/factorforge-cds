# FactorForge v3 Audit Report

> Date: 2026-05-19  
> Scope: STEP 1 audit for v3-alpha repositioning  
> Repository: `C:\Work\eijex\factorforge`

## Summary

The live repository uses package name `factorforge`, not `codonforge`. The v2 engine is a production rule-based optimizer under `src/factorforge/engines/v2/`. The v3 implementation is currently a scaffold plus training scripts: tokenizer, BART decoder skeleton, v3 pipeline fallback path, and `scripts/3_training/train_v3_esm2_bart.py`.

The key failure mode remains the training target: v3 is trained with cross-entropy against native CDS codons. That teaches native codon usage, not expression-optimized synthetic CDS generation. Current loss has a differentiable expected-GC penalty connected to logits, but it is a single-point absolute penalty, has no expected CAI term, and does not apply synonymous codon constraints.

## Files Inspected

- `src/factorforge/engines/v2/optimizer.py`
- `src/factorforge/engines/v2/rules/reverse_translator.py`
- `src/factorforge/engines/v2/rules/rule_engine.py`
- `src/factorforge/engines/v2/rules/domesticator.py`
- `src/factorforge/engines/v2/scoring.py`
- `src/factorforge/engines/v2/utils.py`
- `src/factorforge/engines/v2/validator.py`
- `src/factorforge/engines/v3/tokenizer.py`
- `src/factorforge/engines/v3/metrics.py`
- `src/factorforge/engines/v3/modeling_bart_decoder.py`
- `src/factorforge/engines/v3/pipeline.py`
- `scripts/3_training/dataset.py`
- `scripts/3_training/train_v3_esm2_bart.py`
- `scripts/4_evaluation/compare_baseline.py`
- `configs/v3_training_config_run2.yml`
- `data/nbenthamiana_codons.json`
- `data/nbenthamiana_golden_set.json`
- `tests/engines/v2/*`
- `tests/engines/v3/*`

## Current v2 Flow

`RuleBasedOptimizer.optimize()`:

1. Validates input with `InputValidator`.
2. Normalizes optimization profile.
3. For protein input, calls `ReverseTranslator.generate_candidates(..., n=1)`.
4. For DNA input, reuses the DNA sequence and calculates metrics.
5. Runs `RuleEngine.scan_all()` for rule checks.
6. Returns `OptimizationResult` with DNA sequence, CAI, GC, composite score, violations, and scan metadata.

`ReverseTranslator` builds amino-acid-to-codon maps from bundled codon tables, supports balanced/high-CAI/GC-target/assembly/ramp/viral profiles, calculates CAI using relative adaptiveness weights, and calculates GC as a percentage.

## Current v3 Flow

`V3Pipeline.run()`:

1. Normalizes amino acid input.
2. If encoder embeddings and decoder are supplied, calls greedy decoder generation.
3. If no embeddings are supplied, falls back to v2 `gc_target` through `EngineRegistry`.
4. Applies post-guard checks/fixes for PolyA, rule scan, and domestication.
5. Computes v3 metrics with `compute_cai()` and `compute_gc()`.
6. Builds an explainability/FDA-style report.

The training path is separate in `scripts/3_training/train_v3_esm2_bart.py`. It trains a BART-style decoder using per-token ESM2 embeddings and codon labels from JSONL data.

## CAI Calculation

v2:

- `ReverseTranslator.calculate_cai()` uses golden-set relative adaptiveness weights when available.
- Falls back to the standard codon table and per-amino-acid max frequency.
- Uses geometric mean of codon weights.
- Skips stop codons.
- Returns `0.0` for sequence lengths not divisible by 3.

v3:

- `compute_cai()` in `src/factorforge/engines/v3/metrics.py` uses a `CodonUsageTable`.
- Builds relative adaptiveness weights from the provided codon usage table.
- Uses geometric mean.
- Returns `0.0` if any codon is missing or has non-positive weight.
- Does not explicitly validate that sequence length is divisible by 3; it truncates trailing bases via `len(seq) // 3`.

## GC Calculation

v2:

- `calculate_gc()` in `src/factorforge/engines/v2/utils.py` returns percentage `0-100`.
- `ReverseTranslator.calculate_gc_content()` rounds to two decimals.

v3:

- `compute_gc()` returns percentage `0-100`.
- Training loss builds a codon-level GC vector in fraction units `0.0-1.0`.
- Training config uses `gc_target: 0.425`, also fraction units.

## Codon Vocabulary

`CodonTokenizer.default()` creates:

- 5 special tokens: `[PAD]`, `[BOS]`, `[EOS]`, `[MASK]`, `[UNK]`
- 64 DNA codon tokens in deterministic T/C/A/G base order
- Total vocabulary size: 69

This is a one-token-per-codon vocabulary, so expected GC and expected CAI can be implemented directly over logits/probabilities.

## Special Tokens in Metrics

Runtime v3 metrics operate on decoded DNA strings, so skipped special tokens are excluded when `CodonTokenizer.decode(..., skip_special_tokens=True)` is used. The training GC vector assigns `0.0` GC to special tokens and masks target positions with `labels.ge(5)`, so special-token label positions are excluded from the current GC penalty.

Risk: logits can still allocate probability mass to special tokens at valid codon positions. Because their GC vector value is `0.0`, this can lower expected GC unless special-token logits are masked or their probabilities are excluded/renormalized for expected metrics.

## GC Penalty

Current implementation:

- `_compute_gc_penalty()` in `scripts/3_training/train_v3_esm2_bart.py`
- Computes `softmax(logits)` over vocabulary.
- Multiplies probabilities by a codon GC vector.
- Averages expected GC over positions where `labels != -100` and `labels >= 5`.
- Returns `abs(predicted_gc - gc_target)`.

Findings:

- Differentiable: yes, gradient flows from loss to logits.
- Computed from logits/probabilities: yes.
- Computed from decoded argmax tokens: no.
- Padding labels excluded: yes via `labels.ne(-100)`.
- Special-token labels excluded: yes via `labels.ge(5)`.
- Special-token probabilities at valid positions excluded: no; they contribute `0.0` to expected GC.
- Expected GC implemented: yes, but with the special-probability caveat above.
- Bounded GC penalty implemented: no.

## Expected CAI / Log CAI

Expected CAI or expected log CAI is not implemented in the current v3 training loss. Current training objective is CE plus optional GC penalty.

## Synonymous Codon Masking

No v3 training or decoding synonym mask was found.

- Training does not restrict probability mass to codons synonymous with the input amino acid.
- Greedy decoding uses unrestricted `argmax`.
- v3 can output any codon token, including codons incompatible with the input amino acid.
- Decoding can output stop codons if logits favor `TAA`, `TAG`, or `TGA`.

v2 rule-engine/domestication code performs synonymous substitutions for specific fixes, but that is not a v3 decoder constraint.

## Amino Acid Preservation

v2 reverse translation is designed to preserve amino acid sequence for protein inputs.

v3 neural decoding does not guarantee preservation. The fallback path without embeddings delegates to v2 `gc_target`, but actual decoder generation is unconstrained.

Additional audit note: v2 auto-detection can classify short protein strings made only of DNA/IUPAC-compatible letters, such as `MAST`, as DNA. CLI or diagnostic callers that know an input is protein should force protein handling rather than relying on auto-detection. This also affects the v3 no-embedding fallback because it delegates to v2 through `EngineRegistry`.

## Internal Stop Codon Risk

v2 can include a terminal stop when the input protein contains `*`. For normal proteins without `*`, v2 should not emit internal stops. Existing v2 codon table includes stop codons under `*`, and tests now check that standard amino-acid codon families exclude stop codons.

v3 decoder can emit stop codons because the codon vocabulary includes all 64 DNA triplets and no synonym/stop mask is applied.

## Suspected Causes of GC Collapse to 31-33%

Likely contributors:

1. Native-CDS objective mismatch: CE trains the decoder toward native codon distributions rather than optimized synthetic CDS.
2. Single-point GC penalty conflict: CE to native labels and an absolute GC target can compete without teaching optimized codon selection.
3. No expected CAI reward: nothing directly encourages high-CAI codon probability.
4. No synonym mask: probability mass can move to invalid or context-incompatible codons.
5. Special-token probability leakage: special-token logits at valid positions contribute `0.0` expected GC under the current vector.
6. Narrow global GC target may be infeasible or undesirable for some proteins without a feasibility check.

## Recommended Next Code Changes

1. Add a formal metrics and validator layer with machine-readable validation results.
2. Add Pareto feasibility analysis before judging CAI/GC target failure.
3. Add a v2 adapter and reproducible benchmark harness.
4. Replace single-point GC loss with bounded expected-GC penalty.
5. Add expected log CAI over codon probabilities.
6. Exclude or renormalize special-token probabilities in expected metric losses.
7. Add synonymous codon masks during training and decoding.
8. Add post-decoding validator and fallback-to-v2 policy.

## Tests Added in STEP 1

- GC calculation.
- CAI calculation on a simple known table.
- DNA translation with standard codon mapping.
- Codon-to-amino-acid mapping sanity.
- Reverse translation amino acid preservation.
- v2 output has no internal stop codons for normal protein input.
- Every standard amino acid maps only to valid synonymous non-stop codons.
- Tokenizer special tokens are disjoint from codon tokens and skipped during decode.

## Constraint Feasibility Analysis

STEP 3 adds `src/factorforge/ml/feasibility.py`, which computes exact global GC and CAI feasibility over synonymous codon choices using dynamic programming. The analyzer keeps the best log-CAI sequence for each reachable GC count, so impossible GC/CAI targets can be reported separately from model failure.

The CLI entry point is `scripts/check_constraint_feasibility.py`. It accepts either direct protein input or a FASTA file and reports:

- minimum and maximum possible GC;
- maximum achievable CAI without a GC constraint;
- maximum CAI under GC 41-44%, 40-50%, and 40-55%;
- target feasible/impossible for user-specified CAI and GC bounds;
- local GC summaries for best candidates.

## v2 Baseline Reproducibility

STEP 4 adds `src/factorforge/engines/v3/inference/v2_adapter.py` as the formal v2 adapter for v3-alpha. The adapter treats the input as protein sequence, avoiding the short-protein auto-detection issue found in STEP 1, and returns:

- engine identity;
- protein sequence;
- DNA sequence;
- CAI, GC, and composite score;
- structured validator output;
- v2 rule-scan metadata;
- warnings and errors.

The reproducible benchmark runner is `scripts/run_v2_baseline_benchmark.py`. It uses the deterministic `high_cai` profile over safe synthetic fixtures in `tests/fixtures/benchmark_proteins.py` and writes:

- `experiments/results/metrics/v2_baseline_metrics.csv`
- `experiments/results/metrics/v2_baseline_metrics.json`

## Pre-Run-4 Safety Gate Follow-up

After STEP 6 review, the immediate blocker around synonymous constraints was addressed in code without starting full training:

- `src/factorforge/engines/v3/synonym_mask.py` builds per-position boolean masks from protein sequence and codon tokenizer.
- `scripts/3_training/dataset.py` now requires protein sequence metadata and carries `synonym_mask` through `collate_fn`.
- `scripts/3_training/train_v3_esm2_bart.py` passes `synonym_mask` into expected GC and expected log CAI.
- `src/factorforge/engines/v3/inference/constrained_decoder.py` provides constrained greedy decoding and validator-backed fallback to v2.
- `src/factorforge/engines/v3/pipeline.py` uses constrained decoding when a decoder and encoder embeddings are supplied.

The safety gate does not run training and does not change v2 production behavior.

Remaining Run 4 blockers:

- approved protein FASTA selection and actual train/eval pseudo-label generation;
- tiny training smoke with fixture data after explicit approval;
- final review before user-facing v3-alpha claims.

## Run 4 Prep Follow-up

Added after the safety gate:

- `scripts/1_data_preparation/generate_v2_pseudolabels.py`
- `configs/v3_training_config_alpha_run1.yml`
- `docs/design/run4_experiment_plan.md`
- `docs/design/v3_alpha_product_boundary.md`
- `docs/model_cards/factorforge_v3_alpha_model_card.md`

Verification:

- pseudo-label generator tests pass;
- Run 4 config dry-run initializes the model without training;
- full test suite passes.
