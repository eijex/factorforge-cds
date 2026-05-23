# v3 ML Training Results

Internal R&D record. Not user-facing.

---

## Key Findings

- GFP CAI ≥ 0.92 + GC 41–44% simultaneously is **infeasible** under current *N. benthamiana* codon weights (Pareto max CAI at GC 41–44%: 0.884; at GC 40–55%: 0.969)
- GC collapse (32–33%) root cause: native CDS CE objective mismatch + single-target GC penalty + no expected CAI reward + no synonym mask + special-token probability leakage
- Synonym mask + constrained decoding now prevent non-synonymous codon output and internal stop codons

---

## Training Runs

**Run 1–3** (GFP, early experiments):

| Metric | v2 Baseline | v3 Run 1 | v3 Run 2 | v3 Run 3 |
|--------|-------------|---------|---------|---------|
| CAI | 0.909 | 0.974 | 0.717 | 0.718 |
| GC% | 56.4% | 60.9% | 32.6% | 30.9% |
| Eval loss | — | 6.66 | 2.505 | 2.734 |
| Training data | — | 917 | 38,754 | 38,754 |

**alpha_run1** (pseudo-label high_cai, 34,878 sequences):

| Metric | v2 baseline | v3 alpha_run1 |
|--------|-------------|--------------|
| CAI mean | 0.733 | 0.733 |
| GC% mean | 33.23% | 33.24% |
| Validator pass rate | — | 100% |
| AA identity | — | 100% |

Outcome: **GC target miss** — decoded GC remained below 40% target. Root cause: `gc_weight=0.2` too weak relative to CE loss.

**alpha_run2** (pseudo-label mixed, ~30,000 sequences):

| Metric | v2 baseline | v3 alpha_run2 |
|--------|-------------|--------------|
| CAI mean | 0.8025 | 0.8025 |
| GC% mean | 42.54% | 42.54% |
| GC% range | 40.36–53.81% | 41.99–53.81% |
| Validator pass rate | — | 100% |
| AA identity | — | 100% |
| Eval sequences | — | 3,876 |

Outcome: **✅ GC recovery** — mixed teacher distribution + `gc_weight=2.0` resolved GC collapse.
