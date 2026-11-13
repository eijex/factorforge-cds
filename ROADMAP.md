# FactorForge — Research & Development Roadmap

Baseline: November 2026 · Current public release: **v3.6.4**

> [Changelog](CHANGELOG.md) · [Validation Log](VALIDATION.md) · [Releases](https://github.com/eijex/factorforge-cds/releases) · [Milestones](https://github.com/eijex/factorforge-cds/milestones) · [Web App](https://factorforge.eijex.com)

---

## Mission

FactorForge develops reproducible, constraint-aware coding sequence (CDS) design and pre-synthesis review methods for plant-based protein expression, with primary focus on *Nicotiana benthamiana*.

Our priority is scientific validation and reproducibility across existing design strategies (`1DP`, `3DP`, `4DP`), not expansion into unrelated AI models or biological design domains. Historical release details are maintained in [CHANGELOG.md](CHANGELOG.md).

---

## Phase 1 — Core Reliability *(Current priority)*

Maintain and improve the deterministic CDS optimization engine (`1DP`) and pre-synthesis review workflow:

- Strengthen constraint correctness, computational efficiency, reproducible benchmarks, and host-specific codon references for *N. benthamiana*.
- Verify software-level correctness, hard-constraint satisfaction, and reproducibility across 1DP, experimental 3DP, and experimental 4DP workflows.
- Maintain and calibrate deterministic RNA motif QC (splice-like motifs, poly-A-like motifs, and AU-rich elements surfaced in the 1DP/3DP/4DP comparison matrix) without adding heavy foundation model dependencies.
- Complete independent holdout evaluation before promoting any experimental optimizer to a production default.

**Success criteria:** Verified constraint satisfaction, reproducible results, and documented performance against established baselines.

---

## Phase 2 — Structure & Translation-Aware Research *(Experimental research)*

Evaluate the existing structure-guided (`3DP`) and translation-offset (`4DP`) design approaches already exposed as research workflows:

- Use externally generated protein structures (AlphaFold 3, ESMFold), explicitly reviewed structural annotations, and literature-supported hypotheses.
- Benchmark codon harmonization, local translation features, and `1DP` vs. `3DP` vs. `4DP` sequence differences without introducing additional foundation models by default.
- Pre-existing experimental `learning/` and `sllm/` research adapters in the repository remain isolated prototypes and are not promoted to production defaults or promised features.

**Success criteria:** Reproducible computational comparisons, transparent assumptions, and clearly documented limitations.

---

## Phase 3 — Wet-Lab Validation *(Dependent on experimental collaboration)*

Test whether `3DP` and `4DP` provide measurable biological benefits over baseline and `1DP` CDS designs:

- Conduct matched comparisons of baseline, 1DP, 3DP, and 4DP Adalimumab heavy-chain (HC) and light-chain (LC) CDS designs. Evaluate computational differences independently from experimentally measured expression, yield, and protein quality.
- Record expression, yield, and relevant protein-quality measurements using predefined experimental endpoints (see [VALIDATION.md](VALIDATION.md)).

**Collaboration & Feedback:** Public-safe summaries may be submitted via [GitHub Issues](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml) under the [Wet-lab Validation](https://github.com/eijex/factorforge-cds/milestone/9) milestone. Private or confidential experimental feedback should be sent to `eijex.lab@gmail.com` instead of public issues.

**Success criteria:** Reproducible experimental evidence establishing whether particular design approaches offer measurable benefits under tested conditions.

---

## Phase 4 — Evidence-Driven Improvement *(Conditional future development)*

Pursue model or host expansion only when empirical evidence and data readiness justify it:

- Refine design objectives and host-specific policies using validated experimental evidence.
- Explore learned ranking or constrained machine-learning methods only when suitable training data and independent holdout evaluation datasets are available, without weakening deterministic sequence constraints.
- **Host scope:** Primary research host is *Nicotiana benthamiana*. Broader host expansion (`Sf9`, `Tni`, *Arabidopsis*, tomato, *Lemna*, *Wolffia*) is deferred indefinitely unless supported by curated reference data and validation capacity.

**Success criteria:** Demonstrated improvements supported by independent validation, preserving strict deterministic invariants.

---

## Scope Boundaries

- FactorForge is **not** intended to replace protein structure prediction systems (use AlphaFold / ESMFold directly), train or deploy general-purpose biological foundation models, perform *de novo* therapeutic discovery, support animal/stem-cell or microbial (`E. coli`, yeast) expression systems, or control laboratory robotics.
- External AI tools may assist exploratory research (or separate PRISM research tracks), but are not core FactorForge product dependencies.
- All computational outputs remain in-silico research candidates (`UNVERIFIED_PUBLIC_CORPUS`, `Biological validation: NOT_AVAILABLE`) unless supported by relevant experimental validation.
