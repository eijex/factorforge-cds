# Manuscript and AgentOS integration audit

Date: 2026-07-18  
Status: implementation audit in progress  
Workspace job: eijex-workspace Job 211

## Verified state

- Repository: branch `update-conda-recipe-sha-v3.3.2`, commit `347f4f383e865740b501204e743c75da03eed84f`.
- Package, REST API, and static web version: `3.3.2`.
- Manuscript draft: v0.7.3; archived software boundary: FactorForge v3.3.0.
- Canonical source: `eijex-workspace/papers/manuscripts/_prism/paper1-factorforge/0.7.3/2026-07-14/working/v073-reproducibility-closure-patch/source/`.
- Frozen historical anchor: SGN QLD183 v103, 49,257 filtered unique CDS, FactorForge v3.2.0 artifact layer, `scoring_contract v1.1`, historical packaged legacy codon reference, historical 55–65% scoring window.
- Current production anchor: `nbenthamiana_nbev11_hc_v2`, NbeV1.1 high-confidence CDS-derived reference, default GC policy 40–47%.

## Artifact locations and immutability

- Frozen public benchmark: `reproducibility/benchmark_v0.5.1/`.
- Historical results: `benchmarks/results/v3.2.0/`.
- Sensitivity artifacts: `reproducibility/benchmark_v3.2.2/`.
- Reference manifests: `src/factorforge/data/reference_policy_manifest.json` and related schemas under `schemas/`.

These paths are read-only for Job 211. Table 5 and Figures 2–5 will not be regenerated.

## Discrepancies found

1. REST profile optimization does not accept or propagate a seed even though the profile engine supports it.
2. Omitted profile seeds are generated internally and are unrecoverable.
3. DesignPackage-compatible provenance does not record or hash applicable seeds.
4. The CLI does not expose a profile seed.
5. The web uses generic completion wording when warnings or an unmet GC target may be present.
6. The internal 10-year blueprint retains a v3.3 Algorithm Depth collision already corrected in the public roadmap.
7. Job 210 already defines broad cross-system data architecture, so Job 211 must specialize benchmark governance and preserve ValidationHub wet-lab ownership.

## Reference-verification task

Manually confirm that cited sources collectively support algorithm choice, host context, codon-reference source, sequence-composition treatment, RNA-structure treatment, and evaluation metric. Unsupported subclaims must be narrowed or receive a verified citation; support is not inferred from citation proximity.

## Proposed modifications

- Canonical v0.7.3 LaTeX sections.
- Codon-reference provenance documentation and consistency tests.
- REST handler, profile metadata/pipeline, CLI, internal DesignPackage-compatible model, web result rendering, and focused tests.
- `docs/agentos/` and `schemas/agentos_benchmark_schema.sql`.
- Workspace Job 211 registration and completion report.

## Release-auditor ratings

- Step 0 environment integrity: **4/5** — repositories and `_version/` accessible; `main..develop` has no commits.
- Step 1 software version state: **4/5** — live 3.3.2 surfaces align.
- Step 2 manuscript/code alignment: **4/5** — source correctly separates v3.3.0 software, v3.2.0 evidence, QLD183 corpus, and legacy reference; wording/seed closure remained open at audit time.
- Step 3 roadmap direction: **3/5** — internal blueprint collision remains.
- Step 4 backlog triage: **3/5** — recent ideas are mostly long-horizon or evidence-volume dependent.
- Step 5 ground truth: **4/5 (tests pending)** — Job 210 is indexed; post-change test verification remains open.

Items above are `[Confirmed]` from live files/commands except citation-support mapping and post-change tests, which remain `[Not verified]` until completion.
