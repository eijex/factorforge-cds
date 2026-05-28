# v3 ML Prototype Archive

This directory preserves the historical ML-oriented prototype that was developed
as an attempted successor to the internal v2 rule-based engine.

## Context

FactorForge's development has gone through several generations:

- **v1 (NBent_OptiCodon)** — thesis-derived codon optimization baseline; internal use only
- **v2 (Rule-Based Engine)** — deterministic, constraint-aware redesign; internal; succeeded and became the production foundation
- **v3-alpha (this directory)** — ML-based design attempt; intended as a step up from v2 but did not reach production quality
- **v3.0+ (current public release)** — the v2 engine released as open source under `factorforge.engines.profile`

## Why This Prototype Was Archived

The ML prototype explored sequence-level generative approaches for codon design.
Training and evaluation showed that model performance was insufficient compared
to the deterministic v2 engine on the key metrics (CAI, GC constraint adherence,
PolyA signal elimination). The deterministic approach — constraint-based dynamic
programming over a synonymous codon table — proved more reliable and interpretable
for the target use case.

The prototype was archived rather than deleted to preserve the research record.
The lessons from this track inform the planned ML engine (v3.7+), which will
revisit ML-based design as `--engine ml` once sufficient wet-lab validation data is available.

## What Is Preserved

Training and evaluation code, model card materials, and early decoder experiments.
Large generated experiment outputs (CSV loss logs, decoded comparison tables) are
excluded to keep the repository compact.

The historical Kaggle training notebook is intentionally excluded from the
public archive because it contains external workspace and checkpoint handling
that is not part of the supported release.

## Relationship to Current Package

The archived prototype is not installed, registered, or used by FactorForge
v3.1.x. Current public releases use deterministic, constraint-aware design paths
under `factorforge.engines.profile`.
