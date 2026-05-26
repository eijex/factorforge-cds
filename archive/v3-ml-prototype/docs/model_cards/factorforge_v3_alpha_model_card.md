# FactorForge v3-alpha Model Card

## Model Summary

FactorForge v3-alpha is an experimental ML-assisted CDS candidate generation and evidence layer for plant-oriented codon optimization workflows. It combines a frozen protein-embedding encoder design with a BART-style codon decoder and rule-based v2 fallback components.

## Intended Use

Intended uses:

- in-silico CDS candidate generation;
- candidate scoring and reranking;
- v2 baseline comparison;
- sequence-level validator reporting;
- pre-wet-lab candidate prioritization.

Not intended for:

- guaranteed expression or yield prediction;
- clinical, regulatory, or manufacturing validation;
- use without human review;
- replacement of wet-lab testing.

## Training Data Boundary

Run 4 should use validated v2 pseudo-labels and constrained candidate data, not native CDS labels as the final optimization target. Native CDS may be used as background/reference data.

## Safety Controls

Current controls:

- synonym mask generation from protein sequence;
- constrained codon decoding;
- structured validator;
- v2 fallback on invalid v3 candidates;
- CAI/GC feasibility analysis;
- bounded expected GC loss;
- expected log CAI objective.

## Known Limitations

- No wet-lab validation is currently available.
- ESM2 benefit has not been proven by ablation.
- In-silico metrics are proxies and should not be treated as biological proof.
- Local GC, motif, and synthesis checks reduce known sequence-level risks but do not guarantee expression.

## Evaluation Requirements

Before promotion beyond alpha:

- amino acid preservation must be 100% on benchmark candidates;
- internal stop and invalid codon counts must be zero;
- v2 fallback must remain available;
- v3 must be compared against v2, native CDS, random synonymous baseline, and constrained oracle where feasible;
- ESM2 ablation must be reported;
- claims-control language must be reviewed.

## Approved Claims

- "provides in-silico construct design evidence";
- "supports candidate prioritization";
- "uses validator-backed sequence checks";
- "compares candidates against a v2 baseline."

## Disallowed Claims

- "guarantees expression";
- "increases yield";
- "wet-lab validated";
- "regulatory-ready";
- "replaces v2 production engine."

