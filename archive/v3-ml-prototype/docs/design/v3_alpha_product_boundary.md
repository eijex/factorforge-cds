# FactorForge v3-alpha Product Boundary

> Status: draft  
> Applies to: v3-alpha candidate generation, scoring, reranking, and evidence reports

## What v3-alpha Does

FactorForge v3-alpha is an in-silico CDS design assistant. It can:

- generate candidate CDS sequences;
- compare v2, v3, and future constrained-oracle candidates;
- report CAI, global GC, local GC, first-region GC, motifs, repeats, homopolymers, and validator status;
- provide an in-silico evidence report;
- support human review and candidate prioritization.

## What v3-alpha Does Not Do

FactorForge v3-alpha does not:

- guarantee protein expression or yield;
- replace wet-lab validation;
- optimize promoter, UTR, vector, host growth, agroinfiltration, purification, or assay conditions;
- replace v2 production behavior;
- make regulatory, clinical, safety, or manufacturing validation claims.

## Valid Claims

Allowed wording:

- "designed to reduce known sequence-level risks";
- "provides in-silico construct design evidence";
- "supports prioritization of candidates for future validation";
- "compares candidates against rule-based v2 baseline metrics";
- "uses hard validation checks before recommending candidates."

## Forbidden Claims

Do not use:

- "increases yield";
- "guarantees expression";
- "validated expression optimizer";
- "wet-lab proven";
- "production replacement for v2";
- "regulatory-ready construct";
- "clinically validated";
- "manufacturing-grade yield optimization."

## Human Review Policy

All v3-alpha outputs require human review before synthesis, wet-lab testing, partner sharing, or operational use. The reviewer should inspect:

- amino acid preservation;
- stop codons and invalid codons;
- global and local GC;
- forbidden motifs and repeats;
- feasibility analysis for stated targets;
- v2 fallback comparison.

## Data Logging Requirements

Each v3-alpha output should log:

- input protein sequence hash or approved identifier;
- engine and model version;
- tokenizer hash;
- config file path/hash;
- v2 fallback candidate;
- candidate DNA sequence;
- validator result;
- metrics result;
- fallback decision;
- warnings/errors;
- timestamp and run identifier.

Do not log secrets, unpublished proprietary sequences, or partner-confidential sequence assets unless the deployment has explicit storage approval.

## Future DMTL Loop

Future design-make-test-learn work may incorporate wet-lab expression data. Until such data exists and is reviewed, v3-alpha remains an in-silico evidence layer. Future model updates should separate:

- design features;
- wet-lab assay context;
- measured outcome;
- batch and protocol metadata;
- uncertainty and failure modes.

## BD-Safe Summary

FactorForge v3-alpha is an in-silico CDS design and evidence layer that compares candidate sequences against a rule-based v2 baseline, applies hard sequence validation, and supports prioritization for future experimental validation. It does not claim or guarantee protein yield improvement.

