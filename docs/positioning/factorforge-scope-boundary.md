# FactorForge Positioning: Post-Protein-Design, Pre-Synthesis CDS Engineering

## Decision

FactorForge is not a protein folding, inverse-folding, protein representation,
or generative protein-sequence model.

It operates after the amino-acid sequence has been selected and before DNA
synthesis, cloning, assembly, expression testing, and wet-lab validation.

## Functional boundary

Input:
- a fixed amino-acid sequence
- host and codon-reference context
- configured sequence and assembly constraints

Output:
- synonymous CDS candidates
- computational validation results
- warnings and review metadata
- reproducible design artifacts

## Non-claims

FactorForge does not establish:
- protein structure or folding
- expression success or yield
- cloning or synthesis acceptance
- biological efficacy
- regulatory readiness

A candidate that passes all configured FactorForge checks satisfies those
checks only. It remains unevaluated with respect to protein folding, in
planta expression, yield, and functional performance — those outcomes
require separate protein-level analysis and experimental evidence.

## Relationship to upstream AI

Protein-design and representation models — fold, inverse-fold, read
(representation/embedding), and write (generative sequence design) — may be
used upstream to select or evaluate a protein sequence. Their outputs may
become FactorForge inputs, but those models are not part of FactorForge's
core, and FactorForge does not attempt to replicate or replace them.

An earlier internal research prototype explored an ESM2 (read) + BART
(write) architecture for codon-sequence generation. That prototype is
archived. It differed from generative protein-design models such as ProGen
in an important way even while active: its "write" step generated
conditional synonymous CDS candidates for an already-fixed protein, not new
amino-acid sequences. The current, deterministic constraint-based design and
validation pipeline is unrelated to that archived prototype.

## Strategic rationale

The project prioritizes deterministic CDS design, explicit claim boundaries,
and accumulation of design-test evidence through ValidationHub over
premature development of large learned sequence-generation models. As
design-outcome evidence accumulates through that loop, future versions may add
calibrated, evidence-conditioned review layers on top of this same
computational-CDS-design scope — without changing what FactorForge claims
about downstream biological performance today.
