# Sequence Policy Boundary

FactorForge does not define, approve, or validate a laboratory standard
operating procedure (SOP). It provides a deterministic CDS-design engine,
structured sequence findings, and a versioned mechanism for applying the
sequence-review policy chosen by the user or laboratory.

The public product term is **Sequence Policy Profile**. The existing
`SopProfile` name remains in the Python and API implementation for backward
compatibility; it should not be interpreted as a claim that FactorForge covers
or certifies a complete laboratory SOP.

## Responsibility split

### FactorForge invariants

These are software and CDS-integrity guarantees enforced by FactorForge rather
than optional laboratory preferences:

- the output is a valid CDS;
- the translated amino-acid sequence is preserved;
- no unintended internal stop codon is introduced;
- reading-frame and schema requirements remain valid; and
- deterministic findings and provenance are reported reproducibly.

### User or laboratory policy

The treatment of advisory findings depends on the experimental and assembly
workflow. The user or laboratory decides whether findings such as splice-like
motifs, PolyA-like motifs, AU-rich regions, local-GC warnings, homopolymers, or
specific restriction sites are required, advisory, or ignored.

FactorForge can express and apply that decision consistently. It does not claim
that one setting is biologically correct for every laboratory, host, construct,
or downstream process.

## Bundled example

The bundled FactorForge policy is an **illustrative example configuration**. It
reports supported advisory findings without treating an unvalidated predictor
signal as biological ground truth. It is not:

- a validated laboratory SOP;
- a policy approved by PlantForm or any other laboratory;
- regulatory or biosafety guidance;
- wet-lab evidence; or
- a guarantee of expression, yield, synthesis acceptance, or biological
  performance.

Users should adapt the example to their own reviewed workflow. A laboratory-
specific profile remains the responsibility of that laboratory.

## Provenance contract

When a Sequence Policy Profile is used, FactorForge records its version and
digest with the computational result. This supports reproducibility; it does
not validate the policy itself.

For prospective studies, preserve the generated design separately from the
post-review built construct. Any expert-requested sequence change should be
recorded as a remediation event rather than attributed silently to the design
engine.

## Privacy and contributions

Do not publish confidential laboratory policies, raw sequences, construct
identifiers, or private process parameters. Public contributions must be
sanitized, sequence-free suggestions. Laboratory-specific profiles should be
shared only through an authorized private channel.

