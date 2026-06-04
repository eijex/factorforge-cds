# FactorForge Architecture

FactorForge is a multi-host sequence, codon, and expression-design engine. Its
core responsibility is to turn protein or coding-sequence inputs into candidate
coding sequences that are host-aware, rule-audited, score-calibrated, and ready
for downstream design review.

The architecture is intentionally layered so that sequence design logic, host
configuration, rule scanning, evidence tracking, developer workflows, and future
agent tooling can evolve independently.

## Layer Overview

### 1. Sequence Design Engine

The sequence design engine handles the direct sequence-design workflow:

- Reverse translation from amino-acid sequence to CDS candidates.
- Codon optimization against the selected host codon table.
- CAI scoring against codon adaptiveness weights.
- GC band scoring against profile-specific acceptable GC ranges.
- MFE and mRNA structure scoring when ViennaRNA is available.
- Assembly-friendly translation modes for cloning-oriented workflows.

This layer produces candidate sequences and metrics. It should not own host
profile policy, evidence records, or developer workflow automation.

### 2. Host Profile Registry

The host profile registry will define per-host YAML profiles for plants,
mammalian systems, insect systems, and microbial systems. A profile identifies
the species, strain or cell line, tissue or cell type, expression platform,
codon usage source, confidence level, and validation status.

The same species can have multiple profiles when expression context differs.
For example, whole-plant transient expression and suspension-cell expression
must be tracked as separate profiles even when they share a species-level codon
usage table.

### 3. Rule Engine

The rule engine scans and, where supported, repairs sequence features that can
affect expression, synthesis, cloning, or downstream handling:

- Homopolymer runs.
- Local GC windows.
- Rare codon runs.
- Restriction sites.
- PolyA signal motifs.
- Splice-like motifs.
- Repeat patterns.
- CpG and TpA dinucleotide density.

Rules should report severity and provenance separately from optimization score.
That separation lets users review biological and manufacturing risks without
silently changing the definition of an optimized sequence.

### 4. Evidence Layer

The evidence layer is planned as the source for biological and data-source
traceability. It should support:

- PubMed and NCBI lookup notes.
- Codon usage source audit records.
- Host-specific evidence tables.
- Validation notes for each profile.

Evidence records should distinguish literature support, sequence-database
support, heuristic assumptions, and wet-lab validation. Heuristic or simulated
outputs must not be presented as biological proof.

### 5. Developer Inbox

The developer inbox is planned as a workflow layer for converting operational
signals into reviewable engineering tasks:

- GitHub Actions failure inbox.
- Issue templates for reproducible implementation requests.
- Claude Code handoff documents with relevant files, scope, and commands.

The inbox is not an auto-fix system. Its role is to package context for human
review and implementation.

### 6. MCP / Skills Layer

The MCP and skills layer is planned, not implemented. It will expose selected
FactorForge design tools and evidence workflows to AI agents through explicit,
auditable tool boundaries.

Planned integrations include MCP servers for sequence operations, host profile
operations, rule audits, evidence lookup, and developer-inbox workflows. Planned
agent skills will help scaffold host profiles, review scoring changes, audit
rule semantics, and keep documentation synchronized after code changes.

## Connection Diagram

```text
                        +-----------------------------+
                        |       AI Agent Clients      |
                        | Claude Code, IDE assistants |
                        +--------------+--------------+
                                       |
                         planned MCP / skills boundary
                                       |
        +------------------------------v------------------------------+
        |                    MCP / Skills Layer                       |
        |        planned tool servers and planned agent skills        |
        +------------------------------+------------------------------+
                                       |
                 +---------------------+---------------------+
                 |                                           |
        +--------v---------+                       +---------v--------+
        | Developer Inbox  |                       |  Evidence Layer  |
        | CI, issues,      |                       | PubMed, NCBI,    |
        | handoff docs     |                       | source audits    |
        +--------+---------+                       +---------+--------+
                 |                                           |
                 +---------------------+---------------------+
                                       |
        +------------------------------v------------------------------+
        |                    Host Profile Registry                    |
        | species, cell line, platform, codon source, validation      |
        +------------------------------+------------------------------+
                                       |
        +------------------------------v------------------------------+
        |                         Rule Engine                        |
        | homopolymers, GC windows, rare codons, motifs, sites        |
        +------------------------------+------------------------------+
                                       |
        +------------------------------v------------------------------+
        |                    Sequence Design Engine                   |
        | reverse translation, CAI, GC, MFE, assembly-friendly design |
        +-------------------------------------------------------------+
```
