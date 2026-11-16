---
title: 'FactorForge: Reproducible CDS Candidate Generation and Sequence Review for Plant Molecular Pharming'
tags:
  - Python
  - bioinformatics
  - codon optimization
  - dynamic programming
  - plant synthetic biology
  - reproducibility
authors:
  - name: Mun-Kyu Kim
    orcid: 0009-0004-1413-7726
    corresponding: true
    affiliation: 1
affiliations:
  - name: PlantForm Corporation, 1920 Yonge Street, Suite 200, Toronto, Ontario, M4S 3E2, CANADA
    index: 1
date: 16 November 2026
bibliography: paper.bib
---

# Summary

FactorForge is an open-source Python package for coding-sequence (CDS) design and pre-synthesis review in plant molecular pharming. Synonymous DNA sequences encoding the same protein can differ markedly in codon usage, GC content, and unintended restriction or RNA processing motifs. FactorForge generates synonymous candidates under explicit user constraints, runs multi-category sequence checks, and records structured JSON review artifacts before DNA synthesis.

This paper describes **FactorForge v3.6.8**, verified at commit `8e7473298f07ddb7465267c28a4a2651e84e338b` on `main` (following the version 3.6.8 commit `0d4c0c5411781d89863b7ffb73e664bd9ea33979`). Package releases and algorithmic engines are versioned separately: the stable profile engine is v1.0.0, the GC-constrained `1DP` objective serves as the baseline sequence optimizer, and the stable automaton-constrained dynamic programming (`DP v2`) engine is v2.0.1. Primary host support targets *Nicotiana benthamiana*, with Tobacco BY-2 included as an experimental host profile. All candidates and scores produced by FactorForge are computational outputs; they do not establish in vivo expression, yield, folding, or biological activity without experimental validation.

# Statement of need

In synthetic biology and molecular pharming, teams often need to trace *why* a synthetic CDS was selected—not just store the final FASTA string. Changing the host codon table, target GC window, restriction enzyme set, or optimizer can substantially shift both the output sequence and its QC metrics [@sharp1987; @gustafsson2004]. Researchers also need a clear boundary between deterministic sequence constraints (strict amino-acid preservation and forbidden-motif removal) and untested biological hypotheses.

FactorForge separates synonymous sequence generation from independent sequence review, exposing both through a CLI and Python API that emit machine-readable JSON records [@wilkinson2016].

| Capability / Feature | JCat / OPTIMIZER [@grote2005; @puigbo2007] | GeneDesign [@richardson2010] | DNA Chisel [@zulkower2020] | FactorForge |
| :--- | :--- | :--- | :--- | :--- |
| **Primary Domain** | General microbial / heterologous | General synthetic biology | General sequence optimizer | Plant molecular pharming (*N. benthamiana*, BY-2) |
| **Host References** | Built-in / custom user tables | Built-in organism tables | User-supplied usage tables | Versioned plant host tables (`C1`--`C4`, SGN QLD183, NbeV1.1) |
| **Type IIS Restriction** | Site avoidance / post-hoc checks | Restriction site removal | Local search pattern removal | Automaton-constrained DP state space (Aho--Corasick [@ahocorasick1975]) |
| **Optimization Model** | Frequency / Monte Carlo matching | Rule-based / greedy editing | Local stochastic hill-climbing | Exact DP over configured scores and reachable states |
| **Audit & Provenance** | Plain text / web display | GenBank / FASTA export | Annotated GenBank export | Structured JSON manifests (with optional SHA-256 hashes) |

*Qualitative comparison based on published tool documentation; no cross-tool empirical benchmark or superiority claim is made [@kim2026factorforge].*

General-purpose optimizers such as DNA Chisel [@zulkower2020] resolve sequence constraints via local stochastic search without global optimality or feasibility guarantees. FactorForge instead models codon scoring, global GC bounds, and restriction-site avoidance as an exact dynamic program over reachable states—exact with respect to the configured codon scores and implemented constraints, rather than biological optimality—while bundling curated plant codon tables and structured provenance outputs in a single package.

# State of the field

Synonymous codon choice affects translation elongation, mRNA stability, and recombinant protein accumulation [@gustafsson2004; @sharp1987]. Web utilities such as JCat [@grote2005] and OPTIMIZER [@puigbo2007] compute host-adapted sequences, while toolkits including GeneDesign [@richardson2010], D-Tailor [@guimaraes2014], and DNA Chisel [@zulkower2020] support broader constraint-based editing. For *N. benthamiana* specifically, Tekle [@tekle2025] investigated autoregressive transformer models trained on genomic and transcriptomic features.

FactorForge focuses on reproducible, constraint-exact CDS compilation tailored to Type IIS Golden Gate and Modular Cloning (MoClo) workflows [@engler2008; @weber2011]. Its stable `DP v2` engine (referred to as Type IIS-aware `3DP` in the accompanying preprint [@kim2026factorforge]) tracks codon position, cumulative GC count, and Aho--Corasick automaton state [@ahocorasick1975], pruning forbidden restriction sites (e.g., BsaI, BpiI, BsmBI) across codon junctions and user-supplied assembly flanks during search.

The accompanying preprint (v2) [@kim2026factorforge] evaluates 49,257 SGN QLD183 records [@ranawaka2023], including a $2\times 2$ factorial run (FactorForge v3.6.4) crossing codon reference (legacy composite vs. NbeV1.1 [@kurotani2023]) with target GC window (55--65% vs. 40--47%). For example, the Assembly-friendly preset's combined pass rate rose from 0.00% to 9.39% when the window was aligned with NbeV1.1. These runs used the GC-constrained `1DP` objective (tracking codon position and cumulative GC count only) with Type IIS status scored after generation (no filtering); they are computational pass rates, not expression evidence, and do not benchmark the automaton-constrained `DP v2`, whose proteome-scale evaluation remains future work.

# Software design

FactorForge requires Python $\ge 3.10$, is licensed under `AGPL-3.0-only`, and is distributed on PyPI as `factorforge-cds` (source: `eijex/factorforge-cds`, import namespace: `factorforge`). The package is organized into four components:

1. **Candidate generation.** A rule-based profile engine (v1.0.0) supports seeded stochastic sampling, while `DP v2` (v2.0.1) optimizes additive codon scores subject to integer GC bounds and Aho--Corasick forbidden-motif tracking [@ahocorasick1975], explicitly reporting infeasibility and distinguishing exact feasible solutions from fallback outputs.
2. **Engine scope.** The profile engine (v1.0.0), `1DP` (GC-count state only), and `DP v2` (v2.0.1; adds the Aho--Corasick automaton state) are the stable paths described here. Research-candidate or experimental paths, such as `DP v2.1.1` (initiation-region composition and homopolymer penalties), are outside the scope of this paper.
3. **Sequence review.** Independent validators check amino-acid identity, GC windows, Type IIS sites, homopolymer runs, cryptic splice-like motifs, and plant RNA instability signals, with optional 5' folding energy checks via ViennaRNA [@lorenz2011].
4. **Provenance manifests.** CLI and Python entrypoints export JSON reports capturing configuration parameters, sequence metrics, and software and codon-table versions, with schema fields and SHA-256 hashes varying by entrypoint path.

# Research impact statement

FactorForge provides a reproducible baseline for plant CDS design and pre-synthesis QC. Benchmark outputs, runner scripts, and checksums from the accompanying preprint are archived in the repository and Zenodo deposits; full replay requires pinned software versions and the public SGN input corpus:

- **Automated tests and CI.** The test suite covers profile generation, automaton-constrained DP, sequence validation, and CLI workflows. On `main` commit `8e7473298f07ddb7465267c28a4a2651e84e338b` (v3.6.8 baseline), the 8-job GitHub Actions CI matrix across Python 3.10--3.13 on Ubuntu and Windows passed in [workflow run #43440798315](https://github.com/eijex/factorforge-cds/actions/runs/43440798315) (13 November 2026). Mathematical verification tests in `tests/test_math_v2_dp_optimality.py` (6 brute-force parity and flank-junction checks) and `tests/test_dp_hoisting_equivalence.py` (54 uncached-vs-hoisted equivalence checks) check `DP v2` outputs against brute-force enumeration on small instances.
- **Worked example.** An executable [sfGFP worked example](https://github.com/eijex/factorforge-cds/tree/8e7473298f07ddb7465267c28a4a2651e84e338b/examples/worked_example) demonstrates fixed-seed candidate generation and JSON manifest export.
- **Persistent archival.** FactorForge v3.5.11 is archived on Zenodo at [doi:10.5281/zenodo.23631740](https://doi.org/10.5281/zenodo.23631740), and benchmark datasets are preserved at [doi:10.5281/zenodo.20676276](https://doi.org/10.5281/zenodo.20676276) and [doi:10.5281/zenodo.22074106](https://doi.org/10.5281/zenodo.22074106). Baseline v3.6.8 will receive a dedicated archival DOI upon completion of JOSS review.

While developed within a translational plant-pharming pipeline rather than for broad general-purpose adoption, FactorForge offers a tested, auditable foundation for plant CDS design.

# AI usage disclosure

Generative AI tools, including OpenAI ChatGPT/Codex, Google Gemini, and Anthropic Claude, assisted with research planning, exploration of mathematical and algorithmic formulations, software implementation, test development, documentation, and manuscript drafting and revision. The author retained responsibility for the research objectives and core architectural decisions, reviewed and modified AI-assisted contributions, executed relevant software tests and reproducibility checks, and evaluated mathematical claims and references against available evidence. AI-generated outputs were not treated as independent scientific evidence. The author takes full responsibility for the accuracy, originality, and integrity of the software and manuscript.

# Acknowledgements

The author thanks Doug Cossar (PlantForm Corporation) for discussions on plant molecular pharming workflows.

**Conflict of interest:** The author develops and maintains FactorForge through Eijex and is affiliated with PlantForm Corporation; no PlantForm proprietary data were used in the reported analyses. The author holds no FactorForge patents, equity interests, or paid consulting arrangements, and no external grant funding supported this work.

# References
