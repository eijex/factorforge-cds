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

FactorForge is an open-source Python package for synonymous codon optimization and coding-sequence (CDS) review in plant molecular pharming. During translation, ribosomes read messenger RNA (mRNA) in codons, and the choice among synonymous codons can influence translation efficiency in a given host. FactorForge generates candidate DNA sequences by substituting synonymous codons while preserving the amino acid sequence encoded by the input CDS. It records sequence metrics and checks in structured JSON. The review reports the codon adaptation index (CAI) relative to a host reference, GC content, and Type IIS restriction-site status. An optional RNA folding check estimates minimum free energy (MFE) near the 5' end.

This paper describes the stable dynamic programming engine (`DP v2`) for *Nicotiana benthamiana* in **FactorForge v3.6.8**. All candidates and scores produced by FactorForge are computational outputs; they do not establish in vivo expression, yield, folding, or biological activity without experimental validation.

# Statement of need

In synthetic biology and molecular pharming, teams often need to trace *why* a synthetic CDS was selected—not just store the final FASTA string. Changing the host codon table, target GC window, restriction enzyme set, or optimizer can substantially shift both the output sequence and its QC metrics [@sharp1987; @gustafsson2004]. Researchers also need a clear boundary between deterministic sequence constraints (strict amino-acid preservation and forbidden-motif removal) and untested biological hypotheses.

FactorForge separates synonymous sequence generation from independent sequence review, exposing both through a CLI and Python API that emit machine-readable JSON records [@wilkinson2016].

| Capability / Feature | JCat / OPTIMIZER [@grote2005; @puigbo2007] | GeneDesign [@richardson2010] | DNA Chisel [@zulkower2020] | FactorForge |
| :--- | :--- | :--- | :--- | :--- |
| **Primary Domain** | General microbial / heterologous | General synthetic biology | General sequence optimizer | Plant molecular pharming (*N. benthamiana*) |
| **Host References** | Built-in / custom user tables | Built-in organism tables | User-supplied usage tables | Versioned plant host tables (`C1`--`C4`, SGN QLD183, NbeV1.1) |
| **Type IIS Restriction** | Site avoidance / post-hoc checks | Restriction site removal | Local search pattern removal | Automaton-constrained DP state space (Aho--Corasick [@ahocorasick1975]) |
| **Optimization Model** | Frequency / Monte Carlo matching | Rule-based / greedy editing | Local stochastic hill-climbing | Exact DP over configured scores and reachable states |
| **Audit & Provenance** | Plain text / web display | GenBank / FASTA export | Annotated GenBank export | Structured JSON manifests (with optional SHA-256 hashes) |

*Qualitative comparison based on published tool documentation; no cross-tool empirical benchmark or superiority claim is made [@kim2026factorforge].*

General-purpose optimizers such as DNA Chisel [@zulkower2020] resolve sequence constraints via local stochastic search without global optimality or feasibility guarantees. FactorForge instead models codon scoring, global GC bounds, and restriction-site avoidance as an exact dynamic program over reachable states—exact with respect to the configured codon scores and implemented constraints, rather than biological optimality—while bundling curated plant codon tables and structured provenance outputs in a single package.

# State of the field

Synonymous codon choice affects translation elongation, mRNA stability, and recombinant protein accumulation [@gustafsson2004; @sharp1987]. Web utilities such as JCat [@grote2005] and OPTIMIZER [@puigbo2007] compute host-adapted sequences, while toolkits including GeneDesign [@richardson2010], D-Tailor [@guimaraes2014], and DNA Chisel [@zulkower2020] support broader constraint-based editing. For *N. benthamiana* specifically, Tekle [@tekle2025] investigated autoregressive transformer models trained on genomic and transcriptomic features.

FactorForge focuses on reproducible, constraint-exact CDS compilation tailored to Type IIS Golden Gate and Modular Cloning (MoClo) workflows [@engler2008; @weber2011]. Its stable `DP v2` engine tracks codon position, cumulative GC count, and Aho--Corasick automaton state [@ahocorasick1975], pruning forbidden restriction sites (e.g., BsaI, BpiI, BsmBI) across codon junctions and user-supplied assembly flanks during search.

The accompanying preprint (v2) [@kim2026factorforge] reports a $2\times 2$ factorial benchmark run with FactorForge v3.6.4 on 49,257 SGN QLD183 records [@ranawaka2023]. It compared two codon references (legacy composite and NbeV1.1 [@kurotani2023]) and two target GC windows (55--65% and 40--47%). Combined pass required sequence integrity, GC-window compliance, and absence of the specified Type IIS sites. Mean CAI was reported separately against the reference used in each condition. The factorial comparison did not report 5' RNA folding MFE.

With NbeV1.1 held fixed, the Assembly-friendly preset's combined pass rate increased from 0.00% to 9.39% when the target GC window changed from 55--65% to 40--47%. A separate GC-constrained baseline in the same benchmark was assessed for Type IIS sites after candidate generation. Proteome-scale evaluation of `DP v2` remains future work.

# Software design

FactorForge requires Python $\ge 3.10$, is licensed under `AGPL-3.0-only`, and is distributed on PyPI as `factorforge-cds` (source: `eijex/factorforge-cds`, import namespace: `factorforge`). The package is also available through [Bioconda](https://anaconda.org/bioconda/factorforge-cds), which currently distributes version 3.6.4. The software has three main components:

1. **Candidate generation.** `DP v2` optimizes additive codon scores subject to integer GC bounds and Aho--Corasick forbidden-motif tracking [@ahocorasick1975], explicitly reporting infeasibility and distinguishing exact feasible solutions from fallback outputs.
2. **Sequence review.** The review reports CAI as a measure of codon usage relative to the selected host reference. Independent validators check amino-acid identity, compliance with the specified GC window, Type IIS sites, homopolymer runs, cryptic splice-like motifs, and plant RNA instability signals. An optional ViennaRNA check estimates local 5' RNA folding MFE [@lorenz2011].
3. **Provenance manifests.** CLI and Python entrypoints export JSON reports capturing configuration parameters, sequence metrics, and software and codon-table versions, with schema fields and SHA-256 hashes varying by entrypoint path.

# Research impact statement

FactorForge provides a reproducible baseline for plant CDS design and pre-synthesis QC. Benchmark outputs, runner scripts, and checksums from the accompanying preprint are archived in the repository and Zenodo deposits; full replay requires pinned software versions and the public SGN input corpus:

- **Automated tests and CI.** The test suite covers profile generation, automaton-constrained DP, sequence validation, and CLI workflows. On `main` commit `8e7473298f07ddb7465267c28a4a2651e84e338b` (a development snapshot reporting version 3.6.8), the 8-job GitHub Actions CI matrix across Python 3.10--3.13 on Ubuntu and Windows passed in [workflow run #43440798315](https://github.com/eijex/factorforge-cds/actions/runs/43440798315) (13 November 2026). Tests in `tests/test_math_v2_dp_optimality.py` cover brute-force parity on small instances and flank junctions. Tests in `tests/test_dp_hoisting_equivalence.py` compare uncached and hoisted implementations.
- **Worked example.** An executable [sfGFP worked example](https://github.com/eijex/factorforge-cds/tree/8e7473298f07ddb7465267c28a4a2651e84e338b/examples/worked_example) demonstrates fixed-seed candidate generation and JSON manifest export.
- **Persistent archival.** FactorForge v3.5.11 is archived on Zenodo at [doi:10.5281/zenodo.23631740](https://doi.org/10.5281/zenodo.23631740), and benchmark datasets are preserved at [doi:10.5281/zenodo.20676276](https://doi.org/10.5281/zenodo.20676276) and [doi:10.5281/zenodo.22074106](https://doi.org/10.5281/zenodo.22074106). A dedicated software archive for v3.6.8 is planned after JOSS review.

While developed within a translational plant-pharming pipeline rather than for broad general-purpose adoption, FactorForge offers a tested, auditable foundation for plant CDS design.

# AI usage disclosure

Generative AI tools, including OpenAI ChatGPT/Codex, Google Gemini, and Anthropic Claude, assisted with research planning, exploration of mathematical and algorithmic formulations, software implementation, test development, documentation, and manuscript drafting and revision. The author defined the research objectives and made the core architectural decisions. The author personally reviewed, edited, and validated all AI-assisted contributions to the software, documentation, and manuscript. The author also ran the relevant software tests and reproducibility checks and checked the mathematical claims and references against the available evidence. AI-generated outputs were not treated as independent scientific evidence. The author takes full responsibility for the accuracy, originality, and integrity of the software and manuscript.

# Acknowledgements

The author thanks Doug Cossar, PhD (PlantForm Corporation), for discussions on plant molecular pharming workflows.

**Conflict of interest:** The author develops and maintains FactorForge through Eijex and is affiliated with PlantForm Corporation; no PlantForm proprietary data were used in the reported analyses. The author holds no FactorForge patents, equity interests, or paid consulting arrangements, and no external grant funding supported this work.

# References
