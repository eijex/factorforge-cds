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
  - name: PlantForm Corporation, Toronto, Ontario, Canada
    index: 1
date: 19 November 2026
bibliography: paper.bib
---

# Summary

FactorForge is an open-source Python package for coding-sequence (CDS) candidate generation and pre-synthesis sequence review in plant molecular pharming. Different synonymous DNA sequences can encode the same protein while differing substantially in codon adaptation, nucleotide composition, and restriction-site content. FactorForge helps researchers generate sequence candidates under user-configured constraints, inspect multi-category sequence quality checks, and retain structured, machine-readable computational review artifacts before DNA synthesis and laboratory testing.

This paper describes **FactorForge v3.6.4**, identified by its release commit `ab2972436dfad1aa26cec6d9365308769928f0d9`. Package versions and algorithmic engine versions are strictly separated: the stable profile engine is version 1.0.0, and the stable dynamic programming (DP v2) engine is version 2.0.1. Primary host support focuses on *Nicotiana benthamiana*, with Tobacco BY-2 supported as an experimental host context. All generated sequence candidates and scores are computational outputs; they do not constitute empirical evidence of protein expression, yield, folding, or in vivo biological performance.

# Statement of need

Synthetic biology and molecular pharming workflows require auditable provenance explaining why a specific coding sequence was chosen, rather than merely recording its final nucleotide sequence. Altering the reference codon table, target GC policy, restriction enzyme set, or underlying optimization algorithm can markedly shift both the resulting nucleotide sequence and its computational scores [@sharp1987; @gustafsson2004]. Furthermore, translational researchers must distinguish strict translational preservation and configured motif avoidance from untested biological interpretations.

FactorForge addresses these challenges by decoupling synonymous sequence generation from independent multi-domain sequence review, packaging both within scriptable CLI and Python APIs that produce FAIR-oriented JSON artifacts [@wilkinson2016].

| Capability / Feature | JCat / OPTIMIZER [@grote2005; @puigbo2007] | GeneDesign [@richardson2010] | DNA Chisel [@zulkower2020] | FactorForge |
| :--- | :--- | :--- | :--- | :--- |
| **Primary Domain** | General microbial / heterologous | General synthetic biology | General sequence optimizer | Plant molecular pharming (*N. benthamiana*, BY-2) |
| **Host References** | Built-in / custom user tables | Static built-in tables | User-supplied usage tables | Versioned plant host tables (`C1`--`C4`, SGN QLD183, NbeV1.1) |
| **Type IIS Restriction** | Post-generation filtration | Restriction site removal | Local search pattern removal | Automaton-constrained DP state space (Aho--Corasick [@ahocorasick1975]) |
| **Optimization Model** | Greedy frequency matching | Greedy / rule-based editing | Local stochastic hill-climbing | Exact dynamic programming over reachable states |
| **Audit & Provenance** | Plain text / web display | GenBank / FASTA export | Annotated GenBank export | Machine-readable JSON with SHA-256 provenance hashes |

Whereas general-purpose sequence optimizers such as DNA Chisel [@zulkower2020] resolve sequence conflicts through local stochastic search heuristics without global optimality guarantees, FactorForge formulates codon selection, global GC boundaries, and restriction-site avoidance as an exact dynamic program over reachable state tuples. A standalone, dedicated package is required to enforce this exact state formulation, bundle curated plant host references, and emit structured, reproducible provenance records across automated synthesis pipelines.

# State of the field

Synonymous codon usage influences translation kinetics, transcript stability, and heterologous protein accumulation, though computational adaptation metrics do not by themselves guarantee improved expression [@gustafsson2004; @sharp1987]. Classic web-based utilities such as JCat [@grote2005] and OPTIMIZER [@puigbo2007] calculate host-directed adaptation indexes. Comprehensive synthetic biology toolkits such as GeneDesign [@richardson2010], D-Tailor [@guimaraes2014], and DNA Chisel [@zulkower2020] provide broader sequence design workflows. In plant-specific sequence research, Tekle [@tekle2025] explored transformer-based autoregressive codon prediction integrating genomic and transcriptomic features in *N. benthamiana*.

FactorForge provides a plant-specific, versioned CDS design workflow in which codon-reference provenance, constrained generation, and sequence-review evidence are integrated as explicit research artifacts. This addresses workflow reproducibility requirements that are not provided as a unified, plant-focused application by general-purpose sequence optimization frameworks. The package's core contribution is a multi-constraint framework tailored for Type IIS Golden Gate and Modular Cloning (MoClo) assembly standards [@engler2008; @weber2011]. Its stable DP v2 engine integrates global GC-count bounds with an Aho--Corasick motif automaton [@ahocorasick1975], pruning forbidden restriction sites (such as BsaI, BpiI, and BsmBI) across codon junctions and supplied assembly flanks directly during search. Exactness concerns reachability within the formulated state space; it does not imply biological optimality.

The accompanying preprint (v2) [@kim2026factorforge] reports a v3.2.0 preset benchmark, a matched v3.3.0 DP execution on 49,257 reference-derived SGN QLD183 records [@ranawaka2023], and a v3.6.4 $2\times 2$ factorial evaluation crossing codon reference (legacy composite vs. NbeV1.1 high-confidence [@kurotani2023]) with target GC window (55--65% vs. 40--47%). Outcomes depend jointly on reference and GC policy (e.g., Assembly-friendly combined pass rate shifts from 0.00% to 9.39% when the GC window is aligned with the native NbeV1.1 reference). These are computational pass rates on public sequences, obtained with 1DP (unconstrained Type IIS state space with post-generation filtering) rather than the automaton-constrained DP v2 described here, and are not evidence of expression. Because population-scale factorial benchmarks were executed under 1DP, evaluating the motif-automaton DP v2 at proteome scale remains future work; the algorithmic correctness of the DP v2 state space is independently validated against exhaustive brute-force sequence enumeration on short coding sequences (Preprint v2, Appendix A.3).

# Software design

FactorForge requires Python $\ge 3.10$ and is licensed under the GNU Affero General Public License v3.0 only (`AGPL-3.0-only`). The package is published on PyPI as `factorforge-cds` and maintained on GitHub (`eijex/factorforge-cds`). The primary import namespace is `factorforge`.

The architecture encompasses four functional domains:

1. **Candidate generation.** The stable profile engine provides rule-based sequence generation, including stochastic sampling parameterized by explicit random seeds. The stable DP v2 engine integrates additive codon scores (such as CAI), integer GC bounds, and Aho--Corasick automaton tracking [@ahocorasick1975]. This state representation ensures exact mathematical optimization within the defined model, explicitly reporting when target bounds are unreachable and distinguishing feasible sequences from Eq. (19) fallback outputs.
2. **Research engines vs. stable baseline.** A runtime registry isolates production-stable engines from experimental research candidates. DP v2.1.1 introduces initiation-region composition and homopolymer penalty policies, while newer 3DP (automaton-constrained joint optimization) and 4DP (exit-tunnel translation pause modeling) engines remain experimental research paths.
3. **Multi-domain sequence review.** Independent sequence evaluation modules inspect translation fidelity, GC composition, Type IIS restriction recognition sites, and sequence hazard scans (such as homopolymer runs, splice-site mimics, and plant RNA instability motifs). Optional initiation-region secondary structure folding energy evaluation is supported via ViennaRNA [@lorenz2011].
4. **FAIR provenance and reporting.** CLI and Python interfaces produce structured JSON review manifests recording user configuration, sequence metrics, software and reference table versions, and SHA-256 hashes of input and output artifacts.

# Research impact statement

FactorForge provides reproducible sequence design workflows for plant synthetic biology and molecular pharming. All computational analyses reported in the accompanying preprint can be replayed from public genomic sequence corpora using the provided CLI tools and documentation:

- **Automated test suite and CI.** The full test suite covers profile generation, automaton-constrained dynamic programming, sequence validation, and CLI entrypoints. On release commit `ab2972436dfad1aa26cec6d9365308769928f0d9` (v3.6.4), the complete GitHub Actions CI matrix across Python 3.10, 3.11, 3.12, and 3.13 on Ubuntu and Windows passed all 8 matrix test jobs in [workflow run #43018092835](https://github.com/eijex/factorforge-cds/actions/runs/43018092835) (11 November 2026). In addition, 42 focused mathematical tests verify exact parity between DP v2 reachability and exhaustive brute-force enumeration (Preprint v2, Appendix A.3).
- **Worked examples and reproducibility.** An executable [sfGFP worked example](https://github.com/eijex/factorforge-cds/tree/v3.6.4/examples/worked_example) demonstrates fixed-seed candidate generation and structured JSON review artifact production under documented constraints.
- **Persistent archival.** FactorForge v3.5.11 is archived on Zenodo at [doi:10.5281/zenodo.23631740](https://doi.org/10.5281/zenodo.23631740). Historical benchmark datasets and matched v3.3.0 DP evaluation outputs are preserved at [doi:10.5281/zenodo.20676276](https://doi.org/10.5281/zenodo.20676276) and [doi:10.5281/zenodo.22074106](https://doi.org/10.5281/zenodo.22074106). The audited submission baseline v3.6.4 will receive a dedicated post-review archival DOI upon completion of the JOSS review process.

FactorForge has been developed and evaluated primarily within our translational research pipeline. While broad external community adoption is not claimed, the package provides researchers with an open, rigorously tested foundation for plant CDS generation.

# AI usage disclosure

Generative AI assistance tools (Google Gemini 1.5/2.0, OpenAI GPT-4o, Anthropic Claude 3.5 Sonnet) were used during the software development and manuscript preparation:
1. **Source code and tests:** Generating initial edge-case unit test templates and type annotations in `src/factorforge/` and `tests/`.
2. **Documentation and scripts:** Assisting with CLI documentation drafting and benchmark replay script formatting.
3. **Manuscript preparation:** Language editing, LaTeX/Markdown formatting, and cross-reference checks.

All AI-suggested code, test fixtures, and text were independently reviewed, executed, and verified by the author. No autonomous scientific claims or data fabrication were generated by AI, and the author takes full intellectual responsibility for all software and manuscript contents.

# Acknowledgements

The author conceived the software architecture, developed the implementation, conducted the computational benchmarks, and wrote the manuscript. The author thanks Doug Cossar (PlantForm Corporation) for discussions on plant molecular pharming workflows, and acknowledges exploratory bioinformatics work by Thomas Tekle [@tekle2025] at PlantForm Corporation that provided valuable context for plant codon modeling.

**Conflict of interest:** The author develops and maintains FactorForge through Eijex and is affiliated with PlantForm Corporation; no PlantForm proprietary data were used in the reported analyses. The author reports no FactorForge patent, equity interest, or paid commercial consulting service. No dedicated external grant funding supported this work.

# References
