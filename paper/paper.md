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
date: 30 October 2026
bibliography: paper.bib
---

# Summary

FactorForge is an open-source Python package for coding-sequence (CDS) candidate generation and pre-synthesis sequence review in plant molecular pharming. Different synonymous DNA sequences can encode the same protein while differing in codon adaptation, nucleotide composition, and restriction-site content. FactorForge helps researchers generate candidates, inspect configured sequence checks, and retain computational review artifacts before synthesis and laboratory testing.

This paper describes **FactorForge v3.5.11**, identified by its tagged release and commit `7e54fae4bd3a3e97f6d61d080d6f600ac17c5ece`. Package versions and engine versions are distinct: the stable profile engine is 1.0.0 and the stable DP v2 engine is 2.0.1. Primary host support is *Nicotiana benthamiana*; Tobacco BY-2 is an experimental host context. Generated candidates and scores are computational outputs, not evidence of expression, yield, folding, or other biological performance.

# Statement of need

Sequence-design workflows need to preserve why a candidate was selected, not merely its final nucleotide string. Changing the codon reference, objective, constraints, or software can change both the sequence and its computational score. Researchers also need to distinguish translation preservation and configured motif checks from untested biological interpretations.

FactorForge addresses these needs through scriptable candidate generation, separate sequence review, and structured reporting. Its intended users are plant molecular pharming researchers and computational biologists integrating sequence design into reproducible workflows. Report-producing paths expose configuration, reference and version identifiers, sequence metrics, and available SHA-256 provenance fields. Schemas and enforcement differ across entrypoints, so consumers must verify the recorded contract rather than assume every result contains an identical manifest. These artifacts support FAIR-oriented stewardship without establishing biological validity [@wilkinson2016].

# State of the field

Codon usage is relevant to heterologous expression, but computational adaptation metrics do not by themselves establish improved expression [@gustafsson2004]. JCat [@grote2005] and OPTIMIZER [@puigbo2007] support host-directed codon adaptation. GeneDesign [@richardson2010], D-Tailor [@guimaraes2014], and DNA Chisel [@zulkower2020] provide broader synthetic-sequence design capabilities. DNA Chisel combines local-search-based constraint resolution with objective optimization across diverse applications.

FactorForge's contribution is a plant-focused workflow combining versioned references, constrained synonymous generation, and inspectable review outputs. Rather than replace general-purpose frameworks, it packages these choices for plant CDS preparation and reproducibility. Its stable DP v2 implementation carries GC-count and motif-automaton state, rejecting configured forbidden motifs during search across codon junctions and supplied flanks. Exactness concerns the implemented reachable state space and supplied additive codon scores; it is not biological optimality. Constraint and fallback contracts remain execution-path specific. No cross-tool empirical superiority is claimed.

The published accompanying preprint [@kim2026factorforge] provides historical evidence: a v3.2.0 preset benchmark and a separate v3.3.0 DP execution on 49,257 reference-derived SGN QLD183 CDS/protein records. The historical DP experiment used a GC-constrained objective and post-generation Type IIS scoring. Those results are not a benchmark of v3.5.11 or its automaton-constrained DP v2 engine.

# Software design

The Python package requires Python $\ge 3.10$ and is licensed `AGPL-3.0-only`. PyPI provides the audited 3.5.11 release as `factorforge-cds`; Bioconda is an additional distribution channel whose recipe version must be checked separately. The import namespace is `factorforge`.

Four functional areas explain the architecture; these are responsibilities rather than a claim that the package contains only four modules:

1. **Candidate generation.** The stable profile pathway supplies rule-based generation, including stochastic choices controlled by seeds. Stable DP v2 combines codon scoring, global GC bounds, and configured motif exclusion. Its implemented state representation supports exact optimization within that model, while reported reachability distinguishes the requested GC interval from fallback outcomes. The trade-off is greater state and memory cost than independent per-codon selection.
2. **Explicit research paths.** The registry separates stable engines from research candidates and feature-gated previews. DP v2.1.1 adds initiation-region composition and homopolymer policies. Other local DP candidates, learned-generation previews, and adaptive partial-DP rescue are not interchangeable with the stable baseline. A rescued suffix is conditioned on its retained prefix and constraint state; it does not establish an unrestricted global optimum or biological benefit.
3. **Sequence review and optional folding evidence.** Profile checks and `factorforge.validation` report translation integrity, GC metrics, Type IIS status, and advisory findings according to the selected path. Supplied upstream/downstream context can condition relevant calculations; the DP classes described here vary CDS codons rather than jointly design a 5′ UTR. ViennaRNA [@lorenz2011] is optional. The DP v2.1.1 folding-evidence path uses a 60-nt initiation window comprising 15 upstream nucleotides and 45 CDS nucleotides, with explicit unavailable/incomplete-context status. Separately evaluated folding energy must not be confused with a proxy objective or experimental measurement.
4. **Workflow interfaces and provenance.** CLI and Python interfaces expose scriptable generation and reporting. JSON outputs separate computational metrics, configured checks, warnings, and version/provenance fields. Hashes support artifact comparison when the relevant fields are present; they do not certify scientific correctness. Public web and separate MCP access are independently deployed interfaces whose runtime versions require their own checks.

This separation makes design choices inspectable while retaining a clear distinction between stable functionality, research implementations, and biological evidence. The [tagged source and documentation](https://github.com/eijex/factorforge-cds/tree/v3.5.11) define the audited release; moving package indexes or web deployments do not redefine the historical experiments.

# Research impact statement

The published preprint supplies evidence of research use through its 49,257-record computational analysis and codon-reference sensitivity study [@kim2026factorforge]. These are reference-derived sequence evaluations, not wet-lab or therapeutic outcomes. Historical benchmark artifacts and executable examples provide concrete materials for reproducing that use without implying independent external adoption.

- **Tests and CI.** Tag-specific collection enumerated 1,229 tests in the audited Python 3.11 environment. Collection is not execution: the focused audit ran 20 relevant smoke, DP, rescue, and registry checks successfully. The CI workflow defines Python 3.10–3.13 on Linux and Windows, but its [recorded v3.5.11 run](https://github.com/eijex/factorforge-cds/actions/runs/40154207184) was not successful. A passing release-wide CI result remains a submission-readiness requirement.
- **Worked example.** The archived [sfGFP workflow](https://github.com/eijex/factorforge-cds/tree/v3.3.0/examples/worked_example) illustrates fixed-seed design and JSON review artifacts under its historical software boundary. It is not a v3.5.11 biological-validation result.
- **Persistent archival.** The v3.5.11 software release is archived at [doi:10.5281/zenodo.23631740](https://doi.org/10.5281/zenodo.23631740). Historical benchmark data and the matched v3.3.0 DP execution have separate archives at [doi:10.5281/zenodo.20676276](https://doi.org/10.5281/zenodo.20676276) and [doi:10.5281/zenodo.22074106](https://doi.org/10.5281/zenodo.22074106). Any software revision accepted after JOSS review will require a matching final archival record.

# AI usage disclosure

Generative AI tools, including OpenAI ChatGPT/Codex, Google Gemini, and Anthropic Claude, supported planning, drafting, logical structure, mathematical and algorithmic formulation, software implementation, language revision, and formatting; their role extended beyond language editing. AI-generated content is not independent scientific evidence. All AI-assisted outputs were reviewed and verified by the author. The author takes full intellectual responsibility for the software, manuscript, equations, code, data, references, and final wording.

# Acknowledgements

The author conceived the software architecture, developed the implementation, conducted the computational benchmarks, and wrote the manuscript. The author thanks Doug Cossar (PlantForm Corporation) for discussions of plant molecular pharming workflows. No dedicated external grant or institutional funding supported this work.

# References
