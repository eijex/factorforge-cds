# FactorForge

**Open-source constraint-based CDS design and pre-synthesis sequence review for plant CDS workflows, with primary support for *Nicotiana benthamiana* (Tobacco BY-2: experimental).**

[![License](https://img.shields.io/badge/license-AGPL--3.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![PyPI](https://img.shields.io/pypi/v/factorforge-cds.svg)](https://pypi.org/project/factorforge-cds/)
[![Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/factorforge-cds)
[![CI](https://github.com/eijex/factorforge-cds/actions/workflows/ci.yml/badge.svg)](https://github.com/eijex/factorforge-cds/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/eijex/factorforge-cds/branch/main/graph/badge.svg)](https://codecov.io/gh/eijex/factorforge-cds)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20407330.svg)](https://doi.org/10.5281/zenodo.20407330)
[![Web App](https://img.shields.io/badge/web-factorforge.eijex.com-brightgreen.svg)](https://factorforge.eijex.com)
[![Buy Me A Coffee](https://img.shields.io/badge/Buy%20Me%20a%20Coffee-Donate-yellow?style=flat&logo=buy-me-a-coffee)](https://www.buymeacoffee.com/eijex)


FactorForge performs deterministic CDS design with CAI/GC metrics, PolyA-signal screening, and Golden Gate/MoClo-aware checks. It is positioned as a pre-synthesis review harness: it helps teams generate reproducible CDS candidates, inspect assembly-relevant sequence constraints, and package design metadata before downstream synthesis, cloning, or experimental review. Primary support: *N. benthamiana* (agroinfiltration). Experimental host context: Tobacco BY-2 (`--host by2`).

FactorForge v3.6.8 uses independently versioned engines:

| Generation | Engine | Version | Availability |
|---|---|---:|---|
| Gen 1 | Rule/profile | 1.0.0 | Stable, public |
| Gen 2 | DP v2 | 2.0.1 | Stable, public default feasibility path |
| Gen 2 | DP v2.1.1 | 2.1.1 | Explicit local-guard path; computational evidence only |
| Gen 3 | sLLM Hybrid | 0.2.0-preview.1 | Feature-gated constrained-generation research preview |
| Rescue | Adaptive partial DP | 1.0.0 | Exact suffix rescue conditioned on a verified prefix |

**→ [Full Documentation](https://eijex.github.io/factorforge-cds/)** · **[Roadmap](https://eijex.github.io/factorforge-cds/roadmap/)**

---

## Quick Start

```bash
# Bioconda: Linux/macOS, or Ubuntu inside WSL on Windows
conda create -n factorforge-review --override-channels -c conda-forge -c bioconda --strict-channel-priority python=3.11 factorforge-cds
conda activate factorforge-review

# Via PyPI
pip install factorforge-cds

# Run optimization
factorforge optimize my_protein.fasta -o output.fasta
```

Or use the **[web app](https://factorforge.eijex.com)** — no installation required.

Bioconda does not support native Windows. Use WSL for that installation route;
see [installation and reviewer checks](docs/getting-started.md) for the tested
version, platform limitations and installed-package smoke command. Local DP
optimization needs no Neon database, GCP account or GPU.

---

## Access Options

| Method | Description | Link |
|--------|-------------|------|
| **Bioconda** | Linux/macOS; WSL on Windows | [Isolated installation](docs/getting-started.md#via-bioconda-recommended-for-bioinformatics) |
| **Web App** | No installation, demo & light use | [factorforge.eijex.com](https://factorforge.eijex.com) |
| **CLI / Python** | Local use, batch processing, data privacy | `pip install factorforge-cds` |
| **Docker** | Full web interface locally | `docker pull ghcr.io/eijex/factorforge-cds:latest` |
| **Eijex MCP** | MCP-compatible agent access | [mcp.eijex.com](https://mcp.eijex.com) |

---

## Repository Structure

Experimental codon-distribution and evidence-ledger modules in the development
checkout are scaffolds, not a validated laboratory policy or active learned
recommendation service. Computational checks do not establish synthesis readiness
or biological performance. Private inputs and collaborator packages must remain
outside public repositories and dashboards.

The supported deterministic engines are the profile engine, stable DP v2, and
the explicit DP v2.1.1 development candidate under:

```text
src/factorforge/engines/profile/
src/factorforge/engines/dp_v2.py
src/factorforge/engines/dp_v2_1_1.py
src/factorforge/engines/sllm/
src/factorforge/discovery/
```

DP v2.1.1 adds an exact active-layer 5′ GC guard and Aho-Corasick rejection of
homopolymers of 6 nt or longer to the v2.1 initiation-aware objective. It emits
local-composition metrics and, when ViennaRNA and sufficient transcript context
are available, a separately evaluated 5′ MFE value. It is not the default; its
single-target calibration does not establish holdout generalization or biological
performance.

The v3.5.x discovery-slate surface generates versioned Top-K research candidates,
applies a shared deterministic hard-constraint filter, and records generator and
fallback lineage. The sLLM path is disabled by default. A partial-DP rescue solves
an exact suffix conditioned on the retained prefix; it is not a claim of global
optimality and does not establish biological performance.

Historical implementation tracks are preserved under `archive/` for provenance
and are not imported by the installed package or exposed as supported engines.

---

## ⚠️ Validation Status

FactorForge outputs are **in-silico only** and have not been experimentally validated in wet-lab conditions. These checks support reviewability and reproducibility; they do not guarantee expression, yield, synthesis acceptance, folding, glycosylation, regulatory approval, or downstream biological performance. See [Validation](https://eijex.github.io/factorforge-cds/validation/) and [VALIDATION.md](VALIDATION.md).

---

## Citing

```
FactorForge v3.6.8 (2026). Open-source constraint-based CDS design and sequence review.
Eijex. https://github.com/eijex/factorforge-cds
```

The [v3.5.11 software archive](https://doi.org/10.5281/zenodo.23631740)
and the [accompanying computational preprint (v2)](https://doi.org/10.20944/preprints202610.3128.v2)
are separate citations. The preprint v2 reports a 2×2 factorial benchmark run with FactorForge v3.6.4 on 49,257 SGN QLD183 records, not a benchmark of every later software release.

---

## Maintainer

Mun-Kyu Kim ([@eijex](https://github.com/eijex))

## Optional PostgreSQL support

FactorForge design, CLI, Python API, and file export do not require a database.
PostgreSQL persistence is an explicit Eijex integration path for retaining shared
campaign and candidate identities; it is not enabled by default.

The extra pins `eijex-db-core==0.1.2`. Supply its approved wheel bundle using
`pip install --find-links PATH_TO_BUNDLE "factorforge-cds[postgres]"`; the shared
package is not assumed to have been published on PyPI. Standalone installation
does not require it.
The deployment must also provide a compatible `eijex-db-core` package and set
`FACTORFORGE_DATABASE_URL`. FactorForge contains no default database credentials and
does not silently fall back to another backend. See [Persistence](docs/persistence.md).

Explicit local SQLite research checkpoints in `factorforge.db.connector` remain a
separate local-only utility and are not the shared DBTL system of record.

## License

GNU Affero General Public License v3.0 — see [LICENSE](LICENSE).

**Disclaimer:** FactorForge is provided for research purposes only. Outputs are computational and have not been experimentally validated.

---

## Get in Touch

- **Sponsor** — Support our research via [GitHub Sponsors](https://github.com/sponsors/eijex) or [PayPal](https://paypal.me/munkyukim)
- **Docs** — [eijex.github.io/factorforge-cds](https://eijex.github.io/factorforge-cds/)
- **Wet-lab Feedback** — Public-safe feedback summaries are welcome via [Share Wet-lab Feedback (GitHub)](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml). Do not submit raw sequences, confidential construct details, internal batch IDs, patient data, private contact information, exact process parameters, or confidential partner/customer data. Email `eijex.lab@gmail.com` for private or sensitive summaries. See [VALIDATION.md](VALIDATION.md) before submitting.
- **Sequence Policy Profiles** — FactorForge applies the sequence-review policy
  chosen by the user; it does not define or approve a laboratory SOP. The web
  app provides an illustrative YAML example and active-profile export.
  Sequence-free, non-confidential suggestions can be submitted through the
  public policy-template form; use `eijex.lab@gmail.com` for private or
  laboratory-specific profiles. See the [sequence policy boundary](docs/sequence-policy.md).
- **GitHub Issues** — bugs, features: [github.com/eijex/factorforge-cds/issues](https://github.com/eijex/factorforge-cds/issues)
- **Email** — eijex.lab@gmail.com
- **FactorForge** — [factorforge.eijex.com](https://factorforge.eijex.com)
- **Eijex MCP** — [mcp.eijex.com](https://mcp.eijex.com)
- **Lab** — [www.eijex.com](https://www.eijex.com)


## Testing status

The web UI distinguishes the 1DP sequence-only baseline from experimental structure-conditioned (`3DP`) and translation-offset (`4DP`) workflows (`Biological validation: NOT_AVAILABLE`) and keeps watermark disabled. See [web UI boundaries](web/README.md).

## Public 1DP, 3DP, and 4DP web access

The public web interface (`factorforge.eijex.com`) enables `1DP` (sequence baseline), `3DP` (prepared exact-match PDB/mmCIF upload up to 450 KB or ESMFold auto-prediction with reviewed positions, up to 500 aa), and `4DP` (schema-validated pause sidecar JSON up to 250 KB with reviewed positions and explicit hypothesis acknowledgment) under FactorForge's built-in default SOP (`FACTORFORGE_EXPERIMENTAL_3DP_ENABLED=true` and `FACTORFORGE_PUBLIC_4DP_ENABLED=true` by default; synthetic 4DP fixtures remain restricted to local loopback testing). Both `3DP` and `4DP` are experimental computational hypotheses (`UNVERIFIED_PUBLIC_CORPUS`, `Biological validation: NOT_AVAILABLE`). See [research boundaries](docs/experimental-structure.md) and [web UI boundaries](web/README.md).

This is stateless computational testing, with no external service dependencies or remote database recording. Structure files are not retained in the response or browser history; normal design results may be retained in local browser history. Submit only public-safe inputs. Structural conformity metrics are codon-weight proxies, not folding/expression measurements.

## Current public policy UI

Custom SOP editing and upload remain disabled in the public web form; FactorForge's built-in default SOP is applied automatically across 1DP, 3DP, and 4DP workflows. Saved custom SOPs are ignored without deletion unless reactivated via `ENABLE_CUSTOM_SOP`.

### Scientific web observations

Direct 1DP/3DP browser testing reports sense-CDS RNA motif scans and available initiation-window MFE in results, downloadable reports and local history. Structure confidence observations can suggest positions for explicit acceptance and review; unsupported inputs do not invent positions. These checks do not establish expression, yield or folding benefits. The default policy remains active and custom SOP management remains paused. See `docs/experimental-structure.md` for the API and measurement boundary.

### Local reviewed structure file workflow

`factorforge nbfold annotate` and `factorforge nbfold optimize` provide the experimental local Phase A/B workflow for exact-match structure files and an explicit user review manifest. See [NB-Fold documentation](docs/experimental-structure.md).

Web structure handling now requires exact-match uploaded PDB/mmCIF coordinates and explicit reviewed positions. Automatic external folding, artificial structure/confidence fallback, and auto-approval are disabled. Structure metrics are synthetic reconstruction checks; independent public-structure/DSSP validation remains open.

### Experimental human-reference prototype

`factorforge.engines.stemcell_3dp.StemCell3DPOptimizer` is a local experimental interface using a generic human-frequency reference and provisional Eijex review criteria. It is not an MSC/iPSC-calibrated model or a laboratory-approved SOP. See [prototype boundaries](docs/model_cards/human-reference-prototype.md). It is not a new production web host.
# Local research candidate: DP v2.2

The opt-in Python `DPV220Optimizer` shares the v2.1.1 DP core and adds reviewed
position penalties. The experimental generic-human adapter uses this core;
existing plant web/NB-Fold routes retain their prior objective. See the
[model card](docs/model_cards/dp_v2_2.md) for usage, relaxation reporting and
engineering validation boundaries. This is not an expression or folding claim.
