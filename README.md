# FactorForge

**Open-source constraint-based CDS design and pre-synthesis sequence review for plant CDS workflows, with primary support for *Nicotiana benthamiana* (Tobacco BY-2: experimental).**

[![License](https://img.shields.io/badge/license-AGPL--3.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![PyPI](https://img.shields.io/pypi/v/factorforge-cds.svg)](https://pypi.org/project/factorforge-cds/)
[![CI](https://github.com/eijex/factorforge-cds/actions/workflows/ci.yml/badge.svg)](https://github.com/eijex/factorforge-cds/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/eijex/factorforge-cds/branch/main/graph/badge.svg)](https://codecov.io/gh/eijex/factorforge-cds)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20407330.svg)](https://doi.org/10.5281/zenodo.20407330)
[![Web App](https://img.shields.io/badge/web-factorforge.eijex.com-brightgreen.svg)](https://factorforge.eijex.com)
[![Buy Me A Coffee](https://img.shields.io/badge/Buy%20Me%20a%20Coffee-Donate-yellow?style=flat&logo=buy-me-a-coffee)](https://www.buymeacoffee.com/eijex)


FactorForge performs deterministic CDS design with CAI/GC metrics, PolyA-signal screening, and Golden Gate/MoClo-aware checks. It is positioned as a pre-synthesis review harness: it helps teams generate reproducible CDS candidates, inspect assembly-relevant sequence constraints, and package design metadata before downstream synthesis, cloning, or experimental review. Primary support: *N. benthamiana* (agroinfiltration). Experimental host context: Tobacco BY-2 (`--host by2`).

FactorForge v3.5.11 uses independently versioned engines:

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
pip install factorforge-cds
factorforge optimize my_protein.fasta -o output.fasta
```

Or use the **[web app](https://factorforge.eijex.com)** — no installation required.

---

## Access Options

| Method | Description | Link |
|--------|-------------|------|
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
FactorForge v3.5.11 (2026). Open-source constraint-based CDS design and sequence review.
Eijex. https://github.com/eijex/factorforge-cds
```

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

The web UI now labels watermark as disabled and distinguishes the sequence-only baseline from experimental structure-conditioned testing. The local approved workflow supports explicit 1DP-only or paired 3DP mode. 3DP availability is local research testing; it is not a new public optimization endpoint or a biological-performance claim. See [web UI boundaries](web/README.md).

## Direct 3DP browser testing (unreleased)

The FactorForge form supports explicit 1DP/3DP selection. 3DP uploads a PDB/mmCIF file, selects chain and one-based file model order, and requires reviewed one-based protein positions. The selected observed chain must exactly match the complete protein; partial or homologous structures require the retained local reviewed mapping workflow. Direct requests use `/api/optimize` with `design_mode=paired_3dp` and `structure_context={content,format,chain,model,positions,positions_reviewed,target_adaptiveness,penalty_weight}`.

Limits: 3–300 amino acids, structure content 250 KB, HTTP JSON body 1 MB. Both candidates use the configured codon reference, requested GC and canonical motif exclusions plus custom additions. Invalid inputs and infeasible candidates fail without a baseline fallback. Existing sequence-only settings and API requests remain supported. 3DP uses its explicit DP settings and does not apply Kozak, TpA, template or seed transformations. SOP checks are reported for both candidates; regeneration is manual and never silently changes mode. Watermark stays disabled.

This is stateless computational testing, with no external service dependencies or remote database recording. Those existing integrations remain available for later operational use. The structure file is not retained in the response or browser history; normal design results may be retained in local browser history. Submit only public-safe inputs. Structural conformity metrics are codon-weight proxies, not folding/expression measurements. These source changes are not deployed or a new package release.

## Current public policy UI (2026-10-12)

Custom SOP management is paused. The built-in FactorForge default policy is applied; upload/export/preset/editing controls are hidden and disabled. Their code and previous browser configuration remain available for later reactivation through ENABLE_CUSTOM_SOP. Saved custom SOPs are ignored without deletion. 1DP/3DP mode selection and explicit structure inputs remain active. This supersedes the earlier public SOP upload/export instructions. Public API contracts are retained.

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
