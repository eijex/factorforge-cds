# FactorForge Web Interface

## v3.6.1: current public policy

The public What's New modal contains selected user-facing changes. Historical entries do not override the current availability or data-handling policy.

1DP (sequence baseline), 3DP (PDB/mmCIF upload with reviewed positions), and 4DP (prepared pause sidecar JSON upload with reviewed model positions) are all publicly enabled in the web interface. Prepared-input CPU optimization runs with `FACTORFORGE_EXPERIMENTAL_3DP_ENABLED=true` and `FACTORFORGE_PUBLIC_4DP_ENABLED=true` by default; synthetic 4DP fixtures (`FACTORFORGE_LOCAL_4DP_FIXTURES=true`) remain restricted to explicit local loopback testing. Optimization from prepared structural inputs does not inherently require GCP or a GPU; automatic structure prediction is a separate workload. Structural and exit-tunnel model-offset comparisons are computational hypotheses (`UNVERIFIED_PUBLIC_CORPUS`, `Biological validation: NOT_AVAILABLE`).

FactorForge's built-in default SOP is applied automatically across 1DP, 3DP, and 4DP workflows. These built-in criteria are not an approved laboratory SOP, and custom SOP editing remains disabled in the public form. If you have your own custom SOP and want candidate sequences generated using it, email `eijex.lab@gmail.com`. GNU AGPL-3.0 covers the software. Agree separately on sharing and use before sending custom SOP materials.

The Arg/Ser, splice-remediation and 5-prime boundary comparison controls are
not connected to web execution. They are unchecked, disabled and collapsed
by default. They do not modify API requests; existing engine constraints and
advisory checks remain separate. No translation-efficiency benefit is asserted.

Focused CDS design and pre-synthesis review for N. benthamiana workflows.

Execution Trace shows response-backed engine/profile evidence when available.
It does not certify that every check passed, that platform services ran, or that
artifacts were persisted. Inspect the individual checks and QC matrix.

Sequence policy files can be downloaded as examples, edited, and uploaded. The
active profile persists only in this browser. Optional advanced settings remain
collapsed; the edit button expands them. A separate registry-backed modal builder
is not part of the verified web runtime. Uploaded profiles are examples or user
policies, not approved laboratory SOPs.

The default Design Setup shows the host, applied-settings summary and sequence
design profile. Engine selection, sequence/assembly constraints, review policy,
codon distribution and MoClo controls remain available under Advanced settings.
Watermark and capability-gated ML comparisons are under Experimental features.
Profile upload, example YAML download, active-profile export and editing are
available under Upload, download or edit profile. Collapsing a section does not
disable its settings or change the optimization request.

The baseline example is displayed as **FactorForge Default**, with the description
"Example configuration · Adapt to your workflow". Its internal profile ID and
serialized settings are unchanged; the laboratory-SOP boundary remains in profile details.

Sequence input listeners are connected before server metadata is fetched, so
early paste and slow metadata responses do not leave Generate disabled. The
server rule registry includes BpiI and SapI recognition targets requested by
the bundled profiles; profile-selected enforcement remains separate from detection.

The What's New panel summarizes the versioned CHANGELOG, including v3.5.2's
optional codon balancing, sequence-audit/evidence support and policy-driven review.
These historical notes do not change the current default engine or validation boundary.

## 🌐 Live Demo

**Web App**: https://factorforge.eijex.com
**API**: https://factorforge.eijex.com/api/optimize

## 📁 Files

```
web/
├── index.html          # Main page (3-panel layout)
├── css/
│   └── style.css       # Custom styles
└── js/
    └── app.js          # Main application logic
```

## 🚀 Quick Start

### Local Development

```bash
# Open in browser
open index.html

# Or use a local server
python -m http.server 8000
# Visit: http://localhost:8000
```

### Deploy to GitHub Pages

```bash
# Copy files to gh-pages branch
git checkout -b gh-pages
cp -r web/* .
git add .
git commit -m "Deploy frontend"
git push origin gh-pages

# Enable GitHub Pages in repo settings
# Settings → Pages → Source: gh-pages branch
```

## ⚙️ Configuration

### API Endpoint

Edit `js/app.js` line 3:

```javascript
const API_ENDPOINT = 'https://factorforge.eijex.com/api/optimize';
```

## 🎨 Features

- ✅ File upload (FASTA)
- ✅ Text input (paste sequence)
- ✅ Focused N. benthamiana web workflow
- ✅ Host selector for N. benthamiana and experimental Tobacco BY-2 workflows
- ✅ 4 public CDS design profiles
- 🧪 Deployment-gated ML preview and Rule-vs-ML dual comparison (disabled on public deployments by default)
- ✅ Windowed codon alignment viewer for long CDS sequences
- ✅ Capability-gated canonical DB save control (disabled by default)
- ✅ Real-time results
- ✅ Exact CDS-to-CDS codon shift accounting with explicit terminal-stop handling
- ✅ High-DPI interactive 1D codon tracks and measured-effects inspector; missing host frequencies remain `N/A`
- ✅ Custom restriction site input and removal report
- ✅ Optional reproducibility seed and Type IIS enzyme presets
- ✅ Progressive disclosure for alternative objectives and expert settings
- ✅ Upload-first Sequence Policy Profiles: download a commented example,
  export the active profile, drag and drop YAML or JSON, reset to the bundled
  example, and retain the active profile locally in the browser
- ✅ Privacy-gated public policy suggestions through a dedicated GitHub issue form;
  sensitive or laboratory-specific feedback stays on the private email path
- ✅ Capability-gated DP v2.1.1 2.1.1-dev local-guard development candidate;
  stable DP v2 2.0.1 remains the default feasibility path
- ✅ Domestication, MFE availability, and GC target transparency
- ✅ Researcher Decision Report with authoritative acceptance states, prioritized
  next actions, requested-vs-applied settings, detailed checks, candidate
  comparison, reproducibility provenance, and standalone print-friendly HTML
  that preserves the light or dark theme selected at download time
- ✅ Design Comparison section with reference/candidate metric context, DP v2.1.1
  three-axis evidence classes, and explicit independent-evaluation boundaries
- ✅ Sequence-free machine-readable evidence-record JSON export; sequence-bearing
  FASTA, GenBank, and HTML artifacts remain explicitly separate
- ✅ Download (FASTA, GenBank)
- ✅ Self-contained interactive HTML optimization report generated from the same canonical browser report data
- ✅ Responsive design
- ✅ No login required

## 🔧 Tech Stack

### Researcher design-review workflow

The desktop screen is a three-column workbench: sequence input, a compact Design
Brief, and Design Review. All three columns follow the same document scroll so a
shorter column does not appear frozen while a longer result is reviewed. The brief
exposes the expression host, recommended deterministic method, and applied
requirements. The primary workflow is a single versioned YAML/JSON Sequence
Policy Profile;
optional sequence, assembly, and review-policy controls retain their existing
API semantics inside one collapsed manual-overrides panel. The bundled profile is an illustrative
configuration rather than a complete, recommended, or approved SOP for a named laboratory. Uploaded
profiles are validated and stored only in the local browser; the server does not
retain them. Unavailable
execution modes, disabled objectives, and immutable reference policy are
excluded from the primary path. After a successful run, the same review column
shows computational checks, candidate evidence, sequence output, and downloads.

The Researcher Decision Report uses the API-provided `automated_decision`,
`decision_summary`, and acceptance-criteria rows as its authority; the browser
does not create a separate CAI or GC pass/fail rule. It leads with a decision
brief, translates non-passing checks into bounded review actions, and compares
requested settings with the API-recorded effective settings. Missing and
uncomputed values remain explicit rather than being displayed as zero. New
local-history entries retain the report snapshot and provenance without retaining
the raw input sequence. The evidence-record JSON omits raw input and output
sequences by design; the standalone HTML includes the optimized sequence and
displays a sequence-data handling notice.
The standalone report preserves the app theme selected when the file is
downloaded, while print output remains light for legibility.
For CDS inputs, the result header also shows exact sense-codon totals,
synonymous shifts, unchanged codons, and translation-derived substitutions.
The terminal stop codon is reported separately and excluded from sense-codon
accounting. The dual canvas map uses the same host-frequency scale on both
tracks. If the API does not supply a codon frequency, the track and inspector
show `N/A`; the browser never substitutes a hard-coded frequency or invents an
optimizer rationale. The primary interactive HTML export consumes the same
in-memory report model as the dashboard and embeds its sequences, audit
provenance, metric summary, and responsive canvas inspector without external
runtime dependencies.
For DP v2.1.1 results, both the in-app report and standalone HTML show the three
declared scientific axes and keep RNA folding in a separate independent-
evaluation block. The local guard reports its active 5′ GC layer and 5-nt
homopolymer ceiling; calibration and holdout status remain distinct. CDS inputs
receive a reference-versus-candidate table for GC,
first-30-nt GC, amino-acid identity, and nucleotide changes. Protein-only inputs
receive an explicitly labeled candidate-only view. These metric differences do
not establish expression, yield, or biological superiority.

The experimental comparison panel displays missing or invalid metrics as
`Not evaluated`. AA identity uses an explicit numeric `aa_identity` fraction;
responses without it cannot display a pass. Type IIS clearance requires an
explicit Boolean result, and conflicting site counts do not display as clean.
ML preview and database controls remain gated by deployment capabilities.


- **Frontend**: HTML5 + Tailwind CSS + Vanilla JS
- **Backend**: Vercel Serverless Functions (Python)
- **Hosting**: Vercel

## 📊 Browser Support

- Chrome 90+
- Firefox 88+
- Safari 14+
- Edge 90+

## 📄 License

GNU Affero General Public License v3.0 - see [LICENSE](../LICENSE)

## Current testing UI

The workflow selector offers 1DP baseline and experimental 3DP. Baseline retains the existing public design form. Selecting 3DP opens `http://127.0.0.1:8788/ui?mode=paired_3dp` on the tester's own computer; no entered sequence is copied or sent. An operator must first configure/start the local service. This link is not a hosted public 3DP service.

The local UI offers actual 1DP-only DP v2 generation (one candidate) or 3DP generation with a DP v2 comparator (two candidates). Registered target/structure references and human approval are required. The tester selects the method; no automatic biological winner is chosen.

Watermark is displayed as **Disabled**, remains unchecked even after an uploaded profile requests it, and is not enabled by this UI. Other profile settings and the design objective are preserved. This UI pause does not remove library-level experimental functions.

## Direct 3DP browser testing (unreleased)

The FactorForge form supports explicit 1DP/3DP selection. 3DP uploads a PDB/mmCIF file, selects chain and one-based file model order, and requires reviewed one-based protein positions. The selected observed chain must exactly match the complete protein; partial or homologous structures require the retained local reviewed mapping workflow. Direct requests use `/api/optimize` with `design_mode=paired_3dp` and `structure_context={content,format,chain,model,positions,positions_reviewed,target_adaptiveness,penalty_weight}`.

Limits: 3–300 amino acids, structure content 250 KB, HTTP JSON body 1 MB. Both candidates use the configured codon reference, requested GC and canonical motif exclusions plus custom additions. Invalid inputs and infeasible candidates fail without a baseline fallback. Existing sequence-only settings and API requests remain supported. 3DP uses its explicit DP settings and does not apply Kozak, TpA, template or seed transformations. SOP checks are reported for both candidates; regeneration is manual and never silently changes mode. Watermark stays disabled.

This is stateless computational testing, with no external service dependencies or remote database recording. Those existing integrations remain available for later operational use. The structure file is not retained in the response or browser history; normal design results may be retained in local browser history. Submit only public-safe inputs. Structural conformity metrics are codon-weight proxies, not folding/expression measurements. These source changes are not deployed or a new package release.

## Current public policy UI (2026-10-12)

Custom SOP management is paused. The public form applies the built-in FactorForge default policy and offers standard 1DP. General 3DP/4DP web access is paused.

### Local 4DP preview

For local development, bind the app and API to loopback and set `FACTORFORGE_LOCAL_4DP_PREVIEW=true`. The browser uses the server capability before enabling the control. Upload a prepared sidecar that matches the input protein and review its positions. The comparison uses its declared offset, with a preview limit of 100 amino acids and 100 KB for the sidecar. The optional synthetic test control requires `FACTORFORGE_LOCAL_4DP_FIXTURES=true` and an explicit test selection. Model positions and sequence checks are computational evidence. The 4DP preview has no biological outcome data.

### Scientific web observations

Direct 1DP/3DP browser testing reports sense-CDS RNA motif scans and available initiation-window MFE in results, downloadable reports and local history. Structure confidence observations can suggest positions for explicit acceptance and review; unsupported inputs do not invent positions. These checks do not establish expression, yield or folding benefits. The default policy remains active and custom SOP management remains paused. See `docs/experimental-structure.md` for the API and measurement boundary.

Web structure handling now requires exact-match uploaded PDB/mmCIF coordinates and explicit reviewed positions. Automatic external folding, artificial structure/confidence fallback, and auto-approval are disabled. Structure metrics are synthetic reconstruction checks; independent public-structure/DSSP validation remains open.
