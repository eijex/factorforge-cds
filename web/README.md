# FactorForge Web Interface

Focused CDS design and pre-synthesis review for N. benthamiana workflows.

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
- ✅ Upload-first laboratory SOP profiles: download a human-readable YAML
  template, drag and drop YAML or JSON, reset to the bundled conservative
  template, and retain the active profile locally in the browser
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
requirements. The primary workflow is a single versioned YAML/JSON SOP file;
optional sequence, assembly, and review-policy controls retain their existing
API semantics inside one collapsed manual-overrides panel. The bundled profile is a conservative starting
template rather than a complete or approved SOP for a named laboratory. Uploaded
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
