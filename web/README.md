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
- ✅ Custom restriction site input and removal report
- ✅ Optional reproducibility seed and Type IIS enzyme presets
- ✅ Progressive disclosure for alternative objectives and expert settings
- ✅ Domestication, MFE availability, and GC target transparency
- ✅ Researcher-oriented Design Review Report with authoritative acceptance states,
  detailed checks, candidate comparison, reproducibility provenance, and a
  standalone print-friendly HTML download
- ✅ Sequence-free machine-readable evidence-record JSON export; sequence-bearing
  FASTA, GenBank, and HTML artifacts remain explicitly separate
- ✅ Download (FASTA, GenBank)
- ✅ Responsive design
- ✅ No login required

## 🔧 Tech Stack

### Researcher design-review workflow

The desktop screen is a persistent three-column workbench: sequence input,
a compact Design Brief, and Design Review. The review column keeps the research
workflow and result destination visible before and after generation. The brief
exposes the expression host, recommended deterministic method, and applied
requirements. Optional sequence, assembly, and review-policy controls retain
their existing API semantics through progressive disclosure. Unavailable
execution modes, disabled objectives, and immutable reference policy are
excluded from the primary path. After a successful run, the same review column
shows computational checks, candidate evidence, sequence output, and downloads.

The Design Review Report uses the API-provided `automated_decision`,
`decision_summary`, and acceptance-criteria rows as its authority; the browser
does not create a separate CAI or GC pass/fail rule. Missing and uncomputed values
remain explicit rather than being displayed as zero. New local-history entries
retain the report snapshot and provenance without retaining the raw input
sequence. The evidence-record JSON omits raw input and output sequences by
design; the standalone HTML includes the optimized sequence and displays a
sequence-data handling notice.

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
