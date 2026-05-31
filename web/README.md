# FactorForge Web Interface

Plant CDS Design — N. benthamiana & Tobacco BY-2

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
const API_ENDPOINT = 'https://your-vercel-app.vercel.app/api/optimize';
```

## 🎨 Features

- ✅ File upload (FASTA)
- ✅ Text input (paste sequence)
- ✅ 4 optimization profiles
- ✅ Real-time results
- ✅ Custom restriction site input and removal report
- ✅ Download (FASTA, GenBank)
- ✅ Responsive design
- ✅ No login required

## 🔧 Tech Stack

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
