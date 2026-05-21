# FactorForge Web Interface

3-Panel Dashboard for N. benthamiana Codon Optimization

## 🌐 Live Demo

**Web App**: https://factorforge.vercel.app
**API**: https://factorforge.vercel.app/api/optimize

## 📁 Files

```
web/
├── index.html          # Main page (3-panel layout)
├── css/
│   └── style.css       # Custom styles
├── js/
│   ├── app.js          # Main application logic
│   └── gemini.js       # Gemini AI chatbot (optional)
└── assets/             # Images, logos
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

### Gemini API Key

Edit `js/gemini.js` line 3:

```javascript
const GEMINI_API_KEY = 'your-gemini-api-key';
```

Get your free API key: https://aistudio.google.com/app/apikey

## 📖 Full Documentation

See [DEPLOY_GUIDE.md](../DEPLOY_GUIDE.md) for complete deployment instructions.

## 🎨 Features

- ✅ File upload (FASTA)
- ✅ Text input (paste sequence)
- ✅ 5 optimization profiles
- ✅ Real-time results
- ✅ Download (FASTA, GenBank)
- ✅ Gemini AI chatbot (optional)
- ✅ Responsive design
- ✅ No login required

## 🔧 Tech Stack

- **Frontend**: HTML5 + Tailwind CSS + Vanilla JS
- **Backend**: Vercel Serverless Functions (Python)
- **AI**: Google Gemini 2.0 Flash
- **Hosting**: GitHub Pages (frontend) + Vercel (backend)

## 📊 Browser Support

- Chrome 90+
- Firefox 88+
- Safari 14+
- Edge 90+

## 📄 License

Apache License 2.0 - see [LICENSE](../LICENSE)
