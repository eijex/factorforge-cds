# Contributing to FactorForge

Thank you for your interest in contributing to FactorForge. This project is maintained by [Eijex](https://github.com/eijex) and welcomes contributions from the scientific community.

**Live demo**: [factorforge.vercel.app](https://factorforge.vercel.app) — try FactorForge in your browser without installation.

## Ways to Contribute

- Report bugs and unexpected results
- Suggest new features or biological use cases
- Improve documentation
- Submit pull requests with fixes or enhancements
- Share validation results from wet-lab experiments

## Getting Started

### 1. Fork and Clone

```bash
git clone https://github.com/eijex/factorforge.git
cd factorforge
pip install -e ".[dev]"
```

### 2. Run Tests

```bash
pytest tests/ -v
```

All tests must pass before submitting a pull request.

### 3. Code Style

FactorForge uses [ruff](https://github.com/astral-sh/ruff) for linting and formatting:

```bash
ruff check .
ruff format .
```

## Submitting a Pull Request

1. Create a branch from `main`
2. Make your changes with clear, focused commits
3. Add or update tests as needed
4. Update documentation if behavior changes
5. Open a pull request with a clear description

### PR Title Format

```
feat: add support for Arabidopsis thaliana codon table
fix: correct GC% calculation for short sequences
docs: clarify tAI computation requirements
```

## Reporting Issues

When reporting a bug, please include:
- Python version and OS
- FactorForge version (`pip show factorforge`)
- Minimal reproducible example
- Expected vs actual behavior

## Biological Contributions

If you have wet-lab validation data (expression levels, yield measurements), we especially welcome:
- Comparisons of FactorForge-optimized vs native sequences
- Multi-host validation results
- Edge case sequences (very long proteins, transmembrane regions, etc.)

Please open an issue to discuss before submitting large datasets.

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/). By participating, you agree to uphold these standards.
