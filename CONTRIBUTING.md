# Contributing to FactorForge

Thank you for your interest in contributing to FactorForge. This project is maintained by [Eijex](https://github.com/eijex) and welcomes contributions from the scientific community.

**Live demo**: [factorforge-cds.vercel.app](https://factorforge-cds.vercel.app) — try FactorForge in your browser without installation.

## Ways to Contribute

- Report bugs and unexpected results
- Suggest new features or biological use cases
- Improve documentation
- Submit pull requests with fixes or enhancements
- Share validation results from wet-lab experiments
- Contribute codon usage tables for new plant hosts

## Getting Started

### Prerequisites

- Python 3.10 or higher
- Git

### 1. Fork and Clone

```bash
git clone https://github.com/eijex/factorforge-cds.git
cd factorforge
```

### 2. Install Dependencies

For general development and testing:

```bash
pip install -e ".[dev]"
```

### 3. Run Tests

```bash
pytest tests/ -v
```

All tests must pass before submitting a pull request. Tests are organized under `tests/` mirroring the `src/factorforge/` structure:

```
tests/
├── engines/profile/    # Profile engine registration and contracts
├── test_analysis/       # Sequence metrics and DP feasibility tests
├── api/           # API contract tests
└── validation/    # Design Package schema
```

When adding a new feature, add tests in the corresponding subdirectory.

### 4. Code Style

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
- FactorForge version (`pip show factorforge-cds`)
- Minimal reproducible example
- Expected vs actual behavior

## Biological Contributions

### Wet-lab Validation Data

If you have wet-lab validation data (expression levels, yield measurements), we especially welcome:
- Comparisons of FactorForge-optimized vs native sequences
- Multi-host validation results
- Edge case sequences (very long proteins, transmembrane regions, etc.)

Preferred submission paths:
1. **Google Form** (recommended): [Submit wet-lab result](https://docs.google.com/forms/d/e/1FAIpQLSeSx-wYvF6YwHhSPdLMl-L44frCugdm25X_eDz50OaqTD66qA/viewform?usp=header)
2. **GitHub Issue**: use the `wet_lab_result` issue template
3. See [VALIDATION.md](VALIDATION.md) for the full submission format

Please open an issue to discuss before submitting large datasets.

### Codon Usage Tables

Contributions of codon usage tables for new plant hosts are highly valuable. If you have RSCU data or a curated codon usage table for a plant species not currently supported, please open an issue with:
- Species name and genome reference
- Source of the codon usage data (database, publication, or custom analysis)
- Target expression system (transient, stable, viral)

## Release Process (maintainers)

Releases are managed by the project maintainer using the following workflow:

1. Update `CHANGELOG.md` — move `[Unreleased]` entries to the new version section
2. Run the version bump script:
   ```bash
   python bump_version.py X.Y.Z          # updates 16 version-bearing files
   python bump_version.py X.Y.Z --dry-run # preview changes without writing
   ```
3. Manually update `web/index.html` changelog panel (new Current block)
4. Add summary entry to `docs/changelog.md`
5. Commit, tag, and push — GitHub Actions handles PyPI, Docker, and GitHub Release automatically

See `CHANGELOG.md` for the full release checklist.

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/). By participating, you agree to uphold these standards.
