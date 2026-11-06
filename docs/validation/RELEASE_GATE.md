# FactorForge release-gate status

## v3.5.12 published (2026-11-06)

This release adds consent/policy wording, default-off experimental web
availability and Bioconda installed-package smoke coverage.

- [Exact-source release CI](https://github.com/eijex/factorforge-cds/actions/runs/42181806379): all eight Python/OS combinations passed.
- [PyPI publication and fresh-install smoke](https://github.com/eijex/factorforge-cds/actions/runs/42184356326): passed.
- [Remote browser suite](https://github.com/eijex/factorforge-cds/actions/runs/42184346124): 63 passed.
- Public web v3.5.12 was deployed; synthetic 1DP preserves amino-acid identity and paused structure routes return HTTP 403.
- [Bioconda recipe update](https://github.com/bioconda/bioconda-recipes/pull/71287) uses the verified published source hash; its merge/build is pending. Prior installed-package smoke coverage is for v3.5.11, not a new Bioconda artifact.
- [Zenodo v3.5.12 archive](https://doi.org/10.5281/zenodo.24001252) is bound to the GitHub v3.5.12 tag. Citation metadata was updated after the archive appeared; the release tag remains immutable.

No new biological validation is claimed.

## Historical v3.5.10 gate scope and commands

This gate covers the public Open Bio output and I/O contracts. FactorForge remains
an in-silico CDS design assistant; this gate does not establish wet-lab, yield, or
clinical performance.

The maintenance gate also verifies database-free wheel/source installs, an explicit
PostgreSQL plus pinned DB-core integration, private artifact checksums, transaction
rollback, and concurrent retries. Shared persistence remains experimental and is
not automatically enabled in the public stateless deployment. No legacy-data or
production-schema migration is authorized by the public package release.

## Release Gate Commands

```bash
python -m pytest tests/test_registry_production_sync.py -v -W always
python -m pytest tests/test_design_package_schema.py -v
python -m pytest tests/test_design_package_semantics.py -v
python -m pytest tests/test_design_package_serialization.py -v
python -m pytest tests/test_openbio_missing_metric_contract.py -v
python -m pytest tests/test_iupac_validation.py -v
python -m pytest tests/test_fasta_io.py -v
python -m pytest tests/test_host_profile_metadata.py -v
python -m pytest tests/test_benchmark_regression.py -v
python -m pytest tests/test_no_raw_sequence_logging.py -v
python scripts/parameter_audit.py --mode discovery
python -m pytest tests/ -q
```

`parameter_audit.py --mode discovery` writes all findings and exits `0`.
`parameter_audit.py --mode strict` writes the same findings and exits `1` when
unregistered findings remain; strict mode is a future hard gate.

## Merge Blockers

- Registry-production sync emits a warning or fails.
- `DEFAULT_CAI_TARGET`, `DEFAULT_GC_LOW`, or `DEFAULT_GC_HIGH` cannot be imported
  from configured local source.
- `Domesticator.GOLDEN_GATE_ENZYMES` is absent or differs from
  `("BsaI", "BpiI", "BsmBI")`.
- `BbsI` appears as a registry value, production enzyme list item, active audit
  target, strict sync expected value, or benchmark expected active enzyme.
- Missing MFE is encoded as `0.0` or a string placeholder.
- Raw sequence appears in logs, exceptions, FASTA headers, public Design Package
  metadata, or benchmark CSV/JSON.
- N. benthamiana taxonomy ID differs from `4100`.
- BY-2 is marked stable.
- Wolffia globosa is marked stable.
- Benchmark summary lacks a registry hash.
- Public wording implies expression, yield, wet-lab, or clinical performance.
- A publication-submission artifact is added to the public repository.
- Biopython, DNAChisel, or pySBOL3 is newly added as a required dependency by
  this release gate work.

## Manual Review Checklist

- Confirm the public Design Package contains only a sequence hash, never raw
  input or output sequence.
- Confirm all claim-boundary booleans are `true`.
- Confirm MFE status, value, and used fields agree.
- Confirm FASTA output headers contain only allowlisted metadata.
- Confirm `BbsI` occurrences are alias or provenance explanations only.
- Confirm benchmark results contain no raw-sequence columns.
- Confirm the registry changelog records the canonical label normalization.
- Confirm `docs/host_profiles/HOST_PROFILE_REGISTRY_ROADMAP.md` remains deferred
  until host profile registry migration is done.
- Review `git ls-files --others --exclude-standard` and `git status --short`
  before publishing.

## Non-Scope Reminder

This gate does not add optimizer benchmarking conclusions, ML engines, Kaggle or
Colab adapters, MCP contracts, GenBank/SBOL export, DNAChisel integration,
external optimizer comparisons, new stable host profiles, or external
publication artifacts. The host profile roadmap and Wolffia profile work
remain deferred until the post-070 follow-up.

## Import Environment Note

For direct import checks, use pytest, `PYTHONPATH=src`, or an editable install
instead of plain `python -c`. This avoids accidentally importing an older
site-packages installation.
