# FactorForge v3.2.3 Release Gate

This gate covers the public Open Bio output and I/O contracts. FactorForge remains
an in-silico CDS design assistant; this gate does not establish wet-lab, yield, or
clinical performance.

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
- Confirm the registry changelog records the 091 canonical label normalization.
- Confirm `docs/host_profiles/HOST_PROFILE_REGISTRY_ROADMAP.md` remains deferred
  until Job 070 is Done.
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
site-packages installation, as observed during Job 091.
