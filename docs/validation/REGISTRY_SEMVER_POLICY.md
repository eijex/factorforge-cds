# Registry Semver Policy

The canonical parameter registry (`src/factorforge/registry/current_parameter_registry.yaml`) governs the constraint values used across FactorForge optimization, scoring, and benchmarking. Because benchmark results and reported metrics are derived from registry values, registry changes are versioned explicitly.

## Version impact rules

- **Patch**: Editorial or provenance-only changes that do not alter any resolved parameter value (for example, fixing a citation or comment).
- **Minor**: Adding new parameters or sections that do not change the resolved value of any existing parameter consumed by current benchmarks or the optimizer.
- **Major**: Any change that alters a resolved value consumed by the optimizer, scoring, or benchmark configuration (for example, changing a GC range, CAI target, or forbidden Type IIS site list).

## Required actions on a major registry change

1. The registry `version` field is bumped and the change is documented in the registry's own provenance metadata.
2. Any benchmark results derived under the prior registry version are treated as frozen and are not silently overwritten — a new versioned results directory is created (for example `benchmarks/results/v3.2.1/`), per [`benchmarks/README.md`](../../benchmarks/README.md).
3. The benchmark spec (`benchmarks/benchmark_spec.yaml`) and any dependent documentation are reviewed for consistency with the new registry version.
4. The change is reflected in the project changelog with its semver impact.

## Rationale

Benchmark and validation claims are only meaningful when the parameter values behind them are traceable to a specific registry version. Treating registry value changes as semver-significant ensures that historical results remain interpretable and that public claims can always be traced to the exact configuration that produced them.
