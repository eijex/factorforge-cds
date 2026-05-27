# FactorForge Archive

This directory preserves historical FactorForge implementation tracks for
reference. Archived code is not imported by the installed package and is not
part of the supported runtime API.

Current production code lives under:

```text
src/factorforge/engines/profile/
```

## Archive Layout

| Directory | Generation | Status | Description |
|-----------|-----------|--------|-------------|
| `v1-nbent-opticodon/` | v1 | Internal | Thesis-derived codon optimization baseline (NBent_OptiCodon); not vendored — pointer only |
| `v2-rule-engine/` | v2 | Internal → Production | Deterministic rule-based engine that became `factorforge.engines.profile` |
| `v3-ml-prototype/` | v3-alpha | Archived | ML-based design attempt; insufficient performance vs deterministic baseline; preserved for research provenance |

The archive keeps implementation history visible without mixing old engine
names into the public production package. See the project README for the full
[development history](https://github.com/eijex/factorforge-cds#development-history).
