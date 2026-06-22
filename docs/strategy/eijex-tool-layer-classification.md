# Eijex Tool Layer Classification

**Decision date:** 2026-06-11
**Status:** Approved — governs FactorForge scope and future Eijex module planning

---

## Core Principle

Eijex is a modular bioinformatics workbench for plant molecular farming,
starting with registry-driven CDS design and expanding through optional
protein, reference, and evidence annotation layers.

FactorForge is **not** an all-in-one platform. Its scope is:

> protein sequence → plant-compatible CDS design → sequence-level validation → Design Package

---

## Tool / Resource Classification

| Tool / Resource | FactorForge core | Eijex | Recommended location | Notes |
|---|---|---|---|---|
| Genome sequence (SGN QLD183) | ✅ | ✅ | Host Profile / Reference Manifest | CDS source, codon usage provenance |
| Gene annotation (GFF3) | ✅ | ✅ | Host Profile / Reference Manifest | CDS coordinate provenance |
| Genome Browser (JBrowse/SGN) | ❌ | link-out | External Reference Link | No need to build internally |
| BLAST | limited | ✅ (future) | Reference QC Adapter | Sequence identity / QC only |
| AlphaFold | ❌ | ✅ (future) | Protein Intelligence Layer | Structure annotation, not CDS optimizer |
| ESM / ESMFold | ❌ | ✅ (future) | Protein Intelligence Layer | Protein embedding / risk annotation |
| Proteomics DB | ❌ | ✅ (later) | Protein Evidence Layer | Protein-level expression evidence |
| CRISPR-P | ❌ | future product | EditForge / TargetForge | Genome editing — different problem domain |
| CCTop | ❌ | future product | EditForge / TargetForge | CRISPR guide/off-target — different domain |
| VIGS DB | ❌ | future product | Gene2PhenotypeHub | Functional genomics / phenotype layer |

---

## Immediate Action (now)

- **Genome sequence + Gene annotation** → Reference Provenance Pack
  - `docs/reference/NBENTHAMIANA_REFERENCE.md`
  - `scripts/download_nbenthamiana_reference.py`
  - checksum + filtering rule documentation

---

## Design Package: optional annotation slots

External tools attach as `external_annotations` slots — default is `not_run`.
FactorForge core operates independently with no external tool dependencies.

```json
{
  "external_annotations": {
    "protein_structure": { "status": "not_run", "provider": null },
    "sequence_similarity": { "status": "not_run", "provider": null },
    "proteomics_evidence": { "status": "not_run", "provider": null }
  }
}
```

---

## Claim boundary reminder

- No yield / expression / wet-lab / clinical claims in any layer
- AlphaFold / ESM annotations = structural context only, not performance prediction
- BLAST = sequence identity QC, not functional validation
