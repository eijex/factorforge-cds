# Host Profile Registry

The host profile registry defines explicit expression contexts for sequence
design. The central principle is:

> Same species does not mean same profile. Species, cell line, tissue/cell type,
> expression platform, data source, confidence, and validation status must all
> be specified separately.

A profile should describe the host context that the optimizer is targeting, not
just the taxonomic species. This avoids treating whole plants, suspension cells,
stable mammalian production lines, research cell types, and insect systems as
interchangeable.

## Required Profile Fields

Each profile YAML entry must include:

```yaml
id:
display_name:
species:
strain_or_cell_line:
tissue_or_cell_type:
expression_platform:
profile_status:
default_profile:
codon_usage:
  source_type:
  table_status:
gc_content:
  ideal_band:
  soft_band:
cai:
  target_cai:
  min_feasible_cai:
rules:
  homopolymer_expression_warn_nt:
  homopolymer_synthesis_warn_nt:
recommended_use:
not_recommended_for:
validation_status:
```

`profile_status` should use values such as `supported`, `experimental`, or
`planned`. `validation_status` should separately describe what has and has not
been validated for the expression context.

## Example Profiles

### Nicotiana benthamiana

```yaml
id: nicotiana_benthamiana
display_name: Nicotiana benthamiana transient expression
species: Nicotiana benthamiana
strain_or_cell_line: whole plant
tissue_or_cell_type: leaf tissue
expression_platform: agroinfiltration transient expression
profile_status: supported
default_profile: true
codon_usage:
  source_type: host codon usage table
  table_status: supported
gc_content:
  ideal_band: [55.0, 65.0]
  soft_band: [45.0, 75.0]
cai:
  target_cai: 0.90
  min_feasible_cai: 0.82
rules:
  homopolymer_expression_warn_nt: 6
  homopolymer_synthesis_warn_nt: 8
recommended_use:
  - Whole-plant transient expression design.
  - Default plant expression profile.
not_recommended_for:
  - Mammalian expression design.
  - Insect or microbial expression design.
validation_status: Supported for current plant-first optimization workflows.
```

### Nicotiana tabacum BY-2

```yaml
id: nicotiana_tabacum_by2
display_name: Nicotiana tabacum BY-2 suspension culture
species: Nicotiana tabacum
strain_or_cell_line: BY-2
tissue_or_cell_type: suspension culture cell line
expression_platform: plant cell suspension expression
profile_status: experimental
default_profile: false
codon_usage:
  source_type: species-level codon usage proxy
  table_status: proxy
gc_content:
  ideal_band: [55.0, 65.0]
  soft_band: [45.0, 75.0]
cai:
  target_cai: 0.90
  min_feasible_cai: 0.82
rules:
  homopolymer_expression_warn_nt: 6
  homopolymer_synthesis_warn_nt: 8
recommended_use:
  - Experimental planning for BY-2 suspension expression.
  - Comparison against plant transient-expression defaults.
not_recommended_for:
  - Claiming BY-2 expression performance validation.
  - Production release without profile-specific validation.
validation_status: Uses the N. tabacum codon table as a proxy; not wet-lab validated in FactorForge for BY-2 expression performance.
```

### CHO-K1

```yaml
id: cho_k1
display_name: CHO-K1 stable expression
species: Cricetulus griseus
strain_or_cell_line: CHO-K1
tissue_or_cell_type: mammalian ovary-derived cell line
expression_platform: mammalian stable expression
profile_status: planned
default_profile: false
codon_usage:
  source_type: planned host-specific codon usage table
  table_status: planned
gc_content:
  ideal_band: [45.0, 65.0]
  soft_band: [40.0, 70.0]
cai:
  target_cai: 0.90
  min_feasible_cai: 0.82
rules:
  homopolymer_expression_warn_nt: 6
  homopolymer_synthesis_warn_nt: 8
recommended_use:
  - Future mammalian stable-expression design.
  - Host-profile schema testing.
not_recommended_for:
  - Current production optimization.
  - Plant expression defaults.
validation_status: Planned; not implemented or validated.
```

### Human iPSC

```yaml
id: human_ipsc
display_name: Human iPSC research expression
species: Homo sapiens
strain_or_cell_line: induced pluripotent stem cell
tissue_or_cell_type: pluripotent stem cell
expression_platform: research cell culture expression
profile_status: planned
default_profile: false
codon_usage:
  source_type: planned host- and context-specific codon usage review
  table_status: planned
gc_content:
  ideal_band: [45.0, 65.0]
  soft_band: [40.0, 70.0]
cai:
  target_cai: 0.90
  min_feasible_cai: 0.82
rules:
  homopolymer_expression_warn_nt: 6
  homopolymer_synthesis_warn_nt: 8
recommended_use:
  - Future research-only human iPSC sequence design review.
  - Profile schema planning.
not_recommended_for:
  - Clinical, diagnostic, or regulatory claims.
  - Current production optimization.
validation_status: Research/planned; not implemented or validated.
```

### Apis mellifera

```yaml
id: apis_mellifera
display_name: Apis mellifera expression design
species: Apis mellifera
strain_or_cell_line: unspecified
tissue_or_cell_type: unspecified insect expression context
expression_platform: planned insect expression profile
profile_status: planned
default_profile: false
codon_usage:
  source_type: planned insect codon usage source audit
  table_status: planned
gc_content:
  ideal_band: [35.0, 55.0]
  soft_band: [30.0, 60.0]
cai:
  target_cai: 0.90
  min_feasible_cai: 0.82
rules:
  homopolymer_expression_warn_nt: 6
  homopolymer_synthesis_warn_nt: 8
recommended_use:
  - Future insect host-profile planning.
  - Codon source audit design.
not_recommended_for:
  - Current production optimization.
  - Treating all insect systems as equivalent.
validation_status: Planned; not implemented or validated.
```
