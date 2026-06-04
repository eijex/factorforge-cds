# Rule Engine Roadmap

The rule engine reports sequence features that can affect expression quality,
synthesis, cloning, or downstream handling. Rule findings should remain
separate from composite optimization score so that risk review is explicit.

## Rule Status

| Rule | Current status | Notes |
|------|---------------|-------|
| Homopolymer | Implemented | Two thresholds: expression (>=6 nt) vs synthesis (>=8 nt) |
| GC window | Implemented | 60-nt window, 30-nt step |
| Rare codon run | Implemented | w < 0.3, min_run=3 default |
| Restriction sites | Implemented | Configurable site list |
| PolyA motifs | Implemented | Heuristic; false positive risk noted |
| Splice-like motifs | Implemented | Partial; false positive risk noted |
| Repeat patterns | Planned | Not yet implemented |
| CpG / TpA dinucleotide | Implemented | CpG opt-in for plant expression scoring, TpA active by default |

## Homopolymer Thresholds

FactorForge intentionally uses two homopolymer thresholds because expression
quality and synthesis difficulty are independent risks:

- `HOMOPOLYMER_EXPRESSION_WARN_NT = 6`: expression stability and sequence
  quality warning.
- `HOMOPOLYMER_SYNTHESIS_WARN_NT = 8`: synthesis and manufacturing difficulty
  warning.

The expression threshold is stricter because shorter A/T-rich runs can still
matter for expression-stability review. The synthesis threshold is less strict
and is used by the profile rule engine's synthesis-oriented scan.

## GC Window Rules

GC window scanning is implemented with a 60-nt window and 30-nt step. The
current constants are useful defaults, but the future migration path is to move
GC window configuration into host profile YAML. Host-specific profiles should
define acceptable local GC behavior rather than relying on global defaults.

## Rare Codon Runs

Rare codon run scanning is implemented with a default rare-codon threshold of
`w < 0.3` and `min_run=3`. This rule should continue to report the threshold
used in the finding so snapshots can detect silent semantic changes.

## Restriction Sites

Restriction-site scanning is implemented with a configurable site list. The
assembly-friendly workflow currently focuses on Type IIS avoidance, but the rule
engine should remain capable of scanning arbitrary configured motifs.

## Motif Rules

PolyA motifs and splice-like motifs are implemented as heuristic scanners.
Findings should be treated as review signals rather than definitive biological
validation. Future reports should carry false-positive warnings and context
windows.

## Dinucleotide Rules

CpG and TpA dinucleotide detection and fixing are implemented. Plant scoring
defaults keep CpG penalty opt-in while TpA remains active by default. Mammalian
profiles may need different dinucleotide semantics and should make those choices
explicit in host profile YAML.

## Planned Repeat Pattern Rule

Repeat-pattern penalties are planned, not implemented. The intended scope is to
detect local repeat structures that can affect synthesis, cloning, or sequence
stability without conflating them with homopolymer-only checks.
