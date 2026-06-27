# Rule Engine Roadmap

The rule engine reports sequence features relevant to expression-context review,
synthesis, cloning, or downstream handling. Rule findings should remain
separate from composite optimization score so that risk review is explicit.

## Rule Status

This table lists all 9 advisory sequence-risk scanners run by `RuleEngine.scan_all()`
(`scan_mode="full"`, the default), plus the 2 assembly-review checks.

| Rule | Current status | Notes |
|------|---------------|-------|
| Homopolymer | Implemented | Two thresholds: expression-context review (>=6 nt) vs synthesis (>=8 nt) |
| Local GC extremes | Implemented | 50-nt sliding window, 25-75% threshold (synthesis guard). Distinct from `GC_OPT_MID` (43.5% for *N. benthamiana*, v3.3.0+; a scoring target) and from the 60-nt/30-nt-step GC-window calculation in the archived `analysis/metrics.py` validation path, which is not part of the active default scan. |
| Rare codon run | Implemented | w < 0.3, min_run=3 default |
| Restriction sites (assembly review) | Implemented | Configurable site list; halts the production design pipeline when unresolvable |
| PolyA motifs | Implemented | Heuristic; false positive risk noted |
| ARE (AU-rich elements) | Implemented | Pattern-based heuristic scanner |
| AT-rich runs | Implemented | Minimum run length 6 nt by default |
| Splice-like motifs | Implemented | Partial; false positive risk noted |
| Repeat patterns | Implemented | Perfect/tandem repeats >=15 nt (recombination risk); previously misdocumented in this table as "Planned" — corrected 2026-06-21 |
| CpG / TpA dinucleotide | Implemented | CpG opt-in for plant sequence-context scoring, TpA active by default |
| MoClo overhang validity/collision (assembly review) | Implemented, opt-in | Not run by default (requires `construct_template`); when run, reports warnings only and never halts the pipeline |

## Homopolymer Thresholds

FactorForge intentionally uses two homopolymer thresholds because
expression-context review and synthesis difficulty are independent concerns:

- `HOMOPOLYMER_EXPRESSION_WARN_NT = 6`: expression-context and sequence
  quality warning.
- `HOMOPOLYMER_SYNTHESIS_WARN_NT = 8`: synthesis and manufacturing difficulty
  warning.

The expression-context threshold is stricter because shorter A/T-rich runs can
still matter during design review. The synthesis threshold is less strict and is
used by the profile rule engine's synthesis-oriented scan.

## Local GC Extremes

Local GC extreme scanning is implemented with a 50-nt sliding window and a
25-75% threshold band (a synthesis-hostility guard, not the global GC
optimization target). The current constants are useful defaults, but the
future migration path is to move local-GC threshold configuration into host
profile YAML. Host-specific profiles should define acceptable local GC
behavior rather than relying on global defaults.

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

## Repeat Pattern Rule

Repeat-pattern detection is implemented (`scan_repeats`, minimum length 15 nt)
and runs by default as part of the standard 9-scanner advisory set. It detects
perfect/tandem repeat structures that can affect synthesis, cloning, or
sequence stability, kept separate from homopolymer-only checks.

## AU-Rich Element and AT-Run Rules

AU-rich element (ARE) detection and AT-rich run detection (`min_length=6` nt
default) are both implemented and run by default as part of the standard
9-scanner advisory set, alongside the other rules described above.

## MoClo Overhang Rule

MoClo Level-0 overhang validity and internal-collision checks are implemented
in the construct-builder layer, not the RuleEngine. Unlike the other rules in
this document, this check is opt-in: it only runs when a construct template is
requested, and even then it reports warnings without halting the pipeline. It
is not part of the reported benchmark.
