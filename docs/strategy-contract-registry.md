# Strategy Contract Registry

The strategy contract registry defines, per optimization strategy, **what the
strategy actually does** versus **what its name and exposure surface imply it
does**. It is the strategy-axis counterpart to the host-axis
[Host Profile Registry](host-profile-registry.md). The two compose: a host
profile says *which codon table / reference a host resolves to*, and a strategy
contract says *whether a given strategy actually consumes that host axis or is
pinned to a fixed reference regardless of host*.

The central principle is:

> Accepting a `host` argument does not mean a strategy is host-aware. Being
> exposed at an entry point does not mean it produces a meaningful result there.
> A strategy's contract must state its reference source, its real host-awareness,
> its per-entry-point exposure, and the strongest verification it has actually
> received — separately.

This registry is **documentation only**. It is not wired to runtime and does not
gate or alter any optimization output. It is a draft produced by a strategy/host
semantic contract audit and is intended as the basis for a future, optionally
runtime-enforced contract.

## Verification Level Vocabulary

Each strategy declares the **strongest** verification level it has actually
received (not aspirational). Levels are cumulative in rigor:

| Level | Meaning |
|-------|---------|
| `numeric_verified` | Output exists and named numeric metrics (CAI/GC/etc.) are computed and within expected ranges. |
| `behavioral_verified` | Output is produced and changes as expected across inputs/options (e.g. differs by host when it should). |
| `semantic_verified` | Output has been confirmed to *mean what the strategy name/claim says* — the reference it uses matches the claim, host-awareness matches the claim, and no silent substitution occurs. |
| `externally_validated` | Claims confirmed against an external/wet-lab or third-party reference beyond the codebase's own tests. |

## Required Strategy Fields

```yaml
id:
display_name:
reference_class:        # general_host_table | dedicated_golden_set | dp_hardcoded_table | delegates_to | host_table_plus_fixed_constant
reference_asset:        # concrete file(s) / constant(s) the strategy reads
code_location:          # primary implementing function(s)
host_aware:             # true | false | partial   (does output change by host?)
host_aware_note:        # why, in one line
exposure:               # per-entry-point: supported | rejected | hidden | silently_dropped | not_applicable | unguarded | warned
  python_api_profile:
  python_api_dp:
  cli_profile:
  cli_compare_profiles:
  cli_objective_dp:
  rest_explicit:
  rest_implicit:
  rest_compare:
  rest_batch:
  web_ui_default_host:
  web_ui_nondefault_host:
  get_metadata_listed:
verification_level:
known_gaps:
notes:
```

`exposure` value meanings:

- `supported` — works and returns a meaningful result.
- `rejected` — explicitly refused with an error (HTTP 400 / `UsageError` / `ValueError`).
- `hidden` — present but intentionally non-selectable (e.g. disabled web UI radio).
- `silently_dropped` — the field is accepted at the wire level but ignored, returning HTTP 200 without acting on it. **This is the most dangerous value** and must be called out explicitly.
- `unguarded` — reachable with no host/strategy compatibility guard, so an unsupported combination silently produces fixed-reference output under another name.
- `not_applicable` — the entry point has no concept for this strategy.
- `warned` — reachable with no rejection, but a log warning discloses the host-invariant substitution (chosen over rejection so the library boundary stays consistent with the host-invariant-by-design contract and existing tests).

## Strategy Contracts

### balanced

```yaml
id: balanced
display_name: Balanced (CAI + GC balance)
reference_class: general_host_table
reference_asset: "{host}_codons.json (resolved from the active host)"
code_location: "engines/profile/rules/reverse_translator.py::_balanced_translate"
host_aware: true
host_aware_note: Codon table is selected per host; output differs by host (test-confirmed).
exposure:
  python_api_profile: supported
  python_api_dp: not_applicable
  cli_profile: supported
  cli_compare_profiles: supported
  cli_objective_dp: not_applicable
  rest_explicit: supported
  rest_implicit: supported
  rest_compare: supported
  rest_batch: supported
  web_ui_default_host: supported
  web_ui_nondefault_host: supported
  get_metadata_listed: true
verification_level: behavioral_verified
known_gaps: []
notes: Reference (default) host-aware strategy; the template the others are compared against.
```

### gc_target

```yaml
id: gc_target
display_name: GC Target
reference_class: host_table_plus_fixed_constant
reference_asset: "{host}_codons.json (codon choice) + GC_OPT_MID=60.0 (scoring.py:22, default target when target_gc unset)"
code_location: "engines/profile/rules/reverse_translator.py::_gc_target_translate"
host_aware: partial
host_aware_note: Codon selection is host-aware, but the default GC target is a single N. benthamiana-calibrated constant (already flagged TODO in code at reverse_translator.py:453-454).
exposure:
  python_api_profile: supported
  python_api_dp: not_applicable
  cli_profile: supported
  cli_compare_profiles: supported
  cli_objective_dp: rejected   # see dual-namespace note below
  rest_explicit: supported
  rest_implicit: supported
  rest_compare: supported
  rest_batch: supported
  web_ui_default_host: supported
  web_ui_nondefault_host: supported   # BY-2 fallback target
  get_metadata_listed: true
verification_level: behavioral_verified
known_gaps:
  - Default GC target is host-invariant; per-host GC profiles not yet sourced from the active host.
  - The string "gc_target" is also a declared CLI --objective value for the DP engine, where it ALWAYS raises (see Dual-Namespace Collision section).
notes: Pass an explicit target_gc to override the host-invariant default.
```

### assembly_friendly

```yaml
id: assembly_friendly
display_name: Assembly-Friendly (Golden Gate / MoClo)
reference_class: delegates_to
reference_asset: "{host}_codons.json (via _balanced_translate, preferred_ratio=0.6)"
code_location: "engines/profile/rules/reverse_translator.py::_assembly_friendly_translate -> _balanced_translate"
host_aware: true
host_aware_note: Inherits balanced's host-aware codon selection; adds BsaI/BpiI Type IIS site-avoidance retries.
exposure:
  python_api_profile: supported
  python_api_dp: not_applicable
  cli_profile: supported
  cli_compare_profiles: supported
  cli_objective_dp: not_applicable
  rest_explicit: supported
  rest_implicit: supported
  rest_compare: supported
  rest_batch: supported
  web_ui_default_host: supported
  web_ui_nondefault_host: supported
  get_metadata_listed: true
verification_level: behavioral_verified
known_gaps:
  - Does not yet score local GC-window uniformity or repeat patterns (documented in docstring).
notes: Stochastic retry up to max_attempts; may return a site-containing sequence with a warning if unresolved.
```

### high_cai

```yaml
id: high_cai
display_name: High CAI
reference_class: dedicated_golden_set
reference_asset: "nbenthamiana_golden_set.json (load_golden_set() — no host parameter; always this one file)"
code_location: "engines/profile/rules/reverse_translator.py::_high_cai_translate; engines/profile/utils.py::load_golden_set"
host_aware: false
host_aware_note: Intentionally host-invariant by design — anchored to the N. benthamiana high-expression golden set. No BY-2 golden set exists. This is a designed boundary, not a bug.
exposure:
  python_api_profile: warned      # logs a warning instead of rejecting (host-invariance is by design; rejecting would conflict with existing tests)
  python_api_dp: not_applicable
  cli_profile: supported          # rejected for non-default host
  cli_compare_profiles: supported # rejected for non-default host incl. in --compare-profiles list
  cli_objective_dp: rejected      # always raises (dual-namespace), and not a click.Choice value
  rest_explicit: supported        # HTTP 400 for non-default host
  rest_implicit: not_applicable   # implicit path never resolves to high_cai
  rest_compare: rejected          # any host/host_profile field is rejected (HTTP 400)
  rest_batch: rejected            # same as above
  web_ui_default_host: supported
  web_ui_nondefault_host: rejected # radio disabled
  get_metadata_listed: true        # listed host-independently — discoverability gap (still open)
verification_level: semantic_verified
known_gaps:
  - Listed in GET /api/optimize supported_profiles regardless of host (discoverability gap).
notes: The only profile-engine strategy that is host-invariant by design. It is rejected for non-default hosts at the REST/CLI/web surfaces; the library boundary warns instead of silently staying quiet (python_api_profile), and /api/optimize/compare·/batch reject any host field outright instead of silently dropping it.
```

### ramp

```yaml
id: ramp
display_name: 5' Ramp (N-terminal codon ramp)
reference_class: delegates_to
reference_asset: "{host}_codons.json (via _balanced_translate + _apply_nterminal_ramp, both host-table based)"
code_location: "engines/profile/rules/reverse_translator.py::_ramp_translate -> _balanced_translate + _apply_nterminal_ramp"
host_aware: true
host_aware_note: Implemented and host-aware (output differs by host, confirmed by direct execution). NOT unimplemented — earlier drafts misclassified it.
exposure:
  python_api_profile: supported
  python_api_dp: not_applicable
  cli_profile: unguarded          # --profile has no click.Choice constraint; passes through
  cli_compare_profiles: unguarded # likewise passes through (directly confirmed)
  cli_objective_dp: not_applicable
  rest_explicit: rejected         # not in VALID_PROFILES -> "Invalid profile"
  rest_implicit: not_applicable
  rest_compare: rejected          # not in VALID_PROFILES
  rest_batch: rejected            # not in VALID_PROFILES
  web_ui_default_host: hidden     # radio disabled
  web_ui_nondefault_host: hidden  # radio disabled
  get_metadata_listed: false      # not in supported_profiles
verification_level: behavioral_verified
known_gaps:
  - Three-way exposure asymmetry — the SAME implemented feature is hidden (web), rejected (REST), and passes through unguarded (CLI --profile). No single source of truth for "is ramp public?".
  - ramp_codons=50 default may cover an entire short CDS (TODO in code at reverse_translator.py:557-561).
notes: Genuinely implemented; exposure policy is inconsistent across entry points rather than uniformly "private".
```

### viral_delivery

```yaml
id: viral_delivery
display_name: Viral Delivery (TRV)
reference_class: delegates_to
reference_asset: "{host}_codons.json (codon choice == balanced; NO viral_delivery-specific reverse-translate branch) + scoring.py:89-93 TRV constants gc_opt=47.5/gc_min=37.5/gc_max=57.5 (host-invariant)"
code_location: "engines/profile/rules/reverse_translator.py:358-359 (calls _balanced_translate); engines/profile/scoring.py SCORING_CONFIGS['viral_delivery']"
host_aware: true
host_aware_note: Codon SELECTION is host-aware (identical mechanism to balanced); the 'viral_delivery-ness' lives only in the host-invariant TRV scoring constants. There is no viral_delivery-specific reverse-translation logic.
exposure:
  python_api_profile: supported
  python_api_dp: not_applicable
  cli_profile: unguarded          # passes through (directly confirmed)
  cli_compare_profiles: unguarded
  cli_objective_dp: not_applicable
  rest_explicit: rejected         # not in VALID_PROFILES
  rest_implicit: not_applicable
  rest_compare: rejected
  rest_batch: rejected
  web_ui_default_host: hidden     # radio disabled (web/index.html:316)
  web_ui_nondefault_host: hidden
  get_metadata_listed: false
verification_level: numeric_verified
known_gaps:
  - Same three-way exposure asymmetry as ramp (hidden/rejected/unguarded).
  - 'Strategy' differs from balanced only at the scoring layer, not the sequence-generation layer — the name implies more divergence than the reverse-translator actually applies.
notes: TRV scoring constants are host-invariant; sequence candidates are byte-identical in mechanism to balanced.
```

### feasibility_best

```yaml
id: feasibility_best
display_name: Feasibility Best (DP engine default)
reference_class: dp_hardcoded_table
reference_asset: "nbenthamiana_codons.json (load_codon_usage_table() called with no path argument; always the default file)"
code_location: "api/optimize.py::optimize_feasibility_best (line 771); cli/main.py::_build_dp_result (line 46); analysis/feasibility.py::analyze_feasibility"
host_aware: false
host_aware_note: The host argument is not used for codon-table selection anywhere in the feasibility_best body. Host-invariance is structural, not guarded — the guard lives in a separate layer.
exposure:
  python_api_profile: not_applicable
  python_api_dp: unguarded        # host arg accepted but unused in table selection
  cli_profile: not_applicable
  cli_compare_profiles: not_applicable
  cli_objective_dp: supported     # the ONLY working --objective value
  rest_explicit: supported        # HTTP 400 for non-default host
  rest_implicit: supported        # for non-default host, disclosed resolve to 'balanced'
  rest_compare: not_applicable
  rest_batch: not_applicable
  web_ui_default_host: supported  # default selection
  web_ui_nondefault_host: rejected # radio disabled
  get_metadata_listed: true        # supported_objectives: ["feasibility_best"]
verification_level: semantic_verified
known_gaps:
  - python_api_dp is unguarded (same library-boundary gap as high_cai).
  - The /api/optimize response for feasibility_best bundles a host-blind candidate (feasibility_best itself) with host-aware (gc_target) and host-invariant-by-design (high_cai) comparison candidates in one payload — three different host contracts under one response.
notes: DP default objective. The dual-namespace strings gc_target/high_cai are declared as --objective choices but always raise.
```

## Dual-Namespace Collision

The strings `gc_target` and `high_cai` exist in **two unrelated namespaces**:

| String | Namespace | Behavior |
|--------|-----------|----------|
| `gc_target` | profile engine (`--profile`, REST `profile`) | Works — `_gc_target_translate` |
| `gc_target` | DP `--objective` (declared in `click.Choice`) | **Always raises** `ValueError: DP engine currently supports --objective feasibility_best.` (exit 1) |
| `high_cai` | profile engine (`--profile`, REST `profile`) | Works (host-guarded) |
| `high_cai` | DP `--objective` (declared in `click.Choice`) | **Always raises** (same message) |

The DP `--objective` `click.Choice` advertises three values but only
`feasibility_best` ever succeeds. `docs/how-it-works.md`'s "Design Objectives"
table mirrors this declaration and so describes two never-working values as if
they were functional objectives.

## Cross-Cutting Findings

1. **Library boundary has a soft guard, not a hard one.** The
   `UNSUPPORTED_STRATEGY_HOST_COMBINATION` rejection is implemented only at
   the REST, CLI, and web layers. A direct
   `RuleBasedOptimizer().optimize(profile="high_cai", host="ntabacum")` call
   still returns N. benthamiana golden-set output (this is the host-invariant
   design, not changed) — a `logger.warning` is emitted at this call site so
   the silent substitution is at least observable, without rejecting
   (rejecting would conflict with the host-invariant-by-design contract and
   the existing `test_ntabacum_baseline_unchanged[high_cai]` /
   `test_high_cai_is_host_invariant_by_design` tests).

2. **`/api/optimize/compare` and `/api/optimize/batch` reject `host` outright.**
   Both handlers reject any request body containing `host` or `host_profile`
   with HTTP 400 (`HOST_NOT_SUPPORTED_ON_ENDPOINT`) instead of silently
   ignoring it. This is an API contract change: a request that previously got
   a misleading HTTP 200 now gets an explicit 400.

3. **Exposure policy is not single-sourced.** `ramp` and `viral_delivery` are
   each hidden in the web UI, rejected by REST `VALID_PROFILES`, and passed
   through unguarded by CLI `--profile` (which has no `click.Choice`). The same
   implemented feature has three different public answers depending on the door
   used.

4. **Seven strategies, four host-contract patterns** — not a binary
   host-aware/invariant split:
   - Fully host-aware: `balanced`, `assembly_friendly`, `ramp`
   - Host-aware codon choice + host-invariant constant: `gc_target` (default GC),
     `viral_delivery` (TRV scoring)
   - Host-invariant by design (dedicated reference): `high_cai`
   - Host-invariant structurally (host arg unused): `feasibility_best`
