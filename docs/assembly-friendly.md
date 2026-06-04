# Assembly-Friendly Profile

The `assembly_friendly` profile is a design mode for sequences that must remain
compatible with downstream synthesis and cloning workflows. It is not just "no
restriction sites"; the intended scope is a holistic sequence-design mode that
considers codon adaptation, GC behavior, mRNA structure, synthesis constraints,
and cloning constraints together.

## Currently Supported

- BsaI and BpiI Type IIS restriction site avoidance through synonymous
  substitution.
- Adjusted scoring weights with lower CAI pressure and higher GC/MFE weight than
  the `balanced` profile.

The current `assembly_friendly` composite score uses:

```text
w_cai = 0.3
w_gc  = 0.4
w_mfe = 0.3
```

This lowers the preference for maximum codon adaptation so synonymous repair has
room to preserve cloning compatibility and global sequence quality.

## Not Yet Supported

- Local GC variance scoring.
- Repeat-pattern penalty.
- Synthesis vendor-specific constraints.
- Configurable Type IIS enzyme list in the profile interface; the current Type
  IIS set is hardcoded in the implementation.

These items are planned, not implemented. They should be added with regression
fixtures because each one can change candidate ranking and repair behavior.

## Design Intent

Assembly-friendly design should eventually evaluate both local and global
constraints:

- Avoid cloning-disruptive restriction sites.
- Preserve translation identity through synonymous edits.
- Keep CAI above acceptable guardrails.
- Keep global GC inside the selected profile band.
- Avoid local sequence patterns that are likely to complicate synthesis or
  assembly.

The profile should report what was scanned, what was repaired, what could not be
repaired, and which constraints were not evaluated by the current version.
