# Math Specification: DP v2.1.1 local composition guard

**Module:** `factorforge.engines.dp_v2_1_1`  
**Engine version:** `2.1.1-dev`  
**Status:** development candidate; calibration complete, holdout pending

DP v2.1.1 retains the DP v2.1 state `(i, h, q)`, where `i` is the processed
codon count, `h` is cumulative GC count, and `q` is the configured-motif
Aho-Corasick state. No local-GC state dimension is added.

At `K = 15`, the active layer is filtered by

```text
ceil(45 * 0.20) <= h_15 <= floor(45 * 0.30)
9 <= h_15 <= 13
```

The requested interval is intersected with the prefix's synonymous GC envelope.
If motif or homopolymer vetoes make the resulting interval unreachable, the
nearest reachable GC layer is retained and `initiation_gc_status` reports
`AUTOMATON_REACHABILITY_RELAXED`. Intrinsically high/low-GC prefixes report
`SYNONYMOUS_ENVELOPE_CLAMPED`. Sequences shorter than 15 codons report
`NOT_APPLIED_SHORT_SEQUENCE` because no layer-15 state exists.

The automaton is compiled with all configured forbidden motifs plus the four
6-mers `AAAAAA`, `TTTTTT`, `GGGGGG`, and `CCCCCC`; therefore a returned path
cannot contain a homopolymer longer than 5 nt. Translation preservation remains
exact because every transition selects a codon from the current amino acid's
synonymous set.

With at most six synonymous codons per amino acid, time complexity is
`O(N * |H| * |Q|)` up to that constant branching factor. Sparse predecessor
storage is `O(N * |H| * |Q|)`. Full-prefix path ranks provide deterministic
tie-breaking.

The 30/45-nt and 50-bp GC metrics are post-generation composition measurements.
The 5′ MFE field is produced by a separate ViennaRNA evaluation call and is not
part of the DP objective. Its status, reason, window length, and available
upstream-context length are returned so missing dependency or context is not
misrepresented as a numerical result.
