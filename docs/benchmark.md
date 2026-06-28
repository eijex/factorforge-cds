# FactorForge Benchmark

FactorForge v3.2.7 provides a reproducible in-silico benchmark foundation for CDS design quality.

## What is tested
- amino-acid identity preservation
- internal stop codon absence
- invalid codon absence
- length multiple of 3
- Type IIS / Golden Gate assembly compatibility
- CAI and GC as soft metrics
- baseline comparison against simple synonymous-codon methods

## What is not tested
This benchmark does not demonstrate improved protein expression, yield, folding, solubility, or wet-lab performance.

## Primary benchmark concept
The primary benchmark concept is multi-constraint pass rate, not CAI alone.

## Reproduction
See `benchmarks/README.md`.
