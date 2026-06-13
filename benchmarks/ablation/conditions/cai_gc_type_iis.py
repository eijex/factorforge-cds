# benchmarks/ablation/conditions/cai_gc_type_iis.py
"""L4 ablation condition: CAI×GC-biased stochastic CDS with Type IIS site avoidance.

Combines CAI preference and GC steering (target 60%) in a single composite weight:
    weight = cai_weight × gc_proximity_weight

where:
    cai_weight          = codon frequency from codon usage table (proportional to CAI)
    gc_proximity_weight = 1.0 / (1.0 + |projected_gc - 60.0|)

This correctly implements "CAI + GC + Type IIS" ablation. Using GC weight alone
would make L4 equivalent to "GC + Type IIS", misrepresenting the ablation layer.
"""
from __future__ import annotations
import random
from factorforge.analysis.metrics import translate_dna, load_codon_usage_table

_WEIGHTS = load_codon_usage_table().codon_weights  # shared; codon frequency → CAI proxy; 1e-9 fallback keeps total > 0

_AA_TO_CODONS: dict[str, list[str]] = {}
for _cod in (a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"):
    _aa = translate_dna(_cod)
    if _aa and _aa != "*":
        _AA_TO_CODONS.setdefault(_aa, []).append(_cod)

_FORBIDDEN = ["GGTCTC", "GAGACC", "GAAGAC", "GTCTTC", "CGTCTC", "GAGACG"]
_GC_TARGET = 60.0


_CODON_GC: dict[str, int] = {
    c: sum(1 for b in c if b in "GC")
    for aa_codons in _AA_TO_CODONS.values()
    for c in aa_codons
}


def _cai_gc_biased_cds(protein: str, rng: random.Random) -> str:
    out: list[str] = []
    gc_count = 0  # running GC base count in partial CDS (O(n) instead of O(n²))
    for aa in protein:
        codons = _AA_TO_CODONS[aa]
        current_len = len(out) * 3
        weights = []
        for c in codons:
            cai_w = _WEIGHTS.get(c, 1e-9)
            test_gc = gc_count + _CODON_GC[c]
            test_len = current_len + 3
            gc = test_gc / test_len * 100.0
            gc_w = 1.0 / (1.0 + abs(gc - _GC_TARGET))
            weights.append(cai_w * gc_w)
        total = sum(weights)
        r = rng.random() * total
        acc = 0.0
        chosen = codons[-1]
        for c, w in zip(codons, weights):
            acc += w
            if acc >= r:
                chosen = c
                break
        out.append(chosen)
        gc_count += _CODON_GC[chosen]
    return "".join(out)


def ablation_cai_gc_type_iis_cds(
    protein: str, seed: int = 320, max_attempts: int = 50
) -> str:
    """CAI×GC-biased (target=60%) stochastic CDS with Type IIS site avoidance.

    Weight per codon = cai_weight × gc_proximity_weight.
    Each attempt uses seed + attempt_index so retries are independent.
    Returns last attempt if max_attempts exhausted.
    """
    protein = protein.rstrip("*")
    last = ""
    for attempt in range(max_attempts):
        rng = random.Random(seed + attempt)
        cds = _cai_gc_biased_cds(protein, rng)
        last = cds
        if not any(pat in cds for pat in _FORBIDDEN):
            return cds
    return last
