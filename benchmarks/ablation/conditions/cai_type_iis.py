# benchmarks/ablation/conditions/cai_type_iis.py
"""L3 ablation condition: CAI-weighted stochastic CDS with Type IIS site avoidance.

No GC target constraint. Stochastic retry allows Type IIS-free solutions.
"""
from __future__ import annotations
import random
from factorforge.analysis.metrics import translate_dna, load_codon_usage_table

_WEIGHTS = load_codon_usage_table().codon_weights  # shared; codon frequency → CAI proxy
_AA_TO_CODONS: dict[str, list[str]] = {}
for _cod in (a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"):
    _aa = translate_dna(_cod)
    if _aa and _aa != "*":
        _AA_TO_CODONS.setdefault(_aa, []).append(_cod)

_FORBIDDEN = ["GGTCTC", "GAGACC", "GAAGAC", "GTCTTC", "CGTCTC", "GAGACG"]


def _cai_weighted_cds(protein: str, rng: random.Random) -> str:
    out = []
    for aa in protein:
        codons = _AA_TO_CODONS[aa]
        weights = [_WEIGHTS.get(c, 1e-9) for c in codons]  # 1e-9 fallback keeps total > 0
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
    return "".join(out)


def ablation_cai_type_iis_cds(
    protein: str, seed: int = 320, max_attempts: int = 50
) -> str:
    """CAI-weighted stochastic CDS with Type IIS site avoidance.

    Each attempt uses seed + attempt_index so retries are independent.
    Returns the last attempt if all max_attempts are exhausted.
    """
    protein = protein.rstrip("*")
    last = ""
    for attempt in range(max_attempts):
        rng = random.Random(seed + attempt)
        cds = _cai_weighted_cds(protein, rng)
        last = cds
        if not any(pat in cds for pat in _FORBIDDEN):
            return cds
    return last
