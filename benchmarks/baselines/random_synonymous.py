"""Random synonymous-codon baseline (lowest meaningful baseline)."""
from __future__ import annotations
import random
from factorforge.analysis.metrics import translate_dna

# Build AA -> [codons] from the standard table via translate_dna inverse.
_CODONS = [a+b+c for a in "ACGT" for b in "ACGT" for c in "ACGT"]
_AA_TO_CODONS: dict[str, list[str]] = {}
for _cod in _CODONS:
    _aa = translate_dna(_cod)
    if _aa and _aa != "*":
        _AA_TO_CODONS.setdefault(_aa, []).append(_cod)


def random_synonymous_cds(protein: str, seed: int = 320) -> str:
    rng = random.Random(seed)
    out = []
    for aa in protein.rstrip("*"):
        choices = _AA_TO_CODONS[aa]
        out.append(rng.choice(choices))
    return "".join(out)
