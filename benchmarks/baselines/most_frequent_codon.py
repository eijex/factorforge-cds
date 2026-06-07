"""Most-frequent-codon baseline: pick the highest-weight synonymous codon."""
from __future__ import annotations
from factorforge.analysis.metrics import translate_dna, load_codon_usage_table

_TABLE = load_codon_usage_table()
_WEIGHTS = _TABLE.codon_weights

_AA_TO_CODONS: dict[str, list[str]] = {}
for _cod in (a+b+c for a in "ACGT" for b in "ACGT" for c in "ACGT"):
    _aa = translate_dna(_cod)
    if _aa and _aa != "*":
        _AA_TO_CODONS.setdefault(_aa, []).append(_cod)


def most_frequent_codon_cds(protein: str) -> str:
    out = []
    for aa in protein.rstrip("*"):
        best = max(_AA_TO_CODONS[aa], key=lambda c: _WEIGHTS.get(c, 0.0))
        out.append(best)
    return "".join(out)
