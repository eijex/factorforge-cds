"""Greedy CAI-focused baseline: favors CAI, ignoring GC/assembly constraints.
Expected to score highest on CAI but fail multi-constraint pass rate."""
from __future__ import annotations
from factorforge.analysis.metrics import translate_dna, load_codon_usage_table

_WEIGHTS = load_codon_usage_table().codon_weights
_AA_TO_CODONS: dict[str, list[str]] = {}
for _cod in (a+b+c for a in "ACGT" for b in "ACGT" for c in "ACGT"):
    _aa = translate_dna(_cod)
    if _aa and _aa != "*":
        _AA_TO_CODONS.setdefault(_aa, []).append(_cod)


def greedy_cai_cds(protein: str) -> str:
    out = []
    for aa in protein.rstrip("*"):
        out.append(max(_AA_TO_CODONS[aa], key=lambda c: _WEIGHTS.get(c, 0.0)))
    return "".join(out)
