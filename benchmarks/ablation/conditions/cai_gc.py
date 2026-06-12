# benchmarks/ablation/conditions/cai_gc.py
"""L2 ablation condition: CAI + GC target (60%).

Wraps the production gc_target profile. Deterministic: same protein always
produces the same CDS. No Type IIS constraint active.
"""
from __future__ import annotations
from factorforge.engines.profile.optimizer import RuleBasedOptimizer as _Opt

_opt = _Opt()  # shared instance; construction is expensive (loads codon tables)


def ablation_cai_gc_cds(protein: str, seed: int = 320) -> str:
    """Return a GC-targeted CDS for the given protein sequence.

    Uses the production gc_target profile. seed is accepted for interface
    consistency but gc_target is deterministic — the value has no effect.
    """
    return _opt.optimize(protein.rstrip("*"), profile="gc_target", seed=seed).sequence
