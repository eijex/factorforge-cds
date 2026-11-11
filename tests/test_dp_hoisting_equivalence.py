"""Equivalence, parity, and invariant tests for hoisted DP engine loops."""

from __future__ import annotations

import math

from Bio.Seq import Seq

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE, load_codon_usage_table
from factorforge.constraints.type_iis import (
    get_canonical_forbidden_motifs,
    resolve_restriction_site_motifs,
)
from factorforge.engines.dp_v2 import DPV2Optimizer
from factorforge.engines.dp_v2_4 import DPV240Optimizer
from factorforge.engines.structural_pause import StructuralPauseDPOptimizer

SYNTHETIC_PROTEIN_SHORT = "MASTWQHLEKDPVNRSGYF"
SYNTHETIC_PROTEIN_LONG_472 = (
    "MAKTNLFLFLIFSLLLSLSSAEVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNS"
    "GHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSSASTKGPSVFPL"
    "APSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHK"
    "PSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDG"
    "VEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDEL"
    "TKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNH"
    "YTQKSLSLSPGK"
)[:472]


def _reference_dpv2_uncached(protein: str, weights: dict[str, float], motifs, gc_min=0.40, gc_max=0.47):
    """Uncached reference DPV2 DP matching pre-hoist iteration order."""
    aa_to_codons: dict[str, list[str]] = {}
    for codon, aa in STANDARD_GENETIC_CODE.items():
        if aa != "*":
            aa_to_codons.setdefault(aa, []).append(codon)
    automaton = DPV2Optimizer(forbidden_motifs=motifs).automaton
    n_codons = len(protein)
    total_nt = n_codons * 3
    min_gc = int(math.ceil(gc_min * total_nt))
    max_gc = int(math.floor(gc_max * total_nt))
    layers = [{(0, 0): (0.0, None, None, None)}]
    for aa in protein:
        cur = layers[-1]
        nxt = {}
        for (g, s), (score, _, _, _) in cur.items():
            for codon in aa_to_codons[aa]:
                ns, bad = automaton.step_codon(s, codon)
                if bad:
                    continue
                ng = g + codon.count("G") + codon.count("C")
                w = weights.get(codon, 1e-4)
                nsc = score + math.log(max(w, 1e-6))
                k = (ng, ns)
                if k not in nxt:
                    nxt[k] = (nsc, codon, g, s)
                else:
                    es, ec, _, _ = nxt[k]
                    if nsc > es + 1e-12 or (abs(nsc - es) <= 1e-12 and ec is not None and codon < ec):
                        nxt[k] = (nsc, codon, g, s)
        layers.append(nxt)
    mid = (min_gc + max_gc) / 2.0
    in_band = [item for item in layers[-1].items() if min_gc <= item[0][0] <= max_gc]
    best = max(in_band, key=lambda item: (item[1][0], -abs(item[0][0] - mid), -item[0][0], -item[0][1]))
    codons = []
    cg, cs = best[0]
    for idx in range(n_codons, 0, -1):
        _, c, pg, ps = layers[idx][(cg, cs)]
        codons.append(c)
        cg, cs = pg, ps
    codons.reverse()
    return "".join(codons)


def test_hoisted_dp_v2_matches_uncached_reference_bit_for_bit():
    table = load_codon_usage_table()
    motifs = get_canonical_forbidden_motifs()
    expected = _reference_dpv2_uncached(SYNTHETIC_PROTEIN_SHORT, table.codon_weights, motifs, 0.35, 0.55)
    actual = DPV2Optimizer(forbidden_motifs=motifs).optimize(
        SYNTHETIC_PROTEIN_SHORT, table.codon_weights, target_gc_min=0.35, target_gc_max=0.55
    )["sequence"]
    assert actual == expected


def test_hoisted_structural_pause_matches_dpv2_at_zero_penalty():
    table = load_codon_usage_table()
    motifs = get_canonical_forbidden_motifs()
    base = DPV2Optimizer(forbidden_motifs=motifs).optimize(
        SYNTHETIC_PROTEIN_LONG_472, table.codon_weights, target_gc_min=0.35, target_gc_max=0.55
    )["sequence"]
    zero_3dp = StructuralPauseDPOptimizer(
        positions=[215, 250],
        target_adaptiveness=0.40,
        penalty_weight=0.0,
        structure_reference="synthetic-472",
        forbidden_motifs=motifs,
    ).optimize(
        SYNTHETIC_PROTEIN_LONG_472, table.codon_weights, target_gc_min=0.35, target_gc_max=0.55
    )["sequence"]
    assert zero_3dp == base
    assert str(Seq(zero_3dp).translate()) == SYNTHETIC_PROTEIN_LONG_472


def test_hoisted_dp_v2_4_472aa_deterministic_and_valid():
    table = load_codon_usage_table()
    motifs = resolve_restriction_site_motifs(("BsaI", "BsmBI", "BbsI", "SapI", "NotI"))
    kwargs = dict(
        protein_sequence=SYNTHETIC_PROTEIN_LONG_472,
        codon_weights=table.codon_weights,
        target_gc_min=0.35,
        target_gc_max=0.55,
        pause_positions=[215, 250, 251],
        pause_target_adaptiveness=0.35,
        pause_penalty_weight=3.5,
    )
    run_a = DPV240Optimizer(forbidden_motifs=motifs, **{k: v for k, v in kwargs.items() if k.startswith("pause_")}).optimize(
        protein_sequence=SYNTHETIC_PROTEIN_LONG_472,
        codon_weights=table.codon_weights,
        target_gc_min=0.35,
        target_gc_max=0.55,
    )
    run_b = DPV240Optimizer(forbidden_motifs=motifs, **{k: v for k, v in kwargs.items() if k.startswith("pause_")}).optimize(
        protein_sequence=SYNTHETIC_PROTEIN_LONG_472,
        codon_weights=table.codon_weights,
        target_gc_min=0.35,
        target_gc_max=0.55,
    )
    assert run_a["sequence"] == run_b["sequence"]
    assert len(run_a["sequence"]) == 472 * 3
    assert str(Seq(run_a["sequence"]).translate()) == SYNTHETIC_PROTEIN_LONG_472
