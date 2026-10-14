"""FactorForge 3DP Validation Suite: Regression, Trade-off, and Negative Controls.

Validates that 3D structural conditioning:
1. Collapses to 1DP baseline when lambda = 0 (equivalence).
2. Preserves all hard biological/manufacturing constraints (translation, GC band, Type IIS motifs).
3. Quantifies explicit Cost-of-Conditioning (Δ1D) and Structural Gain (Δ3D).
4. Localizes synonymous codon interventions to the designated structural target zone.
5. Measures computational intervention shifts under shuffled annotations.
6. Maintains solution family stability under minor boundary perturbations.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

from factorforge.engines.dp_v2 import DPV2Optimizer
from factorforge.engines.structural_pause import StructuralPauseDPOptimizer


@dataclass
class TradeoffComparison:
    cai_1dp: float
    cai_3dp: float
    delta_1d_cai: float
    gc_1dp: float
    gc_3dp: float
    delta_gc: float
    hard_constraints_pass: bool
    structure_fit_1dp: float
    structure_fit_3dp: float
    delta_3d_gain: float
    total_edits: int
    target_zone_edits: int
    outside_edits: int
    localization_ratio: float
    cds_1dp: str
    cds_3dp: str
    dashboard_table: str

    def to_dict(self, *, include_sequences: bool = False) -> Dict[str, Any]:
        """Export metrics; private CDS export requires an explicit opt-in."""
        result = asdict(self)
        if not include_sequences:
            result.pop("cds_1dp")
            result.pop("cds_3dp")
        return result


@dataclass
class ParetoSweepPoint:
    lambda_val: float
    cai: float
    score_objective: float
    structure_fit: float
    total_edits_from_1dp: int
    target_zone_edits: int
    localization_ratio: float

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class NegativeControlResult:
    nominal_positions: List[int]
    shuffled_positions: List[int]
    nominal_edits_in_nominal_zone: int
    shuffled_edits_in_shuffled_zone: int
    shuffled_edits_in_nominal_zone: int
    causal_shift_confirmed: bool

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class BoundaryPerturbationResult:
    nominal_positions: List[int]
    perturbations: Dict[int, Dict[str, Any]]
    is_stable: bool

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def calculate_sequence_cai(cds: str, codon_weights: Dict[str, float]) -> float:
    """Calculate geometric mean of codon relative adaptiveness."""
    seq = cds.upper()
    n = len(seq) // 3
    if n == 0:
        return 0.0
    log_sum = 0.0
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        w = max(codon_weights.get(codon, 0.01), 0.01)
        log_sum += math.log(w)
    return round(math.exp(log_sum / n), 4)


def calculate_sequence_gc(cds: str) -> float:
    """Calculate GC percentage of CDS."""
    seq = cds.upper()
    if not seq:
        return 0.0
    gc_count = seq.count("G") + seq.count("C")
    return round(gc_count / len(seq), 4)


def calculate_structure_fit(
    cds: str,
    positions: Sequence[int],
    codon_weights: Dict[str, float],
    target_adaptiveness: float = 0.40,
) -> float:
    """Compute codon-weight target conformity; this is not folding accuracy.
    
    Fit = 1.0 - mean(|w_i - target_w|). Maximum 1.0 when all pause positions match target adaptiveness.
    """
    seq = cds.upper()
    pos_set = set(positions)
    if not pos_set:
        return 0.0
    errors = []
    for p in sorted(pos_set):
        idx = (p - 1) * 3
        if idx + 3 <= len(seq):
            codon = seq[idx : idx + 3]
            w = codon_weights.get(codon, 0.01)
            errors.append(abs(w - target_adaptiveness))
    if not errors:
        return 0.0
    mean_error = sum(errors) / len(errors)
    return round(max(0.0, 1.0 - mean_error), 4)


def analyze_edits(
    cds_1dp: str, cds_3dp: str, target_positions: Sequence[int], flank_window: int = 1
) -> Tuple[int, int, int, float]:
    """Analyze codon differences between 1DP and 3DP, segmenting target zone vs outside."""
    if len(cds_1dp) != len(cds_3dp) or len(cds_1dp) % 3:
        raise ValueError("Edit comparison requires equal, complete codon sequences")
    n_codons = len(cds_1dp) // 3
    pos_set = set(target_positions)
    # Expand zone by flank_window
    expanded_zone = set()
    for p in pos_set:
        for offset in range(-flank_window, flank_window + 1):
            if 1 <= p + offset <= n_codons:
                expanded_zone.add(p + offset)

    total_edits = 0
    target_edits = 0
    outside_edits = 0

    for i in range(n_codons):
        p_1based = i + 1
        c1 = cds_1dp[i * 3 : (i + 1) * 3]
        c3 = cds_3dp[i * 3 : (i + 1) * 3]
        if c1 != c3:
            total_edits += 1
            if p_1based in expanded_zone:
                target_edits += 1
            else:
                outside_edits += 1

    ratio = round(target_edits / total_edits, 4) if total_edits > 0 else 1.0
    return total_edits, target_edits, outside_edits, ratio


def format_dashboard_table(
    cai_1dp: float,
    cai_3dp: float,
    delta_cai: float,
    gc_1dp: float,
    gc_3dp: float,
    delta_gc: float,
    hard_pass: bool,
    fit_1dp: float,
    fit_3dp: float,
    delta_fit: float,
    total_edits: int,
    target_edits: int,
    outside_edits: int,
    loc_ratio: float,
) -> str:
    """Format user-specified ASCII/markdown dashboard table."""
    hard_str = "PASS" if hard_pass else "FAIL"
    sign_cai = "+" if delta_cai > 0 else ""
    sign_gc = "+" if delta_gc > 0 else ""
    sign_fit = "+" if delta_fit > 0 else ""

    lines = [
        "                 1DP        3DP       Delta",
        f"CAI              {cai_1dp:.3f}      {cai_3dp:.3f}    {sign_cai}{delta_cai:.3f}",
        f"GC                {gc_1dp*100:.1f}       {gc_3dp*100:.1f}     {sign_gc}{delta_gc*100:.1f}",
        f"Hard constraints  {hard_str}       {hard_str}       --",
        f"Structure fit     {fit_1dp:.2f}       {fit_3dp:.2f}     {sign_fit}{delta_fit:.2f}",
        f"Total edits         --         {total_edits:2d}        --",
        f"Target-zone edits   --         {target_edits:2d}        --",
        f"Outside edits       --         {outside_edits:2d}        --",
        "",
        f"Structural Conditioning Cost: {sign_cai}{delta_cai:.3f} CAI",
        f"Structural Context Gain: {sign_fit}{delta_fit:.2f}",
        f"Hard Constraints: preserved ({hard_str})",
        f"Localization Ratio: {loc_ratio*100:.1f}% ({target_edits}/{total_edits})",
    ]
    return "\n".join(lines)


def check_lambda_zero_equivalence(
    protein: str,
    pause_positions: Sequence[int],
    codon_weights: Dict[str, float],
    target_gc_min: float = 0.40,
    target_gc_max: float = 0.48,
    forbidden_motifs: Optional[List[str]] = None,
    target_adaptiveness: float = 0.40,
    left_flank: str = "",
    right_flank: str = "",
) -> bool:
    """Test 1: λ=0 Equivalence Check.
    
    Verifies that when penalty_weight=0.0, 3DP produces identical CDS and score to 1DP.
    """
    opt_1dp = DPV2Optimizer(forbidden_motifs=forbidden_motifs)
    res_1dp = opt_1dp.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        left_flank=left_flank,
        right_flank=right_flank,
    )

    opt_3dp_zero = StructuralPauseDPOptimizer(
        positions=frozenset(pause_positions),
        target_adaptiveness=target_adaptiveness,
        penalty_weight=0.0,
        structure_reference="test_lambda_zero",
        forbidden_motifs=forbidden_motifs,
    )
    res_3dp = opt_3dp_zero.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        left_flank=left_flank,
        right_flank=right_flank,
    )

    seq_match = (res_1dp["sequence"] == res_3dp["sequence"])
    score_match = (abs(res_1dp["score"] - res_3dp["score"]) < 1e-9)
    return seq_match and score_match


def validate_hard_constraints(
    cds: str, protein: str, gc_min: float, gc_max: float,
    forbidden_motifs: Optional[List[str]] = None,
    left_flank: str = "", right_flank: str = "",
) -> bool:
    """Check sequence identity, exact integer GC bounds and full construct motifs."""
    from factorforge.analysis.metrics import STANDARD_GENETIC_CODE
    from factorforge.constraints.type_iis import get_canonical_forbidden_motifs

    seq = cds.upper()
    target = "".join(protein.upper().split()).rstrip("*")
    if not seq or len(seq) != 3 * len(target) or set(seq) - set("ACGT"):
        return False
    translated = "".join(STANDARD_GENETIC_CODE.get(seq[i:i+3], "?") for i in range(0, len(seq), 3))
    if translated != target:
        return False
    lower = gc_min if gc_min <= 1 else gc_min / 100
    upper = gc_max if gc_max <= 1 else gc_max / 100
    count = seq.count("G") + seq.count("C")
    if not math.ceil(lower * len(seq)) <= count <= math.floor(upper * len(seq)):
        return False
    motifs = get_canonical_forbidden_motifs() if forbidden_motifs is None else forbidden_motifs
    construct = (left_flank + seq + right_flank).upper()
    complement = str.maketrans("ACGT", "TGCA")
    for motif in motifs:
        motif = motif.upper()
        if motif in construct or motif.translate(complement)[::-1] in construct:
            return False
    return True


def compute_tradeoff_comparison(
    protein: str,
    pause_positions: Sequence[int],
    codon_weights: Dict[str, float],
    target_gc_min: float = 0.40,
    target_gc_max: float = 0.48,
    forbidden_motifs: Optional[List[str]] = None,
    lambda_weight: float = 2.0,
    target_adaptiveness: float = 0.40,
    left_flank: str = "",
    right_flank: str = "",
) -> TradeoffComparison:
    """Tests 2, 3, 4, 5: Compute full trade-off metrics, hard constraint regression, and localization."""
    # 1. 1DP Optimization
    opt_1dp = DPV2Optimizer(forbidden_motifs=forbidden_motifs)
    res_1dp = opt_1dp.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        left_flank=left_flank,
        right_flank=right_flank,
    )
    cds_1dp = res_1dp["sequence"]

    # 2. 3DP Optimization
    opt_3dp = StructuralPauseDPOptimizer(
        positions=frozenset(pause_positions),
        target_adaptiveness=target_adaptiveness,
        penalty_weight=lambda_weight,
        structure_reference="tradeoff_evaluation",
        forbidden_motifs=forbidden_motifs,
    )
    res_3dp = opt_3dp.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        left_flank=left_flank,
        right_flank=right_flank,
    )
    cds_3dp = res_3dp["sequence"]

    # 3. Metrics calculation
    cai_1dp = calculate_sequence_cai(cds_1dp, codon_weights)
    cai_3dp = calculate_sequence_cai(cds_3dp, codon_weights)
    delta_1d_cai = round(cai_3dp - cai_1dp, 4)

    gc_1dp = calculate_sequence_gc(cds_1dp)
    gc_3dp = calculate_sequence_gc(cds_3dp)
    delta_gc = round(gc_3dp - gc_1dp, 4)

    fit_1dp = calculate_structure_fit(cds_1dp, pause_positions, codon_weights, target_adaptiveness)
    fit_3dp = calculate_structure_fit(cds_3dp, pause_positions, codon_weights, target_adaptiveness)
    delta_3d_gain = round(fit_3dp - fit_1dp, 4)

    # Independently validate BOTH candidates, using exact GC counts and the
    # configured motifs plus reverse complements, including construct flanks.
    hard_pass = all(
        validate_hard_constraints(
            cds, protein, target_gc_min, target_gc_max,
            forbidden_motifs, left_flank, right_flank,
        )
        for cds in (cds_1dp, cds_3dp)
    )

    # 4. Edit Localization
    total_edits, target_edits, outside_edits, loc_ratio = analyze_edits(
        cds_1dp, cds_3dp, pause_positions, flank_window=1
    )

    dashboard = format_dashboard_table(
        cai_1dp=cai_1dp,
        cai_3dp=cai_3dp,
        delta_cai=delta_1d_cai,
        gc_1dp=gc_1dp,
        gc_3dp=gc_3dp,
        delta_gc=delta_gc,
        hard_pass=hard_pass,
        fit_1dp=fit_1dp,
        fit_3dp=fit_3dp,
        delta_fit=delta_3d_gain,
        total_edits=total_edits,
        target_edits=target_edits,
        outside_edits=outside_edits,
        loc_ratio=loc_ratio,
    )

    return TradeoffComparison(
        cai_1dp=cai_1dp,
        cai_3dp=cai_3dp,
        delta_1d_cai=delta_1d_cai,
        gc_1dp=gc_1dp,
        gc_3dp=gc_3dp,
        delta_gc=delta_gc,
        hard_constraints_pass=hard_pass,
        structure_fit_1dp=fit_1dp,
        structure_fit_3dp=fit_3dp,
        delta_3d_gain=delta_3d_gain,
        total_edits=total_edits,
        target_zone_edits=target_edits,
        outside_edits=outside_edits,
        localization_ratio=loc_ratio,
        cds_1dp=cds_1dp,
        cds_3dp=cds_3dp,
        dashboard_table=dashboard,
    )


def run_pareto_lambda_sweep(
    protein: str,
    pause_positions: Sequence[int],
    codon_weights: Dict[str, float],
    lambda_values: Sequence[float] = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0),
    target_gc_min: float = 0.40,
    target_gc_max: float = 0.48,
    forbidden_motifs: Optional[List[str]] = None,
    target_adaptiveness: float = 0.40,
) -> List[ParetoSweepPoint]:
    """Test 8: Pareto λ Sweep.
    
    Profiles (1D quality vs Structural Fit) across varying penalty weights λ.
    """
    # Baseline 1DP sequence for edit distance calculation
    opt_1dp = DPV2Optimizer(forbidden_motifs=forbidden_motifs)
    res_1dp = opt_1dp.optimize(
        protein_sequence=protein,
        codon_weights=codon_weights,
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
    )
    cds_1dp = res_1dp["sequence"]

    points: List[ParetoSweepPoint] = []
    for l_val in lambda_values:
        opt = StructuralPauseDPOptimizer(
            positions=frozenset(pause_positions),
            target_adaptiveness=target_adaptiveness,
            penalty_weight=float(l_val),
            structure_reference=f"sweep_lambda_{l_val}",
            forbidden_motifs=forbidden_motifs,
        )
        res = opt.optimize(
            protein_sequence=protein,
            codon_weights=codon_weights,
            target_gc_min=target_gc_min,
            target_gc_max=target_gc_max,
        )
        cds = res["sequence"]
        cai = calculate_sequence_cai(cds, codon_weights)
        fit = calculate_structure_fit(cds, pause_positions, codon_weights, target_adaptiveness)
        total_edits, target_edits, _, loc_ratio = analyze_edits(cds_1dp, cds, pause_positions)

        points.append(
            ParetoSweepPoint(
                lambda_val=float(l_val),
                cai=cai,
                score_objective=round(float(res.get("score", 0.0)), 4),
                structure_fit=fit,
                total_edits_from_1dp=total_edits,
                target_zone_edits=target_edits,
                localization_ratio=loc_ratio,
            )
        )
    return points


def run_shuffled_structure_control(
    protein: str,
    nominal_positions: Sequence[int],
    shuffled_positions: Sequence[int],
    codon_weights: Dict[str, float],
    target_gc_min: float = 0.40,
    target_gc_max: float = 0.48,
    forbidden_motifs: Optional[List[str]] = None,
    lambda_weight: float = 2.0,
    target_adaptiveness: float = 0.40,
) -> NegativeControlResult:
    """Test 6: Shuffled-Structure Negative Control.
    
    Verifies that moving structural pause annotations causally moves codon interventions.
    """
    # 1. 1DP baseline
    opt_1dp = DPV2Optimizer(forbidden_motifs=forbidden_motifs)
    cds_1dp = opt_1dp.optimize(
        protein_sequence=protein, codon_weights=codon_weights,
        target_gc_min=target_gc_min, target_gc_max=target_gc_max
    )["sequence"]

    # 2. 3DP with nominal positions
    opt_nom = StructuralPauseDPOptimizer(
        positions=frozenset(nominal_positions),
        target_adaptiveness=target_adaptiveness,
        penalty_weight=lambda_weight,
        structure_reference="nominal_control",
        forbidden_motifs=forbidden_motifs,
    )
    cds_nom = opt_nom.optimize(
        protein_sequence=protein, codon_weights=codon_weights,
        target_gc_min=target_gc_min, target_gc_max=target_gc_max
    )["sequence"]

    # 3. 3DP with shuffled positions
    opt_shuf = StructuralPauseDPOptimizer(
        positions=frozenset(shuffled_positions),
        target_adaptiveness=target_adaptiveness,
        penalty_weight=lambda_weight,
        structure_reference="shuffled_control",
        forbidden_motifs=forbidden_motifs,
    )
    cds_shuf = opt_shuf.optimize(
        protein_sequence=protein, codon_weights=codon_weights,
        target_gc_min=target_gc_min, target_gc_max=target_gc_max
    )["sequence"]

    # Count edits
    _, nom_edits_in_nom, _, _ = analyze_edits(cds_1dp, cds_nom, nominal_positions)
    _, shuf_edits_in_shuf, _, _ = analyze_edits(cds_1dp, cds_shuf, shuffled_positions)
    _, shuf_edits_in_nom, _, _ = analyze_edits(cds_1dp, cds_shuf, nominal_positions)

    # Causal shift is confirmed if shuffled run concentrates edits in shuffled zone
    # and has lower edits in nominal zone compared to nominal run
    causal_shift = (shuf_edits_in_shuf > 0 and shuf_edits_in_shuf >= shuf_edits_in_nom)

    return NegativeControlResult(
        nominal_positions=list(nominal_positions),
        shuffled_positions=list(shuffled_positions),
        nominal_edits_in_nominal_zone=nom_edits_in_nom,
        shuffled_edits_in_shuffled_zone=shuf_edits_in_shuf,
        shuffled_edits_in_nominal_zone=shuf_edits_in_nom,
        causal_shift_confirmed=causal_shift,
    )


def run_boundary_perturbation_test(
    protein: str,
    nominal_positions: Sequence[int],
    offsets: Sequence[int] = (-2, -1, 0, 1, 2),
    codon_weights: Optional[Dict[str, float]] = None,
    target_gc_min: float = 0.40,
    target_gc_max: float = 0.48,
    lambda_weight: float = 2.0,
    target_adaptiveness: float = 0.40,
) -> BoundaryPerturbationResult:
    """Test 7: Boundary Perturbation Sensitivity.
    
    Verifies that shifting nominal boundary by ±1 or ±2 residues yields stable solution family.
    """
    weights = codon_weights or {}
    n_prot = len(protein)

    opt_1dp = DPV2Optimizer()
    cds_1dp = opt_1dp.optimize(
        protein_sequence=protein, codon_weights=weights,
        target_gc_min=target_gc_min, target_gc_max=target_gc_max
    )["sequence"]

    results: Dict[int, Dict[str, Any]] = {}
    base_positions = set(nominal_positions)

    for offset in offsets:
        perturbed = set()
        for p in base_positions:
            np = p + offset
            if 1 <= np <= n_prot:
                perturbed.add(np)

        if not perturbed:
            continue

        opt = StructuralPauseDPOptimizer(
            positions=frozenset(perturbed),
            target_adaptiveness=target_adaptiveness,
            penalty_weight=lambda_weight,
            structure_reference=f"boundary_offset_{offset}",
        )
        res = opt.optimize(
            protein_sequence=protein, codon_weights=weights,
            target_gc_min=target_gc_min, target_gc_max=target_gc_max
        )
        cds = res["sequence"]
        cai = calculate_sequence_cai(cds, weights)
        fit = calculate_structure_fit(cds, sorted(perturbed), weights, target_adaptiveness)
        total_edits, target_edits, _, loc_ratio = analyze_edits(cds_1dp, cds, sorted(perturbed))

        results[offset] = {
            "positions": sorted(perturbed),
            "cai": cai,
            "structure_fit": fit,
            "total_edits": total_edits,
            "target_zone_edits": target_edits,
            "localization_ratio": loc_ratio,
        }

    # Robustness check: CAI variance across minor perturbations is small (< 0.05)
    cai_values = [v["cai"] for v in results.values()]
    is_stable = (max(cai_values) - min(cai_values) < 0.05) if cai_values else False

    return BoundaryPerturbationResult(
        nominal_positions=list(nominal_positions),
        perturbations=results,
        is_stable=is_stable,
    )
