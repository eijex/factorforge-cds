"""DP v2.3.0 Optimizer: Multi-Class Structural Pause Optimization and Automated Biological Annotation.

Guarantees & Architectural Contracts:
1. Exact Hard Initiation GC Band Pruning (i = K = 15 codons / 45 nt):
   Enforces ceil(45 * gc_min) <= h_15 <= floor(45 * gc_max) clamped to the achievable
   synonymous GC envelope [GC_min_syn, GC_max_syn] with ZERO additional state dimensions.
2. Multi-Class Position-Conditioned Structural Pause Optimization:
   Evaluates w_struct(i, c) = ell_i(c) - sum_{k} lambda_k * I(i+1 in L_k) * |w(c) - tau_k|.
   When lambda_k = 0 or L_k = empty for all k, guarantees exact bit-identical parity with
   DPV220Optimizer and DPV211Optimizer.
3. Automated Structural Pause Profiling (StructuralPauseAnnotator):
   Identifies domain boundaries: Signal Peptides (SP), Disulfide Loops (SS),
   Transmembrane Helices (TM), N-Glycosylation sites (GLYC), and Flexible Hinge transitions.
4. Homopolymer Aho-Corasick Automaton Rejection:
   Rejects >= 6-mer homopolymers via compiled DFA transitions.
5. Configured Multi-Enzyme Type IIS and Assembly Invariants:
   Full construct elimination across upstream flank, CDS, stop codon, and downstream flank.
6. Deterministic LexRank Sequence Tie-Breaking:
   Guarantees canonical byte-identical, reproducible path selection.
"""

from __future__ import annotations
import math
import re
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE, detect_homopolymers
from factorforge.constraints.type_iis import get_canonical_forbidden_motifs
from factorforge.engines.profile.scoring import compute_5p_mfe_evidence
from factorforge.engines.sllm.automaton import AutomatonCompiler, CompiledAutomaton
from factorforge.utils.exceptions import UnsatisfiableDesignError

AA_TO_CODONS: Dict[str, List[str]] = {}
for codon, aa in STANDARD_GENETIC_CODE.items():
    if aa != "*":
        AA_TO_CODONS.setdefault(aa, []).append(codon)


# Kyte-Doolittle hydropathy scale for transmembrane domain detection
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


class StructuralPauseAnnotator:
    """Automated Structural Pause Profiler for Protein Domain Boundaries."""

    @classmethod
    def detect_signal_peptide(cls, protein: str) -> Optional[int]:
        """Detect standard tripartite signal peptide cleavage site (1-based residue position)."""
        clean_prot = protein.upper()
        if len(clean_prot) < 25:
            return None

        h_chars = set("LIVFWAM")
        h_end = None
        for w_start in range(1, 15):
            w_seq = clean_prot[w_start : w_start + 8]
            if sum(1 for c in w_seq if c in h_chars) >= 5:
                curr = w_start + 8
                while curr < min(24, len(clean_prot)) and clean_prot[curr] in h_chars:
                    curr += 1
                h_end = curr
                break

        if h_end is None or h_end < 7:
            return None

        small_minus1 = {"A": 4.0, "S": 3.0, "G": 2.5, "T": 1.5, "C": 0.5}
        small_minus3 = {"A": 2.5, "V": 2.2, "I": 1.8, "L": 1.8, "S": 1.8, "T": 1.5, "G": 1.8}

        best_cleave = None
        best_score = -999.0
        start_search = max(17, h_end + 3)
        for i in range(start_search, min(35, len(clean_prot))):
            c_m1 = clean_prot[i - 1]
            c_m3 = clean_prot[i - 3]
            if c_m1 in small_minus1 and c_m3 in small_minus3:
                len_penalty = abs(i - 21) * 0.2
                score = small_minus1[c_m1] + small_minus3[c_m3] - len_penalty
                if score > best_score:
                    best_score = score
                    best_cleave = i
        return best_cleave

    @classmethod
    def detect_disulfide_loops(cls, protein: str) -> List[int]:
        """Detect intradomain disulfide loop emergence points (1-based codon position of C2)."""
        clean_prot = protein.upper()
        cys_indices = [idx + 1 for idx, char in enumerate(clean_prot) if char == "C"]
        pause_positions: List[int] = []

        # Find pairs of cysteines separated by at least 8 residues
        for idx, c1 in enumerate(cys_indices):
            for c2 in cys_indices[idx + 1 : idx + 3]:
                if 8 <= (c2 - c1) <= 120:
                    pause_positions.append(c2)
                    break
        return sorted(list(set(pause_positions)))

    @classmethod
    def detect_transmembrane_boundaries(cls, protein: str) -> List[int]:
        """Detect transmembrane helix boundaries via Kyte-Doolittle sliding window."""
        clean_prot = protein.upper()
        win = 19
        if len(clean_prot) < win:
            return []

        scores = []
        for i in range(len(clean_prot) - win + 1):
            segment = clean_prot[i : i + win]
            score = sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in segment) / win
            scores.append(score)

        tm_ranges: List[Tuple[int, int]] = []
        in_tm = False
        start_idx = 0
        for i, sc in enumerate(scores):
            if sc >= 1.6 and not in_tm:
                in_tm = True
                start_idx = i
            elif sc < 1.4 and in_tm:
                in_tm = False
                end_idx = i + win - 1
                if end_idx - start_idx >= 16:
                    tm_ranges.append((start_idx + 1, end_idx + 1))
        if in_tm:
            tm_ranges.append((start_idx + 1, len(clean_prot)))

        boundaries: List[int] = []
        for start_pos, end_pos in tm_ranges:
            boundaries.append(start_pos)
            boundaries.append(end_pos)
        return sorted(list(set(boundaries)))

    @classmethod
    def detect_glycosylation_sites(cls, protein: str) -> List[int]:
        """Detect N-linked glycosylation sequons N-X-[S/T] where X != P (1-based)."""
        clean_prot = protein.upper()
        sites: List[int] = []
        for m in re.finditer(r"N[^P][ST][^P]", clean_prot):
            sites.append(m.start() + 1)
        return sites

    @classmethod
    def detect_hinge_transitions(cls, protein: str) -> List[int]:
        """Detect flexible linkers ((GGGGS)n, (EAAAK)n) transition boundaries (1-based)."""
        clean_prot = protein.upper()
        hinges: List[int] = []
        for m in re.finditer(r"(G{2,4}S){2,}|(EAAAK){2,}", clean_prot):
            hinges.append(m.start() + 1)
            hinges.append(m.end())
        return sorted(list(set(hinges)))

    @classmethod
    def annotate(cls, protein_sequence: str) -> Dict[str, List[int]]:
        """Extract multi-class structural pause profile from amino acid sequence."""
        prot = "".join(protein_sequence.upper().split()).rstrip("*")
        sp = cls.detect_signal_peptide(prot)
        return {
            "SP": [sp] if sp is not None else [],
            "SS": cls.detect_disulfide_loops(prot),
            "TM": cls.detect_transmembrane_boundaries(prot),
            "GLYC": cls.detect_glycosylation_sites(prot),
            "HINGE": cls.detect_hinge_transitions(prot),
        }


class DPV230Optimizer:
    """Exact Multi-Constraint Dynamic Programming CDS Compiler (Engine v2.3.0).

    Unifies 5' initiation ramp harmonization, synonymous envelope clamping,
    automaton reachability relaxation, construct tail validation, LexRank tie-breaking,
    and multi-class position-conditioned structural pause optimization.
    """

    DEFAULT_CLASS_PARAMS: Dict[str, Dict[str, float]] = {
        "SP": {"weight": 2.0, "target_tau": 0.35},
        "SS": {"weight": 2.0, "target_tau": 0.45},
        "TM": {"weight": 3.0, "target_tau": 0.25},
        "GLYC": {"weight": 1.0, "target_tau": 0.45},
        "HINGE": {"weight": 0.5, "target_tau": 0.70},
        "CUSTOM": {"weight": 2.0, "target_tau": 0.45},
    }

    def __init__(
        self,
        forbidden_enzymes: Optional[Set[str]] = None,
        forbidden_motifs: Optional[Iterable[str]] = None,
        include_reverse_complement: bool = True,
        ramp_length_codons: int = 15,
        alpha_harmonization: float = 1.0,
        beta_mfe_proxy: float = 1.0,
        gamma_gc_penalty: float = 0.5,
        delta_au_bonus: float = 0.5,
        r_ramp_start: float = 0.50,
        r_ramp_end: float = 0.90,
        homopolymer_max_run: int = 5,
        initiation_gc_min: float = 0.20,
        initiation_gc_max: float = 0.30,
        nominal_ramp_gc: float = 0.25,
        pause_positions: Optional[Iterable[int]] = None,
        pause_target_adaptiveness: float = 0.45,
        pause_penalty_weight: float = 0.0,
        auto_annotate_pauses: bool = False,
        pause_classes: Optional[Dict[str, Dict[str, Any]]] = None,
        structure_reference: Optional[str] = None,
    ) -> None:
        if forbidden_motifs is not None:
            motifs = list(forbidden_motifs)
        else:
            motifs = get_canonical_forbidden_motifs(forbidden_enzymes)

        self.homopolymer_max_run = homopolymer_max_run
        if self.homopolymer_max_run > 0:
            for base in ("A", "T", "G", "C"):
                poly = base * (self.homopolymer_max_run + 1)
                if poly not in motifs:
                    motifs.append(poly)

        self.automaton: CompiledAutomaton = AutomatonCompiler.compile(
            motifs, include_rc=include_reverse_complement
        )
        self.ramp_length_codons = ramp_length_codons
        self.alpha = alpha_harmonization
        self.beta = beta_mfe_proxy
        self.gamma = gamma_gc_penalty
        self.delta = delta_au_bonus
        self.r_start = r_ramp_start
        self.r_end = r_ramp_end
        self.initiation_gc_min = initiation_gc_min
        self.initiation_gc_max = initiation_gc_max
        self.nominal_ramp_gc = nominal_ramp_gc
        self.auto_annotate_pauses = auto_annotate_pauses
        self.structure_reference = structure_reference or "none"

        # Multi-Class Structural Pause Conditioning
        self.configured_classes: Dict[str, Dict[str, Any]] = {}
        for cls_name, defaults in self.DEFAULT_CLASS_PARAMS.items():
            self.configured_classes[cls_name] = {
                "positions": set(),
                "weight": defaults["weight"],
                "target_tau": defaults["target_tau"],
            }

        # Override with pause_classes if provided
        if pause_classes:
            for c_name, c_cfg in pause_classes.items():
                if c_name not in self.configured_classes:
                    self.configured_classes[c_name] = {
                        "positions": set(),
                        "weight": 2.0,
                        "target_tau": 0.45,
                    }
                if "positions" in c_cfg:
                    self.configured_classes[c_name]["positions"] = set(c_cfg["positions"])
                if "weight" in c_cfg:
                    self.configured_classes[c_name]["weight"] = float(c_cfg["weight"])
                if "target_tau" in c_cfg:
                    self.configured_classes[c_name]["target_tau"] = float(c_cfg["target_tau"])

        # Backward-compatible scalar pause_positions -> mapped to CUSTOM class
        if pause_positions is not None:
            c_pos = set(pause_positions)
            self.configured_classes["CUSTOM"]["positions"].update(c_pos)
            self.configured_classes["CUSTOM"]["weight"] = pause_penalty_weight
            self.configured_classes["CUSTOM"]["target_tau"] = pause_target_adaptiveness

        if self.ramp_length_codons < 1:
            raise ValueError("ramp_length_codons must be at least 1")
        if self.homopolymer_max_run < 1:
            raise ValueError("homopolymer_max_run must be at least 1")
        if not 0.0 <= self.initiation_gc_min <= self.initiation_gc_max <= 1.0:
            raise ValueError("initiation GC bounds must satisfy 0 <= min <= max <= 1")

    def _resolve_active_pause_map(
        self, protein: str
    ) -> Dict[int, List[Tuple[str, float, float]]]:
        """Resolves active pause map: position -> list of (class_name, weight, target_tau)."""
        active_map: Dict[int, List[Tuple[str, float, float]]] = {}

        # 1. Base configured classes
        for c_name, c_data in self.configured_classes.items():
            w = c_data["weight"]
            tau = c_data["target_tau"]
            if w > 0.0:
                for pos in c_data["positions"]:
                    active_map.setdefault(pos, []).append((c_name, w, tau))

        # 2. Automated structural pause annotation if enabled
        if self.auto_annotate_pauses:
            detected = StructuralPauseAnnotator.annotate(protein)
            for c_name, positions in detected.items():
                defaults = self.configured_classes.get(
                    c_name, self.DEFAULT_CLASS_PARAMS.get(c_name, {"weight": 2.0, "target_tau": 0.45})
                )
                w = defaults["weight"]
                tau = defaults["target_tau"]
                if w > 0.0:
                    for pos in positions:
                        # Avoid redundant duplicate class registrations
                        existing = [x[0] for x in active_map.get(pos, [])]
                        if c_name not in existing:
                            active_map.setdefault(pos, []).append((c_name, w, tau))

        return active_map

    def compute_codon_log_fitness(
        self,
        position_idx: int,
        codon: str,
        codon_weights: Dict[str, float],
        active_pause_map: Dict[int, List[Tuple[str, float, float]]],
    ) -> float:
        """Computes position-dependent structural log-fitness score w_struct(i, c)."""
        w = codon_weights.get(codon, 1e-4)

        # 1. Base piecewise fitness ell_i(c)
        if position_idx >= self.ramp_length_codons:
            # Stage 2: Elongation Body -> Host Log-CAI
            base_fitness = math.log(max(w, 1e-6))
        else:
            # Stage 1: Initiation Ramp (0 <= position_idx < K)
            denom = max(1, self.ramp_length_codons - 1)
            r_i = self.r_start + (self.r_end - self.r_start) * (position_idx / float(denom))
            h_score = -abs(w - r_i)
            g_count = codon.count("G") + codon.count("C")
            au_count = 3 - g_count
            pen_open = -self.gamma * g_count + self.delta * au_count
            base_fitness = self.alpha * h_score + self.beta * pen_open

        # 2. Multi-class structural pause penalty
        one_based_pos = position_idx + 1
        if one_based_pos in active_pause_map:
            pause_penalty = sum(
                w_pen * abs(w - target_tau)
                for _, w_pen, target_tau in active_pause_map[one_based_pos]
            )
            return base_fitness - pause_penalty

        return base_fitness

    def optimize(
        self,
        protein_sequence: str,
        codon_weights: Dict[str, float],
        target_gc_min: float = 0.40,
        target_gc_max: float = 0.47,
        upstream_context: str = "",
        stop_codon: str = "TAA",
        downstream_context: str = "",
        allow_nearest_infeasible: bool = False,
        left_flank: Optional[str] = None,
        right_flank: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Runs exact position-dependent initiation-aware 3D DP with multi-class pause guards."""
        up_ctx = left_flank if left_flank is not None else upstream_context
        down_ctx = right_flank if right_flank is not None else downstream_context

        protein = "".join(protein_sequence.upper().split()).rstrip("*")
        if not protein:
            raise ValueError("protein_sequence must not be empty")

        for aa in protein:
            if aa not in AA_TO_CODONS:
                raise ValueError(f"Unsupported amino acid residue: {aa}")

        # Initial automaton state from upstream_context
        s0 = 0
        for char in up_ctx.upper():
            s0 = self.automaton.step_nucleotide(s0, char)
            if s0 in self.automaton.terminal_states:
                raise ValueError(f"Upstream context '{up_ctx}' contains a forbidden motif.")

        n_codons = len(protein)
        total_nt = n_codons * 3

        # Resolve active pause positions and penalties
        active_pause_map = self._resolve_active_pause_map(protein)

        # Target GC bounds
        gc_min_frac = target_gc_min if target_gc_min <= 1.0 else target_gc_min / 100.0
        gc_max_frac = target_gc_max if target_gc_max <= 1.0 else target_gc_max / 100.0

        min_gc_count = int(math.ceil(gc_min_frac * total_nt))
        max_gc_count = int(math.floor(gc_max_frac * total_nt))
        if min_gc_count > max_gc_count:
            mid_count = int(round(((gc_min_frac + gc_max_frac) / 2.0) * total_nt))
            min_gc_count = mid_count
            max_gc_count = mid_count

        # 5' Initiation Ramp GC Constraints at codon i = K = ramp_length_codons
        ramp_k = min(self.ramp_length_codons, n_codons)
        ramp_nt = ramp_k * 3
        apply_initiation_gc_guard = n_codons >= self.ramp_length_codons
        nom_init_min = int(math.ceil(self.initiation_gc_min * ramp_nt))
        nom_init_max = int(math.floor(self.initiation_gc_max * ramp_nt))
        if nom_init_min > nom_init_max:
            mid_init = int(round(self.nominal_ramp_gc * ramp_nt))
            nom_init_min = mid_init
            nom_init_max = mid_init

        # Reachable synonymous GC range for the initiation ramp
        ramp_min_possible_gc = sum(
            min(c.count("G") + c.count("C") for c in AA_TO_CODONS[aa]) for aa in protein[:ramp_k]
        )
        ramp_max_possible_gc = sum(
            max(c.count("G") + c.count("C") for c in AA_TO_CODONS[aa]) for aa in protein[:ramp_k]
        )

        init_min_gc_count = max(nom_init_min, ramp_min_possible_gc)
        init_max_gc_count = min(nom_init_max, ramp_max_possible_gc)
        initiation_gc_status = (
            "REQUESTED_BAND" if apply_initiation_gc_guard else "NOT_APPLIED_SHORT_SEQUENCE"
        )
        if apply_initiation_gc_guard and (
            init_min_gc_count != nom_init_min or init_max_gc_count != nom_init_max
        ):
            initiation_gc_status = "SYNONYMOUS_ENVELOPE_CLAMPED"
        if init_min_gc_count > init_max_gc_count:
            if apply_initiation_gc_guard:
                initiation_gc_status = "SYNONYMOUS_ENVELOPE_CLAMPED"
            if ramp_min_possible_gc > nom_init_max:
                init_min_gc_count = ramp_min_possible_gc
                init_max_gc_count = min(ramp_min_possible_gc + 2, ramp_max_possible_gc)
            else:
                init_max_gc_count = ramp_max_possible_gc
                init_min_gc_count = max(ramp_max_possible_gc - 2, ramp_min_possible_gc)

        # dp_layers[i] maps (gc_count, automaton_state) ->
        #   (cumulative_score, codon, prev_gc, prev_state, path_rank)
        dp_layers: List[
            Dict[Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int], int]]
        ] = []
        dp_layers.append({(0, s0): (0.0, None, None, None, 0)})

        for i, aa in enumerate(protein):
            current_layer = dp_layers[-1]
            next_layer: Dict[
                Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int], int]
            ] = {}
            synonymous_codons = AA_TO_CODONS[aa]
            is_ramp_boundary = i == ramp_k - 1

            for (g, s), (score, _, _, _, prev_path_rank) in current_layer.items():
                for codon in synonymous_codons:
                    # 1. Automaton Veto (Hard Invariant: Type IIS + Homopolymers)
                    next_s, is_forbidden = self.automaton.step_codon(s, codon)
                    if is_forbidden:
                        continue

                    # 2. Cumulative GC count
                    codon_gc = codon.count("G") + codon.count("C")
                    next_g = g + codon_gc

                    # 3. Position-dependent multi-class structural fitness score
                    step_log_fitness = self.compute_codon_log_fitness(
                        i, codon, codon_weights, active_pause_map
                    )
                    next_score = score + step_log_fitness

                    state_key = (next_g, next_s)
                    if state_key not in next_layer:
                        next_layer[state_key] = (next_score, codon, g, s, prev_path_rank)
                    else:
                        (
                            existing_score,
                            existing_codon,
                            existing_prev_g,
                            existing_prev_s,
                            existing_prev_rank,
                        ) = next_layer[state_key]
                        if next_score > existing_score + 1e-12:
                            next_layer[state_key] = (next_score, codon, g, s, prev_path_rank)
                        elif abs(next_score - existing_score) <= 1e-12:
                            cand_key = (prev_path_rank, codon)
                            exist_key = (existing_prev_rank, existing_codon)
                            if cand_key < exist_key:
                                next_layer[state_key] = (next_score, codon, g, s, prev_path_rank)

            if not next_layer:
                raise UnsatisfiableDesignError(
                    f"DP v2.3.0: State space exhausted at position {i + 1} ({aa}). "
                    f"No synonymous codon path satisfies the full construct automaton constraints."
                )

            # Exact 5' Initiation GC Band Hard Pruning at i = ramp_k - 1
            if is_ramp_boundary and apply_initiation_gc_guard:
                filtered_layer = {
                    k_st: v
                    for k_st, v in next_layer.items()
                    if init_min_gc_count <= k_st[0] <= init_max_gc_count
                }
                if filtered_layer:
                    next_layer = filtered_layer
                else:
                    initiation_gc_status = "AUTOMATON_REACHABILITY_RELAXED"
                    min_dist_to_band = min(
                        0
                        if init_min_gc_count <= k_st[0] <= init_max_gc_count
                        else min(abs(k_st[0] - init_min_gc_count), abs(k_st[0] - init_max_gc_count))
                        for k_st in next_layer.keys()
                    )
                    next_layer = {
                        k_st: v
                        for k_st, v in next_layer.items()
                        if (
                            0
                            if init_min_gc_count <= k_st[0] <= init_max_gc_count
                            else min(
                                abs(k_st[0] - init_min_gc_count), abs(k_st[0] - init_max_gc_count)
                            )
                        )
                        == min_dist_to_band
                    }
                    reachable_counts = [state[0] for state in next_layer]
                    init_min_gc_count = min(reachable_counts)
                    init_max_gc_count = max(reachable_counts)

            # Assign canonical path_rank to all active states
            sorted_states = sorted(
                next_layer.keys(), key=lambda k_st: (next_layer[k_st][4], next_layer[k_st][1])
            )
            rank_map = {st: idx for idx, st in enumerate(sorted_states)}
            compacted_layer: Dict[
                Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int], int]
            ] = {}
            for st, (sc, cd, pg, ps, pr) in next_layer.items():
                compacted_layer[st] = (sc, cd, pg, ps, rank_map[st])

            dp_layers.append(compacted_layer)

        final_layer = dp_layers[-1]

        # Construct full tail sequence: STOP + downstream_context
        s_tail = (stop_codon + down_ctx).upper()

        # Filter states compatible with full construct tail
        valid_final_states: List[
            Tuple[Tuple[int, int], Tuple[float, Optional[str], Optional[int], Optional[int], int]]
        ] = []
        for (g, s), val in final_layer.items():
            curr = s
            tail_valid = True
            for char in s_tail:
                curr = self.automaton.step_nucleotide(curr, char)
                if curr in self.automaton.terminal_states:
                    tail_valid = False
                    break
            if tail_valid:
                valid_final_states.append(((g, s), val))

        if not valid_final_states:
            raise UnsatisfiableDesignError(
                f"DP v2.3.0: No candidate satisfies full construct tail constraints "
                f"(STOP='{stop_codon}', Downstream='{down_ctx}')."
            )

        # Candidate Selection: Multi-tier tie breaking
        target_mid_gc_count = (min_gc_count + max_gc_count) / 2.0
        in_band_candidates = [
            item for item in valid_final_states if min_gc_count <= item[0][0] <= max_gc_count
        ]

        if in_band_candidates:
            best_item = max(
                in_band_candidates,
                key=lambda item: (
                    item[1][0],  # Primary: Total piecewise DP score J(X)
                    -abs(item[0][0] - target_mid_gc_count),  # Tie 1: Distance to target GC midpoint
                    -item[1][4],  # Tie 2: Full-path canonical LexRank
                ),
            )
            gc_feasible = True
        else:
            if not allow_nearest_infeasible:
                raise UnsatisfiableDesignError(
                    f"DP v2.3.0: No synonymous candidate satisfies the requested GC band "
                    f"[{gc_min_frac:.1%}, {gc_max_frac:.1%}] (requires [{min_gc_count}, {max_gc_count}] GC bases of {total_nt} nt) "
                    f"under full construct Type IIS and multi-class pause constraints."
                )
            best_item = min(
                valid_final_states,
                key=lambda item: (
                    min(abs(item[0][0] - min_gc_count), abs(item[0][0] - max_gc_count)),
                    -item[1][0],
                    abs(item[0][0] - target_mid_gc_count),
                    item[1][4],
                ),
            )
            gc_feasible = False

        # Backtracking
        best_state, (best_score, _, _, _, _) = best_item
        chosen_codons: List[str] = []
        curr_g, curr_s = best_state

        for layer_idx in range(n_codons, 0, -1):
            entry = dp_layers[layer_idx][(curr_g, curr_s)]
            _, codon, prev_g, prev_s, _ = entry
            assert codon is not None
            assert prev_g is not None
            assert prev_s is not None
            chosen_codons.append(codon)
            curr_g, curr_s = prev_g, prev_s

        chosen_codons.reverse()
        optimized_dna = "".join(chosen_codons)

        # Compute segmental and local metrics
        ramp_codons = chosen_codons[:ramp_k]
        body_codons = chosen_codons[ramp_k:]

        ramp_dna = "".join(ramp_codons)
        body_dna = "".join(body_codons)

        ramp_cai = (
            math.exp(
                sum(math.log(max(codon_weights.get(c, 1e-4), 1e-6)) for c in ramp_codons) / ramp_k
            )
            if ramp_k > 0
            else 1.0
        )
        body_cai = (
            math.exp(
                sum(math.log(max(codon_weights.get(c, 1e-4), 1e-6)) for c in body_codons)
                / len(body_codons)
            )
            if body_codons
            else 1.0
        )
        global_cai = math.exp(
            sum(math.log(max(codon_weights.get(c, 1e-4), 1e-6)) for c in chosen_codons) / n_codons
        )

        global_gc = (optimized_dna.count("G") + optimized_dna.count("C")) / float(total_nt)
        ramp_gc = (
            (ramp_dna.count("G") + ramp_dna.count("C")) / float(ramp_nt) if ramp_nt > 0 else 0.0
        )
        body_gc = (
            (body_dna.count("G") + body_dna.count("C")) / float(len(body_dna))
            if len(body_dna) > 0
            else 0.0
        )

        # 50-bp sliding GC window audit
        window_size = 50
        step = 10
        gc_windows = []
        if total_nt >= window_size:
            for w_start in range(0, total_nt - window_size + 1, step):
                w_seq = optimized_dna[w_start : w_start + window_size]
                w_gc = (w_seq.count("G") + w_seq.count("C")) / float(window_size)
                gc_windows.append(w_gc)
        else:
            gc_windows.append(global_gc)

        # Compute 5' MFE proxy evidence
        mfe_evidence = compute_5p_mfe_evidence(
            cds_sequence=optimized_dna,
            upstream_context=up_ctx,
        )

        # Homopolymer runs
        homopolymers = detect_homopolymers(optimized_dna, max_run=self.homopolymer_max_run + 1)

        # Detailed Pause Reporting
        pause_report: Dict[int, Any] = {}
        for pos in sorted(active_pause_map.keys()):
            idx = pos - 1
            if idx < len(chosen_codons):
                cd = chosen_codons[idx]
                w_cd = codon_weights.get(cd, 1e-4)
                classes_info = []
                for c_name, w_pen, target_tau in active_pause_map[pos]:
                    classes_info.append({
                        "class": c_name,
                        "weight": w_pen,
                        "target_tau": target_tau,
                        "delta": abs(w_cd - target_tau),
                    })
                pause_report[pos] = {
                    "codon": cd,
                    "adaptiveness": round(w_cd, 4),
                    "classes": classes_info,
                }

        class_summary: Dict[str, Any] = {}
        for c_name, c_data in self.configured_classes.items():
            class_summary[c_name] = {
                "configured_positions": sorted(list(c_data["positions"])),
                "weight": c_data["weight"],
                "target_tau": c_data["target_tau"],
            }

        return {
            "sequence": optimized_dna,
            "score": best_score,
            "protein_length": n_codons,
            "total_nt": total_nt,
            "gc_percent": round(global_gc * 100.0, 2),
            "gc_feasible": gc_feasible,
            "gc_min_target": gc_min_frac,
            "gc_max_target": gc_max_frac,
            "initiation_gc_status": initiation_gc_status,
            "initiation_gc_requested_min_count": nom_init_min,
            "initiation_gc_requested_max_count": nom_init_max,
            "initiation_gc_active_min_count": init_min_gc_count,
            "initiation_gc_active_max_count": init_max_gc_count,
            "initiation_gc_synonymous_min_count": ramp_min_possible_gc,
            "initiation_gc_synonymous_max_count": ramp_max_possible_gc,
            "gc_5p_45nt_percent": round(ramp_gc * 100.0, 2),
            "gc_body_percent": round(body_gc * 100.0, 2),
            "sliding_window_gc": {
                "window_size_bp": window_size,
                "step_bp": step,
                "min_window_gc_percent": round(min(gc_windows) * 100.0, 2),
                "max_window_gc_percent": round(max(gc_windows) * 100.0, 2),
                "mean_window_gc_percent": round((sum(gc_windows) / len(gc_windows)) * 100.0, 2),
            },
            "cai_global": round(global_cai, 4),
            "cai_ramp_45nt": round(ramp_cai, 4),
            "cai_body": round(body_cai, 4),
            "mfe_5p_evidence": mfe_evidence,
            "max_homopolymer_run": (
                max([hp["length"] for hp in homopolymers], default=0)
                if homopolymers
                else self.homopolymer_max_run
            ),
            "homopolymers_detected": homopolymers,
            "pause_conditioning": {
                "auto_annotated": self.auto_annotate_pauses,
                "classes": class_summary,
                "total_active_positions": len(active_pause_map),
                "details": pause_report,
            },
            "engine_version": "2.3.0",
        }


class DPV230EngineAdapter:
    """FactorForge Profile Engine adapter for DPV230Optimizer."""

    def __init__(self, optimizer: Optional[DPV230Optimizer] = None) -> None:
        self.optimizer = optimizer or DPV230Optimizer()

    def optimize_cds(
        self,
        protein_sequence: str,
        codon_weights: Dict[str, float],
        target_gc_min: float = 0.40,
        target_gc_max: float = 0.47,
        upstream_context: str = "",
        stop_codon: str = "TAA",
        downstream_context: str = "",
        allow_nearest_infeasible: bool = False,
    ) -> Dict[str, Any]:
        return self.optimizer.optimize(
            protein_sequence=protein_sequence,
            codon_weights=codon_weights,
            target_gc_min=target_gc_min,
            target_gc_max=target_gc_max,
            upstream_context=upstream_context,
            stop_codon=stop_codon,
            downstream_context=downstream_context,
            allow_nearest_infeasible=allow_nearest_infeasible,
        )
