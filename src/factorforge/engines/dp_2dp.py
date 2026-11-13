"""Experimental 2DP Optimizer: CpG/TpA Dinucleotide Penalty.

EPISTEMIC NOTICE:
UNVERIFIED_PUBLIC_CORPUS · Biological validation: NOT_AVAILABLE
Implements an experimental 2DP dinucleotide-junction penalty (CpG/TpA) on top
of 1DP finite-state constraints. This module does not implement or validate an
independent codon-pair or hexamer frequency model.
"""

from __future__ import annotations

import math
from typing import Any

from factorforge.analysis.metrics import detect_homopolymers
from factorforge.constraints.type_iis import get_canonical_forbidden_motifs
from factorforge.engines.dp_v2_1_1 import AA_TO_CODONS, DPV211Optimizer
from factorforge.engines.profile.scoring import compute_5p_mfe_evidence
from factorforge.engines.sllm.automaton import AutomatonCompiler
from factorforge.utils.exceptions import UnsatisfiableDesignError


class DP2DPOptimizer:
    """Exact 2DP optimizer adding cross-codon CpG/TpA & codon-pair context to 1DP."""

    def __init__(
        self,
        forbidden_motifs: list[str] | None = None,
        forbidden_enzymes: set[str] | None = None,
        include_reverse_complement: bool = True,
        ramp_length_codons: int = 15,
        cpg_junction_penalty: float = 1.5,
        cpg_intra_penalty: float = 1.0,
        tpa_junction_penalty: float = 1.0,
        tpa_intra_penalty: float = 0.6,
        initiation_gc_min: float = 0.20,
        initiation_gc_max: float = 0.30,
        nominal_ramp_gc: float = 0.25,
    ) -> None:
        self.base_1dp = DPV211Optimizer(
            forbidden_motifs=forbidden_motifs,
            forbidden_enzymes=forbidden_enzymes,
            include_reverse_complement=include_reverse_complement,
            ramp_length_codons=ramp_length_codons,
            initiation_gc_min=initiation_gc_min,
            initiation_gc_max=initiation_gc_max,
            nominal_ramp_gc=nominal_ramp_gc,
        )
        self.cpg_junction_penalty = float(cpg_junction_penalty)
        self.cpg_intra_penalty = float(cpg_intra_penalty)
        self.tpa_junction_penalty = float(tpa_junction_penalty)
        self.tpa_intra_penalty = float(tpa_intra_penalty)
        self.automaton: AutomatonCompiler = self.base_1dp.automaton
        self.forbidden_motifs = (
            list(forbidden_motifs)
            if forbidden_motifs is not None
            else get_canonical_forbidden_motifs(forbidden_enzymes)
        )

    def _codon_2dp_step_score(
        self,
        i: int,
        codon: str,
        prev_last_base: str | None,
        codon_weights: dict[str, float],
    ) -> float:
        base_score = self.base_1dp.compute_codon_log_fitness(i, codon, codon_weights)
        intra_cpg = codon.count("CG")
        intra_tpa = codon.count("TA")
        junc_cpg = 1 if (prev_last_base == "C" and codon[0] == "G") else 0
        junc_tpa = 1 if (prev_last_base == "T" and codon[0] == "A") else 0
        return (
            base_score
            - self.cpg_intra_penalty * intra_cpg
            - self.tpa_intra_penalty * intra_tpa
            - self.cpg_junction_penalty * junc_cpg
            - self.tpa_junction_penalty * junc_tpa
        )

    def optimize(
        self,
        protein_sequence: str,
        codon_weights: dict[str, float],
        target_gc_min: float = 0.45,
        target_gc_max: float = 0.60,
        upstream_context: str = "",
        stop_codon: str = "TGA",
        downstream_context: str = "",
        left_flank: str | None = None,
        right_flank: str | None = None,
    ) -> dict[str, Any]:
        """Run exact 2DP optimization tracking (gc_count, automaton_state, prev_last_base)."""
        up_ctx = left_flank if left_flank is not None else upstream_context
        down_ctx = right_flank if right_flank is not None else downstream_context
        protein = "".join(protein_sequence.upper().split()).rstrip("*")
        if not protein:
            raise ValueError("protein_sequence must not be empty")
        for aa in protein:
            if aa not in AA_TO_CODONS:
                raise ValueError(f"Unsupported amino acid residue: {aa}")

        s0 = 0
        for ch in up_ctx.upper():
            s0 = self.automaton.step_nucleotide(s0, ch)
            if s0 in self.automaton.terminal_states:
                raise ValueError(f"Upstream context '{up_ctx}' contains a forbidden motif.")

        n_codons = len(protein)
        total_nt = n_codons * 3
        gc_min_frac = target_gc_min if target_gc_min <= 1.0 else target_gc_min / 100.0
        gc_max_frac = target_gc_max if target_gc_max <= 1.0 else target_gc_max / 100.0
        min_gc = math.ceil(gc_min_frac * total_nt)
        max_gc = math.floor(gc_max_frac * total_nt)
        if min_gc > max_gc:
            mid = round(((gc_min_frac + gc_max_frac) / 2.0) * total_nt)
            min_gc = max_gc = mid

        init_last = up_ctx[-1].upper() if up_ctx else None
        # dp[i] maps (gc_count, state, last_base) -> (score, codon, prev_key, rank)
        dp_layers: list[
            dict[tuple[int, int, str | None], tuple[float, str | None, tuple[int, int, str | None] | None, int]]
        ] = [{(0, s0, init_last): (0.0, None, None, 0)}]

        for i, aa in enumerate(protein):
            cur = dp_layers[-1]
            nxt: dict[
                tuple[int, int, str | None],
                tuple[float, str | None, tuple[int, int, str | None] | None, int],
            ] = {}
            codons = sorted(AA_TO_CODONS[aa])
            for prev_key, (prev_score, _, _, prev_rank) in cur.items():
                prev_gc, prev_state, prev_last = prev_key
                for c_idx, codon in enumerate(codons):
                    ns, is_bad = self.automaton.step_codon(prev_state, codon)
                    if is_bad:
                        continue
                    ngc = prev_gc + codon.count("G") + codon.count("C")
                    if ngc > max_gc:
                        continue
                    rem_nt = (n_codons - 1 - i) * 3
                    if ngc + rem_nt < min_gc:
                        continue
                    step_sc = self._codon_2dp_step_score(i, codon, prev_last, codon_weights)
                    tot_sc = prev_score + step_sc
                    nkey = (ngc, ns, codon[-1])
                    cand_rank = prev_rank * 8 + c_idx
                    existing = nxt.get(nkey)
                    if existing is None or tot_sc > existing[0] + 1e-12 or (
                        abs(tot_sc - existing[0]) <= 1e-12 and cand_rank < existing[3]
                    ):
                        nxt[nkey] = (tot_sc, codon, prev_key, cand_rank)
            if not nxt:
                raise UnsatisfiableDesignError(
                    f"2DP constraints unsatisfiable at residue {i + 1} ({aa})."
                )
            dp_layers.append(nxt)

        best_key: tuple[int, int, str | None] | None = None
        best_score = -float("inf")
        stop_u = stop_codon.upper()
        for key, (score, _, _, _) in dp_layers[-1].items():
            gc_cnt, st, _ = key
            if not (min_gc <= gc_cnt <= max_gc):
                continue
            ok = True
            cur_st = st
            for ch in stop_u + down_ctx.upper():
                cur_st = self.automaton.step_nucleotide(cur_st, ch)
                if cur_st in self.automaton.terminal_states:
                    ok = False
                    break
            if ok and score > best_score:
                best_score = score
                best_key = key

        if best_key is None:
            raise UnsatisfiableDesignError("No 2DP candidate satisfies terminal GC and motif bounds.")

        chosen: list[str] = []
        k: tuple[int, int, str | None] | None = best_key
        for idx in range(n_codons, 0, -1):
            assert k is not None
            _, codon, prev_k, _ = dp_layers[idx][k]
            assert codon is not None
            chosen.append(codon)
            k = prev_k
        chosen.reverse()
        cds = "".join(chosen)
        full_seq = cds + stop_u

        cpg_count = sum(1 for j in range(len(cds) - 1) if cds[j : j + 2] == "CG")
        tpa_count = sum(1 for j in range(len(cds) - 1) if cds[j : j + 2] == "TA")
        junc_cpg = sum(1 for idx in range(1, len(chosen)) if chosen[idx - 1][-1] == "C" and chosen[idx][0] == "G")
        junc_tpa = sum(1 for idx in range(1, len(chosen)) if chosen[idx - 1][-1] == "T" and chosen[idx][0] == "A")
        mfe_ev = compute_5p_mfe_evidence(cds, upstream_context=up_ctx)

        return {
            "sequence": full_seq,
            "cds_only": cds,
            "score": round(best_score, 4),
            "gc_content": round((cds.count("G") + cds.count("C")) / len(cds) * 100.0, 2),
            "cpg_count": cpg_count,
            "tpa_count": tpa_count,
            "cpg_junction_count": junc_cpg,
            "tpa_junction_count": junc_tpa,
            "homopolymers": len(detect_homopolymers(cds)),
            "mfe_5p": mfe_ev.get("mfe_kcal_mol") if isinstance(mfe_ev, dict) else None,
            "design_tier": "2DP",
            "epistemic_status": "UNVERIFIED_PUBLIC_CORPUS",
            "biological_validation": "NOT_AVAILABLE",
        }
