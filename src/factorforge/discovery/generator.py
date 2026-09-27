"""Multi-Contract Candidate Pool Generator for FactorForge Discovery Slate (Job 283A / 284B).

Generates a diverse set of candidate CDS designs using deterministic mathematical
DP solvers, rule-based optimization contracts, and research-preview sLLM neural priors.
"""

from __future__ import annotations
from dataclasses import dataclass
import random
from typing import Any, Dict, List, Optional

from factorforge.analysis.metrics import load_codon_usage_table
from factorforge.discovery.filter import HardConstraintConfig
from factorforge.engines.dp_v2 import DPV2Optimizer
from factorforge.engines.dp_v2_1_1 import DPV211Optimizer
from factorforge.engines.profile.optimizer import RuleBasedOptimizer
from factorforge.engines.sllm.adapters import SLMModelAdapter


@dataclass
class RawCandidate:
    """A raw candidate sequence produced by a specific generation contract."""

    candidate_id: str
    generation_contract: str
    strategy_cluster: str
    sequence_dna: str
    rationale: str
    metadata: Dict[str, Any]


class MultiContractCandidateGenerator:
    """Generates candidate sequences across orthogonal biological and mathematical contracts."""

    def __init__(
        self,
        host: str = "nbenthamiana",
        enable_sllm_preview: bool = False,
        sllm_model_adapter: Optional[SLMModelAdapter] = None,
        hard_filter_config: Optional[HardConstraintConfig] = None,
    ) -> None:
        self.host = host
        self.enable_sllm_preview = enable_sllm_preview
        self.sllm_model_adapter = sllm_model_adapter
        self.hard_filter_config = hard_filter_config or HardConstraintConfig()
        table = load_codon_usage_table()
        self.codon_weights = table.codon_weights

    def generate_pool(
        self,
        protein_sequence: str,
        stop_codon: str = "TAA",
        seed: int = 42,
        enable_sllm_override: Optional[bool] = None,
    ) -> List[RawCandidate]:
        """Generates candidates across all deterministic contracts and optional sLLM preview."""
        protein = "".join(protein_sequence.upper().split()).rstrip("*")
        candidates: List[RawCandidate] = []
        seen_sequences: set[str] = set()
        sllm_active = (
            self.enable_sllm_preview if enable_sllm_override is None else enable_sllm_override
        )

        # Helper to add unique candidates
        def _add_cand(
            cand_id: str,
            contract: str,
            cluster: str,
            seq: str,
            rationale: str,
            meta: Optional[Dict[str, Any]] = None,
        ):
            seq_clean = seq.strip().upper()
            if not seq_clean.endswith(stop_codon):
                seq_clean = seq_clean + stop_codon
            if seq_clean in seen_sequences:
                return
            seen_sequences.add(seq_clean)
            candidates.append(
                RawCandidate(
                    candidate_id=cand_id,
                    generation_contract=contract,
                    strategy_cluster=cluster,
                    sequence_dna=seq_clean,
                    rationale=rationale,
                    metadata=meta or {},
                )
            )

        # ------------------------------------------------------------------
        # Contract 1: DP v2.1.1 Initiation-Favored (High 5' open topology)
        # ------------------------------------------------------------------
        try:
            dp_init = DPV211Optimizer(
                beta_mfe_proxy=1.5,
                gamma_gc_penalty=0.8,
                delta_au_bonus=0.8,
                initiation_gc_min=0.20,
                initiation_gc_max=0.30,
            )
            res1 = dp_init.optimize(
                protein_sequence=protein,
                codon_weights=self.codon_weights,
                stop_codon=stop_codon,
            )
            _add_cand(
                cand_id="cand_dp_initiation_favored",
                contract="dp_v2_1_1_initiation_favored",
                cluster="Cluster-Initiation: Open 5' Structure",
                seq=res1["sequence"],
                rationale="Exact 3D DP with maximized 5' open initiation proxy and strict BsaI/SapI/NotI guard.",
                meta=res1.get("diagnostics", {}),
            )
        except Exception:
            pass

        # ------------------------------------------------------------------
        # Contract 2: DP v2.1.1 Translation-Speed Balanced (High Golden CAI)
        # ------------------------------------------------------------------
        try:
            dp_speed = DPV211Optimizer(
                alpha_harmonization=1.5,
                beta_mfe_proxy=0.5,
                gamma_gc_penalty=0.5,
                delta_au_bonus=0.2,
            )
            res2 = dp_speed.optimize(
                protein_sequence=protein,
                codon_weights=self.codon_weights,
                stop_codon=stop_codon,
            )
            _add_cand(
                cand_id="cand_dp_speed_balanced",
                contract="dp_v2_1_1_speed_balanced",
                cluster="Cluster-Speed: High Translation Rate",
                seq=res2["sequence"],
                rationale="Exact 3D DP heavily weighting golden-set CAI and host-adapted codon velocities.",
                meta=res2.get("diagnostics", {}),
            )
        except Exception:
            pass

        # ------------------------------------------------------------------
        # Contract 3: DP v2.1.1 Local Composition Guarded (Narrow Local GC Window)
        # ------------------------------------------------------------------
        try:
            dp_comp = DPV211Optimizer(
                gamma_gc_penalty=1.5,
                homopolymer_max_run=4,
                nominal_ramp_gc=0.24,
            )
            res3 = dp_comp.optimize(
                protein_sequence=protein,
                codon_weights=self.codon_weights,
                stop_codon=stop_codon,
            )
            _add_cand(
                cand_id="cand_dp_composition_guarded",
                contract="dp_v2_1_1_composition_guarded",
                cluster="Cluster-Structure: Stable Transcript",
                seq=res3["sequence"],
                rationale="Exact 3D DP with strict sliding GC smoothing and homopolymer run ceiling <= 4.",
                meta=res3.get("diagnostics", {}),
            )
        except Exception:
            pass

        # ------------------------------------------------------------------
        # Contract 4: DP v2.0 Global Optimal Baseline
        # ------------------------------------------------------------------
        try:
            dp_v2 = DPV2Optimizer()
            res4 = dp_v2.optimize(
                protein_sequence=protein,
                codon_weights=self.codon_weights,
                stop_codon=stop_codon,
            )
            _add_cand(
                cand_id="cand_dp_v2_baseline",
                contract="dp_v2_global_optimal",
                cluster="Cluster-Reference: Global DP Baseline",
                seq=res4["sequence"],
                rationale="Production DP v2.0 unweighted global optimal sequence.",
                meta=res4.get("diagnostics", {}),
            )
        except Exception:
            pass

        # ------------------------------------------------------------------
        # Contract 5: Profile Rule-Based Balanced (Gen 1 Benchmark)
        # ------------------------------------------------------------------
        try:
            prof_opt = RuleBasedOptimizer()
            res5 = prof_opt.optimize(
                protein_sequence=protein,
                strategy="balanced",
                stop_codon=stop_codon,
            )
            _add_cand(
                cand_id="cand_rule_balanced",
                contract="profile_rule_balanced",
                cluster="Cluster-Reference: Rule Profile Baseline",
                seq=res5["sequence"],
                rationale="Rule-based heuristic frequency-matching baseline for historical benchmark.",
                meta=res5.get("diagnostics", {}),
            )
        except Exception:
            pass

        # ------------------------------------------------------------------
        # Contract 6: DP v2.1.1 Synonymous Diversity Trajectory A
        # ------------------------------------------------------------------
        try:
            rng_a = random.Random(seed)
            w_div_a = {c: w * rng_a.uniform(0.7, 1.3) for c, w in self.codon_weights.items()}
            dp_div_a = DPV211Optimizer(
                beta_mfe_proxy=1.0,
                initiation_gc_min=0.20,
                initiation_gc_max=0.30,
            )
            res_div_a = dp_div_a.optimize(
                protein_sequence=protein,
                codon_weights=w_div_a,
                stop_codon=stop_codon,
            )
            _add_cand(
                cand_id="cand_dp_diversity_a",
                contract="dp_v2_1_1_diversity_trajectory_a",
                cluster="Cluster-Exploration: Synonymous Trajectory A",
                seq=res_div_a["sequence"],
                rationale="Exact DP with modulated synonymous codon weights exploring alternative Pareto optima.",
                meta=res_div_a.get("diagnostics", {}),
            )
        except Exception:
            pass

        # ------------------------------------------------------------------
        # Contract 7: DP v2.1.1 Synonymous Diversity Trajectory B
        # ------------------------------------------------------------------
        try:
            rng_b = random.Random(seed + 100)
            w_div_b = {c: w * rng_b.uniform(0.5, 1.5) for c, w in self.codon_weights.items()}
            dp_div_b = DPV211Optimizer(
                beta_mfe_proxy=0.6,
                initiation_gc_min=0.18,
                initiation_gc_max=0.32,
            )
            res_div_b = dp_div_b.optimize(
                protein_sequence=protein,
                codon_weights=w_div_b,
                stop_codon=stop_codon,
            )
            _add_cand(
                cand_id="cand_dp_diversity_b",
                contract="dp_v2_1_1_diversity_trajectory_b",
                cluster="Cluster-Exploration: Synonymous Trajectory B",
                seq=res_div_b["sequence"],
                rationale="Exact DP exploring secondary orthogonal synonymous codon manifold.",
                meta=res_div_b.get("diagnostics", {}),
            )
        except Exception:
            pass

        # ------------------------------------------------------------------
        # Contract 8: DP v2.1.1 Synthesis-Friendly Guarded
        # ------------------------------------------------------------------
        try:
            dp_syn = DPV211Optimizer(
                homopolymer_max_run=4,
                gamma_gc_penalty=1.2,
                initiation_gc_min=0.22,
                initiation_gc_max=0.28,
            )
            res_syn = dp_syn.optimize(
                protein_sequence=protein,
                codon_weights=self.codon_weights,
                stop_codon=stop_codon,
            )
            _add_cand(
                cand_id="cand_dp_synthesis_stable",
                contract="dp_v2_1_1_synthesis_stable",
                cluster="Cluster-Assembly: Synthesis De-Risked",
                seq=res_syn["sequence"],
                rationale="Exact DP with strict homopolymer <= 4 and high GC penalty for commercial synthesis yield.",
                meta=res_syn.get("diagnostics", {}),
            )
        except Exception:
            pass

        # ------------------------------------------------------------------
        # Contract 8b: DP v2.1.1 Synonymous Diversity Trajectory C
        # ------------------------------------------------------------------
        try:
            rng_c = random.Random(seed + 200)
            w_div_c = {c: w * rng_c.uniform(0.6, 1.4) for c, w in self.codon_weights.items()}
            dp_div_c = DPV211Optimizer(
                alpha_harmonization=1.2,
                beta_mfe_proxy=0.8,
                gamma_gc_penalty=0.6,
                initiation_gc_min=0.20,
                initiation_gc_max=0.30,
            )
            res_div_c = dp_div_c.optimize(
                protein_sequence=protein,
                codon_weights=w_div_c,
                stop_codon=stop_codon,
            )
            _add_cand(
                cand_id="cand_dp_diversity_c",
                contract="dp_v2_1_1_diversity_trajectory_c",
                cluster="Cluster-Exploration: Synonymous Trajectory C",
                seq=res_div_c["sequence"],
                rationale="Exact DP exploring tertiary synonymous codon Pareto front.",
                meta=res_div_c.get("diagnostics", {}),
            )
        except Exception:
            pass

        # ------------------------------------------------------------------
        # Contract 8c: DP v2.1.1 Synonymous Diversity Trajectory D
        # ------------------------------------------------------------------
        try:
            rng_d = random.Random(seed + 300)
            w_div_d = {c: w * rng_d.uniform(0.4, 1.6) for c, w in self.codon_weights.items()}
            dp_div_d = DPV211Optimizer(
                alpha_harmonization=0.8,
                beta_mfe_proxy=1.2,
                gamma_gc_penalty=0.8,
                initiation_gc_min=0.18,
                initiation_gc_max=0.32,
            )
            res_div_d = dp_div_d.optimize(
                protein_sequence=protein,
                codon_weights=w_div_d,
                stop_codon=stop_codon,
            )
            _add_cand(
                cand_id="cand_dp_diversity_d",
                contract="dp_v2_1_1_diversity_trajectory_d",
                cluster="Cluster-Exploration: Synonymous Trajectory D",
                seq=res_div_d["sequence"],
                rationale="Exact DP exploring quaternary synonymous codon Pareto front.",
                meta=res_div_d.get("diagnostics", {}),
            )
        except Exception:
            pass

        # ------------------------------------------------------------------
        # Contract 9: sLLM Neural Generative Prior (Beam Best - Approx MAP)
        # ------------------------------------------------------------------
        if sllm_active:
            try:
                from factorforge.engines.sllm.compiler import DesignContractCompiler
                from factorforge.engines.sllm.decoder import ConstrainedBeamDecoder

                resolved_contract = DesignContractCompiler.resolve_contract(
                    host=self.host,
                    config=self.hard_filter_config,
                    expected_stop=stop_codon,
                )
                beam_dec = ConstrainedBeamDecoder(
                    contract=resolved_contract,
                    model_adapter=self.sllm_model_adapter,
                    beam_width=8,
                )
                res_beam = beam_dec.decode(protein_sequence=protein, stop_codon=stop_codon)
                is_rescued = res_beam.get("rescue_triggered", False)
                cand_id = res_beam.get("candidate_id", "cand_sllm_beam_best")
                contract = "sllm_dp_rescue" if is_rescued else "sllm_hybrid_v1_beam_best"
                cluster = (
                    "Cluster-Rescue: DP Fallback Under Contract"
                    if is_rescued
                    else "Cluster-LatentPrior: Neural Genomic Context"
                )
                _add_cand(
                    cand_id=cand_id,
                    contract=contract,
                    cluster=cluster,
                    seq=res_beam["sequence"],
                    rationale=res_beam.get(
                        "rationale",
                        "Constrained sLLM beam search under canonical DesignContract veto.",
                    ),
                    meta=res_beam,
                )
            except Exception:
                pass

        # ------------------------------------------------------------------
        # Contract 10: sLLM Neural Generative Prior (Stochastic Diversity Sample)
        # ------------------------------------------------------------------
        if sllm_active:
            try:
                from factorforge.engines.sllm.compiler import DesignContractCompiler
                from factorforge.engines.sllm.decoder import ConstrainedSampleDecoder

                resolved_contract = DesignContractCompiler.resolve_contract(
                    host=self.host,
                    config=self.hard_filter_config,
                    expected_stop=stop_codon,
                )
                sample_dec = ConstrainedSampleDecoder(
                    contract=resolved_contract,
                    model_adapter=self.sllm_model_adapter,
                    temperature=0.7,
                    seed=seed,
                )
                res_sample = sample_dec.decode(protein_sequence=protein, stop_codon=stop_codon)
                is_rescued_s = res_sample.get("rescue_triggered", False)
                cand_id_s = res_sample.get("candidate_id", "cand_sllm_sample")
                contract_s = "sllm_dp_rescue" if is_rescued_s else "sllm_hybrid_v1_sample"
                cluster_s = (
                    "Cluster-Rescue: DP Fallback Under Contract"
                    if is_rescued_s
                    else "Cluster-LatentPrior: Neural Genomic Context"
                )
                _add_cand(
                    cand_id=cand_id_s,
                    contract=contract_s,
                    cluster=cluster_s,
                    seq=res_sample["sequence"],
                    rationale=res_sample.get(
                        "rationale",
                        "Constrained sLLM temperature sampling under canonical DesignContract veto.",
                    ),
                    meta=res_sample,
                )
            except Exception:
                pass

        return candidates
