# factorforge/src/factorforge/engines/sllm/decoder.py
"""Constrained Beam and Sampling Decoders with Exact DP Rescue (Job 284A / 284B).

Ensures that sLLM candidate generation:
1. Operates under unified ResolvedDesignContract (Synonymous mask + Aho-Corasick + Homopolymers).
2. Uses mathematically normalized sequence log-probability (approximate constrained MAP for beam search).
3. Evaluates full construct / stop codon junctions with independent HardConstraintFilter.
4. If dead-ends occur, executes Exact DP Rescue under the exact same DesignContract with allow_nearest_infeasible=False.
5. Accurately tracks requested vs actual generator provenance and fallback chains.
"""

from __future__ import annotations
import math
from typing import Any, Dict, List, Optional, Tuple
import numpy as np

from factorforge.discovery.filter import HardConstraintConfig
from factorforge.engines.dp_v2_1_1 import DPV211Optimizer
from factorforge.engines.sllm.adapters import DeterministicMockSLMAdapter, SLMModelAdapter
from factorforge.engines.sllm.compiler import DesignContractCompiler, ResolvedDesignContract
from factorforge.engines.sllm.interfaces import (
    CODON_VOCAB,
    ConstraintProcessorPipeline,
    ConstraintState,
)


def _log_softmax_masked(logits: np.ndarray) -> np.ndarray:
    """Computes log-softmax probabilities over valid finite logits: log P(c_t | c_<t, valid)."""
    valid_mask = logits > -math.inf
    if not np.any(valid_mask):
        return np.full_like(logits, -math.inf)

    finite_logits = logits[valid_mask]
    max_logit = np.max(finite_logits)
    log_sum_exp = max_logit + np.log(np.sum(np.exp(finite_logits - max_logit)))

    log_probs = np.full_like(logits, -math.inf)
    log_probs[valid_mask] = finite_logits - log_sum_exp
    return log_probs


def _softmax_with_temperature(logits: np.ndarray, temperature: float = 1.0) -> np.ndarray:
    """Computes softmax probabilities over valid finite logits with temperature scaling."""
    valid_mask = logits > -math.inf
    if not np.any(valid_mask):
        return np.zeros_like(logits)

    t = max(1e-3, float(temperature))
    scaled = logits[valid_mask] / t
    shifted = scaled - np.max(scaled)
    exp_scores = np.exp(shifted)
    probs = exp_scores / np.sum(exp_scores)

    full_probs = np.zeros_like(logits)
    full_probs[valid_mask] = probs
    return full_probs


class ConstrainedBeamDecoder:
    """Approximate constrained MAP beam search decoder under canonical ResolvedDesignContract."""

    def __init__(
        self,
        contract: Optional[ResolvedDesignContract] = None,
        model_adapter: Optional[SLMModelAdapter] = None,
        beam_width: int = 8,
        host: str = "nbenthamiana",
        hard_config: Optional[HardConstraintConfig] = None,
        pipeline: Optional[ConstraintProcessorPipeline] = None,
    ) -> None:
        if contract is not None:
            self.contract = contract
        else:
            self.contract = DesignContractCompiler.resolve_contract(
                host=host,
                config=hard_config,
            )

        self.model_adapter = model_adapter or DeterministicMockSLMAdapter(host=self.contract.host)
        self.beam_width = max(1, beam_width)

    def _execute_dp_rescue(
        self,
        protein_sequence: str,
        stop_codon: str,
        requested_generator: str = "sllm_hybrid_beam",
        fallback_chain: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """Executes exact DP solver under the exact same DesignContract when sLLM dead-ends."""
        chain = list(fallback_chain or [])
        chain.append("dp_v2_1_1_rescue")

        try:
            dp_solver = DPV211Optimizer(
                homopolymer_max_run=self.contract.hard_config.homopolymer_max_run,
                forbidden_motifs=list(self.contract.final_validator.forbidden_motifs),
                initiation_gc_min=self.contract.hard_config.initiation_45nt_gc_min / 100.0,
                initiation_gc_max=self.contract.hard_config.initiation_45nt_gc_max / 100.0,
            )
            # Calculate total construct CDS length including stop codon
            stop_gc = (stop_codon.count("G") + stop_codon.count("C")) if stop_codon else 0
            stop_len = len(stop_codon) if stop_codon else 0
            coding_len = len(protein_sequence) * 3
            total_cds_len = coding_len + stop_len

            min_gc_total = int(
                math.ceil(total_cds_len * (self.contract.hard_config.global_gc_min / 100.0))
            )
            max_gc_total = int(
                math.floor(total_cds_len * (self.contract.hard_config.global_gc_max / 100.0))
            )

            gc_min_coding = max(0, min_gc_total - stop_gc)
            gc_max_coding = min(coding_len, max_gc_total - stop_gc)

            target_gc_min = gc_min_coding / float(coding_len)
            target_gc_max = gc_max_coding / float(coding_len)

            dp_result = dp_solver.optimize(
                protein_sequence=protein_sequence,
                codon_weights=self.contract.codon_weights,
                stop_codon=stop_codon,
                target_gc_min=target_gc_min,
                target_gc_max=target_gc_max,
                allow_nearest_infeasible=False,  # Strict: NEVER allow nearest infeasible in hard rescue
            )
            rescued_seq = dp_result["sequence"].strip().upper()
            if stop_codon and not rescued_seq.endswith(stop_codon):
                rescued_seq = rescued_seq + stop_codon

            # Verify rescued candidate with independent validator
            val_res = self.contract.final_validator.evaluate(
                sequence_dna=rescued_seq,
                expected_protein_aa=protein_sequence,
                expected_stop=stop_codon,
                upstream_flank=self.contract.upstream_flank,
                downstream_flank=self.contract.downstream_flank,
            )
            if not val_res.is_feasible:
                raise RuntimeError(
                    f"UNSAT: Rescued sequence violated DesignContract: {', '.join(val_res.violations)}"
                )

            return {
                "sequence": rescued_seq,
                "candidate_id": "cand_sllm_dp_rescue",
                "requested_generator": requested_generator,
                "actual_generator": "dp_v2_1_1",
                "status": "DP_RESCUED",
                "score": None,
                "score_type": "dp_objective",
                "dp_objective_score": round(dp_result.get("score", 0.0), 4),
                "rescue_triggered": True,
                "fallback_chain": chain,
                "rationale": "Exact DP rescue executed following sLLM constraint dead-end under exact DesignContract.",
                "model_name": self.model_adapter.model_name,
                "beam_width": self.beam_width,
                "contract_digest": self.contract.contract_digest,
            }

        except Exception as e:
            raise RuntimeError(
                f"UNSAT: Both sLLM beam search and exact DP rescue failed under DesignContract. Cause: {e}"
            ) from e

    def _execute_partial_dp_rescue(
        self,
        protein_sequence: str,
        stop_codon: str,
        deadend_pos: int,
        prefix_codons: List[str],
        requested_generator: str = "sllm_hybrid_beam",
        fallback_chain: Optional[List[str]] = None,
        rewind_ladder: Optional[List[int]] = None,
        lookback_codons: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Attempts Prefix-Conditioned Adaptive Partial DP Rescue with progressive lookback ladder."""
        chain = list(fallback_chain or [])
        chain.append(f"adaptive_partial_dp_attempt_at_{deadend_pos}")

        if rewind_ladder is not None:
            ladder = rewind_ladder
        elif lookback_codons is not None:
            ladder = [lookback_codons]
        else:
            ladder = [5, 10, 20, 50]
        stop_gc = (stop_codon.count("G") + stop_codon.count("C")) if stop_codon else 0
        stop_len = len(stop_codon) if stop_codon else 0
        coding_len = len(protein_sequence) * 3
        total_cds_len = coding_len + stop_len

        global_gc_min_pct = self.contract.hard_config.global_gc_min / 100.0
        global_gc_max_pct = self.contract.hard_config.global_gc_max / 100.0
        min_gc_total = int(math.ceil(total_cds_len * global_gc_min_pct))
        max_gc_total = int(math.floor(total_cds_len * global_gc_max_pct))

        for lookback_k in ladder:
            cut_idx = max(0, deadend_pos - lookback_k)
            if cut_idx < 3:
                # Dead-end happened near 5' initiation -> jump to global DP rescue
                break

            safe_prefix = prefix_codons[:cut_idx]
            prefix_dna = "".join(safe_prefix)
            suffix_protein = protein_sequence[cut_idx:]
            L_suffix = len(suffix_protein) * 3
            if L_suffix == 0:
                continue

            gc_prefix = prefix_dna.count("G") + prefix_dna.count("C")
            gc_min_rem = max(0, min_gc_total - gc_prefix - stop_gc)
            gc_max_rem = min(L_suffix, max_gc_total - gc_prefix - stop_gc)

            if gc_min_rem > gc_max_rem:
                # Infeasible remaining GC budget for this cut index -> try larger rewind
                chain.append(f"k_{lookback_k}_gc_budget_infeasible")
                continue

            target_gc_min = gc_min_rem / float(L_suffix)
            target_gc_max = gc_max_rem / float(L_suffix)

            try:
                dp_solver = DPV211Optimizer(
                    homopolymer_max_run=self.contract.hard_config.homopolymer_max_run,
                    forbidden_motifs=list(self.contract.final_validator.forbidden_motifs),
                    ramp_length_codons=1,
                )
                full_up_ctx = (self.contract.upstream_flank or "") + prefix_dna

                dp_result = dp_solver.optimize(
                    protein_sequence=suffix_protein,
                    codon_weights=self.contract.codon_weights,
                    stop_codon=stop_codon,
                    target_gc_min=target_gc_min,
                    target_gc_max=target_gc_max,
                    upstream_context=full_up_ctx,
                    downstream_context=self.contract.downstream_flank or "",
                    allow_nearest_infeasible=False,
                )

                suffix_seq = dp_result["sequence"].strip().upper()
                stitched_seq = prefix_dna + suffix_seq
                if stop_codon and not stitched_seq.endswith(stop_codon):
                    stitched_seq = stitched_seq + stop_codon

                val_res = self.contract.final_validator.evaluate(
                    sequence_dna=stitched_seq,
                    expected_protein_aa=protein_sequence,
                    expected_stop=stop_codon,
                    upstream_flank=self.contract.upstream_flank,
                    downstream_flank=self.contract.downstream_flank,
                )

                if not val_res.is_feasible:
                    chain.append(f"k_{lookback_k}_validation_failed")
                    continue

                # Adaptive rescue succeeded
                chain.append(f"partial_dp_success_k_{lookback_k}")
                ppr = round(cut_idx / float(len(protein_sequence)), 4)
                return {
                    "sequence": stitched_seq,
                    "candidate_id": "cand_sllm_partial_dp_rescue",
                    "requested_generator": requested_generator,
                    "actual_generator": "sllm_hybrid_partial_dp_rescue",
                    "status": "PARTIAL_DP_RESCUED",
                    "score": None,
                    "score_type": "dp_objective",
                    "dp_objective_score": round(dp_result.get("score", 0.0), 4),
                    "rescue_triggered": True,
                    "partial_rescue": True,
                    "rewind_k": lookback_k,
                    "preserved_prefix_codons": cut_idx,
                    "prefix_preservation_rate": ppr,
                    "stitch_junction_pos": cut_idx,
                    "fallback_chain": chain,
                    "rationale": f"Prefix-conditioned adaptive DP rescue succeeded at K={lookback_k}, preserving {cut_idx}/{len(protein_sequence)} ({ppr * 100:.1f}%) codons.",
                    "model_name": self.model_adapter.model_name,
                    "beam_width": getattr(self, "beam_width", 1),
                    "contract_digest": self.contract.contract_digest,
                }

            except Exception as e:
                chain.append(f"k_{lookback_k}_exception_{type(e).__name__}")
                continue

        # All adaptive lookback attempts exhausted -> Fallback to global DP rescue
        chain.append("all_adaptive_rewinds_exhausted")
        return self._execute_dp_rescue(
            protein_sequence=protein_sequence,
            stop_codon=stop_codon,
            requested_generator=requested_generator,
            fallback_chain=chain,
        )

    def decode(
        self,
        protein_sequence: str,
        stop_codon: Optional[str] = None,
        requested_generator: str = "sllm_hybrid_beam",
        fallback_chain: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """Generates the highest-scoring constrained candidate (cand_sllm_beam_best)."""
        clean_aa = "".join(protein_sequence.upper().split()).rstrip("*")
        if not clean_aa:
            raise ValueError("Cannot decode empty protein sequence.")

        effective_stop = stop_codon or self.contract.expected_stop

        # Beam state tuple: (cumulative_log_prob, list_of_codons, ConstraintState)
        beams: List[Tuple[float, List[str], ConstraintState]] = [
            (
                0.0,
                [],
                ConstraintState(
                    automaton_node=0,
                    position=0,
                    trailing_nt=self.contract.upstream_flank[-6:]
                    if self.contract.upstream_flank
                    else "",
                ),
            )
        ]

        for pos, aa in enumerate(clean_aa):
            new_candidates: List[Tuple[float, List[str], ConstraintState]] = []

            for cum_log_prob, prefix_codons, state in beams:
                # 1. Query raw unconstrained logits from model adapter
                raw_logits = self.model_adapter.predict_logits(prefix_codons, aa, pos)

                # 2. Apply neuro-symbolic constraint processor pipeline
                masked_logits = self.contract.decode_pipeline(pos, aa, state, raw_logits)

                # 3. Compute normalized step log-probabilities
                log_probs = _log_softmax_masked(masked_logits.scores)

                # 4. Collect all feasible transitions
                valid_indices = np.where(log_probs > -math.inf)[0]
                for idx in valid_indices:
                    codon = CODON_VOCAB[idx]
                    step_log_prob = float(log_probs[idx])
                    next_state = self.contract.decode_pipeline.get_next_state(pos, state, codon)
                    new_candidates.append(
                        (cum_log_prob + step_log_prob, prefix_codons + [codon], next_state)
                    )

            if not new_candidates:
                # All beam branches blocked -> Attempt Partial DP rescue
                cur_chain = list(fallback_chain or [])
                cur_chain.append("sllm_beam_dead_end")
                best_prior_prefix = beams[0][1] if beams else []
                return self._execute_partial_dp_rescue(
                    clean_aa,
                    effective_stop,
                    deadend_pos=pos,
                    prefix_codons=best_prior_prefix,
                    requested_generator=requested_generator,
                    fallback_chain=cur_chain,
                )

            # Prune and retain top beam_width branches by cumulative log probability
            new_candidates.sort(key=lambda x: x[0], reverse=True)
            beams = new_candidates[: self.beam_width]

        # Best beam candidate
        best_log_prob, best_codons, _ = beams[0]
        cds_sequence = "".join(best_codons) + effective_stop

        # Final independent validation check
        val_res = self.contract.final_validator.evaluate(
            sequence_dna=cds_sequence,
            expected_protein_aa=clean_aa,
            expected_stop=effective_stop,
            upstream_flank=self.contract.upstream_flank,
            downstream_flank=self.contract.downstream_flank,
        )

        if not val_res.is_feasible:
            # Stop codon or construct boundary triggered violation -> exact DP rescue
            cur_chain = list(fallback_chain or [])
            cur_chain.append("sllm_final_validation_failure")
            return self._execute_dp_rescue(
                clean_aa,
                effective_stop,
                requested_generator=requested_generator,
                fallback_chain=cur_chain,
            )

        return {
            "sequence": cds_sequence,
            "candidate_id": "cand_sllm_beam_best",
            "requested_generator": requested_generator,
            "actual_generator": "sllm_hybrid_beam",
            "status": "SUCCESS",
            "score": round(best_log_prob, 4),
            "score_type": "constrained_sequence_log_likelihood",
            "rescue_triggered": False,
            "fallback_chain": fallback_chain or [],
            "rationale": "Constrained sLLM beam search (approximate MAP) under canonical DesignContract veto.",
            "model_name": self.model_adapter.model_name,
            "beam_width": self.beam_width,
            "contract_digest": self.contract.contract_digest,
        }


class ConstrainedSampleDecoder:
    """Stochastic temperature-scaled sampling decoder under canonical ResolvedDesignContract."""

    def __init__(
        self,
        contract: Optional[ResolvedDesignContract] = None,
        model_adapter: Optional[SLMModelAdapter] = None,
        temperature: float = 0.7,
        host: str = "nbenthamiana",
        hard_config: Optional[HardConstraintConfig] = None,
        pipeline: Optional[ConstraintProcessorPipeline] = None,
        seed: Optional[int] = None,
    ) -> None:
        if contract is not None:
            self.contract = contract
        else:
            self.contract = DesignContractCompiler.resolve_contract(
                host=host,
                config=hard_config,
            )

        self.model_adapter = model_adapter or DeterministicMockSLMAdapter(host=self.contract.host)
        self.temperature = max(0.1, temperature)
        self.rng = np.random.default_rng(seed)

    def decode(
        self,
        protein_sequence: str,
        stop_codon: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Samples a diverse constrained hypothesis sequence (cand_sllm_sample)."""
        clean_aa = "".join(protein_sequence.upper().split()).rstrip("*")
        if not clean_aa:
            raise ValueError("Cannot decode empty protein sequence.")

        effective_stop = stop_codon or self.contract.expected_stop

        prefix_codons: List[str] = []
        state = ConstraintState(
            automaton_node=0,
            position=0,
            trailing_nt=self.contract.upstream_flank[-6:] if self.contract.upstream_flank else "",
        )
        cum_log_prob = 0.0

        for pos, aa in enumerate(clean_aa):
            raw_logits = self.model_adapter.predict_logits(prefix_codons, aa, pos)
            masked_logits = self.contract.decode_pipeline(pos, aa, state, raw_logits)

            probs = _softmax_with_temperature(masked_logits.scores, temperature=self.temperature)
            if np.all(probs == 0.0):
                # Sampling dead end -> Attempt Local Partial DP Rescue on sampled prefix
                beam_retry = ConstrainedBeamDecoder(
                    contract=self.contract,
                    model_adapter=self.model_adapter,
                )
                return beam_retry._execute_partial_dp_rescue(
                    protein_sequence=clean_aa,
                    stop_codon=effective_stop,
                    deadend_pos=pos,
                    prefix_codons=prefix_codons,
                    requested_generator="sllm_hybrid_sample",
                    fallback_chain=["sllm_sample_dead_end"],
                )

            # Sample from normalized probability distribution
            chosen_idx = int(self.rng.choice(len(probs), p=probs))
            chosen_codon = CODON_VOCAB[chosen_idx]

            # Record normalized step log probability
            log_probs = _log_softmax_masked(masked_logits.scores)
            cum_log_prob += float(log_probs[chosen_idx])

            state = self.contract.decode_pipeline.get_next_state(pos, state, chosen_codon)
            prefix_codons.append(chosen_codon)

        cds_sequence = "".join(prefix_codons) + effective_stop

        # Final independent validation check
        val_res = self.contract.final_validator.evaluate(
            sequence_dna=cds_sequence,
            expected_protein_aa=clean_aa,
            expected_stop=effective_stop,
            upstream_flank=self.contract.upstream_flank,
            downstream_flank=self.contract.downstream_flank,
        )

        if not val_res.is_feasible:
            # Fallback to beam search or DP rescue
            beam_retry = ConstrainedBeamDecoder(
                contract=self.contract,
                model_adapter=self.model_adapter,
            )
            return beam_retry.decode(
                protein_sequence=clean_aa,
                stop_codon=effective_stop,
                requested_generator="sllm_hybrid_sample",
                fallback_chain=["sllm_sample_final_validation_failure"],
            )

        return {
            "sequence": cds_sequence,
            "candidate_id": "cand_sllm_sample",
            "requested_generator": "sllm_hybrid_sample",
            "actual_generator": "sllm_hybrid_sample",
            "status": "SUCCESS",
            "score": round(cum_log_prob, 4),
            "score_type": "constrained_sample_log_likelihood",
            "rescue_triggered": False,
            "fallback_chain": [],
            "rationale": "Constrained sLLM temperature sampling exploring non-obvious codon regularities.",
            "model_name": self.model_adapter.model_name,
            "temperature": self.temperature,
            "contract_digest": self.contract.contract_digest,
        }
