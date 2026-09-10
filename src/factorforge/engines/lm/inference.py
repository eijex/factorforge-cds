"""Experimental FactorForge-LM constrained decoding scaffold."""

from __future__ import annotations

from typing import Any, Optional, Dict
import math

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE
from factorforge.core.interfaces import OptimizationResult, OptimizerEngine
from factorforge.engines.lm.adapter import FactorForgeLogitMasker
from factorforge.engines.lm.models.mbart_codon import TORCH_AVAILABLE
from factorforge.engines.lm.tokenizer.control_tokenizer import FactorForgeControlTokenizer
from factorforge.evaluation.evaluator import SharedEvaluator

# --- [Phase 4] Neuro-symbolic Hybrid Engine Imports ---
from factorforge.engines.sllm.interfaces import (
    CodonLogits, ConstraintProcessorPipeline, ConstraintState,
    SynonymousMaskProcessor, AutomatonConstraintProcessor, CODON_VOCAB
)
from factorforge.engines.sllm.automaton import AutomatonCompiler
# --------------------------------------------------------

if TORCH_AVAILABLE:
    import torch  # noqa: F401

AA_TO_CODONS: dict[str, tuple[str, ...]] = {}
for _codon, _aa in STANDARD_GENETIC_CODE.items():
    if _aa != "*":
        AA_TO_CODONS.setdefault(_aa, tuple())
        AA_TO_CODONS[_aa] = (*AA_TO_CODONS[_aa], _codon)


class ConstrainedBeamSearchEngine:
    def __init__(
        self,
        tokenizer: Optional[FactorForgeControlTokenizer] = None,
        model: Optional[Any] = None,
        masker: Optional[FactorForgeLogitMasker] = None,
        beam_width: int = 5,
        target_gc_min: float = 0.40,
        target_gc_max: float = 0.47,
    ) -> None:
        self.tokenizer = tokenizer or FactorForgeControlTokenizer()
        self.target_gc_min = target_gc_min
        self.target_gc_max = target_gc_max
        self.masker = masker or FactorForgeLogitMasker(
            tokenizer=self.tokenizer,
            target_gc_min=self.target_gc_min,
            target_gc_max=self.target_gc_max,
        )
        self.model = model
        self.beam_width = beam_width
        self.evaluator = SharedEvaluator(version="1.0.0")

    def optimize_cds(self, amino_acids: str, host: str = "nbenthamiana", gc_band: str = "40-47", type2is_clean: bool = True, **kwargs: Any) -> dict[str, Any]:
        protein = "".join(amino_acids.upper().split()).rstrip("*")
        if not protein:
            raise ValueError("amino_acids must contain at least one residue")

        unsupported = sorted({aa for aa in protein if aa not in AA_TO_CODONS})
        if unsupported:
            raise ValueError(f"Unsupported amino-acid residues: {', '.join(unsupported)}")

        self.tokenizer.encode_encoder_input(amino_acids=protein, host=host, gc_band=gc_band, type2is_clean=type2is_clean)

        current_cds = ""
        for aa in protein:
            synonym_codons = AA_TO_CODONS[aa]
            logits = [-float("inf")] * self.tokenizer.vocab_size
            for codon in synonym_codons:
                token_id = self.tokenizer.token_to_id[codon]
                logits[token_id] = 0.0

            masked_logits = self.masker.apply_logit_masks(current_cds, logits, expected_next_aa=aa)
            ranked = sorted(synonym_codons, key=lambda codon: masked_logits[self.tokenizer.token_to_id[codon]], reverse=True)
            best = ranked[0]
            if masked_logits[self.tokenizer.token_to_id[best]] == -float("inf"):
                raise ValueError(f"No valid synonymous codon remained for residue {aa!r}")
            current_cds += best

        terminal_stop_policy = kwargs.get("terminal_stop_policy", "preserve")
        has_stop = amino_acids.strip().endswith("*")
        needs_stop = True if terminal_stop_policy == "append" else (has_stop if terminal_stop_policy == "preserve" else False)
        if needs_stop:
            current_cds += "TAA"
        
        forbidden_sites = {"BsaI", "BsmBI", "BpiI"} if type2is_clean else set()
        
        eval_result = self.evaluator.evaluate_candidate(
            candidate_dna=current_cds,
            expected_protein=protein,
            candidate_id="lm-candidate-01",
            target_gc_min_percent=self.target_gc_min * 100 if self.target_gc_min and self.target_gc_min <= 1.0 else self.target_gc_min,
            target_gc_max_percent=self.target_gc_max * 100 if self.target_gc_max and self.target_gc_max <= 1.0 else self.target_gc_max,
            forbidden_type_iis=forbidden_sites,
        )

        return {
            "generation_engine": "lm",
            "inference_mode": "deterministic_scaffold",
            "trained_model_used": False,
            "optimized_sequence": current_cds,
            "sequence_length": len(current_cds),
            "gc_percent": eval_result.metrics.gc_percent,
            "validator_passed": eval_result.passed,
            "evaluation_report": eval_result.model_dump(),
        }


class ONNXBeamSearchEngine(ConstrainedBeamSearchEngine):
    """Real ONNX-based FactorForge-SLM decoder integrated with Neuro-Symbolic Pipeline."""
    
    def __init__(self, onnx_model_path: str, beam_width: int = 5, strict_proof_mode: bool = False, **kwargs):
        super().__init__(**kwargs)
        self.onnx_model_path = onnx_model_path
        self.beam_width = beam_width
        self.strict_proof_mode = strict_proof_mode
        self.proof_logs = []
        try:
            import onnxruntime as ort
            self.session = ort.InferenceSession(self.onnx_model_path, providers=['CPUExecutionProvider'])
            self.ort_available = True
        except ImportError:
            self.session = None
            self.ort_available = False
        except Exception as e:
            self.session = None
            self.ort_available = False
            print(f"Failed to load ONNX model: {e}")

    def optimize_cds(self, amino_acids: str, host: str = "nbenthamiana", gc_band: str = "40-47", type2is_clean: bool = True, **kwargs) -> dict:
        if not self.ort_available or not self.session:
            raise RuntimeError("ONNX Runtime missing or model failed to load. Fallback to scaffold is disabled for audit.")
            
        import numpy as np
        self.model_forward_calls = 0
        
        protein = "".join(amino_acids.upper().split()).rstrip("*")
        if not protein:
            raise ValueError("amino_acids must contain at least one residue")

        # --- Pipeline Setup (Canonical Registry Parity) ---
        motifs_to_ban = []
        if type2is_clean:
            from factorforge.evaluation.evaluator import TYPE_IIS_MOTIFS
            for enzyme in ["BsaI", "BsmBI", "BpiI"]:
                motifs_to_ban.extend(TYPE_IIS_MOTIFS.get(enzyme, []))
            
        compiled_automaton = AutomatonCompiler.compile(motifs_to_ban, include_rc=True)
        pipeline = ConstraintProcessorPipeline([
            SynonymousMaskProcessor(),
            AutomatonConstraintProcessor(compiled_automaton)
        ])
        # ----------------------

        input_ids = self.tokenizer.encode_encoder_input(protein, host=host, gc_band=gc_band, type2is_clean=type2is_clean)
        src_tensor = np.array([input_ids], dtype=np.int64)
        
        bos_id = getattr(self.tokenizer, "bos_id", 2)
        # 빔 초기화: (score, current_cds, dec_ids, constraint_state)
        beam = [(0.0, "", [bos_id], ConstraintState(automaton_node=0, position=0))]
        
        for i, aa in enumerate(protein):
            synonym_codons = AA_TO_CODONS.get(aa, [])
            next_beam = []
            
            dec_ids_batch = [b[2] for b in beam]
            tgt_tensor = np.array(dec_ids_batch, dtype=np.int64)
            src_tensor_batch = np.repeat(src_tensor, len(beam), axis=0)
            
            ort_inputs = {
                self.session.get_inputs()[0].name: src_tensor_batch,
                self.session.get_inputs()[1].name: tgt_tensor
            }
            ort_outs = self.session.run(None, ort_inputs)
            self.model_forward_calls += 1
            
            logits_batch = ort_outs[0][:, -1, :]
            rejection_log = []
            
            for b_idx, (score, current_cds, dec_ids, state) in enumerate(beam):
                raw_model_logits = logits_batch[b_idx].tolist()
                
                # 1. Adapter: Model Vocab -> CodonLogits(64)
                raw_scores = np.zeros(64, dtype=np.float32)
                for c_idx, c in enumerate(CODON_VOCAB):
                    tok_id = self.tokenizer.token_to_id.get(c, -1)
                    if tok_id != -1 and tok_id < len(raw_model_logits):
                        raw_scores[c_idx] = raw_model_logits[tok_id]
                    else:
                        raw_scores[c_idx] = -float('inf')
                
                codon_logits = CodonLogits(raw_scores)
                
                # 2. Pipeline Execution (Predictive Masking)
                try:
                    masked_codon_logits = pipeline(i, aa, state, codon_logits)
                except RuntimeError as e:
                    if str(e) == "BEAM_DEAD_END":
                        rejection_log.append(f"Beam {b_idx} reached dead end.")
                        continue
                    raise
                
                # 3. Expand Hypotheses
                for c_idx, codon in enumerate(CODON_VOCAB):
                    c_logit = masked_codon_logits.scores[c_idx]
                    if c_logit > -float("inf"):
                        token_id = self.tokenizer.token_to_id.get(codon, -1)
                        if token_id == -1: continue
                        
                        new_score = score + c_logit
                        new_cds = current_cds + codon
                        new_dec_ids = dec_ids + [token_id]
                        # 상태 업데이트 (State Transition)
                        new_state = pipeline.get_next_state(i, state, codon)
                        
                        # Store codon for proof logs
                        next_beam.append((new_score, new_cds, new_dec_ids, new_state, codon))
            
            if not next_beam:
                trace_msg = f"Residue {i+1} {aa}:\n" + "\n".join([f"  {r}" for r in rejection_log]) + "\n  -> ERROR: SEARCH_EXHAUSTED"
                print(f"Constraint Death Trace:\n{trace_msg}")
                raise ValueError(trace_msg)
            
            # Prune to beam width
            next_beam.sort(key=lambda x: x[0], reverse=True)
            beam = [(x[0], x[1], x[2], x[3]) for x in next_beam[:self.beam_width]]
            
        best_score, best_cds, best_dec_ids, _ = beam[0]

        terminal_stop_policy = kwargs.get("terminal_stop_policy", "preserve")
        has_stop = amino_acids.strip().endswith("*")
        needs_stop = True if terminal_stop_policy == "append" else has_stop
        if needs_stop:
            best_cds += "TAA"
            
        forbidden_sites = {"BsaI", "BsmBI", "BpiI"} if type2is_clean else set()
        eval_result = self.evaluator.evaluate_candidate(
            candidate_dna=best_cds,
            expected_protein=protein,
            candidate_id="lm-onnx-candidate-01",
            target_gc_min_percent=self.target_gc_min * 100 if self.target_gc_min <= 1.0 else self.target_gc_min,
            target_gc_max_percent=self.target_gc_max * 100 if self.target_gc_max <= 1.0 else self.target_gc_max,
            forbidden_type_iis=forbidden_sites,
        )

        return {
            "generation_engine": "hybrid_sllm",
            "inference_mode": "constraint_aware_beam_search",
            "trained_model_used": True,
            "beam_width": self.beam_width,
            "final_beam_score": float(best_score),
            "optimized_sequence": best_cds,
            "sequence_length": len(best_cds),
            "model_forward_calls": self.model_forward_calls,
            "decoder_steps": len(protein),
            "gc_percent": eval_result.metrics.gc_percent,
            "validator_passed": eval_result.passed,
            "evaluation_report": eval_result.model_dump(),
            "proof_logs": self.proof_logs if self.strict_proof_mode else []
        }

class LMEngineAdapter(OptimizerEngine):
    def __init__(self, beam_engine: Optional[ConstrainedBeamSearchEngine] = None, onnx_path: str = None) -> None:
        if beam_engine:
            self.beam_engine = beam_engine
        elif onnx_path:
            self.beam_engine = ONNXBeamSearchEngine(onnx_path)
        else:
            self.beam_engine = ConstrainedBeamSearchEngine()

    @property
    def name(self) -> str:
        return "FactorForge-LM-Hybrid"

    @property
    def version(self) -> str:
        return "4.0.0-hybrid"

    def optimize(self, sequence: str, profile: str | None = None, host: str = "nbenthamiana", **kwargs: Any) -> OptimizationResult:
        try:
            res = self.beam_engine.optimize_cds(sequence, host=host, **kwargs)
        except ValueError as e:
            if "SEARCH_EXHAUSTED" in str(e):
                from factorforge.engines.lm.diagnostics import diagnose_global_feasibility
                from factorforge.evaluation.evaluator import TYPE_IIS_MOTIFS
                
                type2is_clean = kwargs.get("type2is_clean", True)
                motifs = []
                if type2is_clean:
                    for enzyme in ["BsaI", "BsmBI", "BpiI"]:
                        motifs.extend(TYPE_IIS_MOTIFS.get(enzyme, []))
                        
                diagnosis = diagnose_global_feasibility(sequence, motifs)
                if not diagnosis["global_feasible"]:
                    raise RuntimeError(f"GLOBAL_INFEASIBLE: {diagnosis['message']}") from e
                else:
                    raise RuntimeError(f"SLM_SEARCH_LIMITATION: {diagnosis['message']}") from e
            raise

        metrics = {
            "cai": res.get("evaluation_report", {}).get("metrics", {}).get("cai", 0.0),
            "gc_percent": res.get("gc_percent", 0.0),
            "score": 0.0,
        }
        return OptimizationResult(
            sequence=res["optimized_sequence"],
            metrics=metrics,
            metadata={
                "engine": "hybrid_sllm", 
                "version": self.version, 
                "host": host,
                "inference_mode": res.get("inference_mode", "constraint_aware_beam_search"),
                "model_forward_calls": res.get("model_forward_calls", 0),
                "decoder_steps": res.get("decoder_steps", 0),
                "validator_passed": res.get("validator_passed", False),
                "evaluation_report": res.get("evaluation_report", {}),
            },
        )

    def validate(self, sequence: str) -> bool:
        return bool(sequence and isinstance(sequence, str))
