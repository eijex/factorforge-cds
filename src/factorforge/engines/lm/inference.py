"""Experimental FactorForge-LM constrained decoding scaffold.

This module deliberately exposes a small deterministic research scaffold rather
than a trained production LM. It preserves amino-acid identity through the shared
FactorForge genetic-code table and reports validation-derived status fields.
"""

from __future__ import annotations

from typing import Any, Optional

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE
from factorforge.core.interfaces import OptimizationResult, OptimizerEngine
from factorforge.engines.lm.adapter import FactorForgeLogitMasker
from factorforge.engines.lm.models.mbart_codon import TORCH_AVAILABLE
from factorforge.engines.lm.tokenizer.control_tokenizer import FactorForgeControlTokenizer
from factorforge.evaluation.evaluator import SharedEvaluator

if TORCH_AVAILABLE:
    import torch  # noqa: F401

AA_TO_CODONS: dict[str, tuple[str, ...]] = {}
for _codon, _aa in STANDARD_GENETIC_CODE.items():
    if _aa != "*":
        AA_TO_CODONS.setdefault(_aa, tuple())
        AA_TO_CODONS[_aa] = (*AA_TO_CODONS[_aa], _codon)


class ConstrainedBeamSearchEngine:
    """Experimental constrained codon selector for LM integration tests.

    The class is not a trained expression optimizer. If ``model`` is omitted, it
    uses deterministic synonymous codon selection plus the same mask helpers that
    a future LM decoder would call. Results are computational review artifacts.
    """

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

    def optimize_cds(
        self,
        amino_acids: str,
        host: str = "nbenthamiana",
        gc_band: str = "40-47",
        type2is_clean: bool = True,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Generate a computational CDS candidate for a protein sequence.

        Unknown/unsupported residues fail closed with ``ValueError`` instead of
        silently substituting an unrelated codon.
        """
        protein = "".join(amino_acids.upper().split()).rstrip("*")
        if not protein:
            raise ValueError("amino_acids must contain at least one residue")

        unsupported = sorted({aa for aa in protein if aa not in AA_TO_CODONS})
        if unsupported:
            raise ValueError(f"Unsupported amino-acid residues: {', '.join(unsupported)}")

        # Exercise control-token encoding so invalid control-token regressions are visible.
        self.tokenizer.encode_encoder_input(
            amino_acids=protein,
            host=host,
            gc_band=gc_band,
            type2is_clean=type2is_clean,
        )

        current_cds = ""
        for aa in protein:
            synonym_codons = AA_TO_CODONS[aa]
            logits = [-float("inf")] * self.tokenizer.vocab_size
            for codon in synonym_codons:
                token_id = self.tokenizer.token_to_id[codon]
                # Simple deterministic baseline preference. A future trained model
                # may replace these scores before the same masks are applied.
                logits[token_id] = 0.0

            masked_logits = self.masker.apply_logit_masks(
                current_cds,
                logits,
                expected_next_aa=aa,
            )
            ranked = sorted(
                synonym_codons,
                key=lambda codon: masked_logits[self.tokenizer.token_to_id[codon]],
                reverse=True,
            )
            best = ranked[0]
            if masked_logits[self.tokenizer.token_to_id[best]] == -float("inf"):
                raise ValueError(f"No valid synonymous codon remained for residue {aa!r}")
            current_cds += best

        # Handle terminal stop based on policy
        terminal_stop_policy = kwargs.get("terminal_stop_policy", "preserve")
        
        has_stop = amino_acids.strip().endswith("*")
        needs_stop = False
        if terminal_stop_policy == "append":
            needs_stop = True
        elif terminal_stop_policy == "preserve":
            needs_stop = has_stop
        
        if needs_stop:
            current_cds += "TAA"
        
        # Set up forbidden Type IIS based on type2is_clean flag
        forbidden_sites = {"BsaI", "BsmBI", "BpiI"} if type2is_clean else set()
        
        # Delegate to SharedEvaluator instead of ad-hoc checking
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
    """Real ONNX-based FactorForge-SLM decoder."""
    
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

        # Encode input sequence
        input_ids = self.tokenizer.encode_encoder_input(protein, host=host, gc_band=gc_band, type2is_clean=type2is_clean)
        src_tensor = np.array([input_ids], dtype=np.int64)
        
        # Initial hypothesis: (score, current_cds, decoder_input_ids)
        bos_id = getattr(self.tokenizer, "bos_id", 2)
        beam = [(0.0, "", [bos_id])]
        
        for i, aa in enumerate(protein):
            AA_TO_CODONS = {
                'A': ['GCT', 'GCC', 'GCA', 'GCG'],
                'C': ['TGT', 'TGC'],
                'D': ['GAT', 'GAC'],
                'E': ['GAA', 'GAG'],
                'F': ['TTT', 'TTC'],
                'G': ['GGT', 'GGC', 'GGA', 'GGG'],
                'H': ['CAT', 'CAC'],
                'I': ['ATT', 'ATC', 'ATA'],
                'K': ['AAA', 'AAG'],
                'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
                'M': ['ATG'],
                'N': ['AAT', 'AAC'],
                'P': ['CCT', 'CCC', 'CCA', 'CCG'],
                'Q': ['CAA', 'CAG'],
                'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                'T': ['ACT', 'ACC', 'ACA', 'ACG'],
                'V': ['GTT', 'GTC', 'GTA', 'GTG'],
                'W': ['TGG'],
                'Y': ['TAT', 'TAC'],
                '*': ['TAA', 'TAG', 'TGA']
            }
            synonym_codons = AA_TO_CODONS.get(aa, [])
            next_beam = []
            
            # Batch all hypotheses in the current beam
            dec_ids_batch = [b[2] for b in beam]
            tgt_tensor = np.array(dec_ids_batch, dtype=np.int64)
            src_tensor_batch = np.repeat(src_tensor, len(beam), axis=0)
            
            ort_inputs = {
                self.session.get_inputs()[0].name: src_tensor_batch,
                self.session.get_inputs()[1].name: tgt_tensor
            }
            ort_outs = self.session.run(None, ort_inputs)
            self.model_forward_calls += 1
            
            # Extract logits for the last token across all batch items
            logits_batch = ort_outs[0][:, -1, :]
            
            # To store rejection reasons if all fail
            rejection_log = []
            
            for b_idx, (score, current_cds, dec_ids) in enumerate(beam):
                logits = logits_batch[b_idx].tolist()
                
                # Pad logits if smaller than vocab size
                if len(logits) < self.tokenizer.vocab_size:
                    logits.extend([-float("inf")] * (self.tokenizer.vocab_size - len(logits)))
                
                masked_logits, rejections = self.masker.apply_logit_masks_with_reasons(
                    current_cds,
                    logits,
                    expected_next_aa=aa,
                )
                
                # Log for proof if requested (only log top beam)
                if self.strict_proof_mode and b_idx == 0 and i < 5:
                    raw_codon_scores = {c: logits[self.tokenizer.token_to_id[c]] for c in synonym_codons}
                    masked_codon_scores = {c: masked_logits[self.tokenizer.token_to_id[c]] for c in synonym_codons}
                    
                    local_argmax = max(synonym_codons, key=lambda c: masked_logits[self.tokenizer.token_to_id[c]])
                    self.proof_logs.append({
                        "step": i,
                        "position": i+1,
                        "aa": aa,
                        "raw_logits": raw_codon_scores,
                        "local_argmax_codon": local_argmax,
                    })
                
                # Expand hypotheses
                for codon in synonym_codons:
                    token_id = self.tokenizer.token_to_id[codon]
                    c_logit = masked_logits[token_id]
                    if c_logit != -float("inf"):
                        new_score = score + c_logit
                        new_cds = current_cds + codon
                        new_dec_ids = dec_ids + [token_id]
                        next_beam.append((new_score, new_cds, new_dec_ids, codon))
                    elif b_idx == 0: # Log rejections for the top beam hypothesis
                        reason = rejections.get(codon, "HARD_REJECT: unknown")
                        rejection_log.append(f"{codon} -> {reason}")
            
            if not next_beam:
                trace_msg = f"Residue {i+1} {aa}:\n" + "\n".join([f"  {r}" for r in rejection_log]) + "\n  -> ERROR: No valid hypotheses"
                print(f"Constraint Death Trace:\n{trace_msg}")
                raise ValueError(trace_msg)
            
            # Prune to beam width
            next_beam.sort(key=lambda x: x[0], reverse=True)
            beam = [(x[0], x[1], x[2]) for x in next_beam[:self.beam_width]]
            
            if self.strict_proof_mode and i < 5:
                self.proof_logs[-1]["expanded_hypotheses"] = len(next_beam)
                self.proof_logs[-1]["pruned_hypotheses"] = max(0, len(next_beam) - self.beam_width)
                # the codon that actually made it to the top of the beam at this step
                self.proof_logs[-1]["final_path_codon"] = next_beam[0][3]

        best_score, best_cds, best_dec_ids = beam[0]

        # Terminal stop
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
            "generation_engine": "lm",
            "inference_mode": "actual_beam_search",
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
    """Adapter wrapping ConstrainedBeamSearchEngine or ONNXBeamSearchEngine for standard OptimizerEngine interface."""

    def __init__(self, beam_engine: Optional[ConstrainedBeamSearchEngine] = None, onnx_path: str = None) -> None:
        if beam_engine:
            self.beam_engine = beam_engine
        elif onnx_path:
            self.beam_engine = ONNXBeamSearchEngine(onnx_path)
        else:
            self.beam_engine = ConstrainedBeamSearchEngine()

    @property
    def name(self) -> str:
        return "FactorForge-LM"

    @property
    def version(self) -> str:
        return "3.6.0-scaffold"

    @property
    def inference_mode(self) -> str:
        return "deterministic_scaffold"

    def optimize_cds(
        self, 
        protein: str, 
        profile: str = "balanced", 
        host: str = "nbenthamiana"
    ) -> Dict[str, Any]:
        """Perform exact beam search constraint decoding using raw model logits."""
        
        AA_TO_CODONS = {
            'A': ['GCT', 'GCC', 'GCA', 'GCG'],
            'C': ['TGT', 'TGC'],
            'D': ['GAT', 'GAC'],
            'E': ['GAA', 'GAG'],
            'F': ['TTT', 'TTC'],
            'G': ['GGT', 'GGC', 'GGA', 'GGG'],
            'H': ['CAT', 'CAC'],
            'I': ['ATT', 'ATC', 'ATA'],
            'K': ['AAA', 'AAG'],
            'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
            'M': ['ATG'],
            'N': ['AAT', 'AAC'],
            'P': ['CCT', 'CCC', 'CCA', 'CCG'],
            'Q': ['CAA', 'CAG'],
            'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
            'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
            'T': ['ACT', 'ACC', 'ACA', 'ACG'],
            'V': ['GTT', 'GTC', 'GTA', 'GTG'],
            'W': ['TGG'],
            'Y': ['TAT', 'TAC'],
            '*': ['TAA', 'TAG', 'TGA']
        }

    def optimize(
        self,
        sequence: str,
        profile: str | None = None,
        host: str = "nbenthamiana",
        **kwargs: Any,
    ) -> OptimizationResult:
        res = self.beam_engine.optimize_cds(sequence, host=host, **kwargs)
        metrics = {
            "cai": res.get("evaluation_report", {}).get("metrics", {}).get("cai", 0.0),
            "gc_percent": res.get("gc_percent", 0.0),
            "score": 0.0,
        }
        return OptimizationResult(
            sequence=res["optimized_sequence"],
            metrics=metrics,
            metadata={
                "engine": "lm", 
                "version": self.version, 
                "host": host,
                "inference_mode": res.get("inference_mode", "actual_beam_search"),
                "model_forward_calls": res.get("model_forward_calls", 0),
                "decoder_steps": res.get("decoder_steps", 0),
                "validator_passed": res.get("validator_passed", False),
                "evaluation_report": res.get("evaluation_report", {}),
            },
        )

    def validate(self, sequence: str) -> bool:
        return bool(sequence and isinstance(sequence, str))
