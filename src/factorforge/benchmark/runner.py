from typing import List, Dict, Any, Type
import time
import hashlib

from factorforge.engines import EngineRegistry, register_builtin_engines
from factorforge.evaluation.evaluator import SharedEvaluator
from factorforge.core.interfaces import OptimizerEngine
from factorforge.benchmark.models import (
    BenchmarkTarget,
    BenchmarkRunConfig,
    EngineBenchmarkResult,
    BenchmarkSuiteReport,
    BenchmarkProvenance
)
from factorforge.analysis.metrics import load_codon_usage_table

# Initialize engines
register_builtin_engines()

class BenchmarkRunner:
    """AgentOS Benchmark Harness for evaluating sequence generation methods."""

    def __init__(self, config: BenchmarkRunConfig):
        self.config = config
        
        # Load canonical host table to enable CAI in the SharedEvaluator
        # Default internal path loads N. benthamiana by default, which is our benchmark host.
        table = load_codon_usage_table() 
        
        self.evaluator = SharedEvaluator(
            version="1.0.1", 
            codon_weights=table.codon_weights
        )

    def run_suite(self, targets: List[BenchmarkTarget], engine_names: List[str]) -> BenchmarkSuiteReport:
        """Execute the benchmark suite across requested engines and targets."""
        
        engines: Dict[str, OptimizerEngine] = {}
        for name in engine_names:
            engine_instance = EngineRegistry.get(name)
            if not engine_instance:
                raise ValueError(f"Engine '{name}' not found in registry.")
            engines[name] = engine_instance

        results: List[EngineBenchmarkResult] = []

        for target in targets:
            input_hash = hashlib.sha256(target.protein_sequence.encode()).hexdigest()

            for engine_name, engine in engines.items():
                
                mode = getattr(engine, "inference_mode", "standard")
                
                provenance = BenchmarkProvenance(
                    input_sequence_hash=input_hash,
                    input_type="protein",
                    terminal_stop_policy=self.config.terminal_stop_policy,
                    codon_reference_id=f"host_{self.config.host}",
                    constraint_profile_id="default_gc",
                    evaluator_version=self.evaluator.version,
                    method_version=engine.version,
                    git_commit="current", # In real usage, this would be git rev-parse HEAD
                    seed=self.config.seed,
                    inference_mode=mode,
                    runtime_environment="default"
                )

                start_time = time.time()
                
                # Execute generation
                opt_result = engine.optimize(
                    sequence=target.protein_sequence,
                    host=self.config.host,
                    target_gc_min=self.config.target_gc_min_percent / 100.0, # Some engines still expect fraction
                    target_gc_max=self.config.target_gc_max_percent / 100.0,
                    forbidden_type_iis=set(self.config.forbidden_type_iis),
                    terminal_stop_policy=self.config.terminal_stop_policy,
                    seed=self.config.seed
                )
                
                runtime = time.time() - start_time
                
                # ONE SharedEvaluator for all generated candidates
                eval_result_obj = self.evaluator.evaluate_candidate(
                    candidate_dna=opt_result.sequence,
                    expected_protein=target.protein_sequence,
                    candidate_id=f"{engine_name}-{target.target_id}",
                    target_gc_min_percent=self.config.target_gc_min_percent,
                    target_gc_max_percent=self.config.target_gc_max_percent,
                    forbidden_type_iis=set(self.config.forbidden_type_iis)
                )

                res = EngineBenchmarkResult(
                    engine_name=engine_name,
                    target_id=target.target_id,
                    provenance=provenance,
                    optimized_sequence=opt_result.sequence,
                    evaluation=eval_result_obj,
                    runtime_seconds=runtime
                )
                results.append(res)
                
        return self._aggregate_results(results)

    def _aggregate_results(self, results: List[EngineBenchmarkResult]) -> BenchmarkSuiteReport:
        """Aggregate results into a full BenchmarkSuiteReport."""
        
        pass_counts = {}
        total_targets = len(set(r.target_id for r in results))
        
        cai_sums = {}
        cai_vals = {}
        gc_sums = {}
        engine_counts = {}
        
        for r in results:
            eng = r.engine_name
            if eng not in pass_counts:
                pass_counts[eng] = 0
                cai_vals[eng] = []
                gc_sums[eng] = 0.0
                engine_counts[eng] = 0
                
            if r.evaluation.passed:
                pass_counts[eng] += 1
                
            if r.evaluation.metrics.cai is not None:
                cai_vals[eng].append(r.evaluation.metrics.cai)
                
            gc_sums[eng] += r.evaluation.metrics.gc_percent
            engine_counts[eng] += 1
            
        pass_counts_str = {e: f"{count}/{total_targets}" for e, count in pass_counts.items()}
        avg_cai = {e: (sum(cai_vals[e]) / len(cai_vals[e]) if cai_vals[e] else None) for e in pass_counts}
        avg_gc = {e: (gc_sums[e] / engine_counts[e]) for e in pass_counts}
        
        return BenchmarkSuiteReport(
            suite_name=self.config.suite_name,
            config=self.config,
            results=results,
            pass_counts=pass_counts_str,
            average_cai=avg_cai,
            average_gc=avg_gc
        )
