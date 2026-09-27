# factorforge/benchmarks/run_paper3_ablation_benchmark.py
"""Paper 3 4-Arm Multi-Seed Ablation Benchmark Suite (Job 286B).

Evaluates 4 experimental arms across the benchmark panel with multi-seed sampling (N = 36 - 108 trials / arm):
1. Stage 1 (Arm 1): Raw Unconstrained sLLM (Measures baseline generative defect rate).
2. Stage 2 (Arm 2): Neuro-Symbolic Veto Only (Measures dead-end frequency without DP rescue).
3. Stage 3A (Arm 3A): Global DP Rescue Control (Job 285 Baseline, PPR = 0% on dead-end).
4. Stage 3B (Arm 3B): Adaptive Partial DP Rescue (Job 286 Innovation, Prefix Preservation & Diversity Retention).

Emits:
- benchmarks/results/paper3_ablation/ablation_metrics.json & .csv
- benchmarks/results/paper3_ablation/figure2_resolution_waterfall.html
"""

import os
import json
import time
import argparse
from typing import Dict, List, Any
import numpy as np
import pandas as pd

from factorforge.discovery.filter import HardConstraintConfig
from factorforge.engines.sllm.adapters import OnnxSLMAdapter, DeterministicMockSLMAdapter
from factorforge.engines.sllm.compiler import DesignContractCompiler
from factorforge.engines.sllm.decoder import (
    ConstrainedBeamDecoder,
    ConstrainedSampleDecoder,
    _softmax_with_temperature,
    CODON_VOCAB,
    ConstraintState,
)


def run_arm1_unconstrained_ai(
    protein: str,
    stop_codon: str,
    model_adapter,
    temperature: float = 0.7,
    seed: int = 42,
) -> Dict[str, Any]:
    """Arm 1: Unconstrained sLLM sampling (synonymous mapping only, zero automaton veto)."""
    from factorforge.engines.sllm.data.tokenizer import CodonSequenceTokenizer, AA_TO_INDEX

    tokenizer = CodonSequenceTokenizer()
    rng = np.random.default_rng(seed)

    clean_aa = protein.rstrip("*")
    prefix_codons = []

    for pos, aa in enumerate(clean_aa):
        raw_logits = model_adapter.predict_logits(prefix_codons, aa, pos)
        aa_idx = AA_TO_INDEX.get(aa, 0)
        syn_mask = tokenizer.get_synonymous_mask(aa_idx)

        scores = np.where(syn_mask, raw_logits.scores, -1e9)
        probs = _softmax_with_temperature(scores, temperature=temperature)
        if np.all(probs == 0.0):
            probs[syn_mask] = 1.0 / np.sum(syn_mask)

        chosen_idx = int(rng.choice(len(probs), p=probs))
        prefix_codons.append(CODON_VOCAB[chosen_idx])

    seq = "".join(prefix_codons) + stop_codon
    return {"sequence": seq, "status": "COMPLETED_UNCONSTRAINED", "prefix_preservation_rate": 1.0}


def run_arm2_veto_only(
    protein: str,
    stop_codon: str,
    contract,
    model_adapter,
    temperature: float = 0.7,
    seed: int = 42,
) -> Dict[str, Any]:
    """Arm 2: Neuro-symbolic automaton veto only (fails on dead-end with NO DP rescue)."""
    rng = np.random.default_rng(seed)
    clean_aa = protein.rstrip("*")
    prefix_codons = []
    state = ConstraintState(
        automaton_node=0,
        position=0,
        trailing_nt=contract.upstream_flank[-6:] if contract.upstream_flank else "",
    )

    for pos, aa in enumerate(clean_aa):
        raw_logits = model_adapter.predict_logits(prefix_codons, aa, pos)
        masked_logits = contract.decode_pipeline(pos, aa, state, raw_logits)
        probs = _softmax_with_temperature(masked_logits.scores, temperature=temperature)

        if np.all(probs == 0.0):
            return {
                "sequence": None,
                "status": "DEAD_END",
                "deadend_pos": pos,
                "prefix_preservation_rate": round(pos / float(len(clean_aa)), 4),
            }

        chosen_idx = int(rng.choice(len(probs), p=probs))
        chosen_codon = CODON_VOCAB[chosen_idx]
        state = contract.decode_pipeline.get_next_state(pos, state, chosen_codon)
        prefix_codons.append(chosen_codon)

    seq = "".join(prefix_codons) + stop_codon
    val = contract.final_validator.evaluate(
        sequence_dna=seq, expected_protein_aa=clean_aa, expected_stop=stop_codon
    )
    if not val.is_feasible:
        return {
            "sequence": seq,
            "status": "VETO_FINAL_INADMISSIBLE",
            "violations": val.violations,
            "prefix_preservation_rate": 1.0,
        }

    return {
        "sequence": seq,
        "status": "SUCCESS_VETO_ONLY",
        "prefix_preservation_rate": 1.0,
    }


def run_arm3a_global_dp_rescue(
    protein: str,
    stop_codon: str,
    contract,
    model_adapter,
    temperature: float = 0.7,
    seed: int = 42,
) -> Dict[str, Any]:
    """Arm 3A: Global DP Rescue control (Job 285 baseline: discard prefix on dead-end, PPR=0.0%)."""
    rng = np.random.default_rng(seed)
    clean_aa = protein.rstrip("*")
    prefix_codons = []
    state = ConstraintState(
        automaton_node=0,
        position=0,
        trailing_nt=contract.upstream_flank[-6:] if contract.upstream_flank else "",
    )

    for pos, aa in enumerate(clean_aa):
        raw_logits = model_adapter.predict_logits(prefix_codons, aa, pos)
        masked_logits = contract.decode_pipeline(pos, aa, state, raw_logits)
        probs = _softmax_with_temperature(masked_logits.scores, temperature=temperature)

        if np.all(probs == 0.0):
            # Dead end -> Global DP fallback (PPR = 0.0%)
            beam_decoder = ConstrainedBeamDecoder(contract=contract, model_adapter=model_adapter)
            res = beam_decoder._execute_dp_rescue(
                protein_sequence=clean_aa,
                stop_codon=stop_codon,
                requested_generator="sllm_hybrid_sample_arm3a",
                fallback_chain=["sllm_sample_deadend_global_fallback"],
            )
            res["prefix_preservation_rate"] = 0.0
            res["preserved_prefix_codons"] = 0
            return res

        chosen_idx = int(rng.choice(len(probs), p=probs))
        chosen_codon = CODON_VOCAB[chosen_idx]
        state = contract.decode_pipeline.get_next_state(pos, state, chosen_codon)
        prefix_codons.append(chosen_codon)

    seq = "".join(prefix_codons) + stop_codon
    val = contract.final_validator.evaluate(
        sequence_dna=seq, expected_protein_aa=clean_aa, expected_stop=stop_codon
    )
    if not val.is_feasible:
        beam_decoder = ConstrainedBeamDecoder(contract=contract, model_adapter=model_adapter)
        res = beam_decoder._execute_dp_rescue(
            protein_sequence=clean_aa,
            stop_codon=stop_codon,
            requested_generator="sllm_hybrid_sample_arm3a",
            fallback_chain=["sllm_sample_validation_failure_global_fallback"],
        )
        res["prefix_preservation_rate"] = 0.0
        res["preserved_prefix_codons"] = 0
        return res

    return {
        "sequence": seq,
        "status": "SUCCESS",
        "rescue_triggered": False,
        "prefix_preservation_rate": 1.0,
        "preserved_prefix_codons": len(clean_aa),
    }


def compute_pairwise_dnt(sequences: List[str]) -> float:
    """Computes mean pairwise normalized Hamming distance (Nucleotide Diversity D_nt)."""
    valid_seqs = [s for s in sequences if s]
    n = len(valid_seqs)
    if n < 2:
        return 0.0

    distances = []
    for i in range(n):
        for j in range(i + 1, n):
            s1 = valid_seqs[i]
            s2 = valid_seqs[j]
            min_len = min(len(s1), len(s2))
            if min_len == 0:
                continue
            diffs = sum(1 for a, b in zip(s1[:min_len], s2[:min_len]) if a != b)
            distances.append(diffs / float(min_len))

    return round(float(np.mean(distances)), 4) if distances else 0.0


def run_paper3_ablation_benchmark(
    output_dir: str,
    seeds: List[int] = (42, 100, 2026, 777),
):
    os.makedirs(output_dir, exist_ok=True)

    # Load targets
    fixture_path = r"c:\Work\eijex\factorforge\benchmarks\fixtures\discovery_novelty_panel.json"
    with open(fixture_path, "r", encoding="utf-8") as f:
        panel_data = json.load(f)

    targets = panel_data.get("targets", [])
    print(
        f"Loaded {len(targets)} benchmark targets across {len(seeds)} seeds ({len(targets) * len(seeds)} total trials/arm)."
    )

    # Resolve canonical contract
    contract = DesignContractCompiler.resolve_contract(
        host="nbenthamiana",
        config=HardConstraintConfig(
            forbidden_type_iis_enzymes={"BsaI", "BsmBI"},
            forbidden_other_motifs=["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"],
            homopolymer_max_run=5,
            global_gc_min=30.0,
            global_gc_max=60.0,
            initiation_45nt_gc_min=15.0,
            initiation_45nt_gc_max=45.0,
        ),
        expected_stop="TAA",
    )

    # Model adapter (try ONNX, fallback to Deterministic Mock)
    try:
        model_adapter = OnnxSLMAdapter()
        print("Using quantized ONNX sLLM adapter.")
    except Exception as e:
        print(f"ONNX adapter unavailable ({e}), using DeterministicMockSLMAdapter.")
        model_adapter = DeterministicMockSLMAdapter(host="nbenthamiana")

    records = []
    arm_sequences: Dict[str, Dict[str, List[str]]] = {
        "arm1": {t["benchmark_id"]: [] for t in targets},
        "arm2": {t["benchmark_id"]: [] for t in targets},
        "arm3a": {t["benchmark_id"]: [] for t in targets},
        "arm3b": {t["benchmark_id"]: [] for t in targets},
    }

    for seed in seeds:
        print(f"\n--- Executing Benchmark Seed: {seed} ---")
        for t in targets:
            tid = t["benchmark_id"]
            pname = t["protein_name"]
            seq_aa = t["sequence"]

            # 1. Arm 1: Unconstrained sLLM
            t0 = time.time()
            res1 = run_arm1_unconstrained_ai(
                seq_aa, "TAA", model_adapter, temperature=0.7, seed=seed
            )
            dur1 = time.time() - t0
            val1 = contract.final_validator.evaluate(
                sequence_dna=res1["sequence"], expected_protein_aa=seq_aa, expected_stop="TAA"
            )
            arm_sequences["arm1"][tid].append(res1["sequence"])

            # 2. Arm 2: Veto Only (No DP)
            t0 = time.time()
            res2 = run_arm2_veto_only(
                seq_aa, "TAA", contract, model_adapter, temperature=0.7, seed=seed
            )
            dur2 = time.time() - t0
            arm2_feasible = res2["status"] == "SUCCESS_VETO_ONLY"
            arm_sequences["arm2"][tid].append(res2["sequence"] or "")

            # 3. Arm 3A: Global DP Rescue (Job 285 Control)
            t0 = time.time()
            res3a = run_arm3a_global_dp_rescue(
                seq_aa, "TAA", contract, model_adapter, temperature=0.7, seed=seed
            )
            dur3a = time.time() - t0
            val3a = contract.final_validator.evaluate(
                sequence_dna=res3a["sequence"], expected_protein_aa=seq_aa, expected_stop="TAA"
            )
            arm_sequences["arm3a"][tid].append(res3a["sequence"])

            # 4. Arm 3B: Adaptive Partial DP Rescue (Job 286 Innovation)
            t0 = time.time()
            decoder_3b = ConstrainedSampleDecoder(
                contract=contract, model_adapter=model_adapter, temperature=0.7, seed=seed
            )
            res3b = decoder_3b.decode(seq_aa, stop_codon="TAA")
            dur3b = time.time() - t0
            val3b = contract.final_validator.evaluate(
                sequence_dna=res3b["sequence"], expected_protein_aa=seq_aa, expected_stop="TAA"
            )
            arm_sequences["arm3b"][tid].append(res3b["sequence"])

            records.append(
                {
                    "seed": seed,
                    "benchmark_id": tid,
                    "protein_name": pname,
                    "length_aa": len(seq_aa),
                    # Arm 1
                    "arm1_pass": val1.is_feasible,
                    "arm1_violations_count": len(val1.violations),
                    "arm1_violations": "; ".join(val1.violations),
                    "arm1_runtime_ms": round(dur1 * 1000, 2),
                    # Arm 2
                    "arm2_status": res2["status"],
                    "arm2_pass": arm2_feasible,
                    "arm2_deadend_pos": res2.get("deadend_pos", None),
                    "arm2_runtime_ms": round(dur2 * 1000, 2),
                    # Arm 3A
                    "arm3a_status": res3a["status"],
                    "arm3a_pass": val3a.is_feasible,
                    "arm3a_rescue_triggered": res3a.get("rescue_triggered", False),
                    "arm3a_ppr": res3a.get("prefix_preservation_rate", 0.0),
                    "arm3a_runtime_ms": round(dur3a * 1000, 2),
                    # Arm 3B
                    "arm3b_status": res3b["status"],
                    "arm3b_pass": val3b.is_feasible,
                    "arm3b_partial_rescue": res3b.get("partial_rescue", False),
                    "arm3b_rescue_triggered": res3b.get("rescue_triggered", False),
                    "arm3b_rewind_k": res3b.get("rewind_k", None),
                    "arm3b_ppr": res3b.get(
                        "prefix_preservation_rate", 1.0 if res3b["status"] == "SUCCESS" else 0.0
                    ),
                    "arm3b_preserved_prefix_codons": res3b.get(
                        "preserved_prefix_codons",
                        len(seq_aa) if res3b["status"] == "SUCCESS" else 0,
                    ),
                    "arm3b_runtime_ms": round(dur3b * 1000, 2),
                }
            )

    df = pd.DataFrame(records)
    csv_path = os.path.join(output_dir, "ablation_metrics.csv")
    df.to_csv(csv_path, index=False)

    # Compute Diversity (D_nt) per target and aggregate
    dnt_arm1 = [compute_pairwise_dnt(arm_sequences["arm1"][t["benchmark_id"]]) for t in targets]
    dnt_arm3a = [compute_pairwise_dnt(arm_sequences["arm3a"][t["benchmark_id"]]) for t in targets]
    dnt_arm3b = [compute_pairwise_dnt(arm_sequences["arm3b"][t["benchmark_id"]]) for t in targets]

    # Compute Aggregated Summary Metrics
    total_n = len(df)
    arm1_pass_rate = round(float(df["arm1_pass"].mean()) * 100, 2)
    arm2_pass_rate = round(float(df["arm2_pass"].mean()) * 100, 2)
    arm2_deadend_rate = round(float((df["arm2_status"] == "DEAD_END").mean()) * 100, 2)
    arm3a_pass_rate = round(float(df["arm3a_pass"].mean()) * 100, 2)
    arm3b_pass_rate = round(float(df["arm3b_pass"].mean()) * 100, 2)

    arm3a_global_fallback_count = int(df["arm3a_rescue_triggered"].sum())
    arm3a_global_fallback_rate = round(arm3a_global_fallback_count / total_n * 100, 2)
    arm3a_mean_ppr = round(float(df["arm3a_ppr"].mean()) * 100, 2)

    arm3b_pure_success_count = int((df["arm3b_status"] == "SUCCESS").sum())
    arm3b_partial_rescue_count = int((df["arm3b_status"] == "PARTIAL_DP_RESCUED").sum())
    arm3b_global_rescue_count = int((df["arm3b_status"] == "DP_RESCUED").sum())
    arm3b_global_fallback_rate = round(arm3b_global_rescue_count / total_n * 100, 2)
    arm3b_mean_ppr = round(float(df["arm3b_ppr"].mean()) * 100, 2)

    # Fallback reduction rate
    if arm3a_global_fallback_count > 0:
        fallback_reduction_rate = round(
            (1.0 - (arm3b_global_rescue_count / float(arm3a_global_fallback_count))) * 100, 2
        )
    else:
        fallback_reduction_rate = 0.0

    summary = {
        "benchmark_suite": "Paper 3 4-Arm Multi-Seed Ablation Suite",
        "total_trials": total_n,
        "num_targets": len(targets),
        "num_seeds": len(seeds),
        "seeds": list(seeds),
        "stage1_arm1_unconstrained_sllm": {
            "pass_rate_pct": arm1_pass_rate,
            "mean_dnt": round(float(np.mean(dnt_arm1)), 4),
            "mean_runtime_ms": round(float(df["arm1_runtime_ms"].mean()), 2),
        },
        "stage2_arm2_veto_only": {
            "pass_rate_pct": arm2_pass_rate,
            "deadend_rate_pct": arm2_deadend_rate,
            "mean_runtime_ms": round(float(df["arm2_runtime_ms"].mean()), 2),
        },
        "stage3a_arm3a_global_dp_rescue_control": {
            "pass_rate_pct": arm3a_pass_rate,
            "global_fallback_rate_pct": arm3a_global_fallback_rate,
            "mean_prefix_preservation_rate_pct": arm3a_mean_ppr,
            "mean_dnt": round(float(np.mean(dnt_arm3a)), 4),
            "mean_runtime_ms": round(float(df["arm3a_runtime_ms"].mean()), 2),
        },
        "stage3b_arm3b_adaptive_partial_dp_rescue": {
            "pass_rate_pct": arm3b_pass_rate,
            "pure_sllm_success_count": arm3b_pure_success_count,
            "pure_sllm_success_pct": round(arm3b_pure_success_count / total_n * 100, 2),
            "partial_dp_rescue_count": arm3b_partial_rescue_count,
            "partial_dp_rescue_pct": round(arm3b_partial_rescue_count / total_n * 100, 2),
            "global_dp_fallback_count": arm3b_global_rescue_count,
            "global_dp_fallback_pct": arm3b_global_fallback_rate,
            "global_fallback_reduction_rate_pct": fallback_reduction_rate,
            "mean_prefix_preservation_rate_pct": arm3b_mean_ppr,
            "mean_dnt": round(float(np.mean(dnt_arm3b)), 4),
            "mean_runtime_ms": round(float(df["arm3b_runtime_ms"].mean()), 2),
        },
        "scientific_invariants_verified": {
            "zero_defect_guarantee": bool(arm3a_pass_rate == 100.0 and arm3b_pass_rate == 100.0),
            "diversity_retention": bool(np.mean(dnt_arm3b) >= np.mean(dnt_arm3a)),
            "prefix_preservation_elevation": bool(arm3b_mean_ppr > arm3a_mean_ppr),
        },
    }

    json_path = os.path.join(output_dir, "ablation_metrics.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # Generate Figure 2 HTML Data Visualization
    waterfall_html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Paper 3 Figure 2: 4-Arm Resolution & Diversity Waterfall</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; padding: 24px; background: #f8fafc; color: #1e293b; }}
  .card {{ background: white; border-radius: 12px; padding: 28px; max-width: 820px; margin: auto; box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.1); border: 1px solid #e2e8f0; }}
  h2 {{ color: #0f172a; margin-top: 0; border-bottom: 2px solid #e2e8f0; padding-bottom: 10px; font-size: 20px; }}
  .subtitle {{ color: #64748b; font-size: 13px; margin-bottom: 20px; }}
  .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; margin-bottom: 24px; }}
  .metric-box {{ background: #f1f5f9; border-radius: 8px; padding: 16px; border-left: 4px solid #3b82f6; }}
  .metric-val {{ font-size: 24px; font-weight: 700; color: #0f172a; }}
  .metric-lbl {{ font-size: 12px; color: #64748b; font-weight: 600; text-transform: uppercase; }}
  .bar-container {{ margin: 14px 0; }}
  .bar-label {{ font-weight: 600; font-size: 13px; margin-bottom: 4px; display: flex; justify-content: space-between; }}
  .bar-track {{ height: 26px; background: #e2e8f0; border-radius: 6px; overflow: hidden; display: flex; }}
  .bar-fill {{ height: 100%; display: flex; align-items: center; justify-content: center; color: white; font-weight: 600; font-size: 12px; }}
  .arm1 {{ background: #ef4444; width: {arm1_pass_rate}%; }}
  .arm2 {{ background: #f59e0b; width: {arm2_pass_rate}%; }}
  .arm3a {{ background: #3b82f6; width: {arm3a_pass_rate}%; }}
  .arm3b {{ background: #10b981; width: {arm3b_pass_rate}%; }}
  .table {{ width: 100%; border-collapse: collapse; margin-top: 24px; }}
  .table th, .table td {{ border: 1px solid #cbd5e1; padding: 10px 12px; text-align: left; font-size: 13px; }}
  .table th {{ background: #f8fafc; font-weight: 600; }}
  .highlight {{ background: #ecfdf5; font-weight: 600; }}
</style>
</head>
<body>
<div class="card">
  <h2>Figure 2: Multi-Stage Resolution, Prefix Preservation & Diversity Waterfall</h2>
  <div class="subtitle">Evaluated on N={total_n} trials ({len(targets)} targets &times; {len(seeds)} seeds) under canonical <i>N. benthamiana</i> hard constraints</div>

  <div class="grid">
    <div class="metric-box" style="border-left-color: #10b981;">
      <div class="metric-lbl">Prefix Preservation Rate (PPR)</div>
      <div class="metric-val">{arm3b_mean_ppr}% <span style="font-size:14px; font-weight:normal; color:#64748b;">(vs 3A: {arm3a_mean_ppr}%)</span></div>
    </div>
    <div class="metric-box" style="border-left-color: #8b5cf6;">
      <div class="metric-lbl">Global Fallback Reduction</div>
      <div class="metric-val">{fallback_reduction_rate}% <span style="font-size:14px; font-weight:normal; color:#64748b;">({arm3a_global_fallback_rate}% &rarr; {arm3b_global_fallback_rate}%)</span></div>
    </div>
  </div>

  <div class="bar-container">
    <div class="bar-label"><span>Stage 1: Raw Unconstrained sLLM (Arm 1)</span> <span>{arm1_pass_rate}% Compliance</span></div>
    <div class="bar-track"><div class="bar-fill arm1">{arm1_pass_rate}%</div></div>
  </div>

  <div class="bar-container">
    <div class="bar-label"><span>Stage 2: Neuro-Symbolic Veto Only (Arm 2)</span> <span>{arm2_pass_rate}% Compliance ({arm2_deadend_rate}% Dead-end)</span></div>
    <div class="bar-track"><div class="bar-fill arm2">{arm2_pass_rate}%</div></div>
  </div>

  <div class="bar-container">
    <div class="bar-label"><span>Stage 3A: Global DP Rescue Baseline (Arm 3A)</span> <span>{arm3a_pass_rate}% Compliance (PPR {arm3a_mean_ppr}%, D_nt {summary["stage3a_arm3a_global_dp_rescue_control"]["mean_dnt"]})</span></div>
    <div class="bar-track"><div class="bar-fill arm3a">{arm3a_pass_rate}%</div></div>
  </div>

  <div class="bar-container">
    <div class="bar-label"><span>Stage 3B: Adaptive Partial DP Rescue Innovation (Arm 3B)</span> <span>{arm3b_pass_rate}% Compliance (PPR {arm3b_mean_ppr}%, D_nt {summary["stage3b_arm3b_adaptive_partial_dp_rescue"]["mean_dnt"]})</span></div>
    <div class="bar-track"><div class="bar-fill arm3b">{arm3b_pass_rate}%</div></div>
  </div>

  <table class="table">
    <thead>
      <tr>
        <th>Experimental Arm</th>
        <th>Compliance</th>
        <th>Prefix Preservation (PPR)</th>
        <th>Nucleotide Diversity (D_nt)</th>
        <th>Latency (ms)</th>
        <th>Architectural Role</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td><strong>Stage 1 (Raw AI)</strong></td>
        <td>{arm1_pass_rate}%</td>
        <td>100.0%</td>
        <td>{summary["stage1_arm1_unconstrained_sllm"]["mean_dnt"]}</td>
        <td>{summary["stage1_arm1_unconstrained_sllm"]["mean_runtime_ms"]} ms</td>
        <td>Baseline unconstrained generation defect rate</td>
      </tr>
      <tr>
        <td><strong>Stage 2 (Veto Only)</strong></td>
        <td>{arm2_pass_rate}%</td>
        <td>N/A (Dead-ends)</td>
        <td>N/A</td>
        <td>{summary["stage2_arm2_veto_only"]["mean_runtime_ms"]} ms</td>
        <td>Neuro-symbolic dead-end barrier</td>
      </tr>
      <tr>
        <td><strong>Stage 3A (Global DP)</strong></td>
        <td>{arm3a_pass_rate}%</td>
        <td>{arm3a_mean_ppr}%</td>
        <td>{summary["stage3a_arm3a_global_dp_rescue_control"]["mean_dnt"]}</td>
        <td>{summary["stage3a_arm3a_global_dp_rescue_control"]["mean_runtime_ms"]} ms</td>
        <td>Job 285 Baseline (Full DP fallback on dead-end)</td>
      </tr>
      <tr class="highlight">
        <td><strong>Stage 3B (Partial DP)</strong></td>
        <td><strong>{arm3b_pass_rate}%</strong></td>
        <td><strong>{arm3b_mean_ppr}%</strong></td>
        <td><strong>{summary["stage3b_arm3b_adaptive_partial_dp_rescue"]["mean_dnt"]}</strong></td>
        <td>{summary["stage3b_arm3b_adaptive_partial_dp_rescue"]["mean_runtime_ms"]} ms</td>
        <td><strong>Job 286 Innovation: Retains AI Context & Diversity</strong></td>
      </tr>
    </tbody>
  </table>
</div>
</body>
</html>
"""
    fig2_path = os.path.join(output_dir, "figure2_resolution_waterfall.html")
    with open(fig2_path, "w", encoding="utf-8") as f:
        f.write(waterfall_html)

    print("\n==================================================")
    print("Paper 3 Multi-Seed Ablation Suite Finished!")
    print(f"Total Trials: {total_n} ({len(targets)} targets x {len(seeds)} seeds)")
    print(f"Arm 1 Pass: {arm1_pass_rate}% | Arm 2 Pass: {arm2_pass_rate}%")
    print(
        f"Arm 3A (Global DP): {arm3a_pass_rate}% (PPR {arm3a_mean_ppr}%, D_nt {summary['stage3a_arm3a_global_dp_rescue_control']['mean_dnt']})"
    )
    print(
        f"Arm 3B (Partial DP): {arm3b_pass_rate}% (PPR {arm3b_mean_ppr}%, D_nt {summary['stage3b_arm3b_adaptive_partial_dp_rescue']['mean_dnt']})"
    )
    print(f"Global Fallback Reduction Rate: {fallback_reduction_rate}%")
    print(f"Results saved to: {output_dir}")
    print("==================================================")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output_dir",
        default=r"c:\Work\eijex\factorforge\benchmarks\results\paper3_ablation",
    )
    args = parser.parse_args()
    run_paper3_ablation_benchmark(args.output_dir)
