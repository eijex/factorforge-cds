"""Herceptin (Trastuzumab) 3-Way Benchmark & Delta Analysis Script.

This script performs ground-truth calibration between:
1. Public Herceptin AA Sequence
2. FactorForge Predicted Optimal CDS (Rule-based & ML-based)
3. Platform Lab Verified CDS (Ground Truth from wet-lab)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Any, List, Tuple

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from factorforge.engines.profile.optimizer import RuleBasedOptimizer  # noqa: E402

DEFAULT_REF = ROOT / "examples" / "herceptin_public_aa.fasta"


def parse_fasta(fasta_path: Path) -> List[Tuple[str, str]]:
    """Parse FASTA file into list of (header, sequence)."""
    if not fasta_path.exists():
        return []
    
    records = []
    current_header = ""
    current_seq = []
    
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    records.append((current_header, "".join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            records.append((current_header, "".join(current_seq)))
            
    return records


def compute_nucleotide_identity(seq1: str, seq2: str) -> float:
    """Compute exact nucleotide match percentage."""
    seq1, seq2 = seq1.upper(), seq2.upper()
    min_len = min(len(seq1), len(seq2))
    if min_len == 0:
        return 0.0
    matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)
    return round((matches / max(len(seq1), len(seq2))) * 100.0, 2)


def compute_codon_match_rate(seq1: str, seq2: str) -> float:
    """Compute codon-by-codon match percentage."""
    seq1, seq2 = seq1.upper(), seq2.upper()
    codons1 = [seq1[i:i+3] for i in range(0, len(seq1) - len(seq1)%3, 3)]
    codons2 = [seq2[i:i+3] for i in range(0, len(seq2) - len(seq2)%3, 3)]
    
    min_len = min(len(codons1), len(codons2))
    if min_len == 0:
        return 0.0
    matches = sum(1 for c1, c2 in zip(codons1[:min_len], codons2[:min_len]) if c1 == c2)
    return round((matches / max(len(codons1), len(codons2))) * 100.0, 2)


def run_herceptin_benchmark(
    ref_fasta: Path, platform_cds_fasta: Path | None = None, host: str = "nbenthamiana"
) -> Dict[str, Any]:
    """Run Herceptin optimization and delta analysis."""
    records = parse_fasta(ref_fasta)
    if not records:
        raise FileNotFoundError(f"Reference FASTA not found at {ref_fasta}")
        
    optimizer = RuleBasedOptimizer()
    benchmark_report = {
        "host_organism": host,
        "results": []
    }
    
    platform_records = parse_fasta(platform_cds_fasta) if platform_cds_fasta else []
    platform_dict = {name.split()[0]: seq for name, seq in platform_records}
    
    print("=================================================================")
    print("          FACTORFORGE HERCEPTIN BENCHMARK & DELTA REPORT         ")
    print("=================================================================\n")
    print(f"Host Organism: {host}")
    print(f"Reference FASTA: {ref_fasta.name}\n")
    
    for header, aa_seq in records:
        header_id = header.split("|")[0].split()[0]
        print(f"--> Processing Target: {header_id} ({len(aa_seq)} AAs)")
        
        # 1. Run FactorForge Rule-based Optimization
        res = optimizer.optimize(
            sequence=aa_seq,
            profile="balanced",
            host=host,
            scan_mode="full",
            seed=350,
        )
        
        # 2. Run FactorForge ML-based Constrained Optimization (LM Engine)
        from factorforge.engines.lm.inference import ConstrainedBeamSearchEngine
        lm_engine = ConstrainedBeamSearchEngine()
        lm_res = lm_engine.optimize_cds(
            amino_acids=aa_seq,
            host=host,
            type2is_clean=True
        )
        
        target_res = {
            "target_id": header_id,
            "aa_length": len(aa_seq),
            "factorforge_rule_cds": res.sequence,
            "factorforge_rule_cai": res.metrics.get("cai", 0.0),
            "factorforge_rule_gc": res.metrics.get("gc_percent", 0.0),
            "factorforge_lm_cds": lm_res["optimized_sequence"],
            "factorforge_lm_gc": lm_res["gc_percent"],
            "factorforge_lm_type2is_clean": lm_res["type2is_clean"],
            "scan_violations": res.metrics.get("violations", 0),
        }
        
        # Compute Rule vs LM agreement
        rule_lm_identity = compute_nucleotide_identity(res.sequence, lm_res["optimized_sequence"])
        rule_lm_codon = compute_codon_match_rate(res.sequence, lm_res["optimized_sequence"])
        
        print("   [FACTORFORGE ENGINE PREDICTIONS]")
        print(f"   Rule-based Engine CDS          : CAI={res.metrics.get('cai', 0.0):.3f}, GC={res.metrics.get('gc_percent', 0.0):.2f}%")
        print(f"   ML Preview Scaffold CDS         : GC={lm_res['gc_percent']:.2f}%, TypeIIS Clean={lm_res['type2is_clean']}")
        print(f"   Rule vs LM Concordance         : {rule_lm_identity}% NT identity ({rule_lm_codon}% codon match)")
        
        # 3. Compare against Platform CDS (if provided)
        if header_id in platform_dict:
            plat_cds = platform_dict[header_id]
            nt_identity_rule = compute_nucleotide_identity(res.sequence, plat_cds)
            nt_identity_lm = compute_nucleotide_identity(lm_res["optimized_sequence"], plat_cds)
            
            codon_match_rule = compute_codon_match_rate(res.sequence, plat_cds)
            codon_match_lm = compute_codon_match_rate(lm_res["optimized_sequence"], plat_cds)
            
            target_res["platform_cds_length"] = len(plat_cds)
            target_res["rule_vs_platform_nt_identity"] = nt_identity_rule
            target_res["lm_vs_platform_nt_identity"] = nt_identity_lm
            
            print("   --------------------------------------------------------------")
            print("   [GROUND TRUTH COMPARISON WITH PLATFORM CDS]")
            print(f"   Platform CDS Length            : {len(plat_cds)} bp")
            print(f"   Rule-based vs Platform Identity : {nt_identity_rule}% (Codon Match: {codon_match_rule}%)")
            print(f"   LM-based vs Platform Identity   : {nt_identity_lm}% (Codon Match: {codon_match_lm}%)")
            print("   --------------------------------------------------------------")
        else:
            print("   [STATUS] Awaiting Platform CDS for ground truth comparison.")
            
        print("")
        benchmark_report["results"].append(target_res)
        
    return benchmark_report


def main() -> int:
    parser = argparse.ArgumentParser(description="Herceptin 3-Way Benchmark Calibration Script")
    parser.add_argument("--ref", type=Path, default=DEFAULT_REF, help="Public Herceptin AA FASTA")
    parser.add_argument("--platform-cds", type=Path, default=None, help="Platform verified CDS FASTA")
    parser.add_argument("--host", type=str, default="nbenthamiana", help="Host organism")
    
    args = parser.parse_args()
    run_herceptin_benchmark(args.ref, args.platform_cds, args.host)
    return 0


if __name__ == "__main__":
    sys.exit(main())
