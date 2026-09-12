import sys, os, hashlib, json
sys.path.insert(0, r'C:\Work\eijex\factorforge\src')
sys.path.insert(0, r'C:\Work\eijex\factorforge')

from factorforge.engines.dp_v2 import DPV2Optimizer
from factorforge.engines.dp_v2_1 import DPV21Optimizer
from factorforge.engines.profile.utils import load_golden_set
from factorforge.engines.profile.rules.reverse_translator import ReverseTranslator
from factorforge.evaluation.evaluator import SharedEvaluator

# PlantForm Humira Constructs (Signal Peptide + Mature Chain)
lc_aa = "MAKTNLFLFLIFSLLLSLSSADIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
hc_aa = "MAKTNLFLFLIFSLLLSLSSAEVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"

table = load_golden_set()
codon_weights = ReverseTranslator._build_ref_weights(table)
evaluator = SharedEvaluator(codon_weights=codon_weights)

# Construct Flanks (PlantForm Vector compatibility)
up_flank = "GGTC"
stop_codon = "TAA"
down_flank = "GAGA"

# --- 1. DP v2.0 Optimization ---
opt_v20 = DPV2Optimizer()
res_lc_v20 = opt_v20.optimize(lc_aa, codon_weights, target_gc_min=0.40, target_gc_max=0.47, left_flank=up_flank, right_flank=stop_codon+down_flank)
res_hc_v20 = opt_v20.optimize(hc_aa, codon_weights, target_gc_min=0.40, target_gc_max=0.47, left_flank=up_flank, right_flank=stop_codon+down_flank)

eval_lc_v20 = evaluator.evaluate_candidate(res_lc_v20["sequence"], lc_aa, upstream_context=up_flank, stop_codon=stop_codon, downstream_context=down_flank)
eval_hc_v20 = evaluator.evaluate_candidate(res_hc_v20["sequence"], hc_aa, upstream_context=up_flank, stop_codon=stop_codon, downstream_context=down_flank)

# --- 2. DP v2.1 Optimization ---
opt_v21 = DPV21Optimizer(ramp_length_codons=15)
res_lc_v21 = opt_v21.optimize(lc_aa, codon_weights, target_gc_min=0.40, target_gc_max=0.47, upstream_context=up_flank, stop_codon=stop_codon, downstream_context=down_flank)
res_hc_v21 = opt_v21.optimize(hc_aa, codon_weights, target_gc_min=0.40, target_gc_max=0.47, upstream_context=up_flank, stop_codon=stop_codon, downstream_context=down_flank)

eval_lc_v21 = evaluator.evaluate_candidate(res_lc_v21["sequence"], lc_aa, upstream_context=up_flank, stop_codon=stop_codon, downstream_context=down_flank)
eval_hc_v21 = evaluator.evaluate_candidate(res_hc_v21["sequence"], hc_aa, upstream_context=up_flank, stop_codon=stop_codon, downstream_context=down_flank)

out_dir = r"C:\Work\eijex\eijex-workspace\_poc\humira\2026-09-11"
os.makedirs(out_dir, exist_ok=True)

# Dump v2.1 FASTA
fasta_path = os.path.join(out_dir, "FF_HUMIRA_DP_v2_1_PLANTFORM.fasta")
with open(fasta_path, "w", encoding="utf-8") as f:
    f.write(f">FF-HUMIRA-LC-DP2.1 | PlantForm Light Chain (SP+Mature) | 235 AA / {len(res_lc_v21['sequence'])} bp | Global_CAI={res_lc_v21['cai']:.3f} | Body_CAI={res_lc_v21['cai_body']:.3f} | 5p_GC={res_lc_v21['gc_5p_ramp_percent']:.1f}% | Global_GC={res_lc_v21['gc_percent']:.1f}% | BsaI=0\n")
    f.write(res_lc_v21["sequence"] + "\n\n")
    f.write(f">FF-HUMIRA-HC-DP2.1 | PlantForm Heavy Chain (SP+Mature) | 472 AA / {len(res_hc_v21['sequence'])} bp | Global_CAI={res_hc_v21['cai']:.3f} | Body_CAI={res_hc_v21['cai_body']:.3f} | 5p_GC={res_hc_v21['gc_5p_ramp_percent']:.1f}% | Global_GC={res_hc_v21['gc_percent']:.1f}% | BsaI=0\n")
    f.write(res_hc_v21["sequence"] + "\n")

print("=== HUMIRA COMPARATIVE BENCHMARK: DP v2.0 vs DP v2.1 ===")
print("\n--- Light Chain (LC: 235 AA, 705 bp CDS) ---")
print(f"DP v2.0: Global CAI={res_lc_v20['cai']:.3f}, GC%={res_lc_v20['gc_percent']:.1f}%, 5' Ramp GC%={eval_lc_v20.metrics.gc_5p_ramp_percent:.1f}%, 5' MFE={eval_lc_v20.metrics.mfe_5p_window} kcal/mol, SHA-256={hashlib.sha256(res_lc_v20['sequence'].encode()).hexdigest()[:16]}")
print(f"DP v2.1: Global CAI={res_lc_v21['cai']:.3f}, Body CAI={res_lc_v21['cai_body']:.3f}, GC%={res_lc_v21['gc_percent']:.1f}%, 5' Ramp GC%={res_lc_v21['gc_5p_ramp_percent']:.1f}%, 5' MFE={eval_lc_v21.metrics.mfe_5p_window} kcal/mol, SHA-256={hashlib.sha256(res_lc_v21['sequence'].encode()).hexdigest()[:16]}")

print("\n--- Heavy Chain (HC: 472 AA, 1,416 bp CDS) ---")
print(f"DP v2.0: Global CAI={res_hc_v20['cai']:.3f}, GC%={res_hc_v20['gc_percent']:.1f}%, 5' Ramp GC%={eval_hc_v20.metrics.gc_5p_ramp_percent:.1f}%, 5' MFE={eval_hc_v20.metrics.mfe_5p_window} kcal/mol, SHA-256={hashlib.sha256(res_hc_v20['sequence'].encode()).hexdigest()[:16]}")
print(f"DP v2.1: Global CAI={res_hc_v21['cai']:.3f}, Body CAI={res_hc_v21['cai_body']:.3f}, GC%={res_hc_v21['gc_percent']:.1f}%, 5' Ramp GC%={res_hc_v21['gc_5p_ramp_percent']:.1f}%, 5' MFE={eval_hc_v21.metrics.mfe_5p_window} kcal/mol, SHA-256={hashlib.sha256(res_hc_v21['sequence'].encode()).hexdigest()[:16]}")

print(f"\nSaved v2.1 PlantForm Humira FASTA: {fasta_path}")
