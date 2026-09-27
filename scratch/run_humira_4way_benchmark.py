import sys, os, itertools, json
sys.path.insert(0, r'C:\Work\eijex\factorforge\src')
sys.path.insert(0, r'C:\Work\eijex\factorforge')

from factorforge.engines.profile.optimizer import RuleBasedOptimizer
from factorforge.engines.dp_adapter import DPEngineAdapter
from factorforge.engines.lm.inference import ONNXBeamSearchEngine
from factorforge.evaluation.evaluator import SharedEvaluator
from factorforge.engines.profile.utils import load_golden_set
from factorforge.engines.profile.rules.reverse_translator import ReverseTranslator

class V3CodonTokenizer:
    def __init__(self):
        self.special_tokens = ['[PAD]', '[BOS]', '[EOS]', '[MASK]', '[UNK]']
        bases = ['T', 'C', 'A', 'G']
        self.codons = [''.join(p) for p in itertools.product(bases, repeat=3)]
        self.vocab = self.special_tokens + self.codons
        self.token_to_id = {t: i for i, t in enumerate(self.vocab)}
        self.id_to_token = {i: t for i, t in enumerate(self.vocab)}
    @property
    def vocab_size(self): return len(self.vocab)
    @property
    def pad_id(self): return 0
    @property
    def bos_id(self): return 1
    @property
    def eos_id(self): return 2
    @property
    def unk_id(self): return 4

    def encode_encoder_input(self, aa_seq, host='nbenthamiana', gc_band='40-47', type2is_clean=True):
        return [self.bos_id] + [self.unk_id for _ in aa_seq] + [self.eos_id]

# Load targets
lc_aa = "MAKTNLFLFLIFSLLLSLSSADIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
hc_aa = "MAKTNLFLFLIFSLLLSLSSAEVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"

# Candidate A: PlantForm Baseline DNA
lc_baseline_dna = "ATGGCAAAAACCAATCTTTTTCTTTTTCTCATTTTCTCTTTGCTTTTGAGTTTATCTTCTGCTGACATTCAAATGACTCAAAGCCCTTCTTCTCTTTCTGCTTCTGTTGGAGACAGAGTCACCATCACTTGCAGAGCTTCTCAAGGAATCAGAAATTACCTTGCTTGGTATCAACAGAAGCCTGGAAAAGCTCCCAAACTTTTGATTTATGCTGCTTCCACTCTTCAATCTGGAGTTCCATCCAGATTTTCTGGATCTGGATCTGGAACTGATTTCACTTTGACCATTTCTTCTTTGCAACCTGAAGATGTTGCCACTTACTATTGTCAAAGATACAACAGAGCACCATACACTTTTGGACAGGGAACCAAAGTTGAAATCAAGAGAACTGTTGCTGCTCCTTCTGTTTTCATTTTTCCTCCTTCTGATGAACAATTGAAATCTGGAACTGCCTCTGTTGTTTGTCTTTTGAACAATTTTTATCCCAGAGAAGCCAAAGTTCAATGGAAGGTTGACAATGCTCTTCAATCTGGAAACAGCCAGGAATCTGTCACTGAACAAGATTCCAAAGATTCCACTTATTCTCTCTCTTCCACTTTGACACTTTCCAAGGCTGATTATGAGAAACACAAAGTTTATGCATGTGAAGTCACTCATCAAGGACTTTCTTCTCCTGTCACCAAATCATTCAACAGAGGAGAATGT"
hc_baseline_dna = "ATGGCCAAAACCAATCTTTTTCTTTTTTTGATTTTTTCTCTTCTTCTTTCTCTTTCTTCTGCTGAAGTGCAATTGGTTGAGTCTGGAGGAGGACTGGTTCAACCTGGAAGATCTTTGAGACTTTCTTGTGCTGCTTCTGGATTCACTTTTGATGATTATGCCATGCATTGGGTCAGGCAGGCTCCTGGAAAGGGACTTGAATGGGTTTCTGCCATCACCTGGAATTCAGGGCACATTGATTATGCTGACTCTGTTGAAGGAAGATTCACCATTTCCAGAGACAATGCCAAAAATTCTCTTTATCTTCAGATGAATTCTTTGAGAGCTGAGGACACTGCTGTTTATTATTGTGCCAAAGTTTCTTATCTTTCCACTGCCAGCTCTTTGGATTATTGGGGACAAGGAACTCTTGTCACTGTTTCTTCTGCTTCCACAAAAGGGCCTTCTGTTTTTCCTCTTGCTCCTTCATCCAAATCCACTTCTGGTGGAACTGCTGCTCTTGGATGTCTTGTCAAGGATTATTTTCCTGAACCTGTCACTGTGTCTTGGAATTCTGGAGCTTTGACTTCTGGAGTTCACACATTTCCTGCAGTGCTTCAGTCCTCTGGCCTTTACAGTCTTTCTTCTGTTGTCACTGTTCCATCTTCTTCTTTGGGAACTCAAACTTACATTTGCAATGTCAATCACAAACCATCCAACACCAAAGTTGACAAAAAAGTTGAACCCAAGTCTTGTGACAAGACTCACACTTGTCCTCCTTGTCCTGCACCTGAACTTCTTGGAGGACCTTCTGTTTTTCTTTTCCCTCCCAAGCCAAAAGACACTTTGATGATTTCCAGAACTCCTGAAGTCACTTGTGTTGTTGTTGATGTTTCTCATGAAGATCCTGAAGTCAAATTCAACTGGTATGTTGATGGAGTGGAAGTTCACAATGCCAAAACCAAACCCAGAGAAGAACAATACAATTCCACATACAGGGTTGTTTCTGTTTTGACTGTTCTTCATCAGGATTGGTTGAATGGAAAAGAGTACAAATGCAAAGTTTCCAACAAAGCTCTTCCTGCTCCCATTGAAAAAACCATTTCCAAAGCCAAAGGACAACCCAGAGAGCCTCAAGTCTACACTCTTCCTCCTTCCAGAGATGAATTGACCAAGAATCAGGTTTCTTTGACCTGTCTTGTCAAAGGATTTTATCCTTCTGACATTGCTGTTGAGTGGGAATCCAATGGACAACCTGAAAACAATTACAAAACCACTCCCCCTGTTCTTGACAGTGATGGATCTTTTTTTCTTTATTCCAAATTGACTGTTGACAAATCCAGATGGCAACAAGGCAATGTTTTTTCTTGCAGTGTCATGCATGAGGCTCTTCACAATCATTACACCCAAAAATCTCTTTCTTTGTCTCCAGGAAAA"

# Initialize Evaluator
table = load_golden_set()
codon_weights = ReverseTranslator._build_ref_weights(table)
evaluator = SharedEvaluator(version="1.0.1", codon_weights=codon_weights)

# Initialize Engines
rule_engine = RuleBasedOptimizer()
dp_engine = DPEngineAdapter()
onnx_path = r'C:\Work\eijex\factorforge\models\factorforge_v3_run3_step15000.onnx'
slm_engine = ONNXBeamSearchEngine(onnx_model_path=onnx_path, beam_width=3)
slm_engine.tokenizer = V3CodonTokenizer()

targets = [
    ("Humira_Light_Chain", lc_aa, lc_baseline_dna),
    ("Humira_Heavy_Chain", hc_aa, hc_baseline_dna)
]

print("==========================================================================================")
print("PLANTFORM TARGET-MAB-A (HUMIRA) 1:1 SIDE-BY-SIDE 4-WAY COMPARISON MATRIX")
print("==========================================================================================\n")

for name, aa, base_dna in targets:
    print(f"--- Target: {name} ({len(aa)} AA / {len(aa)*3} bp) ---")
    
    # Candidate A: PlantForm Baseline
    eval_a = evaluator.evaluate_candidate(base_dna, aa, "PlantForm_Baseline", forbidden_type_iis={"BsaI", "BsmBI", "BpiI"})
    
    # Candidate B: Rule-based Profile (High CAI)
    res_b = rule_engine.optimize(aa, profile="high_cai", host="nbenthamiana")
    eval_b = evaluator.evaluate_candidate(res_b.sequence, aa, "FactorForge_Rule_HighCAI", forbidden_type_iis={"BsaI", "BsmBI", "BpiI"})
    
    # Candidate C: Math DP Engine
    res_c = dp_engine.optimize(aa, host="nbenthamiana")
    eval_c = evaluator.evaluate_candidate(res_c.sequence, aa, "FactorForge_Math_DP", forbidden_type_iis={"BsaI", "BsmBI", "BpiI"})
    
    # Candidate D: Real ONNX sLLM Hybrid
    res_d = slm_engine.optimize_cds(aa, host="nbenthamiana", type2is_clean=True)
    eval_d = evaluator.evaluate_candidate(res_d["optimized_sequence"], aa, "FactorForge_sLLM_Hybrid", forbidden_type_iis={"BsaI", "BsmBI", "BpiI"})

    def summarize(e, seq):
        bsai_hits = sum(1 for c in e.checks if "type_iis" in c.check_name and c.result.value == "fail")
        polya_hits = sum(1 for c in e.checks if "polya" in c.check_name and c.result.value == "warning")
        return {
            "aa_ident": e.sequence_integrity.aa_identity,
            "cai": e.metrics.cai,
            "gc": e.metrics.gc_percent,
            "bsai_hits": bsai_hits,
            "polya_warnings": polya_hits,
            "diff_from_base": sum(1 for x, y in zip(base_dna, seq) if x != y)
        }

    sum_a = summarize(eval_a, base_dna)
    sum_b = summarize(eval_b, res_b.sequence)
    sum_c = summarize(eval_c, res_c.sequence)
    sum_d = summarize(eval_d, res_d["optimized_sequence"])

    b_diff_str = str(sum_b['diff_from_base']) + " bp (" + str(round(sum_b['diff_from_base']/len(base_dna)*100, 1)) + "%)"
    c_diff_str = str(sum_c['diff_from_base']) + " bp (" + str(round(sum_c['diff_from_base']/len(base_dna)*100, 1)) + "%)"
    d_diff_str = str(sum_d['diff_from_base']) + " bp (" + str(round(sum_d['diff_from_base']/len(base_dna)*100, 1)) + "%)"

    print(f"{'Metric':<25} | {'A. PlantForm Base':<18} | {'B. Rule (High CAI)':<18} | {'C. Math DP':<18} | {'D. sLLM Hybrid':<18}")
    print("-" * 105)
    print(f"{'AA Identity':<25} | {sum_a['aa_ident']:<18.3f} | {sum_b['aa_ident']:<18.3f} | {sum_c['aa_ident']:<18.3f} | {sum_d['aa_ident']:<18.3f}")
    print(f"{'CAI (Host Adaptation)':<25} | {sum_a['cai']:<18.3f} | {sum_b['cai']:<18.3f} | {sum_c['cai']:<18.3f} | {sum_d['cai']:<18.3f}")
    print(f"{'GC Content (%)':<25} | {sum_a['gc']:<17.1f}% | {sum_b['gc']:<17.1f}% | {sum_c['gc']:<17.1f}% | {sum_d['gc']:<17.1f}%")
    print(f"{'Type IIS Sites (BsaI)':<25} | {sum_a['bsai_hits']:<18} | {sum_b['bsai_hits']:<18} | {sum_c['bsai_hits']:<18} | {sum_d['bsai_hits']:<18}")
    print(f"{'PolyA-like Warnings':<25} | {sum_a['polya_warnings']:<18} | {sum_b['polya_warnings']:<18} | {sum_c['polya_warnings']:<18} | {sum_d['polya_warnings']:<18}")
    print(f"{'Base Pair Differences':<25} | {'0 (Baseline)':<18} | {b_diff_str:<18} | {c_diff_str:<18} | {d_diff_str:<18}\n")
