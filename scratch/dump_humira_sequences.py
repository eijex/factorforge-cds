import sys, os, hashlib
sys.path.insert(0, r'C:\Work\eijex\factorforge\src')
sys.path.insert(0, r'C:\Work\eijex\factorforge')

from factorforge.engines.dp_adapter import DPEngineAdapter

lc_aa = "MAKTNLFLFLIFSLLLSLSSADIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
hc_aa = "MAKTNLFLFLIFSLLLSLSSAEVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"

dp = DPEngineAdapter()
res_lc = dp.optimize(lc_aa, host="nbenthamiana")
res_hc = dp.optimize(hc_aa, host="nbenthamiana")

lc_dp = res_lc.sequence
hc_dp = res_hc.sequence

out_dir = r"C:\Work\eijex\eijex-workspace\_poc\humira\2026-09-10"
os.makedirs(out_dir, exist_ok=True)

fasta_path = os.path.join(out_dir, "FF_HUMIRA_DP_v2.fasta")
with open(fasta_path, "w", encoding="utf-8") as f:
    f.write(f">FF-HUMIRA-LC-DP2 | Light Chain | 235 AA / {len(lc_dp)} bp | CAI={res_lc.metrics['cai']:.3f} | GC={res_lc.metrics['gc_percent']:.1f}% | BsaI=0\n")
    f.write(lc_dp + "\n\n")
    f.write(f">FF-HUMIRA-HC-DP2 | Heavy Chain | 472 AA / {len(hc_dp)} bp | CAI={res_hc.metrics['cai']:.3f} | GC={res_hc.metrics['gc_percent']:.1f}% | BsaI=0\n")
    f.write(hc_dp + "\n")

print("DUMP_SUCCESS")
print("LC_SHA256:", hashlib.sha256(lc_dp.encode()).hexdigest())
print("HC_SHA256:", hashlib.sha256(hc_dp.encode()).hexdigest())
print("LC_LEN:", len(lc_dp))
print("HC_LEN:", len(hc_dp))
print("LC_SEQ:\n" + lc_dp)
print("HC_SEQ:\n" + hc_dp)
