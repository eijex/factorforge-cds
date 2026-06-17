#!/bin/bash
# Job 130 reopen — manifest/file SHA-256 mismatch fix, B/C1/C2 benchmark rerun.
# Run this from the factorforge repo root. Each profile takes ~35 min (runtime_seconds
# was 2093 for the original run) — run in background / expect ~1h45m total for all three.
#
# Context: src/factorforge/data/profiles/{qld183_v103_derived,nbev11_cds_all_derived,
# nbev11_cds_hc_derived}_manifest.json had a stale `codon_profile_sha256` that did not
# match the committed codon table JSON (confirmed: the JSON files are deterministically
# reproducible from source FASTA via scripts/build_codon_profile.py; only the manifest's
# recorded hash was wrong, likely a STEP 6 rebuild-without-remanifest bug). This was
# already fixed in this session (manifest hash field corrected to match the actual
# committed file) but NOT yet committed — uncommitted manifest fix should still be in
# the working tree when this resumes. Verify with:
#   git status --short src/factorforge/data/profiles/
# If the fix is gone (fresh clone), recompute live before rerunning:
#   python -c "
#   import json, hashlib
#   for j, m in [
#       ('src/factorforge/data/profiles/qld183_v103_derived_codons.json', 'src/factorforge/data/profiles/qld183_v103_derived_manifest.json'),
#       ('src/factorforge/data/profiles/nbev11_cds_all_derived_codons.json', 'src/factorforge/data/profiles/nbev11_cds_all_derived_manifest.json'),
#       ('src/factorforge/data/profiles/nbev11_cds_hc_derived_codons.json', 'src/factorforge/data/profiles/nbev11_cds_hc_derived_manifest.json'),
#   ]:
#       d = json.load(open(m, encoding='utf-8'))
#       d['codon_profile_sha256'] = hashlib.sha256(open(j, 'rb').read()).hexdigest()
#       json.dump(d, open(m, 'w', encoding='utf-8'), indent=4)
#   "
#
# Before rerunning, confirm the hard-error path is now clear for all three:
set -e
cd "$(dirname "$0")/../.."

# `bash script.sh` launched from PowerShell's "(base)" conda env does NOT inherit
# that env's PATH (conda activate only modifies the parent PowerShell session, not
# this bash subshell) -- and even if it did, that conda base env is missing project
# deps (yaml, pandas). The actual interpreter with factorforge's deps installed is
# the hermes-agent venv. Pin to it explicitly so this is robust regardless of which
# shell/session invokes this script.
PYTHON="/c/Users/munky/AppData/Local/hermes/hermes-agent/venv/Scripts/python.exe"
if [ ! -x "$PYTHON" ]; then
  echo "WARNING: hermes venv python not found at $PYTHON -- falling back to PATH python" >&2
  PYTHON="python"
fi
echo "Using: $PYTHON ($("$PYTHON" --version 2>&1))"

for p in qld183_v103_derived nbev11_cds_all_derived nbev11_cds_hc_derived; do
  echo "=== sanity check: $p ==="
  "$PYTHON" benchmarks/run_benchmark.py --dataset synthetic --limit 2 --seed 320 \
    --codon-table-path "src/factorforge/data/profiles/${p}_codons.json" \
    --source-profile-id "$p" \
    --source-profile-manifest "src/factorforge/data/profiles/${p}_manifest.json"
done

echo "=== sanity checks passed. Starting full reruns (B, C1, C2) ==="

"$PYTHON" benchmarks/scripts/run_and_archive_profile.py \
  --profile-id qld183_v103_derived \
  --codon-table-path src/factorforge/data/profiles/qld183_v103_derived_codons.json \
  --manifest-path src/factorforge/data/profiles/qld183_v103_derived_manifest.json \
  --seed 320

"$PYTHON" benchmarks/scripts/run_and_archive_profile.py \
  --profile-id nbev11_cds_all_derived \
  --codon-table-path src/factorforge/data/profiles/nbev11_cds_all_derived_codons.json \
  --manifest-path src/factorforge/data/profiles/nbev11_cds_all_derived_manifest.json \
  --seed 320

"$PYTHON" benchmarks/scripts/run_and_archive_profile.py \
  --profile-id nbev11_cds_hc_derived \
  --codon-table-path src/factorforge/data/profiles/nbev11_cds_hc_derived_codons.json \
  --manifest-path src/factorforge/data/profiles/nbev11_cds_hc_derived_manifest.json \
  --seed 320

echo "=== Reruns complete. Aggregating ==="
"$PYTHON" benchmarks/scripts/aggregate_profiles.py

echo "=== Regenerating SHA256SUMS.txt ==="
"$PYTHON" -c "
import hashlib, pathlib
root = pathlib.Path('.')
files = [
    'src/factorforge/data/nbenthamiana_codons.json',
    'src/factorforge/data/profiles/qld183_v103_derived_codons.json',
    'src/factorforge/data/profiles/nbev11_cds_all_derived_codons.json',
    'src/factorforge/data/profiles/nbev11_cds_hc_derived_codons.json',
    'src/factorforge/data/profiles/qld183_v103_derived_manifest.json',
    'src/factorforge/data/profiles/nbev11_cds_all_derived_manifest.json',
    'src/factorforge/data/profiles/nbev11_cds_hc_derived_manifest.json',
    'src/factorforge/data/profiles/qld183_v103_filtered_stats.json',
    'src/factorforge/data/profiles/nbev11_cds_all_derived_filtered_stats.json',
    'src/factorforge/data/profiles/nbev11_cds_hc_derived_filtered_stats.json',
    'reproducibility/benchmark_v3.2.2/benchmark_summary_profile_legacy_packaged.frozen.json',
    'reproducibility/benchmark_v3.2.2/benchmark_summary_profile_qld183_v103_derived.frozen.json',
    'reproducibility/benchmark_v3.2.2/benchmark_summary_profile_nbev11_cds_all_derived.frozen.json',
    'reproducibility/benchmark_v3.2.2/benchmark_summary_profile_nbev11_cds_hc_derived.frozen.json',
    'reproducibility/benchmark_v3.2.2/combined_profile_comparison.json',
    'reproducibility/benchmark_v3.2.2/combined_profile_comparison.tsv',
    'reproducibility/benchmark_v3.2.2/cai_shift_by_profile.tsv',
    'reproducibility/benchmark_v3.2.2/multi_constraint_pass_by_profile.tsv',
    'reproducibility/benchmark_v3.2.2/codon_weight_correlation_matrix.tsv',
    'reproducibility/benchmark_v3.2.2/preferred_codon_rank_shift.tsv',
    'reproducibility/benchmark_v0.5.1/data/profile_A_reproduction_diff.json',
    'reproducibility/benchmark_v0.5.1/data/profile_A_reproduction_diff.tsv',
    'scripts/build_codon_profile.py',
]
with open('reproducibility/benchmark_v3.2.2/SHA256SUMS.txt', 'w', encoding='utf-8', newline='\n') as out:
    for f in files:
        h = hashlib.sha256(open(f, 'rb').read()).hexdigest()
        out.write(f'{h}  {f}\n')
"

echo "=== DONE. Next steps (manual): ==="
echo "1. Diff new reproducibility/benchmark_v3.2.2/*.frozen.json against the previous"
echo "   committed versions (git diff) -- confirm whether multi_constraint_pass / GC /"
echo "   CAI numbers actually changed vs analysis 020's reported values, or only the"
echo "   hash fields changed."
echo "2. If numbers are unchanged: this confirms the original Job 130 STEP 9/10 findings"
echo "   and analysis 020's verdict stand as-is -- just commit the corrected manifests +"
echo "   re-verified frozen summaries + SHA256SUMS.txt."
echo "3. If numbers DID change: STOP, do not just re-commit. Analysis 020's verdict and"
echo "   the evidence-pack text (papers/evidence/factorforge_master_source_pack/) may"
echo "   need to be revised before Job 130 can close again. Report to user first."
echo "4. Update eijex-workspace Job 130 docs (_meta.yml back to Done, report STEP, "
echo "   130-improv-job.md AC6/AC10/AC13 re-confirmed) and push both repos."
