# Tutorial: Designing a GFP CDS for *N. benthamiana*

This tutorial walks through a complete synonymous CDS design workflow using Green Fluorescent Protein (GFP) as an example target, showing how to go from an amino acid sequence to a *Nicotiana benthamiana*-oriented CDS candidate.

## Prerequisites

```bash
pip install factorforge-cds
```

Verify installation:

```bash
factorforge --version
```

## Input Sequence

We use the *Aequorea victoria* GFP sequence (239 amino acids) as our target protein. Save it as `gfp.fasta`:

```
>GFP|Aequorea_victoria|239aa
MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWP
TLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDT
LVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQL
ADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK
```

## Step 1: Run Optimization (CLI)

Run the `balanced` profile — the recommended starting point for *N. benthamiana* CDS design:

```bash
factorforge optimize gfp.fasta --engine profile --profile balanced -o gfp_optimized.fasta
```

For exploratory Tobacco BY-2 sequence-context review:

```bash
factorforge optimize gfp.fasta --engine profile --profile balanced --host by2
```

Expected output:

```
Optimizing with Profile-based v3.4.1...
Saved to: gfp_optimized.fasta
Metrics:
  - cai: 0.769
  - gc_percent: 58.72
  - score: 0.856
```

The designed CDS begins with `ATGGTGAGCAAGGGCGAGGAA...` and can be reviewed for synthesis or downstream assembly. Reverse translation is stochastic, so exact codons and metrics vary slightly between runs.

**What happened:**

| Metric | Value | Target range |
|--------|-------|-------------|
| CAI | 0.769 | — (higher is better) |
| GC% | 58.72% | 55–65% |
| Internal stop codons | 0 | 0 (hard requirement) |
| Rare codon runs | 0 | 0 (ribosome stalling risk) |

## Step 2: Compare Profiles

Different design goals call for different profiles. Use `--compare-profiles` to evaluate all options at once:

```bash
factorforge optimize gfp.fasta --engine profile \
  --compare-profiles balanced,high_cai,gc_target,assembly_friendly \
  --scan-mode fast
```

Output:

```
Profile comparison results:
─────────────────────────────────────────────
Profile               CAI     GC%    Score
─────────────────────────────────────────────
balanced            0.770   57.74    0.856
high_cai            1.000   31.24    0.889
gc_target           0.774   59.83    0.972
assembly_friendly   0.779   57.18    0.905
─────────────────────────────────────────────
```

**How to choose:**

| Profile | Best for |
|---------|----------|
| `balanced` | General *N. benthamiana* CDS design review; good CAI with GC% in target range |
| `high_cai` | CAI-focused comparison; note GC% may fall outside the active host GC range |
| `gc_target` | When GC% must hit a specific value (defaults to the host midpoint of 43.5% for the current NbeV1.1 software default; pass `--target-gc` for other values, e.g. specific vector requirements) |
| `assembly_friendly` | MoClo / Golden Gate workflows; avoids problematic restriction sites |

For most *N. benthamiana* sequence-design tasks, `balanced` is the recommended starting profile.

## Step 3: Python API

The same optimization is available programmatically:

```python
from factorforge.engines import EngineRegistry

# Load the profile engine
optimizer = EngineRegistry.get("profile")

gfp_aa = (
    "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWP"
    "TLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDT"
    "LVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQL"
    "ADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
)

result = optimizer.optimize(gfp_aa, profile="balanced")

print(f"Optimized CDS: {result.sequence[:60]}...")
print(f"CAI:   {result.metrics['cai']:.3f}")
print(f"GC%:   {result.metrics['gc_percent']:.2f}")
print(f"Score: {result.metrics['score']:.3f}")
```

To compare profiles programmatically:

```python
profiles = ["balanced", "high_cai", "gc_target", "assembly_friendly"]

for p in profiles:
    r = optimizer.optimize(gfp_aa, profile=p, scan_mode="fast")
    cai = r.metrics["cai"]
    gc  = r.metrics["gc_percent"]
    score = r.metrics["score"]
    print(f"{p:<20} CAI={cai:.3f}  GC={gc:.2f}%  Score={score:.3f}")
```

## Step 4: Custom Restriction Site Removal

For MoClo / Golden Gate assembly, remove specific restriction sites:

```bash
factorforge optimize gfp.fasta --engine profile --profile assembly_friendly \
  --scan-include restriction_sites \
  -o gfp_moclo.fasta
```

The engine performs synonymous substitutions to eliminate recognition sequences while preserving the amino acid sequence.

## Step 5: Downstream Use

The output FASTA can be reviewed for:

- **Gene synthesis** — review against vendor and project constraints before ordering
- **MoClo Level 0** — use with the `assembly_friendly` profile; check for BsaI/BpiI site removal
- **Agroinfiltration** — clone into a binary vector (e.g. pEAQ-HT, pK7WG2) for *A. tumefaciens*-mediated delivery

!!! note "Wet-lab validation"
    The `5' Ramp` and `Viral Delivery` profiles are currently **pending wet-lab validation** and are disabled by default. Share public-safe summaries via [Share Wet-lab Results (GitHub)](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml); use `eijex.lab@gmail.com` for private or sensitive summaries.

## Summary

| Step | Command |
|------|---------|
| Install | `pip install factorforge-cds` |
| Optimize (balanced) | `factorforge optimize gfp.fasta --engine profile --profile balanced -o out.fasta` |
| Compare profiles | `factorforge optimize gfp.fasta --engine profile --compare-profiles balanced,high_cai,gc_target` |
| Assembly-oriented | `factorforge optimize gfp.fasta --engine profile --profile assembly_friendly -o out.fasta` |

**GFP optimization results (balanced profile, N=1):**

| Metric | Result |
|--------|--------|
| Input length | 239 aa |
| Output CDS length | 720 bp (239 × 3 + stop) |
| CAI | 0.781 |
| GC% | 57.32% |
| Internal stop codons | 0 |
| Rare codon runs (≥3) | 0 |
| Composite score | 0.843 |
