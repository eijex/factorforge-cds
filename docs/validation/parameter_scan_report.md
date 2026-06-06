# Parameter Scan Report

**Total hits:** 231  **Unregistered:** 224

## Unregistered Public-Biological Numbers

| file | line | category | match | context |
|------|------|----------|-------|---------|
| `src\factorforge\database.py` | 76 | aa_identity | `1.0` | `algorithm_version: str = "2.1.0",` |
| `src\factorforge\analysis\feasibility.py` | 91 | cai | `0.82` | `target_cai: float = 0.82,` |
| `src\factorforge\analysis\feasibility.py` | 104 | cai | `0.82` | `target_cai=0.82 aligns with industry practice (>0.8) and is achievable.` |
| `src\factorforge\analysis\feasibility.py` | 92 | gc_global | `55.0` | `target_gc_low: float = 55.0,` |
| `src\factorforge\analysis\feasibility.py` | 110 | gc_global | `55.0` | `ranges = gc_ranges or [(55.0, 65.0), (50.0, 65.0), (40.0, 65.0)]` |
| `src\factorforge\analysis\feasibility.py` | 93 | gc_global | `65.0` | `target_gc_high: float = 65.0,` |
| `src\factorforge\analysis\feasibility.py` | 110 | gc_global | `65.0` | `ranges = gc_ranges or [(55.0, 65.0), (50.0, 65.0), (40.0, 65.0)]` |
| `src\factorforge\analysis\feasibility.py` | 110 | gc_global | `65.0` | `ranges = gc_ranges or [(55.0, 65.0), (50.0, 65.0), (40.0, 65.0)]` |
| `src\factorforge\analysis\feasibility.py` | 110 | gc_global | `65.0` | `ranges = gc_ranges or [(55.0, 65.0), (50.0, 65.0), (40.0, 65.0)]` |
| `src\factorforge\analysis\feasibility.py` | 26 | aa_identity | `1.0` | `return value * 100.0 if 0.0 <= value <= 1.0 else value` |
| `src\factorforge\analysis\metrics.py` | 283 | aa_identity | `1.0` | `tpa_weight: float = 1.0,` |
| `src\factorforge\analysis\metrics.py` | 288 | aa_identity | `1.0` | `Mammalian opt-in: set cpg_weight=1.0 and tpa_weight=1.0 to penalize both.` |
| `src\factorforge\analysis\metrics.py` | 288 | aa_identity | `1.0` | `Mammalian opt-in: set cpg_weight=1.0 and tpa_weight=1.0 to penalize both.` |
| `src\factorforge\analysis\metrics.py` | 293 | aa_identity | `1.0` | `return 1.0` |
| `src\factorforge\analysis\metrics.py` | 297 | aa_identity | `1.0` | `return 1.0` |
| `src\factorforge\analysis\metrics.py` | 301 | aa_identity | `1.0` | `cpg_score = max(0.0, 1.0 - cpg_ratio / 2.0)` |
| `src\factorforge\analysis\metrics.py` | 302 | aa_identity | `1.0` | `tpa_score = max(0.0, 1.0 - tpa_ratio / 2.0)` |
| `src\factorforge\cli\main.py` | 167 | gc_global | `55.0` | `@click.option("--gc-min", type=float, default=55.0, help="Minimum target GC perc` |
| `src\factorforge\cli\main.py` | 168 | gc_global | `65.0` | `@click.option("--gc-max", type=float, default=65.0, help="Maximum target GC perc` |
| `src\factorforge\engines\profile\codon_table_builder.py` | 33 | aa_identity | `1.0` | `blend_ratio: Weight for high-expression data (0.0-1.0). Default 0.7.` |
| `src\factorforge\engines\profile\codon_table_builder.py` | 39 | aa_identity | `1.0` | `if not 0.0 <= blend_ratio <= 1.0:` |
| `src\factorforge\engines\profile\codon_table_builder.py` | 40 | aa_identity | `1.0` | `raise ValueError(f"blend_ratio must be between 0.0 and 1.0, got {blend_ratio}")` |
| `src\factorforge\engines\profile\codon_table_builder.py` | 70 | aa_identity | `1.0` | `# Normalize per amino acid (frequencies must sum to 1.0)` |
| `src\factorforge\engines\profile\construct_builder.py` | 208 | type_iis | `BsaI` | `"BsaI": ["GGTCTC", "GAGACC"],` |
| `src\factorforge\engines\profile\construct_builder.py` | 210 | type_iis | `BsmBI` | `"BsmBI": ["CGTCTC", "GAGACG"],` |
| `src\factorforge\engines\profile\exporter.py` | 448 | type_iis | `BsaI` | `"assembly_standard": "Golden Gate (BsaI)",` |
| `src\factorforge\engines\profile\exporter.py` | 451 | type_iis | `BsaI` | `"type": "BsaI site",` |
| `src\factorforge\engines\profile\scoring.py` | 20 | gc_global | `55.0` | `GC_OPT_MIN = 55.0` |
| `src\factorforge\engines\profile\scoring.py` | 21 | gc_global | `65.0` | `GC_OPT_MAX = 65.0` |
| `src\factorforge\engines\profile\scoring.py` | 38 | aa_identity | `1.0` | `tpa_weight: float = 1.0  # plant default: TpA active` |
| `src\factorforge\engines\profile\scoring.py` | 48 | aa_identity | `1.0` | `"""Normalize weights to sum to 1.0."""` |
| `src\factorforge\engines\profile\scoring.py` | 52 | aa_identity | `1.0` | `"""Ensure active weights sum to 1.0."""` |
| `src\factorforge\engines\profile\scoring.py` | 74 | aa_identity | `1.0` | `# high_cai: CAI 1.0 mimics naturally high-expression N. benthamiana genes (golde` |
| `src\factorforge\engines\profile\scoring.py` | 158 | aa_identity | `1.0` | `# Map to [0, 1] where 0.0 kcal/mol/nt → 1.0 and -0.5 → 0.0` |
| `src\factorforge\engines\profile\scoring.py` | 159 | aa_identity | `1.0` | `return 1.0 + (clamped / 0.5)` |
| `src\factorforge\engines\profile\scoring.py` | 170 | aa_identity | `1.0` | `Returns 1.0 inside [gc_min, gc_max]; linearly decays to 0.0 after` |
| `src\factorforge\engines\profile\scoring.py` | 187 | aa_identity | `1.0` | `return 1.0` |
| `src\factorforge\engines\profile\scoring.py` | 189 | aa_identity | `1.0` | `return max(0.0, 1.0 - distance / decay_width)` |
| `src\factorforge\engines\profile\scoring.py` | 195 | aa_identity | `1.0` | `tpa_weight: float = 1.0,` |
| `src\factorforge\engines\profile\scoring.py` | 200 | aa_identity | `1.0` | `Mammalian opt-in: set cpg_weight=1.0 and tpa_weight=1.0 to penalize both.` |
| `src\factorforge\engines\profile\scoring.py` | 200 | aa_identity | `1.0` | `Mammalian opt-in: set cpg_weight=1.0 and tpa_weight=1.0 to penalize both.` |
| `src\factorforge\engines\profile\scoring.py` | 221 | aa_identity | `1.0` | `full score (1.0); outside the band the score decays linearly to 0.0 over` |
| `src\factorforge\engines\profile\scoring.py` | 244 | aa_identity | `1.0` | `cai_score = max(0.0, min(1.0, cai))` |
| `src\factorforge\engines\profile\scoring_ml.py` | 112 | aa_identity | `1.0` | `return round(max(0.0, min(1.0, math.exp(mean_log_prob))), 3)` |
| `src\factorforge\engines\profile\utils.py` | 86 | aa_identity | `1.0` | `from mononucleotide composition. Ratio < 1.0 means suppressed,` |
| `src\factorforge\engines\profile\utils.py` | 87 | aa_identity | `1.0` | `> 1.0 means enriched.` |
| `src\factorforge\engines\profile\utils.py` | 98 | aa_identity | `1.0` | `1.0` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 18 | type_iis | `BsaI` | `- Golden Gate (BsaI, BpiI, BsmBI)` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 19 | type_iis | `BsaI` | `- MoClo (BsaI + overhangs)` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 26 | type_iis | `BsaI` | `"enzymes": ["BsaI", "BpiI", "BsmBI"],` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 28 | type_iis | `BsaI` | `"BsaI": ["GGTCTC", "GAGACC"],  # Forward and reverse complement` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 34 | type_iis | `BsaI` | `"enzymes": ["BsaI"],` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 35 | type_iis | `BsaI` | `"sites": {"BsaI": ["GGTCTC", "GAGACC"]},` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 91 | type_iis | `BsaI` | `[{'enzyme': 'BsaI', ...}]` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 377 | type_iis | `BsaI` | `# Test sequence (includes BsaI site)` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 18 | type_iis | `BsmBI` | `- Golden Gate (BsaI, BpiI, BsmBI)` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 26 | type_iis | `BsmBI` | `"enzymes": ["BsaI", "BpiI", "BsmBI"],` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 30 | type_iis | `BsmBI` | `"BsmBI": ["CGTCTC", "GAGACG"],` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 113 | aa_identity | `1.0` | `cumprob_list[-1] = 1.0` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 201 | aa_identity | `1.0` | `CAI value (0.0 ~ 1.0).` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 49 | type_iis | `BsaI` | `4. Assembly-Friendly: avoid BsaI/BpiI` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 121 | type_iis | `BsaI` | `"BsaI": ["GGTCTC", "GAGACC"],` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 492 | type_iis | `BsaI` | `- Retries up to max_attempts times until no BsaI/BpiI Type IIS` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 497 | type_iis | `BsaI` | `- Supported: BsaI/BpiI site avoidance via stochastic retry` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 123 | type_iis | `BsmBI` | `"BsmBI": ["CGTCTC", "GAGACG"],` |
| `src\factorforge\schemas\design_package.py` | 48 | aa_identity | `1.0` | `aa_identity: float = 1.0` |
| `src\factorforge\schemas\design_package.py` | 79 | aa_identity | `1.0` | `design_package_version: str = "1.0"` |
| `src\factorforge\utils\restriction_sites.py` | 15 | type_iis | `BsaI` | `{"name": "BsaI", "sequence": "GGTCTC", "scan_rc": True},` |
| `src\factorforge\utils\restriction_sites.py` | 17 | type_iis | `BsmBI` | `{"name": "BsmBI", "sequence": "CGTCTC", "scan_rc": True},` |
| `src\factorforge\utils\validation.py` | 70 | aa_identity | `1.0` | `if identity < 1.0:` |
| `src\factorforge\utils\validation.py` | 71 | aa_identity | `1.0` | `errors.append("amino acid identity is below 1.0")` |
| `docs\assembly-friendly.md` | 11 | type_iis | `BsaI` | `- BsaI and BpiI Type IIS restriction site avoidance through synonymous` |
| `docs\changelog.md` | 96 | aa_identity | `1.0` | `## v3.1.0 — 2026-05-24` |
| `docs\host-profile-registry.md` | 70 | cai | `0.82` | `min_feasible_cai: 0.82` |
| `docs\host-profile-registry.md` | 102 | cai | `0.82` | `min_feasible_cai: 0.82` |
| `docs\host-profile-registry.md` | 134 | cai | `0.82` | `min_feasible_cai: 0.82` |
| `docs\host-profile-registry.md` | 166 | cai | `0.82` | `min_feasible_cai: 0.82` |
| `docs\host-profile-registry.md` | 198 | cai | `0.82` | `min_feasible_cai: 0.82` |
| `docs\host-profile-registry.md` | 66 | gc_global | `55.0` | `ideal_band: [55.0, 65.0]` |
| `docs\host-profile-registry.md` | 98 | gc_global | `55.0` | `ideal_band: [55.0, 65.0]` |
| `docs\host-profile-registry.md` | 194 | gc_global | `55.0` | `ideal_band: [35.0, 55.0]` |
| `docs\host-profile-registry.md` | 66 | gc_global | `65.0` | `ideal_band: [55.0, 65.0]` |
| `docs\host-profile-registry.md` | 98 | gc_global | `65.0` | `ideal_band: [55.0, 65.0]` |
| `docs\host-profile-registry.md` | 130 | gc_global | `65.0` | `ideal_band: [45.0, 65.0]` |
| `docs\host-profile-registry.md` | 162 | gc_global | `65.0` | `ideal_band: [45.0, 65.0]` |
| `docs\host-profile-registry.md` | 67 | gc_local | `75.0` | `soft_band: [45.0, 75.0]` |
| `docs\host-profile-registry.md` | 99 | gc_local | `75.0` | `soft_band: [45.0, 75.0]` |
| `docs\how-it-works.md` | 24 | type_iis | `BsaI` | `BsaI / BsmBI recognition sites via silent edits` |
| `docs\how-it-works.md` | 24 | type_iis | `BsmBI` | `BsaI / BsmBI recognition sites via silent edits` |
| `docs\output.md` | 11 | type_iis | `BsaI` | `| **Domestication report** | BsaI/BsmBI and custom restriction sites removed, ed` |
| `docs\output.md` | 11 | type_iis | `BsmBI` | `| **Domestication report** | BsaI/BsmBI and custom restriction sites removed, ed` |
| `docs\profiles.md` | 24 | type_iis | `BsaI` | `| `assembly_friendly` | Golden Gate / MoClo workflows — avoids BsaI/BpiI Type II` |
| `docs\scoring-calibration.md` | 11 | cai | `0.82` | `A typical lower bound is `0.82`.` |
| `docs\scoring-calibration.md` | 44 | gc_global | `55.0` | `current plant default, the accepted band is `55.0` to `65.0` GC with a` |
| `docs\scoring-calibration.md` | 44 | gc_global | `65.0` | `current plant default, the accepted band is `55.0` to `65.0` GC with a` |
| `docs\scoring-calibration.md` | 25 | aa_identity | `1.0` | `1.0 - abs(gc - gc_opt) / 50.0` |
| `docs\validation.md` | 17 | type_iis | `BsaI` | `| Forbidden restriction sites | ✅ | BsaI, BsmBI, BpII (Golden Gate) |` |
| `docs\validation.md` | 17 | type_iis | `BsmBI` | `| Forbidden restriction sites | ✅ | BsaI, BsmBI, BpII (Golden Gate) |` |
| `docs\tutorials\gfp-nbenthamiana.md` | 155 | type_iis | `BsaI` | `- **MoClo Level 0** — use with the `assembly_friendly` profile; check for BsaI/B` |
| `tests\test_database.py` | 33 | cai | `0.82` | `"cai": 0.82,` |
| `tests\test_restriction_sites.py` | 193 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\test_restriction_sites.py` | 96 | aa_identity | `1.0` | `"CTC": {"aa": "L", "frequency": 1.0},` |
| `tests\test_restriction_sites.py` | 128 | aa_identity | `1.0` | `"GGT": {"aa": "G", "frequency": 1.0},` |
| `tests\test_restriction_sites.py` | 129 | aa_identity | `1.0` | `"CTC": {"aa": "L", "frequency": 1.0},` |
| `tests\test_restriction_sites.py` | 114 | type_iis | `BsaI` | `{"name": "BsaI", "sequence": "GGTCTC"},` |
| `tests\test_restriction_sites.py` | 182 | type_iis | `BsaI` | `assert result["removed_sites"][0]["enzyme"] == "BsaI"` |
| `tests\test_sequence_validator.py` | 69 | aa_identity | `1.0` | `assert result == {"passed": True, "errors": [], "aa_identity": 1.0}` |
| `tests\api\test_optimize_contract.py` | 55 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 85 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 111 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 138 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 155 | gc_global | `55.0` | `constraints={"gc_min": 45.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 170 | gc_global | `55.0` | `{"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 187 | gc_global | `55.0` | `{"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\engines\profile\test_codon_table_builder.py` | 24 | aa_identity | `1.0` | `"M": {"ATG": 1.0},` |
| `tests\engines\profile\test_codon_table_builder.py` | 48 | aa_identity | `1.0` | `"ATG": {"aa": "M", "frequency": 1.0, "per_thousand": 22.0},` |
| `tests\engines\profile\test_codon_table_builder.py` | 53 | aa_identity | `1.0` | `"TAA": {"aa": "*", "frequency": 0.48, "per_thousand": 1.0},` |
| `tests\engines\profile\test_codon_table_builder.py` | 84 | aa_identity | `1.0` | `"""Blended frequencies per amino acid sum to ~1.0."""` |
| `tests\engines\profile\test_codon_table_builder.py` | 88 | aa_identity | `1.0` | `assert abs(total - 1.0) < 0.02, f"AA {aa}: frequencies sum to {total}"` |
| `tests\engines\profile\test_codon_table_builder.py` | 143 | aa_identity | `1.0` | `assert 0.0 <= w <= 1.0, f"Codon {codon} has weight {w}"` |
| `tests\engines\profile\test_codon_table_builder.py` | 146 | aa_identity | `1.0` | `"""Preferred codon (highest freq per AA) has weight 1.0."""` |
| `tests\engines\profile\test_codon_table_builder.py` | 148 | aa_identity | `1.0` | `assert translator.golden_ref_weights["ATG"] == 1.0` |
| `tests\engines\profile\test_dinucleotide.py` | 76 | aa_identity | `1.0` | `assert ratio > 1.0` |
| `tests\engines\profile\test_dinucleotide.py` | 86 | aa_identity | `1.0` | `assert ratio < 1.0` |
| `tests\engines\profile\test_dinucleotide.py` | 149 | aa_identity | `1.0` | `assert 0.0 <= score <= 1.0` |
| `tests\engines\profile\test_dinucleotide.py` | 152 | aa_identity | `1.0` | `"""Very short sequence scores 1.0 (too short to evaluate)."""` |
| `tests\engines\profile\test_dinucleotide.py` | 153 | aa_identity | `1.0` | `assert calculate_dinucleotide_score("ATG") == 1.0` |
| `tests\engines\profile\test_dinucleotide.py` | 158 | aa_identity | `1.0` | `"CCCCGGGGCCCCGGGG", cpg_weight=1.0, tpa_weight=1.0` |
| `tests\engines\profile\test_dinucleotide.py` | 158 | aa_identity | `1.0` | `"CCCCGGGGCCCCGGGG", cpg_weight=1.0, tpa_weight=1.0` |
| `tests\engines\profile\test_dinucleotide.py` | 161 | aa_identity | `1.0` | `"CGCGCGCGCGCGCGCG", cpg_weight=1.0, tpa_weight=1.0` |
| `tests\engines\profile\test_dinucleotide.py` | 161 | aa_identity | `1.0` | `"CGCGCGCGCGCGCGCG", cpg_weight=1.0, tpa_weight=1.0` |
| `tests\engines\profile\test_dinucleotide.py` | 178 | aa_identity | `1.0` | `assert 0.0 <= score <= 1.0` |
| `tests\engines\profile\test_dinucleotide.py` | 181 | aa_identity | `1.0` | `"""ScoringConfig normalizes all 4 weights to sum to 1.0."""` |
| `tests\engines\profile\test_dinucleotide.py` | 184 | aa_identity | `1.0` | `assert abs(total - 1.0) < 1e-6` |
| `tests\engines\profile\test_dinucleotide_fix.py` | 78 | aa_identity | `1.0` | `"CCG": {"aa": "P", "frequency": 1.0},` |
| `tests\engines\profile\test_domesticator.py` | 111 | aa_identity | `1.0` | `"GGT": {"aa": "G", "frequency": 1.0},` |
| `tests\engines\profile\test_domesticator.py` | 112 | aa_identity | `1.0` | `"CTC": {"aa": "L", "frequency": 1.0},` |
| `tests\engines\profile\test_domesticator.py` | 171 | aa_identity | `1.0` | `"GGT": {"aa": "G", "frequency": 1.0},` |
| `tests\engines\profile\test_domesticator.py` | 172 | aa_identity | `1.0` | `"CTC": {"aa": "L", "frequency": 1.0},` |
| `tests\engines\profile\test_domesticator.py` | 238 | aa_identity | `1.0` | `"GGT": {"aa": "G", "frequency": 1.0},` |
| `tests\engines\profile\test_domesticator.py` | 239 | aa_identity | `1.0` | `"CTC": {"aa": "L", "frequency": 1.0},` |
| `tests\engines\profile\test_domesticator.py` | 272 | aa_identity | `1.0` | `"GAG": {"aa": "E", "frequency": 1.0},` |
| `tests\engines\profile\test_domesticator.py` | 273 | aa_identity | `1.0` | `"ACC": {"aa": "T", "frequency": 1.0},` |
| `tests\engines\profile\test_domesticator.py` | 305 | aa_identity | `1.0` | `"GGT": {"aa": "G", "frequency": 1.0},` |
| `tests\engines\profile\test_domesticator.py` | 306 | aa_identity | `1.0` | `"CTC": {"aa": "L", "frequency": 1.0},` |
| `tests\engines\profile\test_domesticator.py` | 26 | type_iis | `BsaI` | `"""Test BsaI site detection"""` |
| `tests\engines\profile\test_domesticator.py` | 27 | type_iis | `BsaI` | `seq = "ATGGGTCTCGAGGAGCTG"  # Contains GGTCTC (BsaI)` |
| `tests\engines\profile\test_domesticator.py` | 32 | type_iis | `BsaI` | `assert any(s["enzyme"] == "BsaI" for s in sites)` |
| `tests\engines\profile\test_domesticator.py` | 45 | type_iis | `BsaI` | `# BsaI: GGTCTC, BpiI: GAAGAC` |
| `tests\engines\profile\test_domesticator.py` | 68 | type_iis | `BsaI` | `# Create a sequence with BsaI site that can be removed` |
| `tests\engines\profile\test_domesticator.py` | 122 | type_iis | `BsaI` | `assert result["unfixable"][0]["enzyme"] == "BsaI"` |
| `tests\engines\profile\test_domesticator.py` | 156 | type_iis | `BsaI` | `site = {"enzyme": "BsaI", "site": "GGTCTC", "position": 3}` |
| `tests\engines\profile\test_domesticator.py` | 190 | type_iis | `BsaI` | `assert "BsaI" in Domesticator.ASSEMBLY_STANDARDS["golden_gate"]["sites"]` |
| `tests\engines\profile\test_domesticator.py` | 258 | type_iis | `BsaI` | `- GAGACC (BsaI RC, pos:0): GAG+ACC 모두 synonymous 없음 → unfixable` |
| `tests\engines\profile\test_domesticator.py` | 259 | type_iis | `BsaI` | `- GGTCTC (BsaI fwd, pos:6): GGT에 GGC synonymous 있음 → fixable` |
| `tests\engines\profile\test_domesticator.py` | 310 | type_iis | `BsaI` | `seq = "GGTCTC"  # BsaI site, unfixable (no synonymous codons)` |
| `tests\engines\profile\test_exporter.py` | 251 | gc_global | `55.0` | `{"sequence": "ATGGCCTAG", "metadata": {"gene_name": "gene2", "cai": 0.9, "gc": 5` |
| `tests\engines\profile\test_exporter.py` | 264 | gc_global | `55.0` | `{"sequence": "ATGGCCTAG", "metadata": {"gene_name": "gene2", "cai": 0.9, "gc": 5` |
| `tests\engines\profile\test_exporter.py` | 39 | type_iis | `BsaI` | `"assembly_standard": "Golden Gate (BsaI)",` |
| `tests\engines\profile\test_exporter.py` | 42 | type_iis | `BsaI` | `"violations_fixed": [{"type": "BsaI site", "position": 147}],` |
| `tests\engines\profile\test_exporter.py` | 314 | type_iis | `BsaI` | `assert "BsaI site" in report` |
| `tests\engines\profile\test_mfe_evidence.py` | 55 | aa_identity | `1.0` | `assert 0.0 <= score <= 1.0` |
| `tests\engines\profile\test_nterminal_ramp.py` | 57 | aa_identity | `1.0` | `assert 0.0 <= cai_n <= 1.0` |
| `tests\engines\profile\test_nterminal_ramp.py` | 58 | aa_identity | `1.0` | `assert 0.0 <= cai_c <= 1.0` |
| `tests\engines\profile\test_pipeline_integration.py` | 230 | aa_identity | `1.0` | `assert 0.0 <= features["cai_score"] <= 1.0` |
| `tests\engines\profile\test_pipeline_integration.py` | 110 | type_iis | `BsaI` | `"unfixable": [{"enzyme": "BsaI", "site": "GGTCTC", "position": 0}],` |
| `tests\engines\profile\test_rare_codon_runs.py` | 25 | aa_identity | `1.0` | `"ATG": {"aa": "M", "frequency": 1.0},` |
| `tests\engines\profile\test_rare_codon_runs.py` | 27 | aa_identity | `1.0` | `"AAG": {"aa": "K", "frequency": 1.0},` |
| `tests\engines\profile\test_reverse_translator.py` | 182 | gc_global | `55.0` | `assert 45.0 <= gc <= 55.0` |
| `tests\engines\profile\test_reverse_translator.py` | 144 | gc_global | `65.0` | `assert 35.0 <= gc <= 65.0` |
| `tests\engines\profile\test_reverse_translator.py` | 50 | aa_identity | `1.0` | `"""Test that codon frequencies sum to ~1.0 for each AA"""` |
| `tests\engines\profile\test_reverse_translator.py` | 65 | aa_identity | `1.0` | `assert 0.0 <= cai <= 1.0` |
| `tests\engines\profile\test_reverse_translator.py` | 375 | aa_identity | `1.0` | `assert 0.0 <= cand["score"] <= 1.0` |
| `tests\engines\profile\test_reverse_translator.py` | 240 | type_iis | `BsaI` | `"""Test that assembly-friendly avoids BsaI sites"""` |
| `tests\engines\profile\test_reverse_translator.py` | 245 | type_iis | `BsaI` | `# BsaI recognition site` |
| `tests\engines\profile\test_rule_engine.py` | 29 | aa_identity | `1.0` | `"ATG": {"aa": "M", "frequency": 1.0},` |
| `tests\engines\profile\test_rule_engine.py` | 123 | aa_identity | `1.0` | `"ATG": {"aa": "M", "frequency": 1.0},` |
| `tests\engines\profile\test_rule_engine.py` | 124 | aa_identity | `1.0` | `"AAT": {"aa": "N", "frequency": 1.0},` |
| `tests\engines\profile\test_rule_engine.py` | 125 | aa_identity | `1.0` | `"AAA": {"aa": "K", "frequency": 1.0},` |
| `tests\engines\profile\test_rule_engine.py` | 126 | aa_identity | `1.0` | `"TAA": {"aa": "*", "frequency": 1.0},` |
| `tests\engines\profile\test_scoring.py` | 139 | gc_global | `55.0` | `assert gc_band_score(60.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 142 | gc_global | `55.0` | `assert gc_band_score(55.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 142 | gc_global | `55.0` | `assert gc_band_score(55.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 145 | gc_global | `55.0` | `assert gc_band_score(65.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 157 | gc_global | `55.0` | `assert gc_band_score(85.0, 55.0, 65.0, decay_width=20.0) == 0.0` |
| `tests\engines\profile\test_scoring.py` | 160 | gc_global | `55.0` | `assert gc_band_score(10.0, 55.0, 65.0) == 0.0` |
| `tests\engines\profile\test_scoring.py` | 172 | gc_global | `55.0` | `score_at_min = calculate_composite_score(cai=0.8, gc=55.0, profile="balanced")` |
| `tests\engines\profile\test_scoring.py` | 139 | gc_global | `65.0` | `assert gc_band_score(60.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 142 | gc_global | `65.0` | `assert gc_band_score(55.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 145 | gc_global | `65.0` | `assert gc_band_score(65.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 145 | gc_global | `65.0` | `assert gc_band_score(65.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 157 | gc_global | `65.0` | `assert gc_band_score(85.0, 55.0, 65.0, decay_width=20.0) == 0.0` |
| `tests\engines\profile\test_scoring.py` | 160 | gc_global | `65.0` | `assert gc_band_score(10.0, 55.0, 65.0) == 0.0` |
| `tests\engines\profile\test_scoring.py` | 127 | gc_local | `75.0` | `result = normalize_mfe(-75.0, 300)  # -0.25 kcal/mol/nt` |
| `tests\engines\profile\test_scoring.py` | 29 | aa_identity | `1.0` | `"""Active weights should sum to 1.0 after normalization."""` |
| `tests\engines\profile\test_scoring.py` | 32 | aa_identity | `1.0` | `assert abs(total - 1.0) < 1e-6` |
| `tests\engines\profile\test_scoring.py` | 39 | aa_identity | `1.0` | `assert abs(total - 1.0) < 1e-6` |
| `tests\engines\profile\test_scoring.py` | 60 | aa_identity | `1.0` | `score = calculate_composite_score(cai=1.0, gc=60.0, profile="balanced")` |
| `tests\engines\profile\test_scoring.py` | 96 | aa_identity | `1.0` | `cfg = ScoringConfig(w_cai=1.0, w_gc=0.0, w_mfe=0.0, use_mfe=False)` |
| `tests\engines\profile\test_scoring.py` | 102 | aa_identity | `1.0` | `for cai in [0.0, 0.5, 1.0]:` |
| `tests\engines\profile\test_scoring.py` | 105 | aa_identity | `1.0` | `assert 0.0 <= score <= 1.0` |
| `tests\engines\profile\test_scoring.py` | 110 | aa_identity | `1.0` | `assert 0.0 <= score <= 1.0` |
| `tests\engines\profile\test_scoring.py` | 117 | aa_identity | `1.0` | `"""MFE = 0 (no structure) → normalized = 1.0."""` |
| `tests\engines\profile\test_scoring.py` | 118 | aa_identity | `1.0` | `assert normalize_mfe(0.0, 300) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 139 | aa_identity | `1.0` | `assert gc_band_score(60.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 142 | aa_identity | `1.0` | `assert gc_band_score(55.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 145 | aa_identity | `1.0` | `assert gc_band_score(65.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 171 | aa_identity | `1.0` | `# GC=55 is at the lower boundary (score=1.0 with band, was 0.9 with old /50 form` |
| `tests\engines\profile\test_scoring.py` | 205 | aa_identity | `1.0` | `sequence, cpg_weight=0.0, tpa_weight=1.0` |
| `tests\engines\profile\test_scoring.py` | 211 | aa_identity | `1.0` | `"CCCCGGGGCCCCGGGG", cpg_weight=1.0, tpa_weight=1.0` |
| `tests\engines\profile\test_scoring.py` | 211 | aa_identity | `1.0` | `"CCCCGGGGCCCCGGGG", cpg_weight=1.0, tpa_weight=1.0` |
| `tests\engines\profile\test_scoring.py` | 214 | aa_identity | `1.0` | `"CGCGCGCGCGCGCGCG", cpg_weight=1.0, tpa_weight=1.0` |
| `tests\engines\profile\test_scoring.py` | 214 | aa_identity | `1.0` | `"CGCGCGCGCGCGCGCG", cpg_weight=1.0, tpa_weight=1.0` |
| `tests\engines\profile\test_scoring.py` | 245 | aa_identity | `1.0` | `assert abs(total - 1.0) < 1e-6` |
| `tests\engines\profile\test_utils.py` | 57 | aa_identity | `1.0` | `"codons": {"ATG": {"aa": "M", "frequency": 1.0}},` |
| `tests\test_analysis\test_feasibility.py` | 10 | aa_identity | `1.0` | `"ATG": 1.0,` |
| `tests\test_analysis\test_feasibility.py` | 12 | aa_identity | `1.0` | `"GCC": 1.0,` |
| `tests\test_analysis\test_feasibility.py` | 15 | aa_identity | `1.0` | `"AAA": 1.0,` |
| `tests\test_analysis\test_metrics.py` | 43 | aa_identity | `1.0` | `assert amino_acid_identity("MAF", "ATGGCTTTC") == 1.0` |
| `tests\test_analysis\test_metrics.py` | 48 | aa_identity | `1.0` | `weights = {"ATG": 1.0, "GCC": 1.0, "GCT": 0.25}` |
| `tests\test_analysis\test_metrics.py` | 48 | aa_identity | `1.0` | `weights = {"ATG": 1.0, "GCC": 1.0, "GCT": 0.25}` |
| `tests\test_schemas\test_design_package.py` | 72 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\test_schemas\test_design_package.py` | 32 | aa_identity | `1.0` | `assert pkg.design_package_version == "1.0"` |
| `tests\test_utils\test_validation.py` | 12 | aa_identity | `1.0` | `assert result["amino_acid_identity"] == 1.0` |

## Allowlisted (not requiring registry entry)

| file | line | match |
|------|------|-------|
| `tests\engines\profile\test_scoring.py` | 149 | `55.0` |
| `tests\engines\profile\test_scoring.py` | 153 | `55.0` |
| `tests\engines\profile\test_scoring.py` | 149 | `65.0` |
| `tests\engines\profile\test_scoring.py` | 153 | `65.0` |
| `tests\engines\profile\test_scoring.py` | 149 | `75.0` |
| `tests\engines\profile\test_scoring.py` | 194 | `1.0` |
| `tests\test_analysis\test_metrics.py` | 49 | `1.0` |