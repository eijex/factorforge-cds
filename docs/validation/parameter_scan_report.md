# Parameter Scan Report

**Total hits:** 108  **Unregistered:** 103

## Unregistered Public-Biological Numbers

| file | line | category | match | context |
|------|------|----------|-------|---------|
| `src\factorforge\analysis\feasibility.py` | 91 | cai | `0.82` | `target_cai: float = 0.82,` |
| `src\factorforge\analysis\feasibility.py` | 104 | cai | `0.82` | `target_cai=0.82 aligns with industry practice (>0.8) and is achievable.` |
| `src\factorforge\analysis\feasibility.py` | 92 | gc_global | `55.0` | `target_gc_low: float = 55.0,` |
| `src\factorforge\analysis\feasibility.py` | 110 | gc_global | `55.0` | `ranges = gc_ranges or [(55.0, 65.0), (50.0, 65.0), (40.0, 65.0)]` |
| `src\factorforge\analysis\feasibility.py` | 93 | gc_global | `65.0` | `target_gc_high: float = 65.0,` |
| `src\factorforge\analysis\feasibility.py` | 110 | gc_global | `65.0` | `ranges = gc_ranges or [(55.0, 65.0), (50.0, 65.0), (40.0, 65.0)]` |
| `src\factorforge\cli\main.py` | 167 | gc_global | `55.0` | `@click.option("--gc-min", type=float, default=55.0, help="Minimum target GC perc` |
| `src\factorforge\cli\main.py` | 168 | gc_global | `65.0` | `@click.option("--gc-max", type=float, default=65.0, help="Maximum target GC perc` |
| `src\factorforge\engines\profile\construct_builder.py` | 208 | type_iis | `BsaI` | `"BsaI": ["GGTCTC", "GAGACC"],` |
| `src\factorforge\engines\profile\construct_builder.py` | 210 | type_iis | `BsmBI` | `"BsmBI": ["CGTCTC", "GAGACG"],` |
| `src\factorforge\engines\profile\exporter.py` | 448 | type_iis | `BsaI` | `"assembly_standard": "Golden Gate (BsaI)",` |
| `src\factorforge\engines\profile\exporter.py` | 451 | type_iis | `BsaI` | `"type": "BsaI site",` |
| `src\factorforge\engines\profile\scoring.py` | 20 | gc_global | `55.0` | `GC_OPT_MIN = 55.0` |
| `src\factorforge\engines\profile\scoring.py` | 21 | gc_global | `65.0` | `GC_OPT_MAX = 65.0` |
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
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 49 | type_iis | `BsaI` | `4. Assembly-Friendly: avoid BsaI/BpiI` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 121 | type_iis | `BsaI` | `"BsaI": ["GGTCTC", "GAGACC"],` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 492 | type_iis | `BsaI` | `- Retries up to max_attempts times until no BsaI/BpiI Type IIS` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 497 | type_iis | `BsaI` | `- Supported: BsaI/BpiI site avoidance via stochastic retry` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 123 | type_iis | `BsmBI` | `"BsmBI": ["CGTCTC", "GAGACG"],` |
| `src\factorforge\utils\restriction_sites.py` | 15 | type_iis | `BsaI` | `{"name": "BsaI", "sequence": "GGTCTC", "scan_rc": True},` |
| `src\factorforge\utils\restriction_sites.py` | 17 | type_iis | `BsmBI` | `{"name": "BsmBI", "sequence": "CGTCTC", "scan_rc": True},` |
| `docs\assembly-friendly.md` | 11 | type_iis | `BsaI` | `- BsaI and BpiI Type IIS restriction site avoidance through synonymous` |
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
| `docs\validation.md` | 17 | type_iis | `BsaI` | `| Forbidden restriction sites | ✅ | BsaI, BsmBI, BpII (Golden Gate) |` |
| `docs\validation.md` | 17 | type_iis | `BsmBI` | `| Forbidden restriction sites | ✅ | BsaI, BsmBI, BpII (Golden Gate) |` |
| `docs\tutorials\gfp-nbenthamiana.md` | 155 | type_iis | `BsaI` | `- **MoClo Level 0** — use with the `assembly_friendly` profile; check for BsaI/B` |
| `tests\test_database.py` | 33 | cai | `0.82` | `"cai": 0.82,` |
| `tests\test_restriction_sites.py` | 193 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\test_restriction_sites.py` | 114 | type_iis | `BsaI` | `{"name": "BsaI", "sequence": "GGTCTC"},` |
| `tests\test_restriction_sites.py` | 182 | type_iis | `BsaI` | `assert result["removed_sites"][0]["enzyme"] == "BsaI"` |
| `tests\api\test_optimize_contract.py` | 55 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 85 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 111 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 138 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 155 | gc_global | `55.0` | `constraints={"gc_min": 45.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 170 | gc_global | `55.0` | `{"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 187 | gc_global | `55.0` | `{"gc_min": 40.0, "gc_max": 55.0},` |
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
| `tests\engines\profile\test_pipeline_integration.py` | 110 | type_iis | `BsaI` | `"unfixable": [{"enzyme": "BsaI", "site": "GGTCTC", "position": 0}],` |
| `tests\engines\profile\test_reverse_translator.py` | 182 | gc_global | `55.0` | `assert 45.0 <= gc <= 55.0` |
| `tests\engines\profile\test_reverse_translator.py` | 144 | gc_global | `65.0` | `assert 35.0 <= gc <= 65.0` |
| `tests\engines\profile\test_reverse_translator.py` | 240 | type_iis | `BsaI` | `"""Test that assembly-friendly avoids BsaI sites"""` |
| `tests\engines\profile\test_reverse_translator.py` | 245 | type_iis | `BsaI` | `# BsaI recognition site` |
| `tests\engines\profile\test_scoring.py` | 139 | gc_global | `55.0` | `assert gc_band_score(60.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 142 | gc_global | `55.0` | `assert gc_band_score(55.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 145 | gc_global | `55.0` | `assert gc_band_score(65.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 157 | gc_global | `55.0` | `assert gc_band_score(85.0, 55.0, 65.0, decay_width=20.0) == 0.0` |
| `tests\engines\profile\test_scoring.py` | 160 | gc_global | `55.0` | `assert gc_band_score(10.0, 55.0, 65.0) == 0.0` |
| `tests\engines\profile\test_scoring.py` | 172 | gc_global | `55.0` | `score_at_min = calculate_composite_score(cai=0.8, gc=55.0, profile="balanced")` |
| `tests\engines\profile\test_scoring.py` | 139 | gc_global | `65.0` | `assert gc_band_score(60.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 142 | gc_global | `65.0` | `assert gc_band_score(55.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 145 | gc_global | `65.0` | `assert gc_band_score(65.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 157 | gc_global | `65.0` | `assert gc_band_score(85.0, 55.0, 65.0, decay_width=20.0) == 0.0` |
| `tests\engines\profile\test_scoring.py` | 160 | gc_global | `65.0` | `assert gc_band_score(10.0, 55.0, 65.0) == 0.0` |
| `tests\engines\profile\test_scoring.py` | 127 | gc_local | `75.0` | `result = normalize_mfe(-75.0, 300)  # -0.25 kcal/mol/nt` |
| `tests\test_schemas\test_design_package.py` | 72 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |

## Allowlisted (not requiring registry entry)

| file | line | match |
|------|------|-------|
| `tests\engines\profile\test_scoring.py` | 149 | `55.0` |
| `tests\engines\profile\test_scoring.py` | 153 | `55.0` |
| `tests\engines\profile\test_scoring.py` | 149 | `65.0` |
| `tests\engines\profile\test_scoring.py` | 153 | `65.0` |
| `tests\engines\profile\test_scoring.py` | 149 | `75.0` |