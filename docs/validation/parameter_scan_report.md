# Parameter Scan Report

**Total hits:** 160  **Unregistered:** 155

## Unregistered Public-Biological Numbers

| file | line | category | match | context |
|------|------|----------|-------|---------|
| `src\factorforge\analysis\feasibility.py` | 17 | cai | `0.82` | `# DEFAULT_CAI_TARGET=0.82 aligns with industry practice (>0.8) and is achievable` |
| `src\factorforge\analysis\feasibility.py` | 26 | cai | `0.82` | `DEFAULT_CAI_TARGET: float = 0.82` |
| `src\factorforge\analysis\feasibility.py` | 124 | gc_global | `55.0` | `# progressively wider windows. (55.0, 65.0) is retained ONLY as an explicit` |
| `src\factorforge\analysis\feasibility.py` | 127 | gc_global | `55.0` | `ranges = gc_ranges or [(40.0, 47.0), (35.0, 50.0), (55.0, 65.0)]` |
| `src\factorforge\analysis\feasibility.py` | 124 | gc_global | `65.0` | `# progressively wider windows. (55.0, 65.0) is retained ONLY as an explicit` |
| `src\factorforge\analysis\feasibility.py` | 127 | gc_global | `65.0` | `ranges = gc_ranges or [(40.0, 47.0), (35.0, 50.0), (55.0, 65.0)]` |
| `src\factorforge\engines\profile\construct_builder.py` | 211 | type_iis | `BsaI` | `"BsaI": ["GGTCTC", "GAGACC"],` |
| `src\factorforge\engines\profile\construct_builder.py` | 213 | type_iis | `BsmBI` | `"BsmBI": ["CGTCTC", "GAGACG"],` |
| `src\factorforge\engines\profile\construct_builder.py` | 212 | type_iis | `BpiI` | `"BpiI": ["GAAGAC", "GTCTTC"],` |
| `src\factorforge\engines\profile\exporter.py` | 448 | type_iis | `BsaI` | `"assembly_standard": "Golden Gate (BsaI)",` |
| `src\factorforge\engines\profile\exporter.py` | 451 | type_iis | `BsaI` | `"type": "BsaI site",` |
| `src\factorforge\engines\profile\scoring.py` | 31 | gc_global | `55.0` | `GC_RANGE_DEFAULT: tuple[float, float] = (55.0, 65.0)` |
| `src\factorforge\engines\profile\scoring.py` | 31 | gc_global | `65.0` | `GC_RANGE_DEFAULT: tuple[float, float] = (55.0, 65.0)` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 18 | type_iis | `BsaI` | `- Golden Gate (BsaI, BpiI, BsmBI)` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 19 | type_iis | `BsaI` | `- MoClo (BsaI + overhangs)` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 33 | type_iis | `BsaI` | `GOLDEN_GATE_ENZYMES: tuple[str, ...] = ("BsaI", "BpiI", "BsmBI")` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 38 | type_iis | `BsaI` | `"enzymes": ["BsaI", "BpiI", "BsmBI"],` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 40 | type_iis | `BsaI` | `"BsaI": ["GGTCTC", "GAGACC"],  # Forward and reverse complement` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 46 | type_iis | `BsaI` | `"enzymes": ["BsaI"],` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 47 | type_iis | `BsaI` | `"sites": {"BsaI": ["GGTCTC", "GAGACC"]},` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 103 | type_iis | `BsaI` | `[{'enzyme': 'BsaI', ...}]` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 389 | type_iis | `BsaI` | `# Test sequence (includes BsaI site)` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 18 | type_iis | `BsmBI` | `- Golden Gate (BsaI, BpiI, BsmBI)` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 33 | type_iis | `BsmBI` | `GOLDEN_GATE_ENZYMES: tuple[str, ...] = ("BsaI", "BpiI", "BsmBI")` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 38 | type_iis | `BsmBI` | `"enzymes": ["BsaI", "BpiI", "BsmBI"],` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 42 | type_iis | `BsmBI` | `"BsmBI": ["CGTCTC", "GAGACG"],` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 18 | type_iis | `BpiI` | `- Golden Gate (BsaI, BpiI, BsmBI)` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 27 | type_iis | `BpiI` | `# BpiI and BbsI share the same GAAGAC Type IIS recognition/cut behavior in` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 29 | type_iis | `BpiI` | `# production code and documentation consistently use BpiI as the canonical` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 33 | type_iis | `BpiI` | `GOLDEN_GATE_ENZYMES: tuple[str, ...] = ("BsaI", "BpiI", "BsmBI")` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 38 | type_iis | `BpiI` | `"enzymes": ["BsaI", "BpiI", "BsmBI"],` |
| `src\factorforge\engines\profile\rules\domesticator.py` | 41 | type_iis | `BpiI` | `"BpiI": ["GAAGAC", "GTCTTC"],` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 53 | type_iis | `BsaI` | `4. Assembly-Friendly: avoid BsaI/BpiI` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 125 | type_iis | `BsaI` | `"BsaI": ["GGTCTC", "GAGACC"],` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 501 | type_iis | `BsaI` | `- Retries up to max_attempts times until no BsaI/BpiI Type IIS` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 506 | type_iis | `BsaI` | `- Supported: BsaI/BpiI site avoidance via stochastic retry` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 127 | type_iis | `BsmBI` | `"BsmBI": ["CGTCTC", "GAGACG"],` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 53 | type_iis | `BpiI` | `4. Assembly-Friendly: avoid BsaI/BpiI` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 126 | type_iis | `BpiI` | `"BpiI": ["GAAGAC", "GTCTTC"],` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 501 | type_iis | `BpiI` | `- Retries up to max_attempts times until no BsaI/BpiI Type IIS` |
| `src\factorforge\engines\profile\rules\reverse_translator.py` | 506 | type_iis | `BpiI` | `- Supported: BsaI/BpiI site avoidance via stochastic retry` |
| `src\factorforge\registry\CHANGELOG.md` | 8 | type_iis | `BsaI` | `- before: `["BsaI", "BsmBI", "BbsI"]` -> after: `["BsaI", "BpiI", "BsmBI"]`` |
| `src\factorforge\registry\CHANGELOG.md` | 8 | type_iis | `BsmBI` | `- before: `["BsaI", "BsmBI", "BbsI"]` -> after: `["BsaI", "BpiI", "BsmBI"]`` |
| `src\factorforge\registry\CHANGELOG.md` | 8 | type_iis | `BpiI` | `- before: `["BsaI", "BsmBI", "BbsI"]` -> after: `["BsaI", "BpiI", "BsmBI"]`` |
| `src\factorforge\utils\restriction_sites.py` | 15 | type_iis | `BsaI` | `{"name": "BsaI", "sequence": "GGTCTC", "scan_rc": True},` |
| `src\factorforge\utils\restriction_sites.py` | 17 | type_iis | `BsmBI` | `{"name": "BsmBI", "sequence": "CGTCTC", "scan_rc": True},` |
| `src\factorforge\utils\restriction_sites.py` | 16 | type_iis | `BpiI` | `{"name": "BpiI", "sequence": "GAAGAC", "scan_rc": True},` |
| `docs\assembly-friendly.md` | 11 | type_iis | `BsaI` | `- BsaI and BpiI Type IIS restriction site avoidance through synonymous` |
| `docs\assembly-friendly.md` | 11 | type_iis | `BpiI` | `- BsaI and BpiI Type IIS restriction site avoidance through synonymous` |
| `docs\host-profile-registry.md` | 71 | cai | `0.82` | `min_feasible_cai: 0.82` |
| `docs\host-profile-registry.md` | 103 | cai | `0.82` | `min_feasible_cai: 0.82` |
| `docs\host-profile-registry.md` | 135 | cai | `0.82` | `min_feasible_cai: 0.82` |
| `docs\host-profile-registry.md` | 167 | cai | `0.82` | `min_feasible_cai: 0.82` |
| `docs\host-profile-registry.md` | 199 | cai | `0.82` | `min_feasible_cai: 0.82` |
| `docs\host-profile-registry.md` | 67 | gc_global | `55.0` | `ideal_band: [55.0, 65.0]` |
| `docs\host-profile-registry.md` | 99 | gc_global | `55.0` | `ideal_band: [55.0, 65.0]` |
| `docs\host-profile-registry.md` | 195 | gc_global | `55.0` | `ideal_band: [35.0, 55.0]` |
| `docs\host-profile-registry.md` | 67 | gc_global | `65.0` | `ideal_band: [55.0, 65.0]` |
| `docs\host-profile-registry.md` | 99 | gc_global | `65.0` | `ideal_band: [55.0, 65.0]` |
| `docs\host-profile-registry.md` | 131 | gc_global | `65.0` | `ideal_band: [45.0, 65.0]` |
| `docs\host-profile-registry.md` | 163 | gc_global | `65.0` | `ideal_band: [45.0, 65.0]` |
| `docs\host-profile-registry.md` | 68 | gc_local | `75.0` | `soft_band: [45.0, 75.0]` |
| `docs\host-profile-registry.md` | 100 | gc_local | `75.0` | `soft_band: [45.0, 75.0]` |
| `docs\how-it-works.md` | 25 | type_iis | `BsaI` | `BsaI / BsmBI recognition sites via silent edits` |
| `docs\how-it-works.md` | 25 | type_iis | `BsmBI` | `BsaI / BsmBI recognition sites via silent edits` |
| `docs\output.md` | 11 | type_iis | `BsaI` | `| **Domestication report** | BsaI/BsmBI and custom restriction sites removed, ed` |
| `docs\output.md` | 11 | type_iis | `BsmBI` | `| **Domestication report** | BsaI/BsmBI and custom restriction sites removed, ed` |
| `docs\profiles.md` | 26 | type_iis | `BsaI` | `| `assembly_friendly` | Golden Gate / MoClo workflows — avoids BsaI/BpiI Type II` |
| `docs\profiles.md` | 26 | type_iis | `BpiI` | `| `assembly_friendly` | Golden Gate / MoClo workflows — avoids BsaI/BpiI Type II` |
| `docs\scoring-calibration.md` | 11 | cai | `0.82` | `A typical lower bound is `0.82`.` |
| `docs\scoring-calibration.md` | 46 | gc_global | `55.0` | `width. *N. tabacum* (BY-2) keeps the prior `55.0` to `65.0` band, resolved` |
| `docs\scoring-calibration.md` | 46 | gc_global | `65.0` | `width. *N. tabacum* (BY-2) keeps the prior `55.0` to `65.0` band, resolved` |
| `docs\strategy-contract-registry.md` | 143 | type_iis | `BsaI` | `host_aware_note: Inherits balanced's host-aware codon selection; adds BsaI/BpiI ` |
| `docs\strategy-contract-registry.md` | 143 | type_iis | `BpiI` | `host_aware_note: Inherits balanced's host-aware codon selection; adds BsaI/BpiI ` |
| `docs\validation.md` | 74 | type_iis | `BsaI` | `| Forbidden restriction sites | Yes | BsaI, BsmBI, BpiI (Golden Gate); halts the` |
| `docs\validation.md` | 74 | type_iis | `BsmBI` | `| Forbidden restriction sites | Yes | BsaI, BsmBI, BpiI (Golden Gate); halts the` |
| `docs\validation.md` | 74 | type_iis | `BpiI` | `| Forbidden restriction sites | Yes | BsaI, BsmBI, BpiI (Golden Gate); halts the` |
| `docs\tutorials\gfp-nbenthamiana.md` | 162 | type_iis | `BsaI` | `- **MoClo Level 0** — use with the `assembly_friendly` profile; check for BsaI/B` |
| `docs\tutorials\gfp-nbenthamiana.md` | 162 | type_iis | `BpiI` | `- **MoClo Level 0** — use with the `assembly_friendly` profile; check for BsaI/B` |
| `docs\validation\RELEASE_GATE.md` | 34 | type_iis | `BsaI` | ``("BsaI", "BpiI", "BsmBI")`.` |
| `docs\validation\RELEASE_GATE.md` | 34 | type_iis | `BsmBI` | ``("BsaI", "BpiI", "BsmBI")`.` |
| `docs\validation\RELEASE_GATE.md` | 34 | type_iis | `BpiI` | ``("BsaI", "BpiI", "BsmBI")`.` |
| `tests\test_benchmark_regression.py` | 78 | type_iis | `BsaI` | `canonical = {"BsaI", "BpiI", "BsmBI"}` |
| `tests\test_benchmark_regression.py` | 78 | type_iis | `BsmBI` | `canonical = {"BsaI", "BpiI", "BsmBI"}` |
| `tests\test_benchmark_regression.py` | 78 | type_iis | `BpiI` | `canonical = {"BsaI", "BpiI", "BsmBI"}` |
| `tests\test_benchmark_scoring.py` | 52 | gc_global | `55.0` | `result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)` |
| `tests\test_benchmark_scoring.py` | 63 | gc_global | `55.0` | `result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)` |
| `tests\test_benchmark_scoring.py` | 74 | gc_global | `55.0` | `result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)` |
| `tests\test_benchmark_scoring.py` | 85 | gc_global | `55.0` | `result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)` |
| `tests\test_benchmark_scoring.py` | 97 | gc_global | `55.0` | `result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)` |
| `tests\test_benchmark_scoring.py` | 52 | gc_global | `65.0` | `result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)` |
| `tests\test_benchmark_scoring.py` | 63 | gc_global | `65.0` | `result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)` |
| `tests\test_benchmark_scoring.py` | 74 | gc_global | `65.0` | `result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)` |
| `tests\test_benchmark_scoring.py` | 85 | gc_global | `65.0` | `result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)` |
| `tests\test_benchmark_scoring.py` | 97 | gc_global | `65.0` | `result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)` |
| `tests\test_database.py` | 33 | cai | `0.82` | `"cai": 0.82,` |
| `tests\test_registry_production_sync.py` | 55 | type_iis | `BpiI` | `BpiI is the canonical FactorForge label for the GAAGAC target.` |
| `tests\test_registry_production_sync.py` | 61 | type_iis | `BpiI` | `"BbsI must not be in GOLDEN_GATE_ENZYMES — use BpiI (canonical label)"` |
| `tests\test_restriction_sites.py` | 193 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\test_restriction_sites.py` | 114 | type_iis | `BsaI` | `{"name": "BsaI", "sequence": "GGTCTC"},` |
| `tests\test_restriction_sites.py` | 182 | type_iis | `BsaI` | `assert result["removed_sites"][0]["enzyme"] == "BsaI"` |
| `tests\test_restriction_sites.py` | 206 | type_iis | `BsaI` | `hits = d.scan_restriction_sites("AAAGGTCTCAAA", "golden_gate")  # GGTCTC = BsaI` |
| `tests\test_validation_report.py` | 10 | gc_global | `55.0` | `_DEFAULT_CONSTRAINTS = {"gc_min": 55.0, "gc_max": 65.0}` |
| `tests\test_validation_report.py` | 10 | gc_global | `65.0` | `_DEFAULT_CONSTRAINTS = {"gc_min": 55.0, "gc_max": 65.0}` |
| `tests\test_worked_example.py` | 25 | gc_global | `55.0` | `GC_MIN = 55.0` |
| `tests\test_worked_example.py` | 26 | gc_global | `65.0` | `GC_MAX = 65.0` |
| `tests\api\test_optimize_contract.py` | 70 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 100 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 167 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 191 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 231 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 248 | gc_global | `55.0` | `constraints={"gc_min": 45.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 263 | gc_global | `55.0` | `{"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 280 | gc_global | `55.0` | `{"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 357 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
| `tests\api\test_optimize_contract.py` | 382 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |
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
| `tests\engines\profile\test_domesticator.py` | 45 | type_iis | `BpiI` | `# BsaI: GGTCTC, BpiI: GAAGAC` |
| `tests\engines\profile\test_exporter.py` | 251 | gc_global | `55.0` | `{"sequence": "ATGGCCTAG", "metadata": {"gene_name": "gene2", "cai": 0.9, "gc": 5` |
| `tests\engines\profile\test_exporter.py` | 264 | gc_global | `55.0` | `{"sequence": "ATGGCCTAG", "metadata": {"gene_name": "gene2", "cai": 0.9, "gc": 5` |
| `tests\engines\profile\test_exporter.py` | 39 | type_iis | `BsaI` | `"assembly_standard": "Golden Gate (BsaI)",` |
| `tests\engines\profile\test_exporter.py` | 42 | type_iis | `BsaI` | `"violations_fixed": [{"type": "BsaI site", "position": 147}],` |
| `tests\engines\profile\test_exporter.py` | 314 | type_iis | `BsaI` | `assert "BsaI site" in report` |
| `tests\engines\profile\test_pipeline_integration.py` | 110 | type_iis | `BsaI` | `"unfixable": [{"enzyme": "BsaI", "site": "GGTCTC", "position": 0}],` |
| `tests\engines\profile\test_reverse_translator.py` | 198 | gc_global | `55.0` | `assert 45.0 <= gc <= 55.0` |
| `tests\engines\profile\test_reverse_translator.py` | 160 | gc_global | `65.0` | `assert 35.0 <= gc <= 65.0` |
| `tests\engines\profile\test_reverse_translator.py` | 264 | type_iis | `BsaI` | `"""Test that assembly-friendly avoids BsaI sites"""` |
| `tests\engines\profile\test_reverse_translator.py` | 269 | type_iis | `BsaI` | `# BsaI recognition site` |
| `tests\engines\profile\test_scoring.py` | 154 | gc_global | `55.0` | `assert gc_band_score(60.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 157 | gc_global | `55.0` | `assert gc_band_score(55.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 160 | gc_global | `55.0` | `assert gc_band_score(65.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 172 | gc_global | `55.0` | `assert gc_band_score(85.0, 55.0, 65.0, decay_width=20.0) == 0.0` |
| `tests\engines\profile\test_scoring.py` | 175 | gc_global | `55.0` | `assert gc_band_score(10.0, 55.0, 65.0) == 0.0` |
| `tests\engines\profile\test_scoring.py` | 154 | gc_global | `65.0` | `assert gc_band_score(60.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 157 | gc_global | `65.0` | `assert gc_band_score(55.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 160 | gc_global | `65.0` | `assert gc_band_score(65.0, 55.0, 65.0) == 1.0` |
| `tests\engines\profile\test_scoring.py` | 172 | gc_global | `65.0` | `assert gc_band_score(85.0, 55.0, 65.0, decay_width=20.0) == 0.0` |
| `tests\engines\profile\test_scoring.py` | 175 | gc_global | `65.0` | `assert gc_band_score(10.0, 55.0, 65.0) == 0.0` |
| `tests\fixtures\contracts\validation_contract_v1.yaml` | 335 | type_iis | `BsaI` | `display_name: "Type IIS restriction sites (BsaI/BpiI/BsmBI)"` |
| `tests\fixtures\contracts\validation_contract_v1.yaml` | 340 | type_iis | `BsaI` | `threshold: "built-in sites: BsaI(GGTCTC), BpiI(GAAGAC), BsmBI(CGTCTC), scan_rc=T` |
| `tests\fixtures\contracts\validation_contract_v1.yaml` | 335 | type_iis | `BsmBI` | `display_name: "Type IIS restriction sites (BsaI/BpiI/BsmBI)"` |
| `tests\fixtures\contracts\validation_contract_v1.yaml` | 340 | type_iis | `BsmBI` | `threshold: "built-in sites: BsaI(GGTCTC), BpiI(GAAGAC), BsmBI(CGTCTC), scan_rc=T` |
| `tests\fixtures\contracts\validation_contract_v1.yaml` | 335 | type_iis | `BpiI` | `display_name: "Type IIS restriction sites (BsaI/BpiI/BsmBI)"` |
| `tests\fixtures\contracts\validation_contract_v1.yaml` | 340 | type_iis | `BpiI` | `threshold: "built-in sites: BsaI(GGTCTC), BpiI(GAAGAC), BsmBI(CGTCTC), scan_rc=T` |
| `tests\test_schemas\test_design_package.py` | 72 | gc_global | `55.0` | `constraints={"gc_min": 40.0, "gc_max": 55.0},` |

## Allowlisted (not requiring registry entry)

| file | line | match |
|------|------|-------|
| `tests\engines\profile\test_scoring.py` | 164 | `55.0` |
| `tests\engines\profile\test_scoring.py` | 168 | `55.0` |
| `tests\engines\profile\test_scoring.py` | 164 | `65.0` |
| `tests\engines\profile\test_scoring.py` | 168 | `65.0` |
| `tests\engines\profile\test_scoring.py` | 164 | `75.0` |