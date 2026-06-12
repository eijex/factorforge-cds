METHOD_ORDER = [
    "random_synonymous",
    "greedy_cai",
    "native_reference",
    "factorforge_balanced",
    "factorforge_gc_target",
    "factorforge_high_cai",
    "factorforge_assembly_friendly",
]

DISPLAY_NAMES = {
    "random_synonymous": "Random\nsynonymous",
    "greedy_cai": "Greedy\nCAI",
    "native_reference": "Native\nreference",
    "factorforge_balanced": "FF\nbalanced",
    "factorforge_gc_target": "FF\nGC target",
    "factorforge_high_cai": "FF\nhigh CAI",
    "factorforge_assembly_friendly": "FF\nassembly-friendly",
}

BASELINE_METHODS = {"random_synonymous", "greedy_cai", "native_reference"}

# Colors: baselines light gray, FactorForge blue, assembly_friendly dark blue
BAR_COLORS = {
    "random_synonymous": "#BDBDBD",
    "greedy_cai": "#BDBDBD",
    "native_reference": "#BDBDBD",
    "factorforge_balanced": "#1976D2",
    "factorforge_gc_target": "#1976D2",
    "factorforge_high_cai": "#1976D2",
    "factorforge_assembly_friendly": "#0D47A1",
}

REPLICATED_METHOD = "random_synonymous"
FAILURE_CATEGORIES = [
    "passed_all",
    "failed_bio",
    "failed_gc_only",
    "failed_iis_only",
    "failed_gc_and_iis",
    "failed_other",
]
