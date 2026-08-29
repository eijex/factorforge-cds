import re

with open("src/factorforge/engines/profile/optimizer.py", "r") as f:
    code = f.read()

target = """        # 4. Rule checks (PolyA, etc.)
        scan_mode = str(kwargs.get("scan_mode", "full"))
        scan_include = kwargs.get("scan_include")
        scan_exclude = kwargs.get("scan_exclude")
        scan_results = rule_engine.scan_all(
            optimized_dna,
            mode=scan_mode,
            include=scan_include,
            exclude=scan_exclude,
        )

        # 5. Build result
        metrics = {
            "cai": candidates[0]["cai"],
            "cai_authority": dict(translator.cai_authority),
            # Keep both names for compatibility across existing tests/callers.
            "gc_content": candidates[0]["gc"],
            "gc_percent": candidates[0]["gc"],
            "score": candidates[0]["score"],
            "violations": sum(len(v) for v in scan_results.values()),
        }
        if profile_value == "balanced":
            host_gc_min, host_gc_max = resolve_host_gc_range(host)
            requested_gc_min_percent = float(kwargs.get("target_gc_min", host_gc_min))
            requested_gc_max_percent = float(kwargs.get("target_gc_max", host_gc_max))
            achieved_gc_percent = metrics["gc_percent"]
            metrics.update(
                {
                    "gc_target_reached": (
                        requested_gc_min_percent <= achieved_gc_percent <= requested_gc_max_percent
                    ),
                    "requested_gc_min_percent": requested_gc_min_percent,
                    "requested_gc_max_percent": requested_gc_max_percent,
                }
            )
        # MFE provenance: expose whether MFE was actually computed so downstream
        # artifacts (API response, Design Package) never report an uncomputed
        # MFE as a misleading 0.0 (016 audit). Score value is unchanged.
        metrics.update(compute_mfe_evidence(optimized_dna, profile=profile_value))

        return OptimizationResult(
            sequence=optimized_dna,
            metrics=metrics,
            metadata={
                "engine": "profile",
                "profile": profile_value,
                "host": host,
                "scan_mode": scan_mode,
                "scan_results": scan_results,
            },
        )"""

replacement = """        # 4. Shared Evaluator checks
        from factorforge.evaluation.evaluator import SharedEvaluator
        
        terminal_stop_policy = kwargs.get("terminal_stop_policy", "preserve")
        expected_protein = protein if seq_type == "dna" else processed_seq
        
        evaluator = SharedEvaluator(
            version="1.0.1", 
            codon_weights=translator.cai_authority.get("weights")
        )
        
        target_gc_min = kwargs.get("target_gc_min", None)
        target_gc_max = kwargs.get("target_gc_max", None)
        forbidden_type_iis = kwargs.get("forbidden_type_iis", set())
        
        if profile_value == "balanced":
            host_gc_min, host_gc_max = resolve_host_gc_range(host)
            target_gc_min = float(target_gc_min if target_gc_min is not None else host_gc_min)
            target_gc_max = float(target_gc_max if target_gc_max is not None else host_gc_max)
            
        eval_result = evaluator.evaluate_candidate(
            candidate_dna=optimized_dna,
            expected_protein=expected_protein,
            candidate_id="profile-candidate-01",
            target_gc_min=target_gc_min,
            target_gc_max=target_gc_max,
            forbidden_type_iis=forbidden_type_iis,
        )

        scan_mode = str(kwargs.get("scan_mode", "full"))
        scan_include = kwargs.get("scan_include")
        scan_exclude = kwargs.get("scan_exclude")
        scan_results = rule_engine.scan_all(
            optimized_dna,
            mode=scan_mode,
            include=scan_include,
            exclude=scan_exclude,
        )

        # 5. Build result
        metrics = {
            "cai": eval_result.metrics.cai or candidates[0]["cai"],
            "cai_authority": dict(translator.cai_authority),
            "gc_content": eval_result.metrics.gc_percent,
            "gc_percent": eval_result.metrics.gc_percent,
            "score": candidates[0]["score"],
            "violations": sum(len(v) for v in scan_results.values()),
        }
        if profile_value == "balanced":
            achieved_gc_percent = metrics["gc_percent"]
            metrics.update(
                {
                    "gc_target_reached": (
                        target_gc_min <= achieved_gc_percent <= target_gc_max
                    ),
                    "requested_gc_min_percent": target_gc_min,
                    "requested_gc_max_percent": target_gc_max,
                }
            )
        metrics.update(compute_mfe_evidence(optimized_dna, profile=profile_value))

        return OptimizationResult(
            sequence=optimized_dna,
            metrics=metrics,
            metadata={
                "engine": "profile",
                "profile": profile_value,
                "host": host,
                "scan_mode": scan_mode,
                "scan_results": scan_results,
                "validator_passed": eval_result.passed,
                "evaluation_report": eval_result.model_dump(),
            },
        )"""

if target in code:
    with open("src/factorforge/engines/profile/optimizer.py", "w") as f:
        f.write(code.replace(target, replacement))
    print("Patched successfully")
else:
    print("Target not found")
