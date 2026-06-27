"""
FactorForge REST API — /api/optimize endpoint
Product Version: 3.2.6
Default objective: feasibility_best (DP feasibility / constraint-based CDS design)
Profile comparison engine: constraint-aware rule-based profiles
"""

from http.server import BaseHTTPRequestHandler
import hashlib
import json
import sys
import os
import re
import logging
from datetime import datetime, timezone
from typing import Any

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

# Try to import FactorForge
try:
    from factorforge.engines import EngineRegistry
    from factorforge.engines.profile.rules.domesticator import Domesticator
    from factorforge.engines.profile.rules.rule_engine import RuleEngine
    from factorforge.engines.profile.utils import get_data_path, load_codon_table
    from factorforge.analysis.metrics import load_codon_usage_table
    from factorforge.analysis.feasibility import analyze_feasibility
    from factorforge.analysis.metrics import (
        calculate_cai,
        calculate_first_region_gc,
        calculate_gc,
        calculate_gc_windows,
        count_internal_stops,
        detect_forbidden_motifs,
        detect_homopolymers,
        detect_invalid_codons,
        detect_repeats,
    )
    from factorforge.utils.restriction_sites import (
        detect_restriction_sites,
        domesticate_custom_sites,
    )
    from factorforge.utils.sequence_validator import validate_cds_output
    from factorforge.validation_registry import VALIDATION_REGISTRY_VERSION, public_badge_checks
    from factorforge.validation_report import (
        VALIDATION_REPORT_SCHEMA_VERSION,
        build_validation_report,
    )

    FACTORFORGE_AVAILABLE = True
    logger.info("FactorForge v3.x profile engine loaded successfully")
except ImportError as e:
    FACTORFORGE_AVAILABLE = False
    logger.warning(f"FactorForge not available: {e}")

# Constants
MIN_SEQUENCE_LENGTH = 3
# Web API length limits (Vercel serverless timeout constraint). Split by input
# type so the unit matches the message: 5,000 aa ≈ 15,000 bp.
MAX_PROTEIN_LENGTH_AA = 5000
MAX_DNA_LENGTH_BP = MAX_PROTEIN_LENGTH_AA * 3  # 15,000 bp
VALID_PROFILES = [
    "balanced",
    "high_cai",
    "gc_target",
    "assembly_friendly",
]
DEFAULT_COMPARE_PROFILES = [
    "balanced",
    "high_cai",
    "gc_target",
    "assembly_friendly",
]
MAX_COMPARE_PROFILES = 6
MAX_BATCH_SEQUENCES = 20
VALID_OBJECTIVES = ["feasibility_best"]
DEFAULT_OBJECTIVE = "feasibility_best"
DEFAULT_HOST_PROFILE = "nbenthamiana"
VALID_HOSTS = ["nbenthamiana", "by2"]
HOST_MAP = {"nbenthamiana": "nbenthamiana", "by2": "ntabacum"}
# Single source of truth for host display metadata, consumed by web/js/app.js
# via GET /api/optimize so web/index.html never hardcodes host cards.
HOST_METADATA = {
    "nbenthamiana": {
        "display_name": "N. benthamiana",
        "description": "Leaf agroinfiltration design host",
    },
    "by2": {
        "display_name": "Tobacco BY-2",
        "description": "N. tabacum cell culture codon usage",
        "caveat": (
            "Experimental: species-level proxy table, not wet-lab validated "
            "for BY-2 expression performance."
        ),
    },
}
DEFAULT_GC_MIN = 55.0
DEFAULT_GC_MAX = 65.0
ENABLE_MOCK = os.environ.get("FACTORFORGE_ENABLE_MOCK", "false").lower() == "true"
ENGINE_VERSIONS = {
    "product": "3.2.6",
    "rule_engine": "3.2.6",
    "dp_engine": "3.2.6",
}
# Valid characters: ACGT (DNA) or standard 20 Amino Acids (Protein) + * (Stop)
VALID_AA = "ACDEFGHIKLMNPQRSTVWY"
VALID_CHARS_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY*]+$", re.IGNORECASE)


def _generate_construct_id() -> str:
    now = datetime.now()
    return f"CF-{now.strftime('%Y%m%d-%H%M%S')}"


class handler(BaseHTTPRequestHandler):
    """Vercel Serverless Function Handler"""

    def do_POST(self):
        """Handle POST requests"""
        try:
            # Parse request body
            content_length = int(self.headers.get("Content-Length", 0))
            body = self.rfile.read(content_length).decode("utf-8")
            data = json.loads(body)

            logger.info(
                f"Received optimization request: sequence_length={len(data.get('sequence', ''))}"
            )

            request_path = self.path.split("?", 1)[0]
            if request_path == "/api/optimize/compare":
                status_code, result = self.handle_compare_request(data)
                self.send_json_response(status_code, result)
                return
            if request_path == "/api/optimize/batch":
                status_code, result = self.handle_batch_request(data)
                self.send_json_response(status_code, result)
                return

            # Extract parameters
            sequence = data.get("sequence", "")
            profile = data.get("profile", "balanced")
            host = self.validate_host(data.get("host", DEFAULT_HOST_PROFILE))
            internal_host = HOST_MAP[host]
            objective = data.get("objective")
            legacy_profile_request = "profile" in data and "objective" not in data
            if objective is None and not legacy_profile_request:
                if "host" in data and internal_host != DEFAULT_HOST_PROFILE:
                    objective = None
                else:
                    objective = DEFAULT_OBJECTIVE
            host_profile = data.get("host_profile", host)

            # Explicit strategy/host compatibility guard. Reject
            # rather than silently emitting nbenthamiana-table output under
            # a different host's name. feasibility_best (DP engine, hardcoded
            # table) and high_cai (nbenthamiana-only golden-set reference)
            # are both N. benthamiana-only by current design.
            if internal_host != DEFAULT_HOST_PROFILE:
                if data.get("objective") == "feasibility_best":
                    self.send_error_response(
                        400,
                        {
                            "error": (
                                "objective=feasibility_best is only supported "
                                "with host=nbenthamiana"
                            ),
                            "error_code": "UNSUPPORTED_STRATEGY_HOST_COMBINATION",
                            "requested_host": internal_host,
                            "requested_strategy": "feasibility_best",
                        },
                    )
                    return
                if data.get("profile") == "high_cai":
                    self.send_error_response(
                        400,
                        {
                            "error": (
                                "high_cai requires the N. benthamiana golden-set "
                                f"reference and is not available for host={internal_host}"
                            ),
                            "error_code": "UNSUPPORTED_STRATEGY_HOST_COMBINATION",
                            "requested_host": internal_host,
                            "requested_strategy": "high_cai",
                        },
                    )
                    return

            # Implicit case (host-only request, no explicit objective/profile)
            # already resolves to a host-supported strategy via the
            # objective/profile defaulting above — this just discloses that
            # a substitution happened.
            implicit_strategy_disclosure = None
            if (
                internal_host != DEFAULT_HOST_PROFILE
                and "objective" not in data
                and "profile" not in data
            ):
                implicit_strategy_disclosure = {
                    "requested_strategy": "feasibility_best",
                    "resolved_strategy": "balanced",
                    "resolution_reason": (
                        "feasibility_best is not available for this host; "
                        "resolved to a host-supported strategy"
                    ),
                }

            return_candidates = bool(data.get("return_candidates", True))
            constraints = self.parse_constraints(data.get("constraints", {}))
            use_template = data.get("use_template", False)
            kozak = data.get("kozak", False)
            dinuc = data.get("dinuc", False)
            custom_restriction_sites = self.parse_custom_restriction_sites(
                data.get("custom_restriction_sites")
            )

            # Validate input
            validation_error = self.validate_input(
                sequence, profile, objective, host_profile, constraints
            )
            if validation_error:
                logger.warning(f"Validation error: {validation_error}")
                self.send_error_response(400, validation_error)
                return

            # Clean sequence
            sequence = self.clean_sequence(sequence)

            # Check if FactorForge is available
            if not FACTORFORGE_AVAILABLE:
                status_code, result = self.handle_unavailable_engine(
                    sequence,
                    profile,
                    kozak,
                    dinuc,
                    constraints,
                    return_candidates,
                    host_profile,
                )
                self.send_json_response(status_code, result)
                return
            else:
                logger.info(
                    "Running real optimization: "
                    f"profile={profile}, objective={objective}, template={use_template}, "
                    f"kozak={kozak}, dinuc={dinuc}"
                )
                result = self.optimize_sequence(
                    sequence,
                    profile,
                    use_template,
                    kozak,
                    dinuc,
                    objective=objective,
                    host_profile=host_profile,
                    host=internal_host,
                    return_candidates=return_candidates,
                    constraints=constraints,
                    custom_restriction_sites=custom_restriction_sites,
                )

            if implicit_strategy_disclosure and isinstance(result, dict):
                result.update(implicit_strategy_disclosure)

            logger.info("Optimization completed successfully")
            self.send_json_response(200, result)

        except json.JSONDecodeError as e:
            logger.error(f"JSON decode error: {e}")
            self.send_error_response(400, "Invalid JSON format")
        except ValueError as e:
            logger.error(f"Value error: {e}")
            self.send_error_response(400, str(e))
        except Exception as e:
            logger.error(f"Unexpected error: {e}", exc_info=True)
            self.send_error_response(500, "Internal server error")

    def do_GET(self):
        """Handle GET requests (health check)"""
        health_info = {
            "status": "healthy",
            "service": "FactorForge API",
            "version": ENGINE_VERSIONS["product"],
            "factorforge_available": FACTORFORGE_AVAILABLE,
            "endpoints": {
                "POST /api/optimize": "Run codon optimization",
                "POST /api/optimize/compare": "Compare profile optimization results",
                "POST /api/optimize/batch": "Run batch profile optimization",
                "GET /api/optimize": "Health check",
            },
            "supported_profiles": VALID_PROFILES,
            "supported_hosts": VALID_HOSTS,
            "host_metadata": HOST_METADATA,
            "supported_objectives": VALID_OBJECTIVES,
            "mock_enabled": ENABLE_MOCK,
            "engine_versions": ENGINE_VERSIONS,
            "validation_registry_version": VALIDATION_REGISTRY_VERSION,
            "validation_report_schema_version": VALIDATION_REPORT_SCHEMA_VERSION,
            "validation_checks": public_badge_checks(),
        }

        if FACTORFORGE_AVAILABLE:
            try:
                optimizer = EngineRegistry.get("profile")
                health_info["engine"] = {"name": optimizer.name, "version": optimizer.version}
            except Exception:
                pass

        logger.info("Health check requested")
        self.send_json_response(200, health_info)

    def do_OPTIONS(self):
        """Handle OPTIONS requests (CORS preflight)"""
        self.send_response(200)
        self.send_cors_headers()
        self.end_headers()

    def parse_constraints(self, constraints):
        """Parse v1 constraints with defaults."""
        constraints = constraints or {}
        try:
            gc_min = float(constraints.get("gc_min", DEFAULT_GC_MIN))
            gc_max = float(constraints.get("gc_max", DEFAULT_GC_MAX))
        except (TypeError, ValueError):
            raise ValueError("constraints.gc_min and constraints.gc_max must be numeric")
        return {"gc_min": gc_min, "gc_max": gc_max}

    def parse_custom_restriction_sites(self, custom_sites):
        """Parse and validate optional custom restriction-site definitions."""
        if not custom_sites:
            return []

        if not isinstance(custom_sites, list):
            raise ValueError("custom_restriction_sites must be a list")

        parsed = []
        for index, site in enumerate(custom_sites):
            if not isinstance(site, dict):
                raise ValueError(f"custom_restriction_sites[{index}] must be an object")

            name = str(site.get("name", "")).strip()
            sequence = str(site.get("sequence", "")).strip().upper()
            scan_rc = bool(site.get("scan_rc", True))
            parsed.append({"name": name, "sequence": sequence, "scan_rc": scan_rc})

        detect_restriction_sites("", parsed)
        return parsed

    def validate_host(self, host):
        """Validate public host flag and return the normalized public name."""
        host_value = str(host or DEFAULT_HOST_PROFILE).strip().lower()
        if host_value not in VALID_HOSTS:
            raise ValueError(f"Invalid host: {host_value}")
        return host_value

    def validate_input(
        self, sequence, profile, objective=None, host_profile=DEFAULT_HOST_PROFILE, constraints=None
    ):
        """Validate input parameters"""
        # Check sequence exists
        if not sequence or not sequence.strip():
            return "Sequence is required"

        # Clean and check length
        cleaned = self.clean_sequence(sequence)

        if len(cleaned) < MIN_SEQUENCE_LENGTH:
            return f"Sequence must be at least {MIN_SEQUENCE_LENGTH} bp"

        length_error = self._length_limit_error(cleaned)
        if length_error:
            return {
                "error": length_error,
                "cli_install": "pip install factorforge-cds",
                "docker_image": "ghcr.io/eijex/factorforge-cds:latest",
            }

        # Check Valid Characters (DNA or Protein)
        if not VALID_CHARS_PATTERN.match(cleaned):
            invalid_chars = set(cleaned) - set(VALID_AA + "*")
            return f"Sequence contains invalid characters: {', '.join(sorted(invalid_chars))}"

        # Check profile
        if profile not in VALID_PROFILES:
            return f"Invalid profile. Must be one of: {', '.join(VALID_PROFILES)}"

        if objective is not None and objective not in VALID_OBJECTIVES:
            return f"Invalid objective. Must be one of: {', '.join(VALID_OBJECTIVES)}"

        host_profile_value = str(host_profile or DEFAULT_HOST_PROFILE).lower()
        if host_profile_value not in VALID_HOSTS and host_profile_value not in HOST_MAP.values():
            return f"Unsupported host_profile. Must be one of: {', '.join(VALID_HOSTS)}"

        constraints = constraints or {"gc_min": DEFAULT_GC_MIN, "gc_max": DEFAULT_GC_MAX}
        if constraints["gc_min"] > constraints["gc_max"]:
            return "constraints.gc_min must be <= constraints.gc_max"

        return None

    def validate_compare_profiles(self, profiles):
        """Validate and normalize profile comparison list."""
        if profiles is None:
            return list(DEFAULT_COMPARE_PROFILES)
        if not isinstance(profiles, list):
            raise ValueError("profiles must be a list")
        if not profiles:
            raise ValueError("profiles must include at least one profile")
        if len(profiles) > MAX_COMPARE_PROFILES:
            raise ValueError(f"profiles must include at most {MAX_COMPARE_PROFILES} profiles")

        normalized = []
        for profile in profiles:
            profile_name = str(profile).strip()
            if profile_name not in VALID_PROFILES:
                raise ValueError(f"Invalid profile: {profile_name}")
            normalized.append(profile_name)
        return normalized

    def _length_limit_error(self, cleaned):
        """Return a web-API length-limit error message, or None if within limits.

        Branches by input type so the unit matches the message: DNA (only
        ACGTU/N) is capped in bp; everything else is treated as protein and
        capped in aa. An ambiguous short ACGT-only protein falls under the
        larger (DNA) limit and is unaffected in practice.
        """
        is_dna = bool(cleaned) and set(cleaned) <= set("ACGTUN")
        if is_dna:
            if len(cleaned) > MAX_DNA_LENGTH_BP:
                return (
                    f"DNA input exceeds maximum length for web API "
                    f"({MAX_DNA_LENGTH_BP:,} bp, approximately {MAX_PROTEIN_LENGTH_AA:,} aa). "
                    "For longer sequences, use the CLI or Docker."
                )
        elif len(cleaned) > MAX_PROTEIN_LENGTH_AA:
            return (
                f"Protein input exceeds maximum length for web API "
                f"(max {MAX_PROTEIN_LENGTH_AA:,} amino acids). "
                "For longer sequences, use the CLI or Docker."
            )
        return None

    def validate_compare_sequence(self, sequence):
        """Validate and clean a profile comparison input sequence."""
        if not sequence or not str(sequence).strip():
            raise ValueError("Sequence is required")

        cleaned = self.clean_sequence(str(sequence))
        if len(cleaned) < MIN_SEQUENCE_LENGTH:
            raise ValueError(f"Sequence must be at least {MIN_SEQUENCE_LENGTH} bp")
        length_error = self._length_limit_error(cleaned)
        if length_error:
            raise ValueError(length_error)
        if not VALID_CHARS_PATTERN.match(cleaned):
            invalid_chars = set(cleaned) - set(VALID_AA + "*")
            raise ValueError(
                f"Sequence contains invalid characters: {', '.join(sorted(invalid_chars))}"
            )
        return cleaned

    def handle_compare_request(self, data):
        """Handle POST /api/optimize/compare requests."""
        try:
            if "host" in data or "host_profile" in data:
                return 400, {
                    "success": False,
                    "error": (
                        "host selection is not supported on this endpoint; "
                        "use POST /api/optimize for host-specific design"
                    ),
                    "error_code": "HOST_NOT_SUPPORTED_ON_ENDPOINT",
                }
            sequence = self.validate_compare_sequence(data.get("sequence", ""))
            profiles = self.validate_compare_profiles(data.get("profiles"))
            scan_mode = str(data.get("scan_mode", "fast")).lower()
            if scan_mode not in {"fast", "full"}:
                raise ValueError("scan_mode must be one of: fast, full")

            if not FACTORFORGE_AVAILABLE:
                logger.error("FactorForge engine unavailable for profile comparison")
                return 503, {"success": False, "error": "Engine unavailable. Contact support."}

            result = self.optimize_profile_comparison(sequence, profiles, scan_mode)
            return 200, result

        except ValueError as e:
            logger.warning(f"Compare validation error: {e}")
            return 400, {"success": False, "error": str(e)}

    def optimize_profile_comparison(self, sequence, profiles, scan_mode="fast"):
        """Run profile optimization sequentially and return compact comparison rows."""
        optimizer = EngineRegistry.get("profile")
        results = []

        for profile in profiles:
            result = optimizer.optimize(
                sequence=sequence,
                profile=profile,
                scan_mode=scan_mode,
            )
            results.append(
                {
                    "profile": profile,
                    "cai": round(float(result.metrics.get("cai", 0.0)), 3),
                    "gc_percent": round(
                        float(result.metrics.get("gc_percent", result.metrics.get("gc_content", 0.0))),
                        2,
                    ),
                    "score": round(float(result.metrics.get("score", 0.0)), 3),
                    "sequence": result.sequence,
                }
            )

        return {"results": results}

    def validate_batch_sequences(self, sequences):
        """Validate and normalize batch optimization sequence entries."""
        if not isinstance(sequences, list) or not sequences:
            raise ValueError("sequences is required and must be non-empty")
        if len(sequences) > MAX_BATCH_SEQUENCES:
            raise ValueError("Batch limit is 20 sequences")

        normalized = []
        for index, entry in enumerate(sequences, start=1):
            if isinstance(entry, dict):
                sequence_id = str(entry.get("id") or f"seq_{index}").strip() or f"seq_{index}"
                raw_sequence = entry.get("sequence", "")
            else:
                sequence_id = f"seq_{index}"
                raw_sequence = entry

            sequence = self.validate_compare_sequence(raw_sequence)
            normalized.append({"id": sequence_id, "sequence": sequence})
        return normalized

    def handle_batch_request(self, data):
        """Handle POST /api/optimize/batch requests."""
        try:
            if "host" in data or "host_profile" in data:
                return 400, {
                    "success": False,
                    "error": (
                        "host selection is not supported on this endpoint; "
                        "use POST /api/optimize for host-specific design"
                    ),
                    "error_code": "HOST_NOT_SUPPORTED_ON_ENDPOINT",
                }
            profile = str(data.get("profile", "balanced")).strip()
            if profile not in VALID_PROFILES:
                raise ValueError(f"Invalid profile. Must be one of: {', '.join(VALID_PROFILES)}")

            scan_mode = str(data.get("scan_mode", "fast")).lower()
            if scan_mode not in {"fast", "full"}:
                raise ValueError("scan_mode must be one of: fast, full")

            sequences = self.validate_batch_sequences(data.get("sequences"))

            if not FACTORFORGE_AVAILABLE:
                logger.error("FactorForge engine unavailable for batch optimization")
                return 503, {"success": False, "error": "Engine unavailable. Contact support."}

            result = self.optimize_batch_sequences(sequences, profile, scan_mode)
            return 200, result

        except ValueError as e:
            logger.warning(f"Batch validation error: {e}")
            return 400, {"success": False, "error": str(e)}

    def optimize_batch_sequences(self, sequences, profile, scan_mode="fast"):
        """Run profile optimization sequentially for a batch of input sequences."""
        optimizer = EngineRegistry.get("profile")
        results = []

        for entry in sequences:
            result = optimizer.optimize(
                sequence=entry["sequence"],
                profile=profile,
                scan_mode=scan_mode,
            )
            results.append(
                {
                    "id": entry["id"],
                    "sequence": result.sequence,
                    "cai": round(float(result.metrics.get("cai", 0.0)), 3),
                    "gc_percent": round(
                        float(result.metrics.get("gc_percent", result.metrics.get("gc_content", 0.0))),
                        2,
                    ),
                    "score": round(float(result.metrics.get("score", 0.0)), 3),
                    "violations": int(result.metrics.get("violations", 0)),
                }
            )

        return {"results": results, "count": len(results), "profile": profile}

    def clean_sequence(self, sequence):
        """Clean sequence: remove whitespace, FASTA headers, convert to uppercase"""
        # Remove FASTA headers (lines starting with >)
        lines = sequence.split("\n")
        sequence_lines = [line for line in lines if not line.startswith(">")]

        # Join and remove all whitespace
        cleaned = "".join(sequence_lines)
        cleaned = "".join(cleaned.split())

        # Convert to uppercase
        return cleaned.upper()

    def optimize_sequence(
        self,
        sequence,
        profile,
        use_template,
        kozak,
        dinuc,
        objective=None,
        host_profile=DEFAULT_HOST_PROFILE,
        host=DEFAULT_HOST_PROFILE,
        return_candidates=False,
        constraints=None,
        custom_restriction_sites=None,
    ):
        """Run actual FactorForge v3.x profile optimization."""
        try:
            constraints = constraints or {"gc_min": DEFAULT_GC_MIN, "gc_max": DEFAULT_GC_MAX}
            if objective == DEFAULT_OBJECTIVE:
                return self.optimize_feasibility_best(
                    sequence=sequence,
                    profile=profile,
                    host_profile=host_profile,
                    host=host,
                    constraints=constraints,
                    kozak=kozak,
                    dinuc=dinuc,
                    return_candidates=return_candidates,
                    custom_restriction_sites=custom_restriction_sites,
                )

            # Get profile-based optimizer
            optimizer = EngineRegistry.get("profile")
            logger.info(
                f"Using FactorForge v3.x profile engine: {optimizer.name} {optimizer.version}"
            )

            # Run optimization
            result = optimizer.optimize(
                sequence=sequence,
                profile=profile,
                host=host,
                kozak=kozak,
                dinuc=dinuc,
            )

            # Build construct if requested
            if use_template:
                logger.info("Building construct with standard_expression template")
                try:
                    from factorforge.engines.profile import ConstructBuilder

                    builder = ConstructBuilder(template="standard_expression")
                    construct = builder.build(result.sequence)
                    optimized_sequence = str(construct.seq)
                except Exception as e:
                    logger.warning(f"Construct building failed: {e}, using raw sequence")
                    optimized_sequence = result.sequence
            else:
                optimized_sequence = result.sequence

            # Extract metrics safely
            cai = float(result.metrics.get("cai", 0.0))
            gc_percent = float(result.metrics.get("gc_percent", 0.0))
            polya_warnings = int(result.metrics.get("polya_warnings", 0))
            table = load_codon_usage_table()
            general_cai = calculate_cai(result.sequence, table.codon_weights)
            # MFE provenance (016 audit): surface whether MFE was actually
            # computed. ViennaRNA is not installed on Vercel, so MFE is normally
            # not_computed in production — never report it as a misleading 0.0.
            mfe_kcal_mol = result.metrics.get("mfe_kcal_mol")
            mfe_status = result.metrics.get("mfe_status", "not_computed")
            mfe_used = bool(result.metrics.get("mfe_used", False))
            mfe_warning = result.metrics.get("mfe_warning")

            # Validation checks
            polya_check = "PASS" if polya_warnings == 0 else "WARNING"
            gc_check = self.gc_check(gc_percent, constraints)
            # NOTE: this checks Type IIS restriction sites (BsaI/BsmBI/BpiI),
            # not MoClo overhang validity. The "moclo" JSON key is kept for
            # frontend backward compatibility; the displayed UI label was
            # corrected to "Restriction Site Check".
            type_iis_sites = Domesticator().scan_restriction_sites(result.sequence, "golden_gate")
            restriction_site_check = "PASS" if not type_iis_sites else "WARNING"
            moclo_check = restriction_site_check

            logger.info(
                f"Optimization metrics: CAI={cai:.3f}, GC={gc_percent:.1f}%, PolyA={polya_warnings}"
            )

            # Format response
            response = {
                "success": True,
                "optimized_sequence": optimized_sequence,
                "original_length": len(sequence),
                "optimized_length": len(optimized_sequence),
                "metrics": {
                    "cai": round(cai, 3),
                    "cai_reference": "profile_golden_set",
                    "general_cai": round(general_cai, 3),
                    "gc_percent": round(gc_percent, 1),
                    "polya_signals": polya_warnings,
                    "length": len(optimized_sequence),
                    "mfe_kcal_mol": (
                        round(float(mfe_kcal_mol), 2) if mfe_kcal_mol is not None else None
                    ),
                    "mfe_status": mfe_status,
                    "mfe_used": mfe_used,
                    "mfe_warning": mfe_warning,
                },
                "profile": profile,
                "use_template": use_template,
                "validation": {"polya": polya_check, "moclo": moclo_check, "gc": gc_check},
                "engine": {"name": optimizer.name, "version": optimizer.version},
            }
            validation_report = build_validation_report(
                optimized_sequence,
                gc_percent=gc_percent,
                constraints=constraints,
                rule_engine=RuleEngine(host=host),
                moclo_requested=use_template,
            )
            response["validation"]["schema_version"] = validation_report["schema_version"]
            response["validation"]["checks"] = validation_report["checks"]
            response.setdefault("metadata", {})["validation_registry_version"] = (
                VALIDATION_REGISTRY_VERSION
            )
            if return_candidates:
                response["recommended_candidate"] = self.build_candidate(
                    candidate_id=profile,
                    label=self.candidate_label(profile),
                    dna_sequence=result.sequence,
                    codon_weights=table.codon_weights,
                    profile_cai=cai,
                    recommendation_reason=f"Profile engine {profile} result",
                    constraints=constraints,
                    host=host,
                    moclo_requested=use_template,
                )
                response["candidates"] = [response["recommended_candidate"]]
                response["engine_versions"] = ENGINE_VERSIONS
            response = self.apply_custom_restriction_sites(
                response,
                custom_restriction_sites,
                constraints=constraints,
                host=host,
            )
            response = self.add_design_package_fields(
                response=response,
                input_sequence=sequence,
                profile=profile,
                objective=objective,
                host_profile=host_profile,
                kozak=kozak,
                dinuc=dinuc,
                constraints=constraints,
            )
            return response

        except Exception as e:
            logger.error(f"Optimization failed: {e}", exc_info=True)
            raise ValueError(f"Optimization error: {str(e)}")

    def optimize_feasibility_best(
        self,
        sequence,
        profile,
        host_profile,
        constraints,
        kozak,
        dinuc,
        host=DEFAULT_HOST_PROFILE,
        return_candidates=True,
        custom_restriction_sites=None,
    ):
        """Run feasibility_best contract and add profile comparison candidates."""
        table = load_codon_usage_table()
        aa_seq = self.clean_sequence(sequence).rstrip("*")
        optimizer = EngineRegistry.get("profile")

        feasibility = analyze_feasibility(
            protein_sequence=aa_seq,
            codon_weights=table.codon_weights,
            target_gc_low=constraints["gc_min"],
            target_gc_high=constraints["gc_max"],
        )
        best = feasibility["target"]["best_candidate"] or feasibility["best_candidate_without_gc"]
        if not best:
            raise ValueError("No feasibility_best candidate generated")

        candidates = [
            self.build_candidate(
                candidate_id="feasibility_best",
                label="Feasibility Best",
                dna_sequence=best["dna_sequence"],
                codon_weights=table.codon_weights,
                profile_cai=float(best["cai"]),
                recommendation_reason=(
                    f"Maximum CAI under GC {constraints['gc_min']:g}-{constraints['gc_max']:g}%"
                    if feasibility["target"]["best_candidate"]
                    else "Maximum CAI without GC constraint; requested GC range was infeasible"
                ),
                constraints=constraints,
                host=host,
            )
        ]

        for candidate_profile in ("gc_target", "high_cai"):
            result = optimizer.optimize(
                sequence=aa_seq,
                profile=candidate_profile,
                host=host,
                kozak=kozak,
                dinuc=dinuc,
            )
            candidates.append(
                self.build_candidate(
                    candidate_id=candidate_profile,
                    label=self.candidate_label(candidate_profile),
                    dna_sequence=result.sequence,
                    codon_weights=table.codon_weights,
                    profile_cai=float(result.metrics.get("cai", 0.0)),
                    recommendation_reason=f"Profile engine {candidate_profile} comparison candidate",
                    constraints=constraints,
                    host=host,
                )
            )

        response = {
            "success": True,
            "recommended_candidate": candidates[0],
            "candidates": candidates if return_candidates else [],
            "validation": {
                "input_type": "protein",
                "sequence_length": len(aa_seq),
                "host_profile": host_profile,
            },
            "engine_versions": ENGINE_VERSIONS,
        }

        primary_dna = candidates[0]["dna_sequence"]
        validation_report = build_validation_report(
            primary_dna,
            gc_percent=float(candidates[0]["gc_percent"]),
            constraints=constraints,
            rule_engine=RuleEngine(host=host),
        )
        response["validation_report"] = validation_report
        response.setdefault("metadata", {})["validation_registry_version"] = (
            VALIDATION_REGISTRY_VERSION
        )

        response = self.apply_custom_restriction_sites(
            response, custom_restriction_sites, constraints=constraints, host=host
        )
        return self.add_design_package_fields(
            response=response,
            input_sequence=sequence,
            profile=profile,
            objective=DEFAULT_OBJECTIVE,
            host_profile=host_profile,
            kozak=kozak,
            dinuc=dinuc,
            constraints=constraints,
        )

    def add_design_package_fields(
        self,
        response,
        input_sequence,
        profile,
        objective,
        host_profile,
        kozak,
        dinuc,
        constraints,
    ):
        """Add DesignPackage-compatible metadata while preserving existing response keys."""
        output_cds = self.primary_dna_sequence(response)
        selected_profile = self.response_profile(response, profile, objective)
        metrics = response.get("metrics", {})
        recommended = response.get("recommended_candidate") or {}
        cai = metrics.get("cai", recommended.get("cai", 0.0))
        gc_percent = metrics.get("gc_percent", recommended.get("gc_percent", 0.0))
        polya_warnings = int(metrics.get("polya_signals", 0))
        internal_stop_count = int(recommended.get("internal_stop_count", 0))
        cds_validation = self.cds_validation_result(input_sequence, output_cds)
        aa_identity = float(cds_validation.get("aa_identity", 0.0))
        codon_rarity_clusters = self.count_rare_codon_runs(output_cds, host_profile)

        param_payload = {
            "objective": objective or "legacy_profile",
            "profile": selected_profile,
            "host_profile": host_profile,
            "kozak": kozak,
            "dinuc": dinuc,
            "constraints": constraints,
        }
        param_str = json.dumps(param_payload, sort_keys=True, separators=(",", ":"))

        response["construct_id"] = _generate_construct_id()
        response["design_package_version"] = "1.0"
        response["created_at"] = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
        response["host_profile"] = host_profile
        response["profile"] = selected_profile
        response["provenance"] = {
            "input_sequence_hash": self.sha256_prefix(input_sequence),
            "output_cds_hash": self.sha256_prefix(output_cds),
            "parameter_hash": self.sha256_prefix(param_str),
        }
        response["wet_lab_feedback"] = {"status": "pending", "submissions": []}
        response.setdefault("target", None)
        response.setdefault(
            "construct_plan",
            {
                "construct_type": "full_length",
                "tag": None,
                "signal_peptide": None,
                "kozak": bool(kozak),
            },
        )
        response["cds_design"] = {
            "engine": "factorforge_cds",
            "product_version": ENGINE_VERSIONS["product"],
            "host_profile": host_profile,
            "objective": objective or "legacy_profile",
            "profile": selected_profile,
            "input_length_aa": self.input_length_aa(input_sequence),
            "output_length_nt": len(output_cds),
            "cai": float(cai),
            "gc_percent": float(gc_percent),
        }
        response["constraint_report"] = {
            "restriction_sites_removed": response.get("custom_restriction_sites", {}).get(
                "removed", []
            ),
            "restriction_sites_unresolved": response.get("custom_restriction_sites", {}).get(
                "unresolved", []
            ),
            "polya_warnings": polya_warnings,
            "internal_stop_count": internal_stop_count,
            "aa_identity": aa_identity,
            "codon_rarity_clusters": codon_rarity_clusters,
            "cds_validation_errors": cds_validation.get("errors", []),
        }
        response["validation_status"] = self.design_validation_status(response)
        return response

    def cds_validation_result(self, input_sequence, output_cds):
        """Return validate_cds_output() status for protein inputs."""
        cleaned = self.clean_sequence(input_sequence)
        if re.fullmatch(r"[ACGT]+", cleaned):
            return {"passed": True, "errors": [], "aa_identity": 1.0}
        return validate_cds_output(cleaned, output_cds)

    def count_rare_codon_runs(self, output_cds, host_profile):
        """Count rare codon runs using the host-specific rule scanner."""
        internal_host = HOST_MAP.get(str(host_profile or DEFAULT_HOST_PROFILE).lower(), host_profile)
        return len(RuleEngine(host=internal_host).scan_rare_codon_runs(output_cds))

    def response_profile(self, response, profile, objective):
        """Return the selected candidate/profile name for DesignPackage metadata."""
        if objective == DEFAULT_OBJECTIVE:
            recommended = response.get("recommended_candidate")
            if isinstance(recommended, dict) and recommended.get("id"):
                return recommended["id"]
            return DEFAULT_OBJECTIVE
        return profile

    def sha256_prefix(self, value):
        """Return the shortened sha256 digest format used in API provenance."""
        return "sha256:" + hashlib.sha256(str(value).encode()).hexdigest()[:16]

    def input_length_aa(self, sequence):
        """Return amino-acid length when the input does not look like DNA."""
        cleaned = self.clean_sequence(sequence).rstrip("*")
        if re.fullmatch(r"[ACGT]+", cleaned):
            return None
        return len(cleaned)

    def design_validation_status(self, response):
        """Map existing validation fields into DesignPackage validation status."""
        validation = response.get("validation", {})
        constraint_report = response.get("constraint_report", {})
        recommended = response.get("recommended_candidate") or {}
        validator_status = recommended.get("validator_status")
        aa_identity = float(constraint_report.get("aa_identity", 0.0))
        cds_errors = constraint_report.get("cds_validation_errors", [])
        polya = str(validation.get("polya", "UNCHECKED")).lower()
        gc = str(validation.get("gc", "UNCHECKED")).lower()
        moclo = str(validation.get("moclo", "UNCHECKED")).lower()
        aa_identity_check = "pass" if aa_identity == 1.0 and not cds_errors else "fail"
        in_silico = (
            "pass"
            if validator_status in (None, "pass") and gc != "warning" and aa_identity_check == "pass"
            else "warning"
        )
        return {
            "in_silico": in_silico,
            "aa_identity_check": aa_identity_check,
            "gc_check": "pass" if gc == "pass" else gc,
            "polya_check": "pass" if polya == "pass" else polya,
            "moclo_check": "unchecked" if moclo == "unchecked" else moclo,
        }

    def apply_custom_restriction_sites(
        self,
        response,
        custom_restriction_sites,
        constraints=None,
        host=DEFAULT_HOST_PROFILE,
    ):
        """Apply custom restriction-site domestication to the primary response CDS."""
        if not custom_restriction_sites:
            return response

        dna_sequence = self.primary_dna_sequence(response)
        usage_table = load_codon_usage_table()
        codon_table = load_codon_table(DEFAULT_HOST_PROFILE, get_data_path())
        before_metrics = self.custom_site_metrics(dna_sequence, usage_table.codon_weights)

        domestication = domesticate_custom_sites(
            dna_sequence,
            custom_restriction_sites,
            codon_table,
        )
        after_sequence = domestication["sequence"]
        after_metrics = self.custom_site_metrics(after_sequence, usage_table.codon_weights)

        self.update_primary_dna_sequence(
            response,
            dna_sequence,
            after_sequence,
            usage_table.codon_weights,
            constraints=constraints,
            host=host,
        )

        metrics = response.setdefault("metrics", {})
        metrics["before"] = before_metrics
        metrics["after"] = after_metrics
        if "cai" in metrics:
            metrics["cai"] = round(after_metrics["cai"], 3)
        if "gc_percent" in metrics:
            metrics["gc_percent"] = round(after_metrics["gc"], 1)
        if "length" in metrics:
            metrics["length"] = len(after_sequence)

        response["custom_restriction_sites"] = {
            "detected": domestication["detected"],
            "removed": domestication["removed"],
            "unresolved": domestication["unresolved"],
        }
        return response

    def primary_dna_sequence(self, response):
        """Return the primary DNA sequence from either legacy or v1 response shape."""
        if response.get("optimized_sequence"):
            return response["optimized_sequence"]

        candidate = response.get("recommended_candidate")
        if isinstance(candidate, dict) and candidate.get("dna_sequence"):
            return candidate["dna_sequence"]

        raise ValueError("No primary DNA sequence available for custom restriction-site cleanup")

    def update_primary_dna_sequence(
        self,
        response,
        before_sequence,
        after_sequence,
        codon_weights,
        constraints=None,
        host=DEFAULT_HOST_PROFILE,
    ):
        """Update primary DNA fields and candidate evidence after custom domestication."""
        if response.get("optimized_sequence") == before_sequence:
            response["optimized_sequence"] = after_sequence
            response["optimized_length"] = len(after_sequence)

        recommended = response.get("recommended_candidate")
        if isinstance(recommended, dict):
            self.update_candidate_sequence(
                recommended,
                before_sequence,
                after_sequence,
                codon_weights,
                constraints=constraints,
                host=host,
            )

        for candidate in response.get("candidates", []):
            if isinstance(candidate, dict):
                self.update_candidate_sequence(
                    candidate,
                    before_sequence,
                    after_sequence,
                    codon_weights,
                    constraints=constraints,
                    host=host,
                )

    def update_candidate_sequence(
        self,
        candidate,
        before_sequence,
        after_sequence,
        codon_weights,
        constraints=None,
        host=DEFAULT_HOST_PROFILE,
    ):
        """Update one candidate when it points at the primary DNA sequence."""
        if candidate.get("dna_sequence") != before_sequence:
            return

        updated = self.build_candidate(
            candidate_id=candidate.get("id", "custom_domesticated"),
            label=candidate.get("label", "Custom Domesticated"),
            dna_sequence=after_sequence,
            codon_weights=codon_weights,
            recommendation_reason=candidate.get(
                "recommendation_reason",
                "Custom restriction-site domesticated candidate",
            ),
            constraints=constraints or {"gc_min": DEFAULT_GC_MIN, "gc_max": DEFAULT_GC_MAX},
            host=host,
        )
        candidate.update(updated)

    def custom_site_metrics(self, dna_sequence, codon_weights):
        """Return compact before/after metrics for custom site domestication."""
        return {
            "cai": round(calculate_cai(dna_sequence, codon_weights), 3),
            "gc": round(calculate_gc(dna_sequence), 1),
        }

    def build_candidate(
        self,
        candidate_id: str,
        label: str,
        dna_sequence: str,
        codon_weights: dict[str, float],
        recommendation_reason: str,
        profile_cai: float | None = None,
        constraints: dict[str, float] | None = None,
        host: str = DEFAULT_HOST_PROFILE,
        moclo_requested: bool = False,
    ) -> dict[str, Any]:
        """Build a v1 candidate payload with evidence metrics."""
        windows = calculate_gc_windows(dna_sequence)
        window_values = [float(window["gc"]) for window in windows]
        first_region = calculate_first_region_gc(dna_sequence, region_sizes=[30])
        internal_stop_count = count_internal_stops(dna_sequence)
        invalid_codon_count = len(detect_invalid_codons(dna_sequence))
        general_cai = round(calculate_cai(dna_sequence, codon_weights), 3)
        type_iis_sites = Domesticator().scan_restriction_sites(dna_sequence, "golden_gate")
        assembly_pass = len(type_iis_sites) == 0
        gc_percent = calculate_gc(dna_sequence)
        active_constraints = constraints or {"gc_min": DEFAULT_GC_MIN, "gc_max": DEFAULT_GC_MAX}
        validation_report = build_validation_report(
            dna_sequence,
            gc_percent=gc_percent,
            constraints=active_constraints,
            rule_engine=RuleEngine(host=host),
            moclo_requested=moclo_requested,
        )

        return {
            "id": candidate_id,
            "label": label,
            "dna_sequence": dna_sequence,
            "cai": round(profile_cai, 3) if profile_cai is not None else general_cai,
            "cai_reference": (
                "profile_golden_set" if profile_cai is not None else "configured_codon_table"
            ),
            "general_cai": general_cai,
            "gc_percent": round(gc_percent, 1),
            "gc_window_min": round(min(window_values), 1) if window_values else 0.0,
            "gc_window_max": round(max(window_values), 1) if window_values else 0.0,
            "first_region_gc": round(float(first_region["first_30nt_gc"]), 1),
            "internal_stop_count": internal_stop_count,
            "invalid_codon_count": invalid_codon_count,
            "repeat_count": len(detect_repeats(dna_sequence)),
            "homopolymer_count": len(detect_homopolymers(dna_sequence)),
            "forbidden_motif_count": len(detect_forbidden_motifs(dna_sequence, [])),
            "forbidden_type_iis_site_count": len(type_iis_sites),
            "assembly_pass": assembly_pass,
            "validator_status": (
                "pass"
                if internal_stop_count == 0 and invalid_codon_count == 0 and assembly_pass
                else "warning"
            ),
            "recommendation_reason": recommendation_reason,
            "checks": validation_report["checks"],
        }

    def candidate_label(self, candidate_id):
        """Return human-readable candidate label."""
        labels = {
            "feasibility_best": "Feasibility Best",
            "gc_target": "GC Target",
            "high_cai": "High CAI",
            "balanced": "Balanced",
            "assembly_friendly": "Assembly Friendly",
        }
        return labels.get(candidate_id, candidate_id.replace("_", " ").title())

    def gc_check(self, gc_percent, constraints):
        """Return PASS/WARNING for configured global GC constraints."""
        return "PASS" if constraints["gc_min"] <= gc_percent <= constraints["gc_max"] else "WARNING"

    def handle_unavailable_engine(
        self,
        sequence,
        profile,
        kozak,
        dinuc,
        constraints,
        return_candidates,
        host_profile=DEFAULT_HOST_PROFILE,
    ):
        """Return mock only when explicitly enabled; otherwise fail closed."""
        if ENABLE_MOCK:
            logger.info("Using mock optimization (FactorForge not available)")
            return 200, self.generate_mock_result(
                sequence, profile, kozak, dinuc, constraints, return_candidates, host_profile
            )

        logger.error("FactorForge engine unavailable and mock fallback disabled")
        return 503, {
            "success": False,
            "error": "Engine unavailable. Contact support.",
        }

    def generate_mock_result(
        self,
        sequence,
        profile,
        kozak,
        dinuc,
        constraints=None,
        return_candidates=False,
        host_profile=DEFAULT_HOST_PROFILE,
    ):
        """Generate mock result for testing (when FactorForge is not available)"""
        logger.warning("Generating mock result - FactorForge engine not available")
        constraints = constraints or {"gc_min": DEFAULT_GC_MIN, "gc_max": DEFAULT_GC_MAX}

        # Determine if input is DNA or Protein
        is_protein = not re.match(r"^[ACGT]+$", sequence)

        if is_protein:
            # Mock reverse-translation using N. benthamiana preferred codons
            mock_map = {
                "A": "GCT",
                "C": "TGT",
                "D": "GAT",
                "E": "GAA",
                "F": "TTT",
                "G": "GGA",
                "H": "CAT",
                "I": "ATT",
                "K": "AAG",
                "L": "CTT",
                "M": "ATG",
                "N": "AAT",
                "P": "CCA",
                "Q": "CAA",
                "R": "AGA",
                "S": "TCT",
                "T": "ACT",
                "V": "GTT",
                "W": "TGG",
                "Y": "TAT",
                "*": "TAA",
            }
            optimized = "".join([mock_map.get(c, "NNN") for c in sequence])
        else:
            # Is DNA
            optimized = sequence.upper()

        # Calculate real GC content
        gc_count = optimized.count("G") + optimized.count("C")
        gc_percent = (gc_count / len(optimized)) * 100 if len(optimized) > 0 else 0

        # Mock CAI based on profile
        mock_cai_values = {
            "balanced": 0.850,
            "high_cai": 0.920,
            "gc_target": 0.800,
            "assembly_friendly": 0.830,
        }
        mock_cai = mock_cai_values.get(profile, 0.850)

        response = {
            "success": True,
            "optimized_sequence": optimized,
            "original_length": len(sequence),
            "optimized_length": len(optimized),
            "metrics": {
                "cai": mock_cai,
                "gc_percent": round(gc_percent, 1),
                "polya_signals": 0,
                "length": len(optimized),
            },
            "profile": profile,
            "use_template": False,
            "kozak": kozak,
            "dinuc": dinuc,
            "validation": {
                "polya": "PASS",
                "gc": self.gc_check(gc_percent, constraints),
                "moclo": "UNCHECKED",
            },
            "engine": {"name": "Mock Engine", "version": "0.0.0"},
            "note": "⚠️ Mock data - FactorForge engine not available. Deploy with FactorForge for real optimization.",
        }
        if return_candidates:
            response["engine_versions"] = ENGINE_VERSIONS
        return self.add_design_package_fields(
            response=response,
            input_sequence=sequence,
            profile=profile,
            objective=None,
            host_profile=host_profile,
            kozak=kozak,
            dinuc=dinuc,
            constraints=constraints,
        )

    def send_json_response(self, status_code, data):
        """Send JSON response with CORS headers"""
        self.send_response(status_code)
        self.send_header("Content-Type", "application/json")
        self.send_cors_headers()
        self.end_headers()
        self.wfile.write(json.dumps(data).encode("utf-8"))

    def send_error_response(self, status_code, message):
        """Send error response"""
        if isinstance(message, dict):
            payload = {"success": False}
            payload.update(message)
            self.send_json_response(status_code, payload)
            return

        self.send_json_response(status_code, {"success": False, "error": message})

    def send_cors_headers(self):
        """Send CORS headers"""
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
