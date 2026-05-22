"""
FactorForge REST API — /api/optimize endpoint
Version: 3.0.0
Engine: FactorForge v2 (rule-based)
"""

from http.server import BaseHTTPRequestHandler
import json
import sys
import os
import re
import logging
from typing import Any

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

# Try to import FactorForge
try:
    from factorforge.engines import EngineRegistry
    from factorforge.engines.v3.metrics import load_codon_usage_table
    from factorforge.ml.feasibility import analyze_feasibility
    from factorforge.ml.metrics import (
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

    FACTORFORGE_AVAILABLE = True
    logger.info("FactorForge v2 engine loaded successfully")
except ImportError as e:
    FACTORFORGE_AVAILABLE = False
    logger.warning(f"FactorForge not available: {e}")

# Constants
MIN_SEQUENCE_LENGTH = 3
MAX_SEQUENCE_LENGTH = 50000  # 50kb max
VALID_PROFILES = ['balanced', 'high_cai', 'gc_target', 'assembly_friendly', 'ramp', 'viral_delivery']
VALID_OBJECTIVES = ['feasibility_best']
DEFAULT_OBJECTIVE = 'feasibility_best'
DEFAULT_HOST_PROFILE = 'nbenthamiana'
DEFAULT_GC_MIN = 40.0
DEFAULT_GC_MAX = 55.0
ENABLE_MOCK = os.environ.get("FACTORFORGE_ENABLE_MOCK", "false").lower() == "true"
ENGINE_VERSIONS = {
    'product': '3.0.0',
    'rule_engine': '2.5.3',
    'dp_engine': '3.0.0',
}
# Valid characters: ACGT (DNA) or standard 20 Amino Acids (Protein) + * (Stop)
VALID_AA = 'ACDEFGHIKLMNPQRSTVWY'
VALID_CHARS_PATTERN = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY*]+$', re.IGNORECASE)


class handler(BaseHTTPRequestHandler):
    """Vercel Serverless Function Handler"""

    def do_POST(self):
        """Handle POST requests"""
        try:
            # Parse request body
            content_length = int(self.headers.get('Content-Length', 0))
            body = self.rfile.read(content_length).decode('utf-8')
            data = json.loads(body)

            logger.info(f"Received optimization request: {len(data.get('sequence', ''))} bp")

            # Extract parameters
            sequence = data.get('sequence', '')
            profile = data.get('profile', 'balanced')
            objective = data.get('objective')
            legacy_profile_request = 'profile' in data and 'objective' not in data
            if objective is None and not legacy_profile_request:
                objective = DEFAULT_OBJECTIVE
            host_profile = data.get('host_profile', DEFAULT_HOST_PROFILE)
            return_candidates = bool(data.get('return_candidates', True))
            constraints = self.parse_constraints(data.get('constraints', {}))
            use_template = data.get('use_template', False)
            kozak = data.get('kozak', False)
            dinuc = data.get('dinuc', False)

            # Validate input
            validation_error = self.validate_input(sequence, profile, objective, host_profile, constraints)
            if validation_error:
                logger.warning(f"Validation error: {validation_error}")
                self.send_error_response(400, validation_error)
                return

            # Clean sequence
            sequence = self.clean_sequence(sequence)

            # Check if FactorForge is available
            if not FACTORFORGE_AVAILABLE:
                status_code, result = self.handle_unavailable_engine(
                    sequence, profile, kozak, dinuc, constraints, return_candidates
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
                    return_candidates=return_candidates,
                    constraints=constraints,
                )

            logger.info("Optimization completed successfully")
            self.send_json_response(200, result)

        except json.JSONDecodeError as e:
            logger.error(f"JSON decode error: {e}")
            self.send_error_response(400, 'Invalid JSON format')
        except ValueError as e:
            logger.error(f"Value error: {e}")
            self.send_error_response(400, str(e))
        except Exception as e:
            logger.error(f"Unexpected error: {e}", exc_info=True)
            self.send_error_response(500, 'Internal server error')

    def do_GET(self):
        """Handle GET requests (health check)"""
        health_info = {
            'status': 'healthy',
            'service': 'FactorForge API',
            'version': ENGINE_VERSIONS['product'],
            'codonforge_available': FACTORFORGE_AVAILABLE,
            'endpoints': {
                'POST /api/optimize': 'Run codon optimization',
                'GET /api/optimize': 'Health check'
            },
            'supported_profiles': VALID_PROFILES,
            'supported_objectives': VALID_OBJECTIVES,
            'mock_enabled': ENABLE_MOCK,
            'engine_versions': ENGINE_VERSIONS,
        }

        if FACTORFORGE_AVAILABLE:
            try:
                optimizer = EngineRegistry.get('v2')
                health_info['engine'] = {
                    'name': optimizer.name,
                    'version': optimizer.version
                }
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
            gc_min = float(constraints.get('gc_min', DEFAULT_GC_MIN))
            gc_max = float(constraints.get('gc_max', DEFAULT_GC_MAX))
        except (TypeError, ValueError):
            raise ValueError('constraints.gc_min and constraints.gc_max must be numeric')
        return {'gc_min': gc_min, 'gc_max': gc_max}

    def validate_input(self, sequence, profile, objective=None, host_profile=DEFAULT_HOST_PROFILE, constraints=None):
        """Validate input parameters"""
        # Check sequence exists
        if not sequence or not sequence.strip():
            return 'Sequence is required'

        # Clean and check length
        cleaned = self.clean_sequence(sequence)

        if len(cleaned) < MIN_SEQUENCE_LENGTH:
            return f'Sequence must be at least {MIN_SEQUENCE_LENGTH} bp'

        if len(cleaned) > MAX_SEQUENCE_LENGTH:
            return f'Sequence must be less than {MAX_SEQUENCE_LENGTH} bp'

        # Check Valid Characters (DNA or Protein)
        if not VALID_CHARS_PATTERN.match(cleaned):
            invalid_chars = set(cleaned) - set(VALID_AA + '*')
            return f'Sequence contains invalid characters: {", ".join(sorted(invalid_chars))}'

        # Check profile
        if profile not in VALID_PROFILES:
            return f'Invalid profile. Must be one of: {", ".join(VALID_PROFILES)}'

        if objective is not None and objective not in VALID_OBJECTIVES:
            return f'Invalid objective. Must be one of: {", ".join(VALID_OBJECTIVES)}'

        if host_profile != DEFAULT_HOST_PROFILE:
            return f'Unsupported host_profile. Must be: {DEFAULT_HOST_PROFILE}'

        constraints = constraints or {'gc_min': DEFAULT_GC_MIN, 'gc_max': DEFAULT_GC_MAX}
        if constraints['gc_min'] > constraints['gc_max']:
            return 'constraints.gc_min must be <= constraints.gc_max'

        return None

    def clean_sequence(self, sequence):
        """Clean sequence: remove whitespace, FASTA headers, convert to uppercase"""
        # Remove FASTA headers (lines starting with >)
        lines = sequence.split('\n')
        sequence_lines = [line for line in lines if not line.startswith('>')]

        # Join and remove all whitespace
        cleaned = ''.join(sequence_lines)
        cleaned = ''.join(cleaned.split())

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
        return_candidates=False,
        constraints=None,
    ):
        """Run actual FactorForge v2 optimization"""
        try:
            constraints = constraints or {'gc_min': DEFAULT_GC_MIN, 'gc_max': DEFAULT_GC_MAX}
            if objective == DEFAULT_OBJECTIVE:
                return self.optimize_feasibility_best(
                    sequence=sequence,
                    profile=profile,
                    host_profile=host_profile,
                    constraints=constraints,
                    kozak=kozak,
                    dinuc=dinuc,
                    return_candidates=return_candidates,
                )

            # Get v2 optimizer
            optimizer = EngineRegistry.get('v2')
            logger.info(f"Using FactorForge v2 engine: {optimizer.name} {optimizer.version}")

            # Run optimization
            result = optimizer.optimize(
                sequence=sequence,
                profile=profile,
                kozak=kozak,
                dinuc=dinuc
            )

            # Build construct if requested
            if use_template:
                logger.info("Building construct with standard_expression template")
                try:
                    from factorforge.engines.v2 import ConstructBuilder
                    builder = ConstructBuilder(template='standard_expression')
                    construct = builder.build(result.sequence)
                    optimized_sequence = str(construct.seq)
                except Exception as e:
                    logger.warning(f"Construct building failed: {e}, using raw sequence")
                    optimized_sequence = result.sequence
            else:
                optimized_sequence = result.sequence

            # Extract metrics safely
            cai = float(result.metrics.get('cai', 0.0))
            gc_percent = float(result.metrics.get('gc_percent', 0.0))
            polya_warnings = int(result.metrics.get('polya_warnings', 0))

            # Validation checks
            polya_check = 'PASS' if polya_warnings == 0 else 'WARNING'
            gc_check = self.gc_check(gc_percent, constraints)
            moclo_check = 'PASS'  # From v2 domesticator

            logger.info(f"Optimization metrics: CAI={cai:.3f}, GC={gc_percent:.1f}%, PolyA={polya_warnings}")

            # Format response
            response = {
                'success': True,
                'optimized_sequence': optimized_sequence,
                'original_length': len(sequence),
                'optimized_length': len(optimized_sequence),
                'metrics': {
                    'cai': round(cai, 3),
                    'gc_percent': round(gc_percent, 1),
                    'polya_signals': polya_warnings,
                    'length': len(optimized_sequence)
                },
                'profile': profile,
                'use_template': use_template,
                'validation': {
                    'polya': polya_check,
                    'moclo': moclo_check,
                    'gc': gc_check
                },
                'engine': {
                    'name': optimizer.name,
                    'version': optimizer.version
                }
            }
            if return_candidates:
                table = load_codon_usage_table()
                response['recommended_candidate'] = self.build_candidate(
                    candidate_id=profile,
                    label=self.candidate_label(profile),
                    dna_sequence=result.sequence,
                    codon_weights=table.codon_weights,
                    recommendation_reason=f'Backward-compatible v2 {profile} profile result',
                )
                response['candidates'] = [response['recommended_candidate']]
                response['engine_versions'] = ENGINE_VERSIONS
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
        return_candidates=True,
    ):
        """Run v1 feasibility_best contract and add v2 comparison candidates."""
        table = load_codon_usage_table()
        aa_seq = self.clean_sequence(sequence).rstrip('*')
        optimizer = EngineRegistry.get('v2')

        feasibility = analyze_feasibility(
            protein_sequence=aa_seq,
            codon_weights=table.codon_weights,
            target_gc_low=constraints['gc_min'],
            target_gc_high=constraints['gc_max'],
        )
        best = feasibility['target']['best_candidate'] or feasibility['best_candidate_without_gc']
        if not best:
            raise ValueError('No feasibility_best candidate generated')

        candidates = [
            self.build_candidate(
                candidate_id='feasibility_best',
                label='Feasibility Best',
                dna_sequence=best['dna_sequence'],
                codon_weights=table.codon_weights,
                recommendation_reason=(
                    f"Maximum CAI under GC {constraints['gc_min']:g}-{constraints['gc_max']:g}%"
                    if feasibility['target']['best_candidate']
                    else 'Maximum CAI without GC constraint; requested GC range was infeasible'
                ),
            )
        ]

        for candidate_profile in ('gc_target', 'high_cai'):
            result = optimizer.optimize(
                sequence=aa_seq,
                profile=candidate_profile,
                kozak=kozak,
                dinuc=dinuc,
            )
            candidates.append(
                self.build_candidate(
                    candidate_id=candidate_profile,
                    label=self.candidate_label(candidate_profile),
                    dna_sequence=result.sequence,
                    codon_weights=table.codon_weights,
                    recommendation_reason=f'v2 {candidate_profile} comparison candidate',
                )
            )

        response = {
            'success': True,
            'recommended_candidate': candidates[0],
            'candidates': candidates if return_candidates else [],
            'validation': {
                'input_type': 'protein',
                'sequence_length': len(aa_seq),
                'host_profile': host_profile,
            },
            'engine_versions': ENGINE_VERSIONS,
        }

        return response

    def build_candidate(
        self,
        candidate_id: str,
        label: str,
        dna_sequence: str,
        codon_weights: dict[str, float],
        recommendation_reason: str,
    ) -> dict[str, Any]:
        """Build a v1 candidate payload with evidence metrics."""
        windows = calculate_gc_windows(dna_sequence)
        window_values = [float(window['gc']) for window in windows]
        first_region = calculate_first_region_gc(dna_sequence, region_sizes=[30])
        internal_stop_count = count_internal_stops(dna_sequence)
        invalid_codon_count = len(detect_invalid_codons(dna_sequence))

        return {
            'id': candidate_id,
            'label': label,
            'dna_sequence': dna_sequence,
            'cai': round(calculate_cai(dna_sequence, codon_weights), 3),
            'gc_percent': round(calculate_gc(dna_sequence), 1),
            'gc_window_min': round(min(window_values), 1) if window_values else 0.0,
            'gc_window_max': round(max(window_values), 1) if window_values else 0.0,
            'first_region_gc': round(float(first_region['first_30nt_gc']), 1),
            'internal_stop_count': internal_stop_count,
            'invalid_codon_count': invalid_codon_count,
            'repeat_count': len(detect_repeats(dna_sequence)),
            'homopolymer_count': len(detect_homopolymers(dna_sequence)),
            'forbidden_motif_count': len(detect_forbidden_motifs(dna_sequence, [])),
            'validator_status': (
                'pass' if internal_stop_count == 0 and invalid_codon_count == 0 else 'fail'
            ),
            'recommendation_reason': recommendation_reason,
        }

    def candidate_label(self, candidate_id):
        """Return human-readable candidate label."""
        labels = {
            'feasibility_best': 'Feasibility Best',
            'gc_target': 'GC Target',
            'high_cai': 'High CAI',
            'balanced': 'Balanced',
            'assembly_friendly': 'Assembly Friendly',
            'ramp': 'Ramp',
            'viral_delivery': 'Viral Delivery',
        }
        return labels.get(candidate_id, candidate_id.replace('_', ' ').title())

    def gc_check(self, gc_percent, constraints):
        """Return PASS/WARNING for configured global GC constraints."""
        return 'PASS' if constraints['gc_min'] <= gc_percent <= constraints['gc_max'] else 'WARNING'

    def handle_unavailable_engine(
        self,
        sequence,
        profile,
        kozak,
        dinuc,
        constraints,
        return_candidates,
    ):
        """Return mock only when explicitly enabled; otherwise fail closed."""
        if ENABLE_MOCK:
            logger.info("Using mock optimization (FactorForge not available)")
            return 200, self.generate_mock_result(
                sequence, profile, kozak, dinuc, constraints, return_candidates
            )

        logger.error("FactorForge engine unavailable and mock fallback disabled")
        return 503, {
            'success': False,
            'error': 'Engine unavailable. Contact support.',
        }

    def generate_mock_result(self, sequence, profile, kozak, dinuc, constraints=None, return_candidates=False):
        """Generate mock result for testing (when FactorForge is not available)"""
        logger.warning("Generating mock result - CodonForge engine not available")
        constraints = constraints or {'gc_min': DEFAULT_GC_MIN, 'gc_max': DEFAULT_GC_MAX}

        # Determine if input is DNA or Protein
        is_protein = not re.match(r'^[ACGT]+$', sequence)

        if is_protein:
            # Mock reverse-translation using N. benthamiana preferred codons
            mock_map = {
                'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA',
                'F': 'TTT', 'G': 'GGA', 'H': 'CAT', 'I': 'ATT',
                'K': 'AAG', 'L': 'CTT', 'M': 'ATG', 'N': 'AAT',
                'P': 'CCA', 'Q': 'CAA', 'R': 'AGA', 'S': 'TCT',
                'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT',
                '*': 'TAA',
            }
            optimized = ''.join([mock_map.get(c, 'NNN') for c in sequence])
        else:
            # Is DNA
            optimized = sequence.upper()

        # Calculate real GC content
        gc_count = optimized.count('G') + optimized.count('C')
        gc_percent = (gc_count / len(optimized)) * 100 if len(optimized) > 0 else 0

        # Mock CAI based on profile
        mock_cai_values = {
            'balanced': 0.850,
            'high_cai': 0.920,
            'gc_target': 0.800,
            'assembly_friendly': 0.830,
            'ramp': 0.870
        }
        mock_cai = mock_cai_values.get(profile, 0.850)

        response = {
            'success': True,
            'optimized_sequence': optimized,
            'original_length': len(sequence),
            'optimized_length': len(optimized),
            'metrics': {
                'cai': mock_cai,
                'gc_percent': round(gc_percent, 1),
                'polya_signals': 0,
                'length': len(optimized)
            },
            'profile': profile,
            'use_template': False,
            'kozak': kozak,
            'dinuc': dinuc,
            'validation': {
                'polya': 'PASS',
                'gc': self.gc_check(gc_percent, constraints),
                'moclo': 'PASS'
            },
            'engine': {
                'name': 'Mock Engine',
                'version': '0.0.0'
            },
            'note': '⚠️ Mock data - FactorForge engine not available. Deploy with FactorForge for real optimization.'
        }
        if return_candidates:
            response['engine_versions'] = ENGINE_VERSIONS
        return response

    def send_json_response(self, status_code, data):
        """Send JSON response with CORS headers"""
        self.send_response(status_code)
        self.send_header('Content-Type', 'application/json')
        self.send_cors_headers()
        self.end_headers()
        self.wfile.write(json.dumps(data).encode('utf-8'))

    def send_error_response(self, status_code, message):
        """Send error response"""
        self.send_json_response(status_code, {
            'success': False,
            'error': message
        })

    def send_cors_headers(self):
        """Send CORS headers"""
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
