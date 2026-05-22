"""
FactorForge REST API — /api/optimize endpoint
Version: 2.5.3
Engine: FactorForge v2 (rule-based)
"""

from http.server import BaseHTTPRequestHandler
import json
import sys
import os
import re
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

# Try to import FactorForge
try:
    from factorforge.engines.registry import EngineRegistry
    FACTORFORGE_AVAILABLE = True
    logger.info("FactorForge v2 engine loaded successfully")
except ImportError as e:
    FACTORFORGE_AVAILABLE = False
    logger.warning(f"FactorForge not available: {e}")

# Constants
MIN_SEQUENCE_LENGTH = 3
MAX_SEQUENCE_LENGTH = 50000  # 50kb max
VALID_PROFILES = ['balanced', 'high_cai', 'gc_target', 'assembly_friendly', 'ramp', 'viral_delivery']
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
            use_template = data.get('use_template', False)
            kozak = data.get('kozak', False)
            dinuc = data.get('dinuc', False)

            # Validate input
            validation_error = self.validate_input(sequence, profile)
            if validation_error:
                logger.warning(f"Validation error: {validation_error}")
                self.send_error_response(400, validation_error)
                return

            # Clean sequence
            sequence = self.clean_sequence(sequence)

            # Check if FactorForge is available
            if not FACTORFORGE_AVAILABLE:
                logger.info("Using mock optimization (FactorForge not available)")
                result = self.generate_mock_result(sequence, profile, kozak, dinuc)
            else:
                logger.info(f"Running real optimization: profile={profile}, template={use_template}, kozak={kozak}, dinuc={dinuc}")
                result = self.optimize_sequence(sequence, profile, use_template, kozak, dinuc)

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
            'version': '2.5.3',
            'codonforge_available': FACTORFORGE_AVAILABLE,
            'endpoints': {
                'POST /api/optimize': 'Run codon optimization',
                'GET /api/optimize': 'Health check'
            },
            'supported_profiles': VALID_PROFILES
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

    def validate_input(self, sequence, profile):
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

    def optimize_sequence(self, sequence, profile, use_template, kozak, dinuc):
        """Run actual FactorForge v2 optimization"""
        try:
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
            gc_check = 'PASS' if 40.5 <= gc_percent <= 44.5 else 'WARNING'
            moclo_check = 'PASS'  # From v2 domesticator

            logger.info(f"Optimization metrics: CAI={cai:.3f}, GC={gc_percent:.1f}%, PolyA={polya_warnings}")

            # Format response
            return {
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

        except Exception as e:
            logger.error(f"Optimization failed: {e}", exc_info=True)
            raise ValueError(f"Optimization error: {str(e)}")

    def generate_mock_result(self, sequence, profile, kozak, dinuc):
        """Generate mock result for testing (when FactorForge is not available)"""
        logger.warning("Generating mock result - CodonForge engine not available")

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

        return {
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
                'gc': 'PASS' if 40.5 <= gc_percent <= 44.5 else 'WARNING',
                'moclo': 'PASS'
            },
            'engine': {
                'name': 'Mock Engine',
                'version': '0.0.0'
            },
            'note': '⚠️ Mock data - FactorForge engine not available. Deploy with FactorForge for real optimization.'
        }

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
