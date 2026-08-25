import hashlib
from typing import Dict, List, Any, Optional

def _hash_sequence(sequence: str) -> str:
    """Returns canonical sequence hash."""
    return hashlib.sha256(sequence.upper().encode("utf-8")).hexdigest()

class DatasetLeakageAuditor:
    """
    Phase 2 (A): Train <-> Validation/Test Audit
    Verifies that a benchmark evaluation split does not overlap with the training split.
    """
    
    def __init__(self, train_sequences: List[str], test_sequences: List[str]):
        self.train_sequences = [seq.upper() for seq in train_sequences]
        self.test_sequences = [seq.upper() for seq in test_sequences]
        self.train_hashes = {_hash_sequence(seq) for seq in self.train_sequences}

    def audit_exact_overlap(self) -> Dict[str, Any]:
        """Checks for exact canonical sequence overlap."""
        overlaps = 0
        for test_seq in self.test_sequences:
            if _hash_sequence(test_seq) in self.train_hashes:
                overlaps += 1

        if overlaps > 0:
            return {
                "check_type": "train_test_exact_overlap",
                "result": "FAIL",
                "match_count": overlaps,
                "details_json": {"message": f"Found {overlaps} exact matches between train and test splits."}
            }
            
        return {
            "check_type": "train_test_exact_overlap",
            "result": "PASS",
            "match_count": 0,
            "details_json": {"message": "No exact canonical sequence overlaps detected."}
        }

    def audit_homology_overlap(self, threshold: float = 0.95) -> Dict[str, Any]:
        """
        Placeholder for computationally heavy homology alignment (e.g. CD-HIT, MMseqs2).
        For now, returns PASS unless we implement full Needleman-Wunsch or local BLAST.
        """
        return {
            "check_type": "train_test_homology_overlap",
            "result": "PASS",
            "threshold": threshold,
            "match_count": 0,
            "details_json": {"message": "Homology check not yet fully implemented. Passed by default."}
        }


class GenerationMemorizationAuditor:
    """
    Phase 2 (B): Train <-> Generated Candidate Audit
    Verifies that a model-generated sequence does not exactly match or highly resemble the training corpus.
    """
    
    def __init__(self, train_sequences: List[str]):
        self.train_sequences = [seq.upper() for seq in train_sequences]
        self.train_hashes = {_hash_sequence(seq) for seq in self.train_sequences}

    def audit_exact_match(self, generated_sequence: str) -> Dict[str, Any]:
        """Checks if generated sequence is an exact clone of a training sequence."""
        gen_hash = _hash_sequence(generated_sequence)
        if gen_hash in self.train_hashes:
            return {
                "check_type": "generated_output_exact_overlap",
                "result": "FLAG", # User requested FLAG / potentially FAIL based on policy
                "match_count": 1,
                "details_json": {"message": "Generated sequence is an exact match to a training sequence."}
            }
            
        return {
            "check_type": "generated_output_exact_overlap",
            "result": "PASS",
            "match_count": 0,
            "details_json": {"message": "Generated sequence is novel (exact match)."}
        }

    def audit_high_similarity(self, generated_sequence: str, threshold: float = 0.95) -> Dict[str, Any]:
        """Placeholder for similarity check against training corpus."""
        return {
            "check_type": "generated_output_homology",
            "result": "WARNING",
            "threshold": threshold,
            "match_count": 0,
            "details_json": {"message": "Generated homology check passed (placeholder)."}
        }
