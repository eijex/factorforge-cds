"""
Unit tests for CpG/TpA dinucleotide fixing in RuleEngine.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.v2.pipeline import OptimizationPipeline
from factorforge.engines.v2.rules.rule_engine import RuleEngine
from factorforge.engines.v2.utils import count_dinucleotides


@pytest.fixture
def engine():
    """Create RuleEngine instance."""
    return RuleEngine()


def translate_with_engine(seq: str, engine: RuleEngine) -> str:
    """Translate a CDS using the engine's active codon table."""
    codons = engine.codon_table["codons"]
    return "".join(codons[seq[i : i + 3]]["aa"] for i in range(0, len(seq), 3))


class TestDinucleotideFix:
    """Tests for RuleEngine.fix_dinucleotides."""

    def test_basic_cpg_reduction(self, engine):
        """ACG repeats have internal CpGs that can be reduced by Thr synonyms."""
        seq = "ACG" * 10
        result = engine.fix_dinucleotides(seq)

        assert result["success"] is True
        assert result["final_count"] < result["initial_count"]
        assert count_dinucleotides(result["modified_seq"], "CG") < count_dinucleotides(seq, "CG")

    def test_amino_acid_sequence_is_preserved(self, engine):
        """Synonymous substitutions preserve the translated amino acid sequence."""
        seq = "ACG" * 10
        result = engine.fix_dinucleotides(seq)

        assert translate_with_engine(result["modified_seq"], engine) == translate_with_engine(
            seq, engine
        )

    def test_non_triplet_length_returns_original_sequence(self, engine):
        """A non-CDS-length sequence is returned unchanged."""
        result = engine.fix_dinucleotides("ACGA")

        assert result["success"] is False
        assert result["modified_seq"] == "ACGA"
        assert result["rounds"] == 0

    def test_sequence_without_target_dinucleotides_is_not_changed(self, engine):
        """Sequence with no CG or TA target dinucleotides reports no fix."""
        seq = "GCT" * 10
        result = engine.fix_dinucleotides(seq)

        assert result["success"] is False
        assert result["initial_count"] == 0
        assert result["final_count"] == 0
        assert result["modified_seq"] == seq

    def test_reduction_pct_exact_for_partial_reduction(self):
        """A 10 to 6 reduction reports 40.0%."""
        codon_table = {
            "amino_acids": {
                "T": {"codons": ["ACG", "ACC"]},
                "P": {"codons": ["CCG"]},
            },
            "codons": {
                "ACG": {"aa": "T", "frequency": 0.5},
                "ACC": {"aa": "T", "frequency": 0.5},
                "CCG": {"aa": "P", "frequency": 1.0},
            },
        }
        local_engine = RuleEngine(codon_table)
        seq = ("ACG" * 4) + ("CCG" * 6)
        result = local_engine.fix_dinucleotides(seq)

        assert result["initial_count"] == 10
        assert result["final_count"] == 6
        assert result["reduction_pct"] == 40.0

    def test_pipeline_integration_runs_without_error(self):
        """Pipeline invokes the dinucleotide fix path and completes."""
        pipeline = OptimizationPipeline(profile="balanced")
        result = pipeline.run("ACG" * 10)

        assert result.sequence
        assert len(result.sequence) == 30
        assert "scan_results" in result.metadata

    def test_round_count_does_not_exceed_max_rounds(self, engine):
        """Reported rounds stay within the requested max_rounds."""
        result = engine.fix_dinucleotides("ACG" * 10, max_rounds=3)

        assert result["rounds"] <= 3
