"""
Unit tests for CpG/TpA dinucleotide fixing in RuleEngine.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.profile.pipeline import OptimizationPipeline
from factorforge.engines.profile.rules.rule_engine import RuleEngine
from factorforge.engines.profile.utils import count_dinucleotides


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
        result = local_engine.fix_dinucleotides(seq, mode="aggressive")

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


class TestDinucleotideModes:
    """Tests for fix_dinucleotides() mode parameter."""

    def test_aggressive_matches_default_call_for_cpg_sequence(self, engine):
        """Aggressive mode matches the default call for a simple reducible sequence."""
        seq = "ACG" * 10
        default_result = engine.fix_dinucleotides(seq)
        aggressive_result = engine.fix_dinucleotides(seq, mode="aggressive")

        assert aggressive_result["mode"] == "aggressive"
        assert aggressive_result["initial_count"] == default_result["initial_count"]
        assert aggressive_result["final_count"] == default_result["final_count"]

    def test_aggressive_returns_cai_fields(self, engine):
        """Aggressive mode reports CAI before and after fixing."""
        result = engine.fix_dinucleotides("ACG" * 10, mode="aggressive")

        assert result["cai_before"] > 0.0
        assert result["cai_after"] >= 0.0
        assert result["mode"] == "aggressive"

    def test_balanced_respects_low_cai_floor(self, engine):
        """Balanced mode reports CAI above a permissive floor."""
        result = engine.fix_dinucleotides("ACG" * 20, mode="balanced", cai_floor=0.5)

        assert result["cai_after"] >= 0.5
        assert result["mode"] == "balanced"

    def test_balanced_high_cai_floor_limits_reduction(self, engine):
        """An unreachable high floor rolls back or keeps CAI above the floor."""
        result = engine.fix_dinucleotides("ACG" * 20, mode="balanced", cai_floor=0.999)

        assert (
            result["final_count"] >= result["initial_count"] * 0.5 or result["cai_after"] >= 0.999
        )
        assert result["mode"] == "balanced"

    def test_cai_preserving_limits_cai_drop(self, engine):
        """CAI-preserving mode keeps CAI loss within the requested budget."""
        result = engine.fix_dinucleotides(
            "ACG" * 20,
            mode="cai_preserving",
            max_cai_drop=0.001,
        )

        assert result["cai_before"] - result["cai_after"] <= 0.002
        assert result["mode"] == "cai_preserving"

    def test_mode_field_returned_for_non_triplet_sequence(self, engine):
        """Mode is reported even when the input sequence cannot be fixed."""
        result = engine.fix_dinucleotides("ACGA", mode="cai_preserving")

        assert result["success"] is False
        assert result["mode"] == "cai_preserving"

    def test_pipeline_uses_balanced_mode_without_error(self):
        """Pipeline default dinucleotide handling completes with balanced mode."""
        pipeline = OptimizationPipeline(profile="balanced")
        result = pipeline.run("ACG" * 10)

        assert result.sequence
        assert "scan_results" in result.metadata
