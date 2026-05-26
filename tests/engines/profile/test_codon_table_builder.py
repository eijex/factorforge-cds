"""
Unit tests for CodonTableBuilder (golden set blending).
"""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.profile.codon_table_builder import build_golden_set


@pytest.fixture
def high_expression_data(tmp_path):
    """Create mock high-expression reference JSON."""
    data = {
        "species": "Nicotiana benthamiana",
        "codon_usage": {
            "F": {"TTC": 0.65, "TTT": 0.35},
            "L": {"TTG": 0.30, "CTT": 0.20, "CTG": 0.20, "CTC": 0.15, "TTA": 0.05, "CTA": 0.10},
            "M": {"ATG": 1.0},
            "A": {"GCT": 0.40, "GCC": 0.30, "GCA": 0.20, "GCG": 0.10},
            "*": {"TAA": 0.50, "TGA": 0.30, "TAG": 0.20},
        },
    }
    path = tmp_path / "high_expression.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


@pytest.fixture
def empirical_data(tmp_path):
    """Create mock empirical codon table JSON."""
    data = {
        "organism": "Nicotiana benthamiana",
        "codons": {
            "TTC": {"aa": "F", "frequency": 0.55, "per_thousand": 22.0},
            "TTT": {"aa": "F", "frequency": 0.45, "per_thousand": 18.0},
            "TTG": {"aa": "L", "frequency": 0.22, "per_thousand": 17.6},
            "CTT": {"aa": "L", "frequency": 0.22, "per_thousand": 17.6},
            "CTG": {"aa": "L", "frequency": 0.21, "per_thousand": 16.8},
            "CTC": {"aa": "L", "frequency": 0.16, "per_thousand": 12.8},
            "TTA": {"aa": "L", "frequency": 0.11, "per_thousand": 8.8},
            "CTA": {"aa": "L", "frequency": 0.08, "per_thousand": 6.4},
            "ATG": {"aa": "M", "frequency": 1.0, "per_thousand": 22.0},
            "GCT": {"aa": "A", "frequency": 0.34, "per_thousand": 18.7},
            "GCC": {"aa": "A", "frequency": 0.28, "per_thousand": 15.4},
            "GCA": {"aa": "A", "frequency": 0.24, "per_thousand": 13.2},
            "GCG": {"aa": "A", "frequency": 0.14, "per_thousand": 7.7},
            "TAA": {"aa": "*", "frequency": 0.48, "per_thousand": 1.0},
            "TGA": {"aa": "*", "frequency": 0.28, "per_thousand": 0.6},
            "TAG": {"aa": "*", "frequency": 0.24, "per_thousand": 0.5},
        },
        "amino_acids": {
            "F": {"codons": ["TTC", "TTT"], "preferred": "TTC"},
            "L": {"codons": ["TTG", "CTT", "CTG", "CTC", "TTA", "CTA"], "preferred": "TTG"},
            "M": {"codons": ["ATG"], "preferred": "ATG"},
            "A": {"codons": ["GCT", "GCC", "GCA", "GCG"], "preferred": "GCT"},
            "*": {"codons": ["TAA", "TGA", "TAG"], "preferred": "TAA"},
        },
        "gc_content": {"overall": 0.44},
    }
    path = tmp_path / "empirical.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


class TestBuildGoldenSet:
    """Test golden set blending logic."""

    def test_blending_produces_valid_schema(self, high_expression_data, empirical_data):
        """Blended table has required sections: codons, amino_acids, organism."""
        result = build_golden_set(high_expression_data, empirical_data, blend_ratio=0.7)
        assert "codons" in result
        assert "amino_acids" in result
        assert "organism" in result
        assert "source" in result
        assert "blend_ratio" in result

    def test_frequencies_normalized_per_aa(self, high_expression_data, empirical_data):
        """Blended frequencies per amino acid sum to ~1.0."""
        result = build_golden_set(high_expression_data, empirical_data, blend_ratio=0.7)
        for aa, info in result["amino_acids"].items():
            total = sum(result["codons"][c]["frequency"] for c in info["codons"])
            assert abs(total - 1.0) < 0.02, f"AA {aa}: frequencies sum to {total}"

    def test_blend_ratio_effect(self, high_expression_data, empirical_data):
        """Higher blend_ratio shifts frequencies toward high-expression data."""
        result_low = build_golden_set(high_expression_data, empirical_data, blend_ratio=0.3)
        result_high = build_golden_set(high_expression_data, empirical_data, blend_ratio=0.9)
        # For Phe: high-expression favors TTC (0.65), empirical also favors TTC (0.55)
        # With higher ratio, TTC frequency should be closer to 0.65
        freq_low = result_low["codons"]["TTC"]["frequency"]
        freq_high = result_high["codons"]["TTC"]["frequency"]
        assert freq_high >= freq_low

    def test_invalid_blend_ratio(self, high_expression_data, empirical_data):
        """Blend ratio outside 0-1 raises ValueError."""
        with pytest.raises(ValueError, match="blend_ratio"):
            build_golden_set(high_expression_data, empirical_data, blend_ratio=1.5)
        with pytest.raises(ValueError, match="blend_ratio"):
            build_golden_set(high_expression_data, empirical_data, blend_ratio=-0.1)

    def test_output_file_written(self, high_expression_data, empirical_data, tmp_path):
        """If output_path specified, writes valid JSON to disk."""
        out = tmp_path / "golden_set_output.json"
        build_golden_set(high_expression_data, empirical_data, output_path=out)
        assert out.exists()
        data = json.loads(out.read_text(encoding="utf-8"))
        assert "codons" in data

    def test_preferred_codon_is_highest_frequency(self, high_expression_data, empirical_data):
        """Preferred codon for each AA has the highest blended frequency."""
        result = build_golden_set(high_expression_data, empirical_data, blend_ratio=0.7)
        for aa, info in result["amino_acids"].items():
            preferred = info["preferred"]
            preferred_freq = result["codons"][preferred]["frequency"]
            for codon in info["codons"]:
                assert result["codons"][codon]["frequency"] <= preferred_freq + 0.001


class TestCAIGoldenSet:
    """Test CAI calculation using golden set reference weights."""

    @pytest.fixture
    def translator(self):
        from factorforge.engines.profile.rules.reverse_translator import ReverseTranslator
        return ReverseTranslator()

    def test_golden_ref_weights_loaded(self, translator):
        """Golden ref weights dict is populated."""
        assert len(translator.golden_ref_weights) > 0
        # Should include common codons
        assert "ATG" in translator.golden_ref_weights
        assert "GCT" in translator.golden_ref_weights

    def test_golden_ref_weights_range(self, translator):
        """All ref weights are in [0, 1]."""
        for codon, w in translator.golden_ref_weights.items():
            assert 0.0 <= w <= 1.0, f"Codon {codon} has weight {w}"

    def test_preferred_codon_has_weight_1(self, translator):
        """Preferred codon (highest freq per AA) has weight 1.0."""
        # ATG is the only Met codon
        assert translator.golden_ref_weights["ATG"] == 1.0

    def test_cai_uses_golden_set(self, translator):
        """CAI with golden set differs from naive working-table calculation."""
        # All GCT codons (preferred Ala in golden set)
        seq_preferred = "ATGGCTGCTGCTGCT"
        cai_preferred = translator.calculate_cai(seq_preferred)
        # All GCG codons (least-preferred Ala in golden set)
        seq_rare = "ATGGCGGCGGCGGCG"
        cai_rare = translator.calculate_cai(seq_rare)
        # Preferred should have higher CAI
        assert cai_preferred > cai_rare

    def test_stop_codons_excluded_from_cai(self, translator):
        """Stop codons should not affect CAI calculation."""
        seq_with_stop = "ATGGCTTAA"
        seq_no_stop = "ATGGCT"
        # Both should have the same CAI (stop excluded)
        assert translator.calculate_cai(seq_with_stop) == translator.calculate_cai(seq_no_stop)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
