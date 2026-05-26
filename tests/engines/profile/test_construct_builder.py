"""
Unit tests for ConstructBuilder
"""

import sys
from pathlib import Path

import pytest

Bio = pytest.importorskip("Bio")
from Bio.SeqRecord import SeqRecord  # noqa: E402

# Add project src to path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.profile.construct_builder import ConstructBuilder  # noqa: E402


@pytest.fixture
def builder() -> ConstructBuilder:
    """ConstructBuilder instance"""
    template_dir = Path(__file__).parent.parent.parent.parent / "data" / "templates"
    return ConstructBuilder(template_dir)


@pytest.fixture
def sample_gene() -> str:
    """Test GFP CDS sequence (simplified)"""
    return (
        "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGC"
        "ACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCAC"
        "TACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCGGTTATGGTGTTCAATGCTTTGCGAGA"
        "TACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTA"
        "TATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAG"
        "AATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAAC"
        "TCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACA"
        "TTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTT"
        "ACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTC"
        "CTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAA"
    )


def test_load_template_standard(builder: ConstructBuilder) -> None:
    """Load standard_expression template"""
    template = builder.load_template("standard_expression")

    assert template["name"]
    assert "components" in template
    assert "overhangs" in template


def test_load_template_high(builder: ConstructBuilder) -> None:
    """Load high_expression template"""
    template = builder.load_template("high_expression")

    assert template["name"]
    assert "components" in template


def test_load_template_not_found(builder: ConstructBuilder) -> None:
    """Fail on missing template"""
    with pytest.raises(FileNotFoundError):
        builder.load_template("missing_template")


def test_assemble_parts(builder: ConstructBuilder, sample_gene: str) -> None:
    """Assemble parts with gene insertion"""
    template = builder.load_template("standard_expression")
    construct_seq = builder.assemble_parts(sample_gene, template)

    assert sample_gene in construct_seq
    assert "USER_INPUT" not in construct_seq

    expected_length = 0
    for component in template["components"]:
        if component["sequence"] == "USER_INPUT":
            expected_length += len(sample_gene)
        else:
            expected_length += len(component["sequence"])
    assert len(construct_seq) == expected_length


def test_add_features(builder: ConstructBuilder, sample_gene: str) -> None:
    """Add features to construct"""
    template = builder.load_template("standard_expression")
    construct_seq = builder.assemble_parts(sample_gene, template)
    record = builder.add_features(construct_seq, template)

    assert isinstance(record, SeqRecord)
    assert len(record.features) == 4
    feature_types = [feature.type for feature in record.features]
    assert feature_types == ["promoter", "5'UTR", "CDS", "terminator"]


def test_validate_construct_success(builder: ConstructBuilder, sample_gene: str) -> None:
    """Validate a normal construct"""
    template = builder.load_template("standard_expression")
    construct_seq = builder.assemble_parts(sample_gene, template)
    record = builder.add_features(construct_seq, template)
    valid, warnings = builder.validate_construct(record, template)

    assert valid is True
    assert all("outside expected range" not in warning for warning in warnings)


def test_validate_construct_warnings(builder: ConstructBuilder) -> None:
    """Detect warnings for short constructs"""
    template = builder.load_template("standard_expression")
    construct_seq = builder.assemble_parts("ATG", template)
    record = builder.add_features(construct_seq, template)
    valid, warnings = builder.validate_construct(record, template)

    assert valid is False
    assert warnings


def test_generate_construct_end_to_end(builder: ConstructBuilder, sample_gene: str) -> None:
    """End-to-end construct generation"""
    record = builder.generate_construct(sample_gene, "standard_expression")

    assert isinstance(record, SeqRecord)
    assert len(record.features) == 4
    assert len(record.seq) > len(sample_gene)


# --- Phase 1: Overhang Validation Tests ---


class TestOverhangValidation:
    """Test MoClo overhang validation and collision detection."""

    def test_validate_overhangs_correct_moclo(self, builder: ConstructBuilder) -> None:
        """Valid MoClo Level 0 chain passes validation."""
        parts = [
            {"overhang_5": "AATG", "overhang_3": "AGGT"},
            {"overhang_5": "AGGT", "overhang_3": "GCTT"},
        ]
        valid, warnings = builder.validate_overhangs(parts)
        assert valid is True
        assert warnings == []

    def test_validate_overhangs_incorrect_5prime(self, builder: ConstructBuilder) -> None:
        """Wrong 5' overhang on first part is detected."""
        parts = [
            {"overhang_5": "TTTT", "overhang_3": "AGGT"},
            {"overhang_5": "AGGT", "overhang_3": "GCTT"},
        ]
        valid, warnings = builder.validate_overhangs(parts)
        assert valid is False
        assert any("5' overhang" in w for w in warnings)

    def test_validate_overhangs_incorrect_3prime(self, builder: ConstructBuilder) -> None:
        """Wrong 3' overhang on last part is detected."""
        parts = [
            {"overhang_5": "AATG", "overhang_3": "AGGT"},
            {"overhang_5": "AGGT", "overhang_3": "CCCC"},
        ]
        valid, warnings = builder.validate_overhangs(parts)
        assert valid is False
        assert any("3' overhang" in w for w in warnings)

    def test_validate_overhangs_chain_mismatch(self, builder: ConstructBuilder) -> None:
        """Adjacent parts with non-matching overhangs are detected."""
        parts = [
            {"overhang_5": "AATG", "overhang_3": "AGGT"},
            {"overhang_5": "XXXX", "overhang_3": "GCTT"},
        ]
        valid, warnings = builder.validate_overhangs(parts)
        assert valid is False
        assert any("mismatch" in w.lower() for w in warnings)

    def test_validate_overhangs_empty_parts(self, builder: ConstructBuilder) -> None:
        """Empty parts list returns invalid."""
        valid, warnings = builder.validate_overhangs([])
        assert valid is False

    def test_internal_overhang_collision_detected(self, builder: ConstructBuilder) -> None:
        """CDS with internal AATG is detected as collision."""
        # Construct a CDS containing AATG internally
        cds_with_aatg = "ATGCCGAATGGCCTAA"  # Contains AATG at position 6
        collisions = builder.check_internal_overhang_collisions(cds_with_aatg)
        assert len(collisions) > 0
        assert any(c["overhang"] == "AATG" for c in collisions)

    def test_internal_overhang_collision_reverse_complement(
        self, builder: ConstructBuilder
    ) -> None:
        """CDS with reverse complement CATT (rc of AATG) is detected."""
        # CATT is the reverse complement of AATG
        cds_with_rc = "ATGCCGCATTGCCTAA"  # Contains CATT at position 6
        collisions = builder.check_internal_overhang_collisions(cds_with_rc)
        assert len(collisions) > 0
        assert any(c["strand"] == "reverse_complement" for c in collisions)

    def test_internal_overhang_collision_clean(self, builder: ConstructBuilder) -> None:
        """CDS without overhang sequences returns empty list."""
        clean_cds = "ATGCCGCCGCCGTAA"
        collisions = builder.check_internal_overhang_collisions(clean_cds)
        # Filter out any hits that aren't AATG or GCTT
        relevant = [c for c in collisions if c["overhang"] in ("AATG", "GCTT")]
        assert len(relevant) == 0
