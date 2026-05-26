"""
Integration tests for OptimizationPipeline
"""

import re
import sys
from pathlib import Path

import pytest

try:
    from Bio import SeqIO
except ImportError:
    pytest.skip("Biopython not installed in CI yet", allow_module_level=True)

# Add project src to path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.profile.pipeline import OptimizationPipeline


@pytest.fixture
def short_protein() -> str:
    """Short test protein sequence"""
    return "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLK"


@pytest.fixture
def gfp_protein() -> str:
    """Full GFP protein sequence"""
    return (
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSR"
        "YPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYN"
        "SHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMV"
        "LLEFVTAAGITHGMDELYK"
    )


def test_pipeline_without_construct(short_protein: str) -> None:
    """Basic optimization without construct"""
    pipeline = OptimizationPipeline(profile="balanced")
    result = pipeline.run(short_protein)

    assert result.sequence
    assert len(result.sequence) > 0
    assert result.construct is None
    assert "profile" in result.metadata


def test_pipeline_with_standard_template(short_protein: str) -> None:
    """standard_expression template usage"""
    pipeline = OptimizationPipeline(profile="balanced", construct_template="standard_expression")
    result = pipeline.run(short_protein)

    assert result.sequence
    assert result.construct is not None
    assert len(result.construct.features) == 4
    assert result.metadata.get("construct_template") == "standard_expression"


def test_pipeline_with_high_expression(short_protein: str) -> None:
    """high_expression template usage"""
    pipeline = OptimizationPipeline(profile="balanced", construct_template="high_expression")
    result = pipeline.run(short_protein)

    assert result.construct is not None
    assert len(result.construct.features) == 4


def test_pipeline_invalid_profile_raises(short_protein: str) -> None:
    """Unknown profile should raise instead of silently falling back."""
    pipeline = OptimizationPipeline(profile="balanced")
    with pytest.raises(ValueError, match="Unknown profile"):
        pipeline.run(short_protein, profile="expression")


def test_pipeline_fast_scan_mode(short_protein: str) -> None:
    """Fast scan mode should be reflected in metadata."""
    pipeline = OptimizationPipeline(profile="balanced")
    result = pipeline.run(short_protein, scan_mode="fast")

    assert result.metadata["scan_mode"] == "fast"
    assert "scan_results" in result.metadata
    assert "dinucleotides" not in result.metadata["scan_results"]


def test_pipeline_run_batch(short_protein: str) -> None:
    """Batch pipeline execution should return ordered results with ids."""
    pipeline = OptimizationPipeline(profile="balanced")
    payload = [
        {"id": "s1", "sequence": short_protein},
        {"id": "s2", "sequence": short_protein},
    ]
    results = pipeline.run_batch(payload, scan_mode="fast")

    assert len(results) == 2
    assert results[0].metadata["input_id"] == "s1"
    assert results[1].metadata["input_id"] == "s2"


def test_pipeline_save_genbank(short_protein: str, tmp_path: Path) -> None:
    """GenBank save for constructs"""
    pipeline = OptimizationPipeline(profile="balanced", construct_template="standard_expression")
    result = pipeline.run(short_protein)

    output = tmp_path / "result.gb"
    result.save(output, format="genbank")

    assert output.exists()
    record = SeqIO.read(output, "genbank")
    assert len(record.features) == 4
    assert len(record.seq) > 0


def test_pipeline_save_fasta_without_construct(short_protein: str, tmp_path: Path) -> None:
    """FASTA save without construct"""
    pipeline = OptimizationPipeline(profile="balanced")
    result = pipeline.run(short_protein)

    output = tmp_path / "result.fasta"
    result.save(output, format="fasta")

    assert output.exists()


def test_viral_delivery_pipeline_integration(short_protein: str) -> None:
    """pipeline.run()에서 viral_delivery 프로필 성공 확인 (001-fix)"""
    pipeline = OptimizationPipeline(profile="viral_delivery")
    result = pipeline.run(short_protein)

    assert result.sequence
    assert len(result.sequence) == len(short_protein) * 3
    assert result.metadata.get("profile") == "viral_delivery"


def test_pipeline_end_to_end_gfp(gfp_protein: str) -> None:
    """Full end-to-end GFP construct"""
    pipeline = OptimizationPipeline(profile="balanced", construct_template="standard_expression")
    result = pipeline.run(gfp_protein)

    assert result.sequence
    assert result.construct
    assert len(result.construct.features) == 4

    cds_length = len(gfp_protein) * 3
    assert len(result.construct.seq) > cds_length

    assert "construct_template" in result.metadata
    assert "construct_features" in result.metadata
    assert result.metadata["construct_features"] == 4


class TestConstructIdAndExportFeatures:
    """014-small: construct_id + export_features() 검증"""

    def test_construct_id_in_metadata(self, short_protein: str) -> None:
        """run() 결과 metadata에 CF-YYYYMMDD-HHMMSS 형식 construct_id 포함"""
        pipeline = OptimizationPipeline(profile="balanced")
        result = pipeline.run(short_protein)

        assert "construct_id" in result.metadata
        assert re.match(r"CF-\d{8}-\d{6}", result.metadata["construct_id"])

    def test_export_features_required_keys(self, short_protein: str) -> None:
        """export_features() 반환 dict에 schema.md 필수 8개 컬럼 존재"""
        pipeline = OptimizationPipeline(profile="balanced")
        result = pipeline.run(short_protein)
        features = result.export_features()

        required = {
            "construct_id", "protein_name", "optimization_profile",
            "cai_score", "gc_content_pct", "mfe_kcal_mol",
            "polya_signal_count", "domestication_edits",
        }
        assert required.issubset(features.keys())

    def test_export_features_purity_is_none(self, short_protein: str) -> None:
        """purity_pct는 None — 실험 후 수동 입력 필드"""
        pipeline = OptimizationPipeline(profile="balanced")
        result = pipeline.run(short_protein)
        features = result.export_features()

        assert features["purity_pct"] is None

    def test_export_features_value_ranges(self, short_protein: str) -> None:
        """cai_score 0~1, gc_content_pct 0~100 범위 확인"""
        pipeline = OptimizationPipeline(profile="balanced")
        result = pipeline.run(short_protein)
        features = result.export_features()

        assert 0.0 <= features["cai_score"] <= 1.0
        assert 0.0 <= features["gc_content_pct"] <= 100.0

    def test_export_features_construct_id_matches_metadata(self, short_protein: str) -> None:
        """export_features()의 construct_id가 metadata와 동일"""
        pipeline = OptimizationPipeline(profile="balanced")
        result = pipeline.run(short_protein)
        features = result.export_features()

        assert features["construct_id"] == result.metadata["construct_id"]
