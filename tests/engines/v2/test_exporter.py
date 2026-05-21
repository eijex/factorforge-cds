"""
Unit tests for SequenceExporter
"""

import os
import sys
import tempfile
from pathlib import Path

import pytest

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.v2.exporter import SequenceExporter


@pytest.fixture
def exporter():
    """Create exporter instance"""
    return SequenceExporter()


@pytest.fixture
def sample_sequence():
    """Sample DNA sequence (GFP partial)"""
    return "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"


@pytest.fixture
def sample_metadata():
    """Sample metadata"""
    return {
        "gene_name": "GFP",
        "protein_seq": "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHKVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK*",
        "profile": "Balanced",
        "cai": 0.87,
        "gc": 51.2,
        "assembly_standard": "Golden Gate (BsaI)",
        "run_id": "abc12345",
        "timestamp": "2026-01-22T12:00:00",
        "violations_fixed": [{"type": "BsaI site", "position": 147}],
        "warnings": [{"message": "Test warning", "suggestion": "Test suggestion"}],
    }


class TestRunIDGeneration:
    """Test run_id generation"""

    def test_generate_run_id(self, exporter):
        """Test run_id generation"""
        seq = "ATGGCCTAA"
        params = {"profile": "balanced"}

        run_id = exporter.generate_run_id(seq, params)

        assert isinstance(run_id, str)
        assert len(run_id) == 8
        assert run_id.isalnum()

    def test_run_id_reproducibility(self, exporter):
        """Test run_id is reproducible"""
        seq = "ATGGCCTAA"
        params = {"profile": "balanced"}

        run_id1 = exporter.generate_run_id(seq, params)
        run_id2 = exporter.generate_run_id(seq, params)

        assert run_id1 == run_id2

    def test_run_id_uniqueness(self, exporter):
        """Test different inputs produce different run_ids"""
        seq1 = "ATGGCCTAA"
        seq2 = "ATGGCCTAG"
        params = {"profile": "balanced"}

        run_id1 = exporter.generate_run_id(seq1, params)
        run_id2 = exporter.generate_run_id(seq2, params)

        assert run_id1 != run_id2


class TestFASTAExport:
    """Test FASTA export functionality"""

    def test_export_fasta_basic(self, exporter, sample_sequence, sample_metadata):
        """Test basic FASTA export"""
        fasta = exporter.export_fasta(sample_sequence, sample_metadata)

        assert fasta.startswith(">")
        assert "PFORM_" in fasta
        assert "GFP" in fasta
        assert "CAI=0.870" in fasta
        assert "GC=51.2" in fasta
        # Strip newlines and check the sequence
        fasta_no_newlines = fasta.replace("\n", "").replace("\r", "")
        assert sample_sequence in fasta_no_newlines

    def test_fasta_header_format(self, exporter, sample_sequence, sample_metadata):
        """Test FASTA header format"""
        fasta = exporter.export_fasta(sample_sequence, sample_metadata)

        lines = fasta.split("\n")
        header = lines[0]

        assert header.startswith(">PFORM_")
        assert "gene=GFP" in header
        assert "profile=Balanced" in header

    def test_fasta_line_width(self, exporter, sample_sequence, sample_metadata):
        """Test FASTA line width formatting"""
        fasta = exporter.export_fasta(sample_sequence, sample_metadata, line_width=60)

        lines = fasta.split("\n")[1:]  # Skip header
        lines = [l for l in lines if l]  # Remove empty lines

        # Check that lines (except possibly last) are 60 chars
        for line in lines[:-1]:
            assert len(line) == 60

    def test_fasta_no_line_breaks(self, exporter, sample_sequence, sample_metadata):
        """Test FASTA without line breaks"""
        fasta = exporter.export_fasta(sample_sequence, sample_metadata, line_width=0)

        lines = fasta.strip().split("\n")
        assert len(lines) == 2  # Header + sequence
        assert lines[1] == sample_sequence

    def test_fasta_file_output(self, exporter, sample_sequence, sample_metadata):
        """Test FASTA file output"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as f:
            temp_file = f.name

        try:
            result = exporter.export_fasta(sample_sequence, sample_metadata, output_file=temp_file)

            assert "written to" in result
            assert os.path.exists(temp_file)

            with open(temp_file, "r") as f:
                content = f.read()
                assert content.startswith(">")
                # Strip newlines and check the sequence
                content_no_newlines = content.replace("\n", "").replace("\r", "")
                assert sample_sequence in content_no_newlines
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)


class TestGenBankExport:
    """Test GenBank export functionality"""

    def test_export_genbank_basic(self, exporter, sample_sequence, sample_metadata):
        """Test basic GenBank export"""
        try:
            genbank = exporter.export_genbank(sample_sequence, sample_metadata)

            assert "LOCUS" in genbank
            assert "PFORM_" in genbank
            assert "FactorForge v2.0" in genbank
            assert "CDS" in genbank
        except ImportError:
            pytest.skip("Biopython not installed")

    def test_genbank_metadata(self, exporter, sample_sequence, sample_metadata):
        """Test GenBank metadata inclusion"""
        try:
            genbank = exporter.export_genbank(sample_sequence, sample_metadata)

            assert "Run ID: abc12345" in genbank
            assert "CAI: 0.870" in genbank
            assert "GC%: 51.2" in genbank
            assert "Profile: Balanced" in genbank
        except ImportError:
            pytest.skip("Biopython not installed")

    def test_genbank_cds_feature(self, exporter, sample_sequence, sample_metadata):
        """Test GenBank CDS feature"""
        try:
            genbank = exporter.export_genbank(sample_sequence, sample_metadata)

            assert "CDS" in genbank
            assert "/translation=" in genbank
            assert "/codon_opt=" in genbank
        except ImportError:
            pytest.skip("Biopython not installed")

    def test_genbank_feature_annotations(self, exporter, sample_sequence, sample_metadata):
        """Test GenBank feature annotations"""
        try:
            metadata = dict(sample_metadata)
            metadata["features"] = [
                {
                    "start": 0,
                    "end": 12,
                    "type": "promoter",
                    "qualifiers": {"note": ["synthetic promoter"]},
                }
            ]
            genbank = exporter.export_genbank(sample_sequence, metadata)

            assert "promoter" in genbank
            assert "synthetic promoter" in genbank
        except ImportError:
            pytest.skip("Biopython not installed")

    def test_genbank_empty_features(self, exporter, sample_sequence, sample_metadata):
        """Test GenBank export with empty features"""
        try:
            metadata = dict(sample_metadata)
            metadata["features"] = []
            genbank = exporter.export_genbank(sample_sequence, metadata)

            assert "FEATURES" in genbank
            assert "CDS" in genbank
        except ImportError:
            pytest.skip("Biopython not installed")

    def test_genbank_file_output(self, exporter, sample_sequence, sample_metadata):
        """Test GenBank file output"""
        try:
            with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".gb") as f:
                temp_file = f.name

            try:
                result = exporter.export_genbank(
                    sample_sequence, sample_metadata, output_file=temp_file
                )

                assert "written to" in result
                assert os.path.exists(temp_file)

                with open(temp_file, "r") as f:
                    content = f.read()
                    assert "LOCUS" in content
            finally:
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
        except ImportError:
            pytest.skip("Biopython not installed")


class TestBatchExport:
    """Test batch export functionality"""

    def test_batch_fasta_export(self, exporter):
        """Test batch FASTA export"""
        sequences = [
            {"sequence": "ATGGCCTAA", "metadata": {"gene_name": "gene1", "cai": 0.8, "gc": 50.0}},
            {"sequence": "ATGGCCTAG", "metadata": {"gene_name": "gene2", "cai": 0.9, "gc": 55.0}},
        ]

        result = exporter.export_batch(sequences, output_format="fasta")

        assert result.count(">") == 2
        assert "gene1" in result
        assert "gene2" in result

    def test_batch_fasta_file_output(self, exporter):
        """Test batch FASTA file output"""
        sequences = [
            {"sequence": "ATGGCCTAA", "metadata": {"gene_name": "gene1", "cai": 0.8, "gc": 50.0}},
            {"sequence": "ATGGCCTAG", "metadata": {"gene_name": "gene2", "cai": 0.9, "gc": 55.0}},
        ]

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as f:
            temp_file = f.name

        try:
            result = exporter.export_batch(sequences, output_format="fasta", output_file=temp_file)

            assert "2 sequences" in result
            assert os.path.exists(temp_file)

            with open(temp_file, "r") as f:
                content = f.read()
                assert content.count(">") == 2
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)

    def test_batch_genbank_requires_output_file(self, exporter):
        """Test GenBank batch export requires output_file"""
        sequences = [
            {
                "sequence": "ATGGCCTAA",
                "metadata": {"gene_name": "gene1", "protein_seq": "MA*"},
            }
        ]

        with pytest.raises(ValueError):
            exporter.export_batch(sequences, output_format="genbank")


class TestReportExport:
    """Test report export functionality"""

    def test_export_report_basic(self, exporter, sample_sequence, sample_metadata):
        """Test basic report export"""
        report = exporter.export_report(sample_sequence, sample_metadata)

        assert "FactorForge v2.0" in report
        assert "Run ID:" in report
        assert "GFP" in report
        assert "CAI: 0.870" in report
        assert "GC Content: 51.2%" in report

    def test_report_violations_section(self, exporter, sample_sequence, sample_metadata):
        """Test report violations section"""
        report = exporter.export_report(sample_sequence, sample_metadata)

        assert "Violations Fixed" in report
        assert "BsaI site" in report

    def test_report_warnings_section(self, exporter, sample_sequence, sample_metadata):
        """Test report warnings section"""
        report = exporter.export_report(sample_sequence, sample_metadata)

        assert "Warnings" in report
        assert "Test warning" in report

    def test_report_sequence_preview(self, exporter, sample_sequence, sample_metadata):
        """Test report sequence preview"""
        report = exporter.export_report(sample_sequence, sample_metadata)

        assert "Sequence Preview" in report
        assert sample_sequence[:100] in report

    def test_report_file_output(self, exporter, sample_sequence, sample_metadata):
        """Test report file output"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            temp_file = f.name

        try:
            result = exporter.export_report(sample_sequence, sample_metadata, output_file=temp_file)

            assert "written to" in result
            assert os.path.exists(temp_file)

            with open(temp_file, "r", encoding="utf-8") as f:
                content = f.read()
                assert "FactorForge v2.0" in content
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)


class TestEdgeCases:
    """Test edge cases"""

    def test_minimal_metadata(self, exporter):
        """Test export with minimal metadata"""
        seq = "ATGGCCTAA"
        metadata = {}

        fasta = exporter.export_fasta(seq, metadata)
        assert fasta.startswith(">")

        report = exporter.export_report(seq, metadata)
        assert "FactorForge v2.0" in report

    def test_empty_sequence(self, exporter):
        """Test export with empty sequence"""
        seq = ""
        metadata = {"gene_name": "test"}

        fasta = exporter.export_fasta(seq, metadata)
        assert ">PFORM_" in fasta

    def test_special_characters_in_gene_name(self, exporter):
        """Test gene name with special characters"""
        seq = "ATGGCCTAA"
        metadata = {"gene_name": "test-gene_v2.1"}

        fasta = exporter.export_fasta(seq, metadata)
        assert "test-gene_v2.1" in fasta


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
