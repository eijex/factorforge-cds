"""
Unit tests for v2 utils
"""

import json
import sys
from pathlib import Path

import pytest

# Add project src to path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.v2.utils import (
    calculate_gc,
    load_codon_table,
    parse_fasta_records,
    translate_codon,
)


class TestCalculateGC:
    """Tests for calculate_gc"""

    def test_calculate_gc_normal(self):
        """GC for mixed sequence"""
        assert calculate_gc("ATGC") == 50.0

    def test_calculate_gc_empty(self):
        """GC for empty sequence"""
        assert calculate_gc("") == 0.0

    def test_calculate_gc_all_gc(self):
        """GC for all GC sequence"""
        assert calculate_gc("GGCC") == 100.0

    def test_calculate_gc_no_gc(self):
        """GC for no GC sequence"""
        assert calculate_gc("AAATTT") == 0.0

    def test_calculate_gc_lowercase(self):
        """GC for lowercase sequence"""
        assert calculate_gc("atgc") == 50.0

    def test_calculate_gc_mixed_case(self):
        """GC for mixed-case sequence"""
        assert calculate_gc("AtGc") == 50.0


class TestLoadCodonTable:
    """Tests for load_codon_table"""

    def test_load_codon_table_success(self, tmp_path: Path):
        """Load codon table from file"""
        codon_table = {
            "amino_acids": {"M": {"codons": ["ATG"]}},
            "codons": {"ATG": {"aa": "M", "frequency": 1.0}},
        }
        path = tmp_path / "testorg_codons.json"
        path.write_text(json.dumps(codon_table), encoding="utf-8")

        result = load_codon_table("testorg", tmp_path)

        assert result["codons"]["ATG"]["aa"] == "M"

    def test_load_codon_table_not_found(self, tmp_path: Path):
        """Raise when codon table file missing"""
        with pytest.raises(FileNotFoundError):
            load_codon_table("missing", tmp_path)


class TestTranslateCodon:
    """Tests for translate_codon"""

    def test_translate_codon_success(self):
        """Translate known codon"""
        codon_table = {"ATG": "M"}
        assert translate_codon("atg", codon_table) == "M"

    def test_translate_codon_invalid(self):
        """Raise for unknown codon"""
        codon_table = {"ATG": "M"}
        with pytest.raises(KeyError):
            translate_codon("TAA", codon_table)


class TestParseFastaRecords:
    """Tests for parse_fasta_records."""

    def test_parse_single_record(self):
        records = parse_fasta_records(">seq1\nATGC\n")
        assert records == [("seq1", "ATGC")]

    def test_parse_multiple_records(self):
        records = parse_fasta_records(">a\nATGC\n>b desc\nTTAA\n")
        assert records == [("a", "ATGC"), ("b", "TTAA")]

    def test_parse_invalid_fasta_raises(self):
        with pytest.raises(ValueError, match="before header"):
            parse_fasta_records("ATGC\n>seq\nTTAA\n")
