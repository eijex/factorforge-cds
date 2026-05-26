"""
Exporter for FactorForge v3.x profile outputs.
GenBank and FASTA export module (P0-5)
"""

from __future__ import annotations

import hashlib
from datetime import datetime
from io import StringIO
from typing import Any


class SequenceExporter:
    """
    Export optimized sequences in GenBank and FASTA formats

    Features:
    - GenBank format (with metadata)
    - FASTA format (key info in header)
    - Reproducibility via run ID
    """

    def __init__(self) -> None:
        """Initialize"""
        pass

    def generate_run_id(self, sequence: str, params: dict[str, Any]) -> str:
        """
        Generate a reproducible run_id

        Args:
            sequence: DNA sequence
            params: Optimization parameters

        Returns:
            8-character hash string

        Raises:
            None.

        Examples:
            >>> exporter = SequenceExporter()
            >>> run_id = exporter.generate_run_id("ATG", {"profile": "balanced"})
            >>> len(run_id) == 8
            True
        """
        # Create a hash from sequence + parameters
        content = f"{sequence}_{params.get('profile', 'balanced')}_{params.get('assembly_standard', 'none')}"
        hash_obj = hashlib.md5(content.encode())
        return hash_obj.hexdigest()[:8]

    def export_genbank(
        self,
        sequence: str,
        metadata: dict[str, Any],
        output_file: str | None = None,
    ) -> str:
        """
        Export in GenBank format

        Args:
            sequence: Optimized DNA sequence
            metadata: {
                "protein_seq": "MAKLFG...",
                "profile": "Balanced",
                "cai": 0.87,
                "gc": 51.2,
                "run_id": "abc12345",
                "timestamp": "2026-01-22T12:00:00",
                "organism": "Nicotiana benthamiana",
                "gene_name": "GFP",
                "violations_fixed": [...],
                "warnings": [...]
            }
            output_file: Output file path (returns string if None)

        Returns:
            GenBank-formatted string

        Raises:
            ImportError: If Biopython is not installed.

        Examples:
            >>> exporter = SequenceExporter()
            >>> gb = exporter.export_genbank("ATG", {"protein_seq": "M"})
            >>> "LOCUS" in gb
            True
        """
        try:
            from Bio import SeqIO
            from Bio.Seq import Seq
            from Bio.SeqFeature import FeatureLocation, SeqFeature
            from Bio.SeqRecord import SeqRecord
        except ImportError:
            raise ImportError("Biopython is required: pip install biopython")

        # Set defaults
        run_id = metadata.get("run_id", self.generate_run_id(sequence, metadata))
        timestamp = metadata.get("timestamp", datetime.now().strftime("%Y%m%d"))
        gene_name = metadata.get("gene_name", "optimized_gene")
        organism = metadata.get("organism", "Nicotiana benthamiana")

        # Build locus ID
        locus_id = f"PFORM_{run_id}_{timestamp}"

        # Build SeqRecord
        record = SeqRecord(
            Seq(sequence),
            id=locus_id,
            name=gene_name[:16],  # GenBank name is limited to 16 chars
            description=f"Codon-optimized for {organism}",
        )

        # Add annotations
        record.annotations["molecule_type"] = "DNA"
        record.annotations["topology"] = "linear"
        record.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()
        record.annotations["organism"] = organism

        # Build COMMENT section
        comment_lines = [
            "FactorForge v3.x - Plant Codon Optimization Tool",
            f"Run ID: {run_id}",
            f"Timestamp: {metadata.get('timestamp', datetime.now().isoformat())}",
            f"Profile: {metadata.get('profile', 'N/A')}",
            f"CAI: {metadata.get('cai', 0.0):.3f}",
            f"GC%: {metadata.get('gc', 0.0):.1f}",
        ]

        # Assembly standard info
        if metadata.get("assembly_standard"):
            comment_lines.append(f"Assembly Standard: {metadata['assembly_standard']}")

        # Fixed violations
        violations_fixed = metadata.get("violations_fixed", [])
        if violations_fixed:
            comment_lines.append(f"Violations Fixed: {len(violations_fixed)}")
            for v in violations_fixed[:5]:  # Show at most 5
                comment_lines.append(
                    f"  - {v.get('type', 'unknown')} at position {v.get('position', 'N/A')}"
                )

        # Warnings
        warnings = metadata.get("warnings", [])
        if warnings:
            comment_lines.append(f"Warnings: {len(warnings)}")
            for w in warnings[:3]:  # Show at most 3
                comment_lines.append(f"  - {w.get('message', 'N/A')}")

        record.annotations["comment"] = "\n".join(comment_lines)

        # Add CDS feature
        if metadata.get("protein_seq"):
            cds_feature = SeqFeature(  # type: ignore[no-untyped-call]
                FeatureLocation(0, len(sequence)),  # type: ignore[no-untyped-call]
                type="CDS",
                qualifiers={
                    "codon_opt": ["Nicotiana benthamiana"],
                    "translation": [metadata["protein_seq"]],
                    "note": [
                        f"CAI={metadata.get('cai', 0.0):.3f}, GC={metadata.get('gc', 0.0):.1f}%"
                    ],
                    "gene": [gene_name],
                },
            )
            record.features.append(cds_feature)

        # Additional feature annotations (promoter, terminator, etc.)
        if metadata.get("features"):
            for feat in metadata["features"]:
                feature = SeqFeature(  # type: ignore[no-untyped-call]
                    FeatureLocation(feat["start"], feat["end"]),  # type: ignore[no-untyped-call]
                    type=feat["type"],
                    qualifiers=feat.get("qualifiers", {}),
                )
                record.features.append(feature)

        # Write file or return string
        if output_file:
            SeqIO.write(record, output_file, "genbank")
            return f"GenBank file written to {output_file}"
        else:
            output = StringIO()
            SeqIO.write(record, output, "genbank")
            return output.getvalue()

    def export_fasta(
        self,
        sequence: str,
        metadata: dict[str, Any],
        output_file: str | None = None,
        line_width: int = 80,
    ) -> str:
        """
        Export in FASTA format

        Args:
            sequence: Optimized DNA sequence
            metadata: Metadata (same as export_genbank)
            output_file: Output file path (returns string if None)
            line_width: Line wrap width (0 for no wrapping)

        Returns:
            FASTA-formatted string

        Raises:
            None.

        Examples:
            >>> exporter = SequenceExporter()
            >>> fasta = exporter.export_fasta("ATG", {"gene_name": "x"})
            >>> fasta.startswith(">")
            True
        """
        # Set defaults
        run_id = metadata.get("run_id", self.generate_run_id(sequence, metadata))
        gene_name = metadata.get("gene_name", "optimized_gene")

        # Build header
        header_parts = [
            f"PFORM_{run_id}",
            f"gene={gene_name}",
            f"CAI={metadata.get('cai', 0.0):.3f}",
            f"GC={metadata.get('gc', 0.0):.1f}",
            f"profile={metadata.get('profile', 'N/A')}",
        ]

        if metadata.get("assembly_standard"):
            header_parts.append(f"assembly={metadata['assembly_standard']}")

        header = ">{}".format("|".join(header_parts))

        # Wrap sequence
        if line_width > 0:
            seq_lines = [sequence[i : i + line_width] for i in range(0, len(sequence), line_width)]
            seq_formatted = "\n".join(seq_lines)
        else:
            seq_formatted = sequence

        fasta_content = f"{header}\n{seq_formatted}\n"

        # Write file or return string
        if output_file:
            with open(output_file, "w") as f:
                f.write(fasta_content)
            return f"FASTA file written to {output_file}"
        else:
            return fasta_content

    def export_batch(
        self,
        sequences: list[dict[str, Any]],
        output_format: str = "fasta",
        output_file: str | None = None,
    ) -> str:
        """
        Export batch sequences

        Args:
            sequences: [{"sequence": "ATG...", "metadata": {...}}, ...]
            output_format: "fasta" or "genbank"
            output_file: Output file path

        Returns:
            Output message

        Raises:
            ValueError: Unsupported format or missing output_file for GenBank batch.

        Examples:
            >>> exporter = SequenceExporter()
            >>> msg = exporter.export_batch([{"sequence": "ATG", "metadata": {}}])
            >>> "sequence" in msg or "FASTA" in msg
            True
        """
        if output_format.lower() == "fasta":
            # FASTA allows multiple sequences in one file
            all_fasta = []
            for seq_data in sequences:
                fasta = self.export_fasta(seq_data["sequence"], seq_data.get("metadata", {}))
                all_fasta.append(fasta.strip())

            combined = "\n".join(all_fasta) + "\n"

            if output_file:
                with open(output_file, "w") as f:
                    f.write(combined)
                return f"Batch FASTA written to {output_file} ({len(sequences)} sequences)"
            else:
                return combined

        elif output_format.lower() == "genbank":
            # GenBank stores each sequence in a separate file
            if not output_file:
                raise ValueError("GenBank batch export requires output_file")

            import os

            base_name, ext = os.path.splitext(output_file)

            for i, seq_data in enumerate(sequences):
                file_name = f"{base_name}_{i + 1:03d}{ext}"
                self.export_genbank(
                    seq_data["sequence"], seq_data.get("metadata", {}), output_file=file_name
                )

            return f"Batch GenBank written ({len(sequences)} files)"

        else:
            raise ValueError(f"Unsupported format: {output_format}")

    def export_report(
        self,
        sequence: str,
        metadata: dict[str, Any],
        output_file: str | None = None,
    ) -> str:
        """
        Create a human-readable report

        Args:
            sequence: Optimized DNA sequence
            metadata: Metadata
            output_file: Output file path

        Returns:
            Report string

        Raises:
            None.

        Examples:
            >>> exporter = SequenceExporter()
            >>> report = exporter.export_report("ATG", {"gene_name": "x"})
            >>> "Optimization Report" in report
            True
        """
        run_id = metadata.get("run_id", self.generate_run_id(sequence, metadata))

        report_lines = [
            "=" * 70,
            "FactorForge v3.x - Optimization Report",
            "=" * 70,
            "",
            f"Run ID: {run_id}",
            f"Timestamp: {metadata.get('timestamp', datetime.now().isoformat())}",
            f"Gene: {metadata.get('gene_name', 'N/A')}",
            "",
            "--- Sequence Information ---",
            f"Length: {len(sequence)} bp",
            f"GC Content: {metadata.get('gc', 0.0):.1f}%",
            f"CAI: {metadata.get('cai', 0.0):.3f}",
            "",
            "--- Optimization Settings ---",
            f"Profile: {metadata.get('profile', 'N/A')}",
            f"Assembly Standard: {metadata.get('assembly_standard', 'None')}",
            f"Organism: {metadata.get('organism', 'Nicotiana benthamiana')}",
            "",
        ]

        # Violations fixed
        violations_fixed = metadata.get("violations_fixed", [])
        if violations_fixed:
            report_lines.append("--- Violations Fixed ---")
            for v in violations_fixed:
                report_lines.append(
                    f"  • {v.get('type', 'Unknown')} at position {v.get('position', 'N/A')}"
                )
                if v.get("fix_description"):
                    report_lines.append(f"    → {v['fix_description']}")
            report_lines.append("")

        # Warnings
        warnings = metadata.get("warnings", [])
        if warnings:
            report_lines.append("--- Warnings ---")
            for w in warnings:
                report_lines.append(f"  ⚠ {w.get('message', 'N/A')}")
                if w.get("suggestion"):
                    report_lines.append(f"    → {w['suggestion']}")
            report_lines.append("")

        # Quality score
        if metadata.get("quality_score"):
            report_lines.append("--- Quality Assessment ---")
            score = metadata["quality_score"]
            stars = "⭐" * int(score)
            report_lines.append(f"Overall Quality: {stars} ({score}/5)")
            report_lines.append("")

        # Sequence preview
        report_lines.append("--- Sequence Preview ---")
        preview_len = min(120, len(sequence))
        report_lines.append(sequence[:preview_len])
        if len(sequence) > preview_len:
            report_lines.append(f"... ({len(sequence) - preview_len} more bp)")
        report_lines.append("")

        report_lines.append("=" * 70)

        report_content = "\n".join(report_lines)

        if output_file:
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(report_content)
            return f"Report written to {output_file}"
        else:
            return report_content


# --- Usage example ---
if __name__ == "__main__":
    exporter = SequenceExporter()

    # Test data
    test_sequence = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"

    test_metadata = {
        "gene_name": "GFP",
        "protein_seq": "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHKVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK*",
        "profile": "Balanced",
        "cai": 0.87,
        "gc": 51.2,
        "assembly_standard": "Golden Gate (BsaI)",
        "violations_fixed": [
            {
                "type": "BsaI site",
                "position": 147,
                "fix_description": "Synonymous substitution R→R (CGT→AGA)",
            }
        ],
        "warnings": [
            {
                "message": "High local GC content at position 450-500",
                "suggestion": "Consider manual review",
            }
        ],
        "quality_score": 5,
    }

    print("=== FASTA Export ===")
    fasta = exporter.export_fasta(test_sequence, test_metadata)
    print(fasta[:200])

    print("\n=== Report Export ===")
    report = exporter.export_report(test_sequence, test_metadata)
    print(report)

    print("\n=== GenBank Export ===")
    try:
        genbank = exporter.export_genbank(test_sequence, test_metadata)
        print(genbank[:500])
    except ImportError as e:
        print(f"Biopython not installed: {e}")
