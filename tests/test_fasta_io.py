"""Native FASTA I/O and public-header privacy tests."""

import pytest

from factorforge.io.fasta import FastaRecord, format_fasta, parse_fasta, read_fasta, write_fasta


def test_parse_multirecord_fasta() -> None:
    records = parse_fasta(">a\nATG\nAAA\n>b\nCCC\n", validation_mode="dna_strict")
    assert [(record.identifier, record.sequence) for record in records] == [
        ("a", "ATGAAA"),
        ("b", "CCC"),
    ]


def test_write_and_read_round_trip(tmp_path) -> None:
    path = tmp_path / "designs.fasta"
    records = [
        FastaRecord(
            "CF-001",
            "ATGAAACCC",
            {"host_profile": "nbenthamiana", "profile": "balanced"},
        )
    ]
    write_fasta(path, records, line_width=6)
    assert read_fasta(path, validation_mode="dna_strict")[0].sequence == "ATGAAACCC"


def test_header_contains_only_allowlisted_metadata() -> None:
    output = format_fasta(
        [
            FastaRecord(
                "CF-001",
                "ATGAAACCC",
                {
                    "host_profile": "nbenthamiana",
                    "profile": "balanced",
                    "long_private_notes": "must stay in Design Package JSON",
                },
            )
        ]
    )
    header = output.splitlines()[0]
    assert "host_profile=nbenthamiana" in header
    assert "profile=balanced" in header
    assert "long_private_notes" not in header
    assert "must stay" not in header


@pytest.mark.parametrize(
    "identifier",
    [
        "PlantForm-confidential-001",
        "private-partner-id",
        "clinical-yield-claim",
        "ATGAAACCCATGAAACCCATG",
    ],
)
def test_sensitive_or_raw_sequence_header_is_rejected(identifier: str) -> None:
    with pytest.raises(ValueError):
        format_fasta([FastaRecord(identifier, "ATGAAACCC")])


def test_raw_sequence_is_not_written_to_header() -> None:
    sequence = "ATGAAACCC" * 4
    output = format_fasta([FastaRecord("CF-001", sequence)])
    assert sequence not in output.splitlines()[0]
