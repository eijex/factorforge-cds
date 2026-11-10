"""Unit tests for offline AlphaFold Server (alphafoldserver.com) job JSON export."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from factorforge.cli.main import cli
from factorforge.structure.alphafold_server_export import (
    AlphaFoldExportError,
    build_alphafold_server_job_payload,
    export_alphafold_server_job,
    parse_protein_fasta,
    sanitize_job_name,
    validate_protein_sequence,
)

SAMPLE_PROTEIN_SEQ_1 = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSE"
SAMPLE_PROTEIN_SEQ_2 = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWV"


def test_validate_protein_sequence_success():
    """Verify normal sequence validation and whitespace/casing normalization."""
    raw = "  mkwvt fisllflfss aysrgvfrrd thkse\n "
    clean = validate_protein_sequence(raw, "TestChain")
    assert clean == SAMPLE_PROTEIN_SEQ_1
    assert all(c in "ACDEFGHIKLMNPQRSTVWY" for c in clean)


def test_validate_protein_sequence_failures():
    """Verify boundaries and non-standard amino acid rejection."""
    # Empty
    with pytest.raises(AlphaFoldExportError, match="empty"):
        validate_protein_sequence("", "EmptyChain")

    with pytest.raises(AlphaFoldExportError, match="below minimum 4 aa"):
        validate_protein_sequence("MKW", "ShortChain")

    # Non-standard amino acids (X, B, Z, J, U)
    with pytest.raises(AlphaFoldExportError, match="invalid or non-standard amino acid"):
        validate_protein_sequence("MKWVTXISLL", "BadChain")

    with pytest.raises(AlphaFoldExportError, match="invalid or non-standard amino acid"):
        validate_protein_sequence("MKW123ISLL", "NumChain")


def test_sanitize_job_name():
    """Verify safe job naming without leaking paths or special characters."""
    assert sanitize_job_name(None) == "factorforge_job"
    assert sanitize_job_name("") == "factorforge_job"
    assert sanitize_job_name("C:/Users/secret/path/my_mAb_candidate.fasta") == "my_mAb_candidate_fasta"
    assert sanitize_job_name("complex @ test #1!") == "complex_test_1"


def test_parse_protein_fasta(tmp_path: Path):
    """Verify parsing single and multi-record FASTA files, stripping headers."""
    fasta_path = tmp_path / "antibodies.fasta"
    fasta_path.write_text(
        ">Heavy_Chain_Proprietary_Clone_42 | internal_id=999\n"
        f"{SAMPLE_PROTEIN_SEQ_1}\n"
        ">Light_Chain_Clone_42\n"
        f"{SAMPLE_PROTEIN_SEQ_2}\n",
        encoding="utf-8",
    )

    records = parse_protein_fasta(fasta_path)
    assert len(records) == 2
    assert records[0] == SAMPLE_PROTEIN_SEQ_1
    assert records[1] == SAMPLE_PROTEIN_SEQ_2
    # Ensure header text is nowhere in the sequence output
    assert "Proprietary" not in records[0]
    assert "internal_id" not in records[0]


def test_build_alphafold_server_job_payload_schema():
    """Verify payload complies strictly with AlphaFold Server official schema."""
    chains = [
        (SAMPLE_PROTEIN_SEQ_1, 2),  # (sequence, count)
        {"sequence": SAMPLE_PROTEIN_SEQ_2, "count": 1},  # dict
    ]
    payload = build_alphafold_server_job_payload(
        chains=chains,
        job_name="antibody_complex",
        model_seeds=42,
    )

    assert isinstance(payload, list)
    assert len(payload) == 1
    job = payload[0]
    assert job["name"] == "antibody_complex"
    assert job["modelSeeds"] == ["42"]
    assert job["dialect"] == "alphafoldserver"
    assert job["version"] == 1
    assert "sequences" in job
    assert len(job["sequences"]) == 2

    # Check chain 1
    c1 = job["sequences"][0]
    assert "proteinChain" in c1
    assert c1["proteinChain"]["sequence"] == SAMPLE_PROTEIN_SEQ_1
    assert c1["proteinChain"]["count"] == 2

    # Check chain 2
    c2 = job["sequences"][1]
    assert "proteinChain" in c2
    assert c2["proteinChain"]["sequence"] == SAMPLE_PROTEIN_SEQ_2
    assert c2["proteinChain"]["count"] == 1


def test_build_alphafold_server_job_payload_errors():
    """Verify invalid inputs raise descriptive errors."""
    with pytest.raises(AlphaFoldExportError, match="At least one protein chain is required"):
        build_alphafold_server_job_payload([])

    with pytest.raises(AlphaFoldExportError, match="count must be >= 1"):
        build_alphafold_server_job_payload([(SAMPLE_PROTEIN_SEQ_1, 0)])

    with pytest.raises(AlphaFoldExportError, match="missing required 'sequence' key"):
        build_alphafold_server_job_payload([{"count": 1}])


def test_export_alphafold_server_job_file_io(tmp_path: Path):
    """Verify file writing and overwrite protections."""
    out_file = tmp_path / "af3_job.json"
    chains = [SAMPLE_PROTEIN_SEQ_1]

    # First write
    payload, _ = export_alphafold_server_job(
        chains=chains,
        output_path=out_file,
        job_name="test_io",
        overwrite=False,
    )
    assert out_file.exists()
    loaded = json.loads(out_file.read_text(encoding="utf-8"))
    assert loaded == payload

    # Re-write without overwrite flag -> FileExistsError
    with pytest.raises(FileExistsError, match="already exists"):
        export_alphafold_server_job(
            chains=chains,
            output_path=out_file,
            job_name="test_io",
            overwrite=False,
        )

    # Re-write with overwrite=True -> success
    export_alphafold_server_job(
        chains=chains,
        output_path=out_file,
        job_name="test_io_overwritten",
        overwrite=True,
    )
    loaded_new = json.loads(out_file.read_text(encoding="utf-8"))
    assert loaded_new[0]["name"] == "test_io_overwritten"


def test_cli_export_af3_single_sequence():
    """Test CLI export with direct --sequence option to stdout."""
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "structure",
            "export-af3",
            "-s",
            SAMPLE_PROTEIN_SEQ_1,
            "-n",
            "cli_test_job",
        ],
    )
    assert result.exit_code == 0
    parsed = json.loads(result.output)
    assert isinstance(parsed, list)
    assert parsed[0]["name"] == "cli_test_job"
    assert parsed[0]["sequences"][0]["proteinChain"]["sequence"] == SAMPLE_PROTEIN_SEQ_1
    assert parsed[0]["sequences"][0]["proteinChain"]["count"] == 1


def test_cli_export_af3_multi_chain_fasta(tmp_path: Path):
    """Test CLI export with multi-chain FASTA and explicit copy counts to file."""
    fasta_file = tmp_path / "mab.fasta"
    fasta_file.write_text(
        f">HeavyChain\n{SAMPLE_PROTEIN_SEQ_1}\n>LightChain\n{SAMPLE_PROTEIN_SEQ_2}\n",
        encoding="utf-8",
    )
    out_file = tmp_path / "af3_request.json"

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "structure",
            "export-af3",
            "-f",
            str(fasta_file),
            "-c",
            "2,2",
            "-n",
            "mab_heterotetramer",
            "-o",
            str(out_file),
        ],
    )
    assert result.exit_code == 0
    assert "saved to" in result.output
    assert "UNVERIFIED_PUBLIC_CORPUS" in result.output

    data = json.loads(out_file.read_text(encoding="utf-8"))
    assert len(data[0]["sequences"]) == 2
    assert data[0]["sequences"][0]["proteinChain"]["count"] == 2
    assert data[0]["sequences"][1]["proteinChain"]["count"] == 2


def test_cli_export_af3_overwrite_protection(tmp_path: Path):
    """Test CLI overwrite rejection and explicit flag bypass."""
    out_file = tmp_path / "existing.json"
    out_file.write_text("{}", encoding="utf-8")

    runner = CliRunner()
    # Without --overwrite -> fails
    res_fail = runner.invoke(
        cli,
        ["structure", "export-af3", "-s", SAMPLE_PROTEIN_SEQ_1, "-o", str(out_file)],
    )
    assert res_fail.exit_code != 0
    assert "already exists" in res_fail.output

    # With --overwrite -> succeeds
    res_ok = runner.invoke(
        cli,
        ["structure", "export-af3", "-s", SAMPLE_PROTEIN_SEQ_1, "-o", str(out_file), "--overwrite"],
    )
    assert res_ok.exit_code == 0
    loaded = json.loads(out_file.read_text(encoding="utf-8"))
    assert loaded[0]["sequences"][0]["proteinChain"]["sequence"] == SAMPLE_PROTEIN_SEQ_1


def test_cli_export_af3_invalid_input():
    """Test CLI error handling for missing options and bad sequences."""
    runner = CliRunner()

    # Neither fasta nor sequence
    res_empty = runner.invoke(cli, ["structure", "export-af3"])
    assert res_empty.exit_code != 0
    assert "exactly one" in res_empty.output

    # Bad sequence
    res_bad_seq = runner.invoke(cli, ["structure", "export-af3", "-s", "MKW123X"])
    assert res_bad_seq.exit_code != 0
    assert "validation error" in res_bad_seq.output

    # Bad counts
    res_bad_count = runner.invoke(
        cli, ["structure", "export-af3", "-s", SAMPLE_PROTEIN_SEQ_1, "-c", "0"]
    )
    assert res_bad_count.exit_code != 0
    assert "Invalid chain count" in res_bad_count.output


def test_provider_token_and_seed_boundaries():
    chains = [("A" * 2500, 2)]
    job = build_alphafold_server_job_payload(chains, model_seeds=[0, 4294967295])[0]
    assert job["modelSeeds"] == ["0", "4294967295"]
    assert build_alphafold_server_job_payload(chains, model_seeds=[])[0]["modelSeeds"] == []
    with pytest.raises(AlphaFoldExportError, match="including copies"):
        build_alphafold_server_job_payload([*chains, "ACDE"])
    for seed in (-1, 4294967296, True, [0.5]):
        with pytest.raises(AlphaFoldExportError, match="Model seeds"):
            build_alphafold_server_job_payload(chains, model_seeds=seed)
    for count in (1.5, True):
        with pytest.raises(AlphaFoldExportError, match="integer"):
            build_alphafold_server_job_payload([("ACDE", count)])


@pytest.mark.parametrize("text", [">empty\n>next\nACDE", ">seq\nACDE\n>empty", "ACDE\n>seq\nACDE"])
def test_fasta_rejects_missing_sequence_records(text):
    with pytest.raises(AlphaFoldExportError):
        parse_protein_fasta(text)


def test_cli_export_af3_does_not_include_fasta_name_or_headers(tmp_path: Path):
    fasta = tmp_path / "private-filename.fasta"
    fasta.write_text(">private-header\nACDE\n>second\nFGHI", encoding="utf-8-sig")
    runner = CliRunner()
    base = ["structure", "export-af3", "-f", str(fasta)]
    result = runner.invoke(cli, base)
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)[0]["name"] == "factorforge_af3_job"
    assert "private-" not in result.stdout
    for counts in ("1", "1,2,3", "1,,2", ""):
        result = runner.invoke(cli, [*base, "-c", counts])
        assert result.exit_code == 2
        assert "Invalid chain count" in result.output
    result = runner.invoke(cli, [*base, "-s", "ACDE"])
    assert result.exit_code == 2
    original = fasta.read_bytes()
    result = runner.invoke(cli, [*base, "-o", str(fasta), "--overwrite"])
    assert result.exit_code == 2
    assert fasta.read_bytes() == original


def test_cli_export_af3_reports_invalid_fasta_encoding(tmp_path: Path):
    fasta = tmp_path / "bad.fasta"
    fasta.write_bytes(b">seq\n\xff")
    result = CliRunner().invoke(cli, ["structure", "export-af3", "-f", str(fasta)])
    assert result.exit_code == 1
    assert "Could not read the FASTA" in result.output
