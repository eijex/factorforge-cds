"""Exercise experimental 4DP coordinate binding with synthetic inputs."""

import json
from pathlib import Path

import pytest
from Bio.Seq import Seq
from click.testing import CliRunner

from factorforge.cli.main import cli
from factorforge.structure.alphafold_server_ingest import parse_af3_cif

PROTEIN = "M" + "KIF" * 12 + "LMM"
COMPONENT = {"M": "MET", "K": "LYS", "I": "ILE", "F": "PHE", "L": "LEU"}


@pytest.fixture
def model_inputs(tmp_path: Path):
    fasta = tmp_path / "synthetic.fasta"
    fasta.write_text(f">synthetic\n{PROTEIN}\n", encoding="utf-8")
    cif = tmp_path / "synthetic.cif"
    cif.write_text(
        "data_test\nloop_\n_atom_site.label_atom_id\n_atom_site.label_asym_id\n"
        "_atom_site.label_seq_id\n_atom_site.label_comp_id\n_atom_site.B_iso_or_equiv\n"
        + "".join(
            f"CA A {i} {COMPONENT[aa]} {45 if i <= 3 else 90}\n"
            for i, aa in enumerate(PROTEIN, start=1)
        ),
        encoding="utf-8",
    )
    return fasta, cif, tmp_path / "model.json", tmp_path / "result.fasta"


def test_cli_ingest_and_optimize_4dp(model_inputs):
    fasta, cif, sidecar, output = model_inputs
    runner = CliRunner()
    ingest = runner.invoke(cli, ["structure", "ingest-af3", "-c", str(cif), "-o", str(sidecar)])
    assert ingest.exit_code == 0, ingest.output
    args = ["optimize", str(fasta), "--pause-sidecar", str(sidecar),
            "--gc-min", "0", "--gc-max", "100", "-o", str(output)]
    result = runner.invoke(cli, args)
    assert result.exit_code == 0, result.output
    assert "sequence-matched model offset positions: [38]" in result.output
    header, *dna = output.read_text(encoding="utf-8").splitlines()
    assert "engine=dp_4dp" in header and "pauses=38" in header
    assert str(Seq("".join(dna)).translate()).removesuffix("*") == PROTEIN
    original = output.read_bytes()
    assert runner.invoke(cli, args).exit_code != 0
    assert output.read_bytes() == original


@pytest.mark.parametrize(
    ("fault", "expected"),
    [
        ("different_protein", "does not match"),
        ("discontinuous_positions", "consecutive one-based"),
        ("tampered_position", "disagree with the residue data"),
        ("multiple_chains", "exactly one selected chain"),
        ("empty_candidates", "No in-bounds model offset"),
        ("non_dp", "only supported with --engine dp"),
        ("compare_profiles", "only supported with --engine dp"),
    ],
)
def test_cli_sidecar_rejects_unbound_and_unsupported_input(model_inputs, fault, expected):
    fasta, cif, sidecar, output = model_inputs
    data = parse_af3_cif(cif, tunnel_offset=100 if fault == "empty_candidates" else 35)
    chain = data["chains"]["A"]
    flags = []
    if fault == "different_protein":
        fasta.write_text(f">other\n{PROTEIN[:-1]}A\n", encoding="utf-8")
    elif fault == "discontinuous_positions":
        chain["residues"][1]["seq_id"] = 3
    elif fault == "tampered_position":
        chain["pause_candidates"][0]["tunnel_pause_position"] = 39
    elif fault == "multiple_chains":
        data["chains"]["B"] = dict(chain, chain_id="B")
    elif fault == "non_dp":
        flags = ["--engine", "profile"]
    elif fault == "compare_profiles":
        flags = ["--compare-profiles", "balanced"]
    sidecar.write_text(json.dumps(data), encoding="utf-8")
    result = CliRunner().invoke(
        cli, ["optimize", str(fasta), "--pause-sidecar", str(sidecar), "-o", str(output), *flags],
    )
    assert result.exit_code != 0
    assert expected in result.output
    assert not output.exists()
