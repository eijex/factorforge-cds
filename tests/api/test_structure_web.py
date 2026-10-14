"""Synthetic direct 3DP tests with actual engines."""

from io import BytesIO
import json
import pytest
from api.optimize import handler

PROTEIN = "MASTWQHLEKDPVNRSGYF"
NAMES = dict(
    zip(
        "MASTWQHLEKDPVNRGYF",
        "MET ALA SER THR TRP GLN HIS LEU GLU LYS ASP PRO VAL ASN ARG GLY TYR PHE".split(),
    )
)


def pdb(protein=PROTEIN):
    return (
        "".join(
            f"ATOM  {i:5d}  CA  {NAMES[a]:3s} A{i:4d}    {float(i):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C  \n"
            for i, a in enumerate(protein, 1)
        )
        + "END\n"
    )


def request(**changes):
    data = {
        "sequence": PROTEIN,
        "host": "nbenthamiana",
        "design_mode": "paired_3dp",
        "constraints": {"gc_min": 0, "gc_max": 100},
        "structure_context": {
            "content": pdb(),
            "format": "pdb",
            "chain": "A",
            "model": 1,
            "positions": [3, 4, 5],
            "positions_reviewed": True,
            "target_adaptiveness": 0.4,
            "penalty_weight": 2.0,
        },
    }
    data.update(changes)
    return data


def post(data):
    h = object.__new__(handler)
    body = json.dumps(data).encode()
    h.path = "/api/optimize"
    h.headers = {"Content-Length": str(len(body))}
    h.rfile = BytesIO(body)
    h.wfile = BytesIO()
    codes = []
    h.send_response = codes.append
    h.send_header = lambda *a: None
    h.end_headers = lambda: None
    h.do_POST()
    return codes[0], json.loads(h.wfile.getvalue())


def test_actual_paired_generation():
    status, r = post(request())
    assert status == 200, r
    assert [c["id"] for c in r["candidates"]] == ["baseline_1dp", "structure_3dp"]
    assert r["structure_comparison"]["hard_constraints_pass"]
    assert r["metadata"]["approval_status"] == "not_collected"
    assert r["metadata"]["evidence_recording"] == "not_connected"
    assert "ATOM" not in json.dumps(r)


@pytest.mark.parametrize(
    "field,value",
    [
        ("positions_reviewed", False),
        ("positions", [True]),
        ("positions", [999]),
        ("positions", [3, 3]),
        ("model", True),
        ("model", 2),
        ("chain", "B"),
        ("penalty_weight", float("nan")),
        ("target_adaptiveness", 0),
        ("format", "url"),
        ("content", "invalid"),
        ("content", pdb(PROTEIN[:-1])),
    ],
)
def test_invalid_structure_no_fallback(field, value):
    d = request()
    d["structure_context"][field] = value
    status, r = post(d)
    assert status == 400, r
    assert "candidates" not in r


@pytest.mark.parametrize(
    "changes",
    [
        {"structure_context": None},
        {"save_db": True},
        {"host": "by2"},
        {"kozak": True},
        {"mode": "slm"},
        {"design_mode": "automatic_best"},
    ],
)
def test_reject_unsupported(changes):
    assert post(request(**changes))[0] == 400


def test_zero_penalty_equivalent():
    d = request()
    d["structure_context"]["penalty_weight"] = 0
    status, r = post(d)
    assert status == 200, r
    assert r["candidates"][0]["dna_sequence"] == r["candidates"][1]["dna_sequence"]


def test_custom_motif():
    assert (
        post(
            request(
                custom_restriction_sites=[
                    {"name": "Blocked start", "sequence": "ATG", "scan_rc": True}
                ]
            )
        )[0]
        == 400
    )


def test_cds_stop_preserved():
    status, r = post(request())
    assert status == 200, r
    status, r = post(request(sequence=r["candidates"][0]["dna_sequence"] + "TAA"))
    assert status == 200, r
    assert all(c["dna_sequence"].endswith("TAA") for c in r["candidates"])


def test_mmcif_input_actual_execution():
    from Bio.PDB import PDBParser, MMCIFIO
    from io import StringIO

    structure = PDBParser(QUIET=True).get_structure("synthetic", StringIO(pdb()))
    writer = MMCIFIO()
    writer.set_structure(structure)
    output = StringIO()
    writer.save(output)
    data = request()
    data["structure_context"].update(content=output.getvalue(), format="mmcif")
    status, result = post(data)
    assert status == 200, result
    assert result["metadata"]["structure_mapping"]["format"] == "mmcif"


@pytest.mark.parametrize(
    "changes", [{"sequence": "M" * 301}, {"constraints": {"gc_min": float("nan"), "gc_max": 100}}]
)
def test_resource_and_nonfinite_gc_rejected(changes):
    assert post(request(**changes))[0] == 400


def test_structure_content_size_limit():
    data = request()
    data["structure_context"]["content"] = " " * 250001
    assert post(data)[0] == 400


def test_request_body_size_limit():
    data = request()
    data["structure_context"]["content"] = " " * 1000001
    assert post(data)[0] == 413


def test_missing_structure_never_auto_approves_or_calls_network(monkeypatch):
    import urllib.request
    def forbidden(*args, **kwargs):
        raise AssertionError("Unexpected external sequence transfer")
    monkeypatch.setattr(urllib.request, "urlopen", forbidden)
    data = request()
    data["structure_context"] = {"chain": "A", "positions": [3, 4, 5]}
    status, result = post(data)
    assert status != 200
    assert "NOT_AVAILABLE" in str(result)


def test_local_uploaded_resolver_keeps_review_false_and_checks_mapping(monkeypatch):
    import urllib.request
    from factorforge.structure.auto_resolver import auto_resolve_structure_for_protein
    def forbidden(*args, **kwargs):
        raise AssertionError("Unexpected external sequence transfer")
    monkeypatch.setattr(urllib.request, "urlopen", forbidden)
    resolved = auto_resolve_structure_for_protein(PROTEIN, pdb(), "pdb")
    assert resolved["positions_reviewed"] is False
    assert resolved["protein_sha256"]
    assert "provider unverified" in resolved["source_label"]
    with pytest.raises(ValueError, match="NOT_AVAILABLE"):
        auto_resolve_structure_for_protein(PROTEIN)
    with pytest.raises(ValueError, match="match"):
        auto_resolve_structure_for_protein("M" * len(PROTEIN), pdb(), "pdb")
