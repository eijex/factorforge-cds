"""Synthetic direct 3DP tests with actual engines."""

import json
from io import BytesIO
from unittest.mock import patch

import pytest

from api.optimize import handler

PROTEIN = "MASTWQHLEKDPVNRSGYF"
NAMES = dict(
    zip(
        "MASTWQHLEKDPVNRGYF",
        ["MET", "ALA", "SER", "THR", "TRP", "GLN", "HIS", "LEU", "GLU", "LYS", "ASP", "PRO", "VAL", "ASN", "ARG", "GLY", "TYR", "PHE"],
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


def post(data, *, research_enabled=True):
    """Research fixtures explicitly opt in; gate tests pass research_enabled=False."""
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
    with patch("api.optimize.EXPERIMENTAL_3DP_ENABLED", research_enabled):
        h.do_POST()
    return codes[0], json.loads(h.wfile.getvalue())


@pytest.mark.parametrize("data", [
    {"design_mode": "paired_3dp", "sequence": PROTEIN},
    {"action": "recommend_structure", "sequence": PROTEIN},
    {"action": "auto_resolve_structure", "sequence": PROTEIN},
])
def test_public_structure_gate_blocks_before_work(monkeypatch, data):
    def forbidden(*args, **kwargs):
        raise AssertionError("Research optimizer must not run while disabled")

    monkeypatch.setattr(handler, "optimize_structure_web", forbidden)
    status, result = post(data, research_enabled=False)
    assert status == 403
    assert result["error_code"] == "EXPERIMENTAL_STRUCTURE_DISABLED"
    assert "candidates" not in result


def test_public_gate_keeps_baseline_available():
    status, result = post({
        "sequence": PROTEIN, "design_mode": "baseline_1dp",
        "objective": "feasibility_best", "constraints": {"gc_min": 0, "gc_max": 100},
    }, research_enabled=False)
    assert status == 200, result
    assert result["success"] is True


def test_public_capabilities_report_paused():
    h = object.__new__(handler)
    h.path = "/api/optimize"
    h.wfile = BytesIO()
    h.send_response = lambda *a: None
    h.send_header = lambda *a: None
    h.end_headers = lambda: None
    with patch("api.optimize.EXPERIMENTAL_3DP_ENABLED", False):
        h.do_GET()
    health = json.loads(h.wfile.getvalue())
    assert health["web_deployment_revision"] == f"v{health['version']}"
    capabilities = health["capabilities"]
    assert capabilities["direct_structure"]["available"] is False
    assert capabilities["direct_structure"]["status"] == "paused"
    assert capabilities["structure_recommendation"]["available"] is False


def test_public_capabilities_default_enable_3dp_and_4dp():
    h = object.__new__(handler)
    h.path = "/api/optimize"
    h.headers = {"Host": "factorforge.eijex.com"}
    h.wfile = BytesIO()
    h.send_response = lambda *a: None
    h.send_header = lambda *a: None
    h.end_headers = lambda: None
    h.do_GET()
    health = json.loads(h.wfile.getvalue())
    capabilities = health["capabilities"]
    assert capabilities["direct_structure"]["available"] is True
    assert capabilities["direct_structure"]["status"] == "experimental"
    assert capabilities["local_4dp"]["available"] is True
    assert capabilities["local_4dp"]["fixtures"] is False


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
    from io import StringIO

    from Bio.PDB import MMCIFIO, PDBParser

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
    "changes", [{"sequence": "M" * 501}, {"constraints": {"gc_min": float("nan"), "gc_max": 100}}]
)
def test_resource_and_nonfinite_gc_rejected(changes):
    assert post(request(**changes))[0] == 400


def test_structure_content_size_limit():
    data = request()
    data["structure_context"]["content"] = " " * 450001
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
