"""End-to-end API -> RNA QC alignment tests across 1DP, 3DP, and 4DP (Job 640)."""

import json
from io import BytesIO
from unittest.mock import patch

from Bio.Seq import Seq

from api.optimize import handler
from factorforge.guards.plant_rna_guards import DETECTOR_VERSION, scan_plant_rna_risks
from tests.api.test_4dp_web import META_PROTEIN, meta_sample_sidecar
from tests.api.test_structure_web import PROTEIN as STRUCT_PROTEIN
from tests.api.test_structure_web import pdb as make_pdb


def _post_optimize(data: dict, *, research_enabled: bool = True) -> tuple[int, dict]:
    h = object.__new__(handler)
    body = json.dumps(data).encode()
    h.path = "/api/optimize"
    h.client_address = ("127.0.0.1", 1234)
    h.headers = {"Content-Length": str(len(body)), "Host": "127.0.0.1"}
    h.rfile = BytesIO(body)
    h.wfile = BytesIO()
    codes: list[int] = []
    h.send_response = codes.append
    h.send_header = lambda *a: None
    h.end_headers = lambda: None
    with patch("api.optimize.EXPERIMENTAL_3DP_ENABLED", research_enabled):
        h.do_POST()
    return codes[0], json.loads(h.wfile.getvalue())


def _coord_triples(scan_dict: dict) -> list[tuple[str, int, int]]:
    return [(f["rule_id"], f["start_0based"], f["end_exclusive"]) for f in scan_dict["findings"]]


def _assert_candidate_qc_alignment(
    response: dict,
    expected_protein: str,
    *,
    expected_primary_id: str | None = None,
) -> None:
    assert response["success"] is True
    candidates = response["candidates"]
    assert len(candidates) >= 1

    for cand in candidates:
        dna = cand["dna_sequence"]
        assert len(dna) % 3 == 0
        assert str(Seq(dna).translate()).removesuffix("*") == expected_protein.rstrip("*")

        rna_risks = cand["scientific_review"]["rna_risks"]
        expected = scan_plant_rna_risks(dna).to_dict()

        assert rna_risks["detector_version"] == DETECTOR_VERSION == "1.1.0"
        assert rna_risks["summary_status"] in {"CLEAN_SCAN", "ADVISORY_FINDINGS"}
        assert rna_risks["summary_status"] == expected["summary_status"]
        assert rna_risks["donor_count"] == expected["donor_count"]
        assert rna_risks["polya_count"] == expected["polya_count"]
        assert rna_risks["are_count"] == expected["are_count"]
        assert _coord_triples(rna_risks) == _coord_triples(expected)

    primary = response["recommended_candidate"]
    if expected_primary_id is not None:
        assert primary["id"] == expected_primary_id
    assert primary["dna_sequence"] == response["optimized_sequence"]

    primary_expected = scan_plant_rna_risks(primary["dna_sequence"]).to_dict()
    top_rna = response["scientific_review"]["rna_risks"]
    assert top_rna["donor_count"] == primary_expected["donor_count"]
    assert top_rna["polya_count"] == primary_expected["polya_count"]
    assert top_rna["are_count"] == primary_expected["are_count"]
    assert _coord_triples(top_rna) == _coord_triples(primary_expected)


def test_1dp_api_qc_matches_independent_rescan() -> None:
    status, res = _post_optimize(
        {
            "sequence": STRUCT_PROTEIN,
            "host": "nbenthamiana",
            "design_mode": "baseline_1dp",
            "objective": "feasibility_best",
            "constraints": {"gc_min": 0, "gc_max": 100},
        },
        research_enabled=False,
    )
    assert status == 200, res
    _assert_candidate_qc_alignment(res, STRUCT_PROTEIN)


def test_3dp_api_qc_matches_recommended_candidate_not_baseline() -> None:
    status, res = _post_optimize(
        {
            "sequence": STRUCT_PROTEIN,
            "host": "nbenthamiana",
            "design_mode": "paired_3dp",
            "constraints": {"gc_min": 0, "gc_max": 100},
            "structure_context": {
                "content": make_pdb(STRUCT_PROTEIN),
                "format": "pdb",
                "chain": "A",
                "model": 1,
                "positions": [3, 4, 5],
                "positions_reviewed": True,
                "target_adaptiveness": 0.4,
                "penalty_weight": 2.0,
            },
        },
        research_enabled=True,
    )
    assert status == 200, res
    assert [c["id"] for c in res["candidates"]] == ["baseline_1dp", "structure_3dp"]
    _assert_candidate_qc_alignment(res, STRUCT_PROTEIN, expected_primary_id="structure_3dp")


def test_4dp_api_qc_matches_recommended_candidate_not_baseline(monkeypatch) -> None:
    monkeypatch.delenv("VERCEL", raising=False)
    monkeypatch.setenv("FACTORFORGE_LOCAL_4DP_FIXTURES", "true")
    sidecar = meta_sample_sidecar()
    status, res = _post_optimize(
        {
            "sequence": META_PROTEIN,
            "design_mode": "paired_4dp",
            "pause_sidecar": sidecar,
            "positions_reviewed": True,
            "allow_mock": True,
            "constraints": {"gc_min": 35, "gc_max": 55},
        },
        research_enabled=True,
    )
    assert status == 200, res
    assert [c["id"] for c in res["candidates"]] == ["baseline_1dp", "structure_4dp"]
    _assert_candidate_qc_alignment(res, META_PROTEIN, expected_primary_id="structure_4dp")


def test_flanked_construct_marks_rna_qc_not_available() -> None:
    h = object.__new__(handler)
    response = {
        "optimized_sequence": "GGG" + "ATGGCATCAACTTGA" + "CCC",
        "candidates": [
            {"id": "flanked", "dna_sequence": "GGG" + "ATGGCATCAACTTGA" + "CCC"},
        ],
    }
    observed = h.attach_scientific_observations(response, expected_protein="MAST")
    assert observed["scientific_review"]["rna_risks"]["summary_status"] == "NOT_AVAILABLE"
    assert observed["candidates"][0]["scientific_review"]["rna_risks"]["summary_status"] == "NOT_AVAILABLE"
