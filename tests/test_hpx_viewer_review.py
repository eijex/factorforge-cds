"""Portable evidence checks; do not certify biological annotations or speed curves."""

import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "examples/hpx_evidence_viewer"


def test_local_structure_binding_and_confidence():
    data = json.loads((ROOT / "reviewed-data.json").read_text(encoding="utf-8"))
    raw = (ROOT / "reference.pdb").read_bytes()
    assert hashlib.sha256(raw).hexdigest() == data["meta"]["reference_pdb_sha256"]
    ca = [
        line
        for line in raw.decode().splitlines()
        if line.startswith("ATOM") and line[12:16].strip() == "CA"
    ]
    assert (
        len(ca)
        == len(data["sequence"])
        == len(data["structural_confidence"])
        == len(data["ca_atoms"])
    )
    for line, atom, confidence in zip(ca, data["ca_atoms"], data["structural_confidence"]):
        assert int(line[22:26]) == atom["num"]
        assert abs(float(line[60:66]) - confidence) < 1e-5
        for start, end, key in [(30, 38, "x"), (38, 46, "y"), (46, 54, "z")]:
            assert abs(float(line[start:end]) - atom[key]) < 1e-3


def test_unbound_profiles_are_unavailable():
    data = json.loads((ROOT / "reviewed-data.json").read_text(encoding="utf-8"))
    assert data["meta"]["annotation_provenance"].endswith("UNVERIFIED")
    for profile in data["candidate_profiles"].values():
        assert profile["status"] == "NOT_AVAILABLE"
        assert all(value is None for value in profile["wi"])


def test_embedded_data_matches_reviewed_snapshot_and_no_external_assets():
    html = (ROOT / "reviewed-index.html").read_text(encoding="utf-8")
    embedded, _ = json.JSONDecoder().raw_decode(html.split("const HPX_DATA = ", 1)[1])
    assert embedded == json.loads((ROOT / "reviewed-data.json").read_text(encoding="utf-8"))
    assert not re.search(r'<(?:script|img|iframe)[^>]+src=["\'](?:https?:)?//', html, re.IGNORECASE)
    assert not re.search(r"@import|fetch\(|XMLHttpRequest|WebSocket", html)
    assert "not Mol*" in html


def test_shift_intervals_are_bound_and_clipped():
    data = json.loads((ROOT / "reviewed-data.json").read_text(encoding="utf-8"))
    spec = data["design_hypotheses"]["shifted_4dp"]
    n = len(data["sequence"])
    for interval in spec["intervals"]:
        center = interval["boundary"] + spec["shift_aa"]
        assert interval["shifted_window"] == [
            max(1, center - spec["half_width"]),
            min(n, center + spec["half_width"]),
        ]
