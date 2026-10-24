"""Unit tests for Hemopexin (HPX, P02790) Evidence and Design Assumption Viewer."""

import json
from pathlib import Path
import re
import pytest

VIEWER_DIR = Path(__file__).resolve().parents[1] / "examples" / "hpx_evidence_viewer"
JSON_PATH = VIEWER_DIR / "viewer-data.json"
HTML_PATH = VIEWER_DIR / "index.html"


@pytest.fixture(scope="module")
def viewer_data():
    assert JSON_PATH.is_file(), f"Viewer data JSON missing at {JSON_PATH}"
    with open(JSON_PATH, "r", encoding="utf-8") as f:
        return json.load(f)


def test_hpx_sequence_and_coordinate_mapping(viewer_data):
    seq = viewer_data["sequence"]
    assert len(seq) == 462, f"Expected 462 residues for full precursor HPX, got {len(seq)}"
    assert seq.startswith("MARVLGAPVALGLWSLCWSLAIA"), "Signal peptide sequence mismatch"

    # Construct mapping check
    mapping = viewer_data["construct_mapping"]
    assert mapping["uniprot_range"] == [1, 462]
    assert mapping["precursor_range"] == [1, 462]
    assert mapping["mature_range"] == [24, 462]
    assert mapping["mature_length"] == 439
    assert mapping["signal_peptide"] == [1, 23]


def test_curated_features_and_provenance(viewer_data):
    features = viewer_data["curated_features"]

    # Domains
    domains = features["domains"]
    assert len(domains) == 10
    d1 = next(d for d in domains if d["name"] == "Hemopexin Repeat 4")
    assert d1["start"] == 185 and d1["end"] == 231
    hinge = next(d for d in domains if d["type"] == "hinge")
    assert hinge["start"] == 232 and hinge["end"] == 258

    # Disulfides
    disulfides = features["disulfides"]
    assert any(ds["start"] == 50 and ds["end"] == 231 for ds in disulfides), "Missing Domain 1 closure (50-231)"
    assert any(ds["start"] == 257 and ds["end"] == 460 for ds in disulfides), "Missing Domain 2 closure (257-460)"

    # Glycosylation
    glyco_positions = {g["pos"] for g in features["glycosylations"]}
    assert {24, 29, 64, 187, 240, 246, 453}.issubset(glyco_positions)


def test_alphafold_plddt_confidence_integrity(viewer_data):
    plddts = viewer_data["structural_confidence"]
    assert len(plddts) == 462
    assert all(0.0 <= p <= 100.0 for p in plddts)

    # Metadata must contain EMBL-EBI disclaimer
    meta = viewer_data["meta"]
    assert "plddt_disclaimer" in meta
    assert "local prediction confidence" in meta["plddt_disclaimer"].lower()
    assert "EMBL-EBI" in meta["plddt_disclaimer"]


def test_4dp_exit_tunnel_hypothesis_and_out_of_bounds_exclusion(viewer_data):
    design = viewer_data["design_hypotheses"]["shifted_4dp"]
    assert design["shift_aa"] == 35
    assert design["half_width"] == 2

    intervals = design["intervals"]
    assert len(intervals) == 2
    # Domain 1 clearance: 231 + 35 +- 2 -> [264, 268]
    assert intervals[0]["boundary"] == 231
    assert intervals[0]["shifted_window"] == [264, 268]
    # Linker clearance: 258 + 35 +- 2 -> [291, 295]
    assert intervals[1]["boundary"] == 258
    assert intervals[1]["shifted_window"] == [291, 295]

    # Out of bounds excluded
    excluded = design["out_of_bounds_excluded"]
    assert len(excluded) >= 1
    assert excluded[0]["theoretical_window"][0] > 462


def test_candidate_profiles_represent_relative_adaptability(viewer_data):
    profiles = viewer_data["candidate_profiles"]
    assert "1dp" in profiles
    assert "3dp" in profiles
    assert "4dp" in profiles

    w1 = profiles["1dp"]["wi"]
    w3 = profiles["3dp"]["wi"]
    w4 = profiles["4dp"]["wi"]
    assert len(w1) == len(w3) == len(w4) == 462

    # In 1DP, the mature chain (24..462) has high average adaptation (mean > 0.80)
    assert sum(w1[23:]) / len(w1[23:]) >= 0.80
    # In 1DP, the hinge region (232..258) has no pause penalty (wi >= 0.68)
    assert min(w1[231:258]) >= 0.68
    assert min(w1[231:258]) > min(w3[231:258]) + 0.15

    # In 3DP, wi is lowered in static hinge 232..258 (wi <= 0.65)
    assert min(w3[231:258]) <= 0.65

    # In 4DP, wi is lowered in shifted window [264..268] (wi <= 0.45)
    assert min(w4[263:268]) <= 0.45


def test_html_viewer_standalone_and_offline_resilient():
    assert HTML_PATH.is_file(), f"HTML viewer missing at {HTML_PATH}"
    html = HTML_PATH.read_text(encoding="utf-8")

    # Must contain embedded data
    assert "const HPX_DATA = {" in html
    # Must NOT contain external script tags
    external_scripts = re.findall(r'<script\s+src=["\']http', html, re.IGNORECASE)
    assert not external_scripts, f"Found external scripts in HTML: {external_scripts}"

    # Must contain EMBL-EBI disclaimer and 4DP hypothesis statement
    assert "EMBL-EBI Methodological Caution" in html
    assert "Design Assumption Scope" in html
    assert "UNVERIFIED_PUBLIC_CORPUS" in html
    assert "Biological validation: NOT_AVAILABLE" in html
