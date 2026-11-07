"""Public release summaries must not advertise internal or unsupported results."""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def modal():
    html = (ROOT / "web/index.html").read_text(encoding="utf-8")
    return html.split('id="changelogModal"', 1)[1].split("<!-- Third-Party Linkout", 1)[0]


@pytest.mark.parametrize(
    "phrase",
    [
        "Internal Milestone",
        "Stage 009",
        "Training Telemetry",
        "_machine-learning",
        "_improvement-jobs",
        "leak-free",
        "67/67",
        "100% bit-identical",
        "AlphaFold-grade",
        "resolving structural exit-tunnel latency",
        "Submitted sequences are not logged or stored",
        "Synchronized structural domain emergence",
    ],
)
def test_modal_excludes_internal_records_and_unbound_guarantees(phrase):
    assert phrase.lower() not in modal().lower()


def test_current_release_and_historical_boundary():
    text = modal()
    assert re.search(r"<h3[^>]*>v3\.5\.12</h3>", text)
    assert len(re.findall(r">Current</span>", text)) == 1
    assert "Historical entries do not override" in text
    assert "General 3DP/4DP web access remains paused" in text
    assert "not a guarantee of infrastructure-wide zero retention" in text
    assert "No infrastructure retention" not in text
    assert "https://github.com/eijex/factorforge-cds/blob/main/CHANGELOG.md" in text


@pytest.mark.parametrize("name", ["CHANGELOG.md", "docs/changelog.md"])
def test_linked_history_withdraws_zero_retention_promise(name):
    text = (ROOT / name).read_text(encoding="utf-8")
    assert "Submitted sequences are not logged or stored" not in text
    assert "Historical privacy wording corrected 2026-11-07" in text


def test_historical_research_notes_do_not_claim_calibrated_kinetics():
    text = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    section = text.split("## [3.5.11]", 1)[1].split("## 3.5.7", 1)[0]
    assert "not calibrated translation speeds" in section
    assert "leak-free" not in section
    assert "Stage009" not in section and "Stage 009" not in section
    assert "resolving structural exit-tunnel latency" not in section
