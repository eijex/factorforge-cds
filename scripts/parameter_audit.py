"""Scan factorforge codebase for hardcoded biological/benchmark numbers.
Discovery scan mode: identifies candidates for registry review, does not require registry to run.
Outputs docs/validation/parameter_scan_report.md.
Run: python scripts/parameter_audit.py
Exit 0 = clean or only allowlisted; Exit 1 = unregistered public-biological numbers found."""
from __future__ import annotations
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs" / "validation" / "parameter_scan_report.md"

# Numbers that are registry-required (public-biological)
REGISTRY_TARGETS = {
    "cai": [r"\b0\.82\b"],
    "gc_global": [r"\b55\.0\b", r"\b65\.0\b"],
    "gc_local": [r"\b25\.0\b", r"\b75\.0\b"],
    "aa_identity": [r"\b1\.0\b"],  # context-dependent — see allowlist
    "type_iis": [r"\bBsaI\b", r"\bBsmBI\b", r"\bBbsI\b"],
    "sgn_n": [r"\b3876\b"],
}

# Numbers allowed without registry entry (formatting, test tolerances, etc.)
ALLOWLIST_PATTERNS = [
    r"FASTA.*60",       # line wrap
    r"round\(.*,\s*[12]\)",  # display rounding
    r"pytest",          # test tolerance
    r"timeout",         # timeout values
    r"#.*allowlist",    # explicit inline allowlist comment
]

SCAN_DIRS = ["src/factorforge", "benchmarks", "docs", "tests"]
SKIP_FILES = {"parameter_audit.py", "current_parameter_registry.yaml", "parameter_scan_report.md", "conftest.py"}
SCAN_EXTENSIONS = {".py", ".md", ".yaml", ".yml"}


def _scan() -> list[dict]:
    hits = []
    for d in SCAN_DIRS:
        for f in (ROOT / d).rglob("*"):
            if f.suffix not in SCAN_EXTENSIONS or f.name in SKIP_FILES:
                continue
            text = f.read_text(encoding="utf-8", errors="replace")
            for category, patterns in REGISTRY_TARGETS.items():
                for pat in patterns:
                    for m in re.finditer(pat, text):
                        line_no = text[:m.start()].count("\n") + 1
                        line = text.splitlines()[line_no - 1].strip()
                        if any(re.search(a, line) for a in ALLOWLIST_PATTERNS):
                            status = "allowlisted"
                        else:
                            status = "UNREGISTERED"
                        hits.append({
                            "file": str(f.relative_to(ROOT)),
                            "line": line_no,
                            "category": category,
                            "match": m.group(),
                            "context": line[:120],
                            "status": status,
                        })
    return hits


def _write_report(hits: list[dict]) -> int:
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    unregistered = [h for h in hits if h["status"] == "UNREGISTERED"]
    lines = [
        "# Parameter Scan Report",
        "",
        f"**Total hits:** {len(hits)}  **Unregistered:** {len(unregistered)}",
        "",
        "## Unregistered Public-Biological Numbers",
        "",
        "| file | line | category | match | context |",
        "|------|------|----------|-------|---------|",
    ]
    for h in unregistered:
        lines.append(f"| `{h['file']}` | {h['line']} | {h['category']} | `{h['match']}` | `{h['context'][:80]}` |")
    lines += ["", "## Allowlisted (not requiring registry entry)", "",
              "| file | line | match |", "|------|------|-------|"]
    for h in hits:
        if h["status"] == "allowlisted":
            lines.append(f"| `{h['file']}` | {h['line']} | `{h['match']}` |")
    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print(f"Report written: {REPORT}")
    return len(unregistered)


if __name__ == "__main__":
    hits = _scan()
    n_bad = _write_report(hits)
    if n_bad > 0:
        print(f"WARNING: {n_bad} unregistered public-biological numbers found. See {REPORT}")
        sys.exit(1)
    print("Clean: all detected numbers are registry-registered or allowlisted.")
