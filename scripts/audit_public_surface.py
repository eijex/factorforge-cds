"""Audit FactorForge public release surfaces for leaked internal info and
unsupported claims. Run with no args from the repo root for a local scan
(used by CI); pass --workspace/--live/--external for the full multi-repo,
multi-surface maintainer audit.
"""

from __future__ import annotations

import argparse
import re
import sys
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from html.parser import HTMLParser
from pathlib import Path


ROOT = Path.cwd()
DEFAULT_LIVE_BASE = "https://eijex.github.io/factorforge-cds/"
DEFAULT_EXTERNAL_URLS = [
    "https://www.eijex.com/",
    "https://factorforge.eijex.com/",
    "https://mcp.eijex.com/",
    "https://pypi.org/project/factorforge-cds/",
    "https://zenodo.org/records/20407330",
    "https://github.com/eijex",
    "https://github.com/eijex/factorforge-cds",
    "https://raw.githubusercontent.com/eijex/factorforge-cds/main/.github/ISSUE_TEMPLATE/wet_lab_result.yml",
]

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

INCLUDE_SUFFIXES = {".md", ".cff", ".toml", ".yml", ".yaml", ".html", ".js", ".ts", ".tsx", ".json"}
FACTORFORGE_INCLUDE_DIRS = {
    ".github",
    "benchmarks",
    "docs",
    "recipes",
    "web",
}
FACTORFORGE_INCLUDE_ROOT_FILES = {
    "CHANGELOG.md",
    "CITATION.cff",
    "CONTRIBUTING.md",
    "README.md",
    "ROADMAP.md",
    "SECURITY.md",
    "VALIDATION.md",
    "mkdocs.yml",
    "pyproject.toml",
}
WEB_INCLUDE_DIRS = {"docs", "src"}
WEB_INCLUDE_ROOT_FILES = {"README.md", "CHANGELOG.md", "SECURITY.md", "package.json"}
MCP_INCLUDE_DIRS = {"skills", "src"}
MCP_INCLUDE_ROOT_FILES = {"README.md", "CHANGELOG.md", "SECURITY.md", "package.json"}
PUBLIC_WORKSPACE_REPOS = {
    "factorforge": (FACTORFORGE_INCLUDE_DIRS, FACTORFORGE_INCLUDE_ROOT_FILES),
    "eijex-web": (WEB_INCLUDE_DIRS, WEB_INCLUDE_ROOT_FILES),
    "eijex-mcp": (MCP_INCLUDE_DIRS, MCP_INCLUDE_ROOT_FILES),
}
EXCLUDE_PARTS = {
    ".git",
    ".pytest_cache",
    ".ruff_cache",
    "archive",
    "build",
    "dist",
    "node_modules",
    "package-lock.json",
    "pnpm-lock.yaml",
    "test-results",
    "RELEASE_GATE.md",
}

PATTERNS = {
    "unsupported_claim": re.compile(
        r"increases yield|guarantees expression|validated expression optimizer|"
        r"wet-lab proven|regulatory-ready|ready for synthesis|"
        r"superior|outperform|promote efficient|delivery optimized|"
        r"plant expression workflows|expression workflows|"
        r"predicts expression|predicts yield|wet-lab success|"
        r"guarantees cloning|predicts synthesis success|guarantees synthesis",
        re.IGNORECASE,
    ),
    "submission_meta": re.compile(
        r"\bJOSS\b|Bioinformatics|Application Note|manuscript|peer review",
        re.IGNORECASE,
    ),
    "internal_reference": re.compile(
        r"\bClaude\b|\bCodex\b|\bAnthropic\b|\bOpenAI\b|CLAUDE\.md|"
        r"C:\\Work|PlantFormOrg|PlantForm-AI|"
        r"eijex-web\.vercel\.app|"
        r"\bJob\s*\d+[A-Za-z-]*\b",
        re.IGNORECASE,
    ),
    "ml_public_exposure": re.compile(
        r"SynCodonLM|language model|HuggingFace|BERT|ml_enhanced|v3-alpha",
        re.IGNORECASE,
    ),
    "sequence_hash_publication": re.compile(
        r"sequence hash|sequence hashes|per-sequence hash",
        re.IGNORECASE,
    ),
    "dangerous_feedback_field": re.compile(
        r"^\s*\d+\.?\s*(Construct ID|Protein name|Promoter|Institution)\b|"
        r"FactorForge\s+better|Detected\s+by\s+Western\s+blot|Detected\s+by\s+ELISA|"
        r"Functional\s+activity\s+confirmed|Relative\s+to\s+native\s+control",
        re.IGNORECASE,
    ),
}

ALLOW = {
    "internal stop",
    "internal stop codon",
    "internal stop codons",
    "internal_stop",
    "non-public",
    "private reporting",
    "private contact",
    "private or sensitive",
    "confidential",
    "non-confidential",
    "unpublished constructs",
    "bioinformatics",
    "bioconda",
    "agents.md",
    "do not include",
    "contains no",
    "removed",
    "coarse",
    "non-identifying",
    "v3-alpha",
    "v3-ml-prototype",
    "not by gene id",
    "never raw input or output sequence",
}


def iter_public_files(root: Path, include_dirs: set[str], include_root_files: set[str]) -> list[Path]:
    files: list[Path] = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(root)
        if any(part in EXCLUDE_PARTS for part in rel.parts):
            continue
        if path.suffix.lower() not in INCLUDE_SUFFIXES:
            continue
        if rel.name in include_root_files or rel.parts[0] in include_dirs:
            files.append(path)
    return sorted(files)


def local_scan_roots(args: argparse.Namespace) -> list[tuple[str, Path, set[str], set[str]]]:
    if args.workspace:
        workspace = Path(args.workspace).resolve()
        scans: list[tuple[str, Path, set[str], set[str]]] = []
        for repo_name, (include_dirs, include_root_files) in PUBLIC_WORKSPACE_REPOS.items():
            repo_root = workspace / repo_name
            if repo_root.exists():
                scans.append((repo_name, repo_root, include_dirs, include_root_files))
        return scans

    root = ROOT.resolve()
    if (root / "src" / "factorforge").exists():
        return [("factorforge", root, FACTORFORGE_INCLUDE_DIRS, FACTORFORGE_INCLUDE_ROOT_FILES)]
    if (root / "src" / "app").exists() and (root / "src" / "app" / "api" / "mcp").exists():
        return [("eijex-mcp", root, MCP_INCLUDE_DIRS, MCP_INCLUDE_ROOT_FILES)]
    if (root / "src" / "app").exists():
        return [("eijex-web", root, WEB_INCLUDE_DIRS, WEB_INCLUDE_ROOT_FILES)]
    return []


def allowed(line: str) -> bool:
    lower = line.lower()
    return any(term in lower for term in ALLOW)


class TextExtractor(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.parts: list[str] = []
        self.skip = False

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag in {"script", "style", "noscript"}:
            self.skip = True

    def handle_endtag(self, tag: str) -> None:
        if tag in {"script", "style", "noscript"}:
            self.skip = False

    def handle_data(self, data: str) -> None:
        if not self.skip:
            text = " ".join(data.split())
            if text:
                self.parts.append(text)


def fetch_text(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "factorforge-public-surface-audit/1.0"})
    with urllib.request.urlopen(request, timeout=20) as response:
        raw = response.read()
    content = raw.decode("utf-8", errors="replace")
    if "<html" in content[:500].lower():
        parser = TextExtractor()
        parser.feed(content)
        return "\n".join(parser.parts)
    return content


def discover_live_urls(base_url: str) -> list[str]:
    base_url = base_url.rstrip("/") + "/"
    sitemap_url = base_url + "sitemap.xml"
    try:
        sitemap = fetch_text(sitemap_url)
        root = ET.fromstring(sitemap)
        urls = [
            elem.text.strip()
            for elem in root.iter()
            if elem.tag.endswith("loc") and elem.text and elem.text.strip().startswith(base_url)
        ]
        return sorted(set(urls)) or [base_url]
    except (ET.ParseError, urllib.error.URLError, TimeoutError, OSError):
        return [base_url]


def scan_lines(source: str, text: str) -> list[tuple[str, str, int, str]]:
    findings: list[tuple[str, str, int, str]] = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        if allowed(line):
            continue
        for name, pattern in PATTERNS.items():
            if pattern.search(line):
                findings.append((name, source, lineno, line.strip()))
    return findings


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit FactorForge public release surfaces.")
    parser.add_argument(
        "--live",
        nargs="?",
        const=DEFAULT_LIVE_BASE,
        help="Also audit deployed docs. Defaults to the FactorForge GitHub Pages URL.",
    )
    parser.add_argument(
        "--external",
        action="store_true",
        help="Audit default external public surfaces: PyPI, Zenodo, GitHub org, and GitHub repo.",
    )
    parser.add_argument(
        "--url",
        action="append",
        default=[],
        help="Additional public URL to audit. May be repeated.",
    )
    parser.add_argument(
        "--workspace",
        help="Audit known public repos under this workspace: factorforge, eijex-web, and eijex-mcp.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    scan_roots = local_scan_roots(args)
    if not scan_roots:
        print("Run from factorforge, eijex-web, or eijex-mcp, or pass --workspace <path-to-eijex-monorepo>.")
        return 2

    findings: list[tuple[str, str, int, str]] = []
    for label, root, include_dirs, include_root_files in scan_roots:
        prefix = "" if len(scan_roots) == 1 else f"{label}/"
        for path in iter_public_files(root, include_dirs, include_root_files):
            text = path.read_text(encoding="utf-8", errors="replace")
            findings.extend(scan_lines(prefix + str(path.relative_to(root)), text))

    if args.live:
        for url in discover_live_urls(args.live):
            try:
                findings.extend(scan_lines(url, fetch_text(url)))
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                findings.append(("live_fetch_error", url, 0, str(exc)))

    external_urls: list[str] = []
    if args.external:
        external_urls.extend(DEFAULT_EXTERNAL_URLS)
    external_urls.extend(args.url)
    for url in sorted(set(external_urls)):
        try:
            findings.extend(scan_lines(url, fetch_text(url)))
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            findings.append(("external_fetch_error", url, 0, str(exc)))

    if not findings:
        print("No public-surface findings.")
        return 0

    for name, source, lineno, line in findings:
        location = f"{source}:{lineno}" if lineno else source
        print(f"{name}: {location}: {line}")
    print(f"\nTotal findings: {len(findings)}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
