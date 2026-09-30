import importlib.util
from pathlib import Path


def load_scanner():
    path = Path(__file__).resolve().parents[1] / "scripts" / "audit_public_surface.py"
    spec = importlib.util.spec_from_file_location("public_surface_scanner", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_explicit_negative_claim_is_not_affirmative():
    scanner = load_scanner()
    for text in (
        "Metrics do not establish expression, yield, or biological superiority.",
        "Evaluation boundary without implying biological superiority.",
        "not establish expression, yield, or biological superiority.",
    ):
        assert not scanner.scan_lines("fixture", text)


def test_second_affirmative_claim_on_same_line_is_not_hidden():
    text = "Metrics do not establish biological superiority. This guarantees expression."
    assert load_scanner().scan_lines("fixture", text)


def test_generated_results_are_excluded_but_curated_docs_are_scanned(tmp_path):
    results = tmp_path / "benchmarks" / "results"
    results.mkdir(parents=True)
    (results / "figure.html").write_text("frozen artifact", encoding="utf-8")
    docs = tmp_path / "docs"
    docs.mkdir()
    public_doc = docs / "guide.md"
    public_doc.write_text("curated documentation", encoding="utf-8")
    files = load_scanner().iter_public_files(tmp_path, {"docs", "benchmarks"}, set())
    assert files == [public_doc]
