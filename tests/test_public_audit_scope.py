"""Root filename allowlists must not accidentally scan excluded manuscript trees."""

import importlib.util
from pathlib import Path

path = Path(__file__).resolve().parents[1] / "scripts/audit_public_surface.py"
spec = importlib.util.spec_from_file_location("surface_audit", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_root_name_allowlist_does_not_include_nested_readmes(tmp_path):
    for relative in ["README.md", "paper/README.md", "web/README.md", "docs/usage.md"]:
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("Synthetic scope fixture", encoding="utf-8")
    actual = module.iter_public_files(
        tmp_path, module.FACTORFORGE_INCLUDE_DIRS, module.FACTORFORGE_INCLUDE_ROOT_FILES
    )
    assert {item.relative_to(tmp_path).as_posix() for item in actual} == {
        "README.md", "web/README.md", "docs/usage.md"
    }
