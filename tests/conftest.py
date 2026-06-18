import sys
import importlib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _prefer_local_factorforge() -> None:
    for name in list(sys.modules):
        if name != "factorforge" and not name.startswith("factorforge."):
            continue
        module = sys.modules.get(name)
        module_file = getattr(module, "__file__", None)
        if module_file is None:
            sys.modules.pop(name, None)
            continue
        loaded_path = Path(module_file).resolve()
        if SRC not in loaded_path.parents:
            sys.modules.pop(name, None)


_prefer_local_factorforge()
local_factorforge = importlib.import_module("factorforge")
local_factorforge.__path__ = [str(SRC / "factorforge")]
_prefer_local_factorforge()


def pytest_configure(config):
    _prefer_local_factorforge()


@pytest.fixture
def design_package() -> dict:
    """Valid public Open Bio Design Package used by contract tests."""
    return {
        "design_package_version": "1.0",
        "design_id": "CF-OPENBIO-001",
        "sequence_type": "cds",
        "input_type": "protein",
        "host_profile": {
            "id": "nbenthamiana",
            "display_name": "N. benthamiana",
            "scientific_name": "Nicotiana benthamiana",
            "ncbi_taxonomy_id": 4100,
            "status": "stable",
        },
        "optimization": {"engine": "factorforge_cds", "profile": "balanced"},
        "metrics": {
            "cai": 0.91,
            "gc_percent": 60.0,
            "mfe_kcal_mol": None,
            "mfe_status": "not_computed",
            "mfe_used": False,
        },
        "validation": {
            "biological_pass": True,
            "assembly_pass": True,
            "aa_identity": 1.0,
            "forbidden_type_iis_site_count": 0,
        },
        "evidence": {
            "sequence_hash": f"sha256:{'a' * 64}",
            "registry_version": "3.2.2",
            "registry_hash": f"sha256:{'b' * 64}",
        },
        "claim_boundary": {
            "in_silico_only": True,
            "no_yield_claim": True,
            "no_wet_lab_claim": True,
            "no_clinical_claim": True,
        },
    }
