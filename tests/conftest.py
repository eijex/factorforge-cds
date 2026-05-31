import sys
import importlib
from pathlib import Path

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
