"""The core SQLite path must not require optional PostgreSQL packages."""

import subprocess
import sys
from pathlib import Path


def test_sqlite_works_with_postgres_import_blocked(tmp_path):
    script = """
import builtins
import sys
original_import = builtins.__import__
def without_postgres(name, *args, **kwargs):
    if name == "psycopg2" or name.startswith("psycopg2."):
        raise ModuleNotFoundError("PostgreSQL dependency absent", name="psycopg2")
    return original_import(name, *args, **kwargs)
builtins.__import__ = without_postgres
sys.path.insert(0, "src")
from factorforge.db.connector import FactorForgeDBConnector
db = FactorForgeDBConnector(db_path=sys.argv[1])
with db._sqlite_connection() as conn:
    assert conn.execute("SELECT count(*) FROM factorforge_candidates").fetchone()[0] == 0
try:
    with FactorForgeDBConnector().get_connection():
        raise AssertionError("PostgreSQL connected without its driver")
except RuntimeError as exc:
    assert "factorforge-cds[postgres]" in str(exc)
"""
    result = subprocess.run(
        [sys.executable, "-c", script, str(tmp_path / "checkpoint.db")],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr
