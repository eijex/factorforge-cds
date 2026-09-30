import builtins
import importlib
import subprocess
import sys
from pathlib import Path

import pytest

from factorforge import database


def test_persistence_is_disabled_without_explicit_configuration():
    assert database.persistence_status({}) == {
        "enabled": False,
        "backend": None,
        "database_url": None,
    }
    with pytest.raises(database.PersistenceNotConfiguredError, match="persistence is disabled"):
        database.configured_database_url({})


def test_non_postgresql_url_is_rejected():
    with pytest.raises(database.PersistenceNotConfiguredError, match="requires.*PostgreSQL"):
        database.configured_database_url({"FACTORFORGE_DATABASE_URL": "sqlite:///local.db"})


def test_explicit_session_factory_url_is_validated():
    with pytest.raises(database.PersistenceNotConfiguredError):
        database.create_session_factory("sqlite:///local.db")


def test_postgresql_urls_pin_the_declared_driver():
    for scheme in ("postgres", "postgresql", "postgresql+psycopg2"):
        assert database.configured_database_url({
            "FACTORFORGE_DATABASE_URL": f"{scheme}://user:password@localhost/db"
        }) == "postgresql+psycopg2://user:password@localhost/db"


def test_unconfigured_write_fails_before_optional_dependencies(monkeypatch):
    monkeypatch.delenv("DATABASE_URL", raising=False)
    monkeypatch.delenv("FACTORFORGE_DATABASE_URL", raising=False)
    with pytest.raises(database.PersistenceNotConfiguredError):
        database.save_optimization("test", "test", "MK", "ATGAAA", {})


def test_url_query_credentials_are_not_returned():
    safe = database.redact_database_url(
        "postgresql://alice:secret@localhost/db?password=query-secret#fragment-secret"
    )
    assert safe == "postgresql://alice@localhost/db"


def test_status_redacts_password():
    status = database.persistence_status(
        {"FACTORFORGE_DATABASE_URL": "postgresql://alice:secret@db.example:5432/eijex"}
    )
    assert status == {
        "enabled": True,
        "backend": "postgresql",
        "database_url": "postgresql+psycopg2://alice@db.example:5432/eijex",
    }
    assert "secret" not in repr(status)


def test_database_module_import_does_not_load_optional_dependencies(monkeypatch):
    original_import = builtins.__import__

    def blocked_import(name, *args, **kwargs):
        if name == "sqlalchemy" or name.startswith("sqlalchemy."):
            raise ModuleNotFoundError("blocked", name="sqlalchemy")
        if name == "eijex_db_core" or name.startswith("eijex_db_core."):
            raise ModuleNotFoundError("blocked", name="eijex_db_core")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked_import)
    importlib.reload(database)
    assert database.persistence_status({})["enabled"] is False


def test_clean_standalone_import_with_postgres_dependencies_blocked():
    script = """
import builtins
import sys
original_import = builtins.__import__
def blocked(name, *args, **kwargs):
    if name in {'psycopg2', 'sqlalchemy', 'eijex_db_core'} or name.startswith(('psycopg2.', 'sqlalchemy.', 'eijex_db_core.')):
        raise ModuleNotFoundError('optional persistence dependency blocked', name=name)
    return original_import(name, *args, **kwargs)
builtins.__import__ = blocked
sys.path.insert(0, 'src')
import factorforge
from factorforge import database
assert factorforge.__version__
assert database.persistence_status({})['enabled'] is False
"""
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr


def test_artifact_store_verifies_actual_content(tmp_path):
    store = database.LocalArtifactStore(tmp_path)
    uri = store.put(b"ATGAAA")
    digest = uri.rsplit("/", 1)[1]
    assert store.get(digest) == b"ATGAAA"
    assert store.put(b"ATGAAA") == uri
    (tmp_path / f"{digest}.seq").write_bytes(b"CORRUPTED")
    with pytest.raises(RuntimeError, match="checksum"):
        store.get(digest)


def test_artifact_store_rejects_path_traversal(tmp_path):
    with pytest.raises(ValueError, match="digest"):
        database.LocalArtifactStore(tmp_path).get("../../private")
