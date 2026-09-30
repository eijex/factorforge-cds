import os
import uuid
import hashlib

import pytest

if not os.environ.get("DATABASE_URL", "").strip():
    pytest.skip("DATABASE_URL is not configured", allow_module_level=True)

pytest.importorskip("sqlalchemy")
pytest.importorskip("psycopg2")

from sqlalchemy import create_engine, text
from sqlalchemy.orm import sessionmaker

from factorforge.database import get_batch, save_optimization


@pytest.fixture(autouse=True)
def private_artifact_store(tmp_path, monkeypatch):
    monkeypatch.setenv("FACTORFORGE_ARTIFACT_DIR", str(tmp_path))


DATABASE_URL = os.environ["DATABASE_URL"]
engine = create_engine(DATABASE_URL, pool_pre_ping=True)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


def _db_available() -> bool:
    try:
        with engine.connect() as conn:
            conn.execute(text("SELECT 1"))
        return True
    except Exception:
        return False


@pytest.mark.skipif(not _db_available(), reason="Database not available")
def test_save_and_retrieve():
    study_number = f"TEST-{uuid.uuid4().hex[:8]}"

    batch_id = save_optimization(
        study_number=study_number,
        protein_name="Test Protein",
        input_sequence="MKLLVV",
        optimized_sequence="ATGAAACTGCTGGTGGTG",
        metrics={
            "gc_content": 0.389,
            "cai": 0.82,
            "execution_time": 0.5,
        },
        session_factory=SessionLocal,
    )

    assert batch_id

    batch = get_batch(study_number, session_factory=SessionLocal)
    assert batch is not None
    assert batch["protein"] == "Test Protein"
    assert len(batch["sequences"]) == 2
    assert [row["type"] for row in batch["sequences"]] == ["input", "optimized"]
    assert batch["sequences"][0]["sha256"] == hashlib.sha256(b"MKLLVV").hexdigest()
    assert batch["sequences"][1]["sha256"] == hashlib.sha256(
        b"ATGAAACTGCTGGTGGTG"
    ).hexdigest()
    assert all(not row["data"].startswith("Sequence<") for row in batch["sequences"])
    assert batch["sequences"][0]["data"] == "MKLLVV..."
    assert batch["sequences"][1]["data"] == "ATGAAACTGCTGGTGGTG..."

    repeated_batch_id = save_optimization(
        study_number=study_number,
        protein_name="Test Protein",
        input_sequence="MKLLVV",
        optimized_sequence="ATGAAACTGCTGGTGGTG",
        metrics={"gc_content": 0.389, "cai": 0.82, "execution_time": 0.5},
        session_factory=SessionLocal,
    )
    assert repeated_batch_id == batch_id
    with pytest.raises(ValueError, match="provenance"):
        save_optimization(study_number, "Test Protein", "MKLLVV",
                          "ATGAAACTGCTGGTGGTG", {"cai": 0.9},
                          session_factory=SessionLocal)


@pytest.mark.skipif(not _db_available(), reason="Database not available")
def test_failed_metric_write_rolls_back_campaign_and_sequences():
    study_number = f"ROLLBACK-{uuid.uuid4().hex}"
    with pytest.raises(ValueError):
        save_optimization(
            study_number, "Rollback fixture", "MK", "ATGAAA",
            {"cai": "invalid-number"}, session_factory=SessionLocal,
        )
    assert get_batch(study_number, session_factory=SessionLocal) is None


@pytest.mark.skipif(not _db_available(), reason="Database not available")
def test_missing_campaign_returns_none():
    assert get_batch(f"MISSING-{uuid.uuid4().hex}", session_factory=SessionLocal) is None


def test_concurrent_retry_has_one_campaign():
    from concurrent.futures import ThreadPoolExecutor
    study = f"CONCURRENT-{uuid.uuid4().hex}"

    def save():
        return save_optimization(study, "Fixture", "MK", "ATGAAA", {"cai": 0.8},
                                 session_factory=SessionLocal)

    with ThreadPoolExecutor(max_workers=4) as pool:
        ids = list(pool.map(lambda _: save(), range(4)))
    assert len(set(ids)) == 1
    assert get_batch(study, session_factory=SessionLocal)["sequences"][0]["data"] == "MK..."


def test_legacy_campaign_without_sequence_roles_is_not_completed():
    from eijex_db_core.models import Campaign
    study = f"LEGACY-{uuid.uuid4().hex}"
    with SessionLocal() as session, session.begin():
        session.add(Campaign(name=study, description="Legacy fixture"))
    with pytest.raises(RuntimeError, match="legacy-record migration"):
        get_batch(study, session_factory=SessionLocal)
