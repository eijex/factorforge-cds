import uuid

import pytest

pytest.importorskip("sqlalchemy")
pytest.importorskip("psycopg2")

from sqlalchemy import text

from src.factorforge.database import engine, get_batch, save_optimization


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
    )

    assert batch_id

    batch = get_batch(study_number)
    assert batch is not None
    assert batch["protein"] == "Test Protein"
    assert len(batch["sequences"]) == 2
