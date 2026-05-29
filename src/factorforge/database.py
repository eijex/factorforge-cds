"""Database models and CRUD operations for FactorForge."""

from __future__ import annotations

import os
import uuid
from typing import Dict, Optional

from sqlalchemy import ARRAY, DECIMAL, TIMESTAMP, Column, ForeignKey, String, Text, create_engine, func
from sqlalchemy.dialects.postgresql import JSONB, UUID
from sqlalchemy.orm import DeclarativeBase, sessionmaker

DATABASE_URL = os.getenv("DATABASE_URL", "")

engine = create_engine(DATABASE_URL, pool_pre_ping=True)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


class Base(DeclarativeBase):
    pass


class Batch(Base):
    __tablename__ = "batches"
    __table_args__ = {"schema": "factorforge"}

    batch_id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    study_number = Column(String(50), unique=True, nullable=False)
    organism = Column(String(100), nullable=False)
    target_protein = Column(String(255), nullable=False)
    created_at = Column(TIMESTAMP, server_default=func.now())
    updated_at = Column(TIMESTAMP, server_default=func.now(), onupdate=func.now())
    status = Column(String(20), default="pending")
    created_by = Column(String(100))


class Sequence(Base):
    __tablename__ = "sequences"
    __table_args__ = {"schema": "factorforge"}

    sequence_id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    batch_id = Column(
        UUID(as_uuid=True),
        ForeignKey("factorforge.batches.batch_id", ondelete="CASCADE"),
        nullable=False,
    )
    sequence_type = Column(String(20), nullable=False)
    sequence_data = Column(Text, nullable=False)
    gc_content = Column(DECIMAL(5, 4))
    cai = Column(DECIMAL(5, 4))
    tm = Column(DECIMAL(5, 2))
    created_at = Column(TIMESTAMP, server_default=func.now())
    metadata_ = Column("metadata", JSONB)


class OptimizationResult(Base):
    __tablename__ = "optimization_results"
    __table_args__ = {"schema": "factorforge"}

    result_id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    batch_id = Column(UUID(as_uuid=True), ForeignKey("factorforge.batches.batch_id"))
    sequence_id = Column(UUID(as_uuid=True), ForeignKey("factorforge.sequences.sequence_id"))
    algorithm_version = Column(String(20), nullable=False)
    execution_time_sec = Column(DECIMAL(8, 3))
    avoided_motifs = Column(ARRAY(Text))
    warnings = Column(ARRAY(Text))
    created_at = Column(TIMESTAMP, server_default=func.now())


def save_optimization(
    study_number: str,
    protein_name: str,
    input_sequence: str,
    optimized_sequence: str,
    metrics: Dict,
    algorithm_version: str = "2.1.0",
) -> str:
    """Save optimization result to database."""
    with SessionLocal() as session:
        batch = Batch(
            study_number=study_number,
            organism="nicotiana_benthamiana",
            target_protein=protein_name,
            status="completed",
        )
        session.add(batch)
        session.flush()

        input_seq = Sequence(
            batch_id=batch.batch_id,
            sequence_type="input",
            sequence_data=input_sequence,
        )
        output_seq = Sequence(
            batch_id=batch.batch_id,
            sequence_type="optimized",
            sequence_data=optimized_sequence,
            gc_content=metrics.get("gc_content"),
            cai=metrics.get("cai"),
            tm=metrics.get("tm"),
            metadata_=metrics,
        )
        session.add_all([input_seq, output_seq])
        session.flush()

        result = OptimizationResult(
            batch_id=batch.batch_id,
            sequence_id=output_seq.sequence_id,
            algorithm_version=algorithm_version,
            execution_time_sec=metrics.get("execution_time"),
            avoided_motifs=metrics.get("avoided_motifs", []),
            warnings=metrics.get("warnings", []),
        )
        session.add(result)
        session.commit()

        return str(batch.batch_id)


def get_batch(study_number: str) -> Optional[Dict]:
    """Retrieve batch by study number."""
    with SessionLocal() as session:
        batch = (
            session.query(Batch)
            .filter(Batch.study_number == study_number)
            .first()
        )
        if not batch:
            return None

        sequences = session.query(Sequence).filter(Sequence.batch_id == batch.batch_id).all()

        return {
            "batch_id": str(batch.batch_id),
            "study_number": batch.study_number,
            "protein": batch.target_protein,
            "status": batch.status,
            "sequences": [
                {
                    "type": seq.sequence_type,
                    "data": f"{seq.sequence_data[:50]}...",
                    "gc": float(seq.gc_content) if seq.gc_content is not None else None,
                    "cai": float(seq.cai) if seq.cai is not None else None,
                }
                for seq in sequences
            ],
        }
