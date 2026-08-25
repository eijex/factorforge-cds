"""FactorForge Database Connector Module (v0.3 Postgres PostgreSQL Track).

Implements a Context-Aware, Policy-Driven Rule Engine connector that connects to the
canonical v0.3 Tri-Tier architecture:
1. common: Canonical Sequence Registry & Biological Context
2. factorforge: Rule Engine & Execution Provenance
3. validationhub: Physical Evidence & Claims
4. agentops: Policy Enforcement
5. ml: Research Scaffolds
"""

import hashlib
import json
import os
import uuid
from contextlib import contextmanager
from typing import Any, Dict, List, Mapping, Optional

import psycopg2
from psycopg2.extras import DictCursor


class FactorForgeDBConnector:
    """Python DB connector for FactorForge Context-Aware Rule Engine."""

    def __init__(self, dsn: Optional[str] = None):
        """Initialize the PostgreSQL connector.
        
        Reads standard Postgres environment variables if dsn is not provided.
        """
        self.dsn = dsn or os.environ.get(
            "FACTORFORGE_DB_DSN", 
            "host=127.0.0.1 port=5432 user=postgres password=postgres dbname=postgres"
        )
        self._init_db_extensions()

    @contextmanager
    def get_connection(self):
        """Context manager for DB connections."""
        conn = psycopg2.connect(self.dsn)
        try:
            yield conn
        finally:
            conn.close()

    def _init_db_extensions(self):
        """Ensure necessary PostgreSQL extensions/structures exist if needed."""
        # Our DDL (schema.sql) is now the authority. We do not CREATE TABLE here.
        # This connector merely operates on the existing v0.3 schema.
        pass

    @staticmethod
    def _manifest_hash(manifest: Mapping[str, Any]) -> str:
        """Return the canonical SHA-256 hash for a manifest payload."""
        payload = dict(manifest)
        payload.pop("manifest_hash", None)
        canonical = json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        )
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    def register_dataset_snapshot(self, manifest: Mapping[str, Any]) -> str:
        """Register verified manifest metadata in the ml.dataset_snapshots research checkpoint."""
        required = {"snapshot_name", "version", "status", "manifest_hash"}
        missing = sorted(required.difference(manifest))
        if missing:
            raise ValueError(f"Dataset snapshot manifest missing fields: {missing}")

        expected_hash = self._manifest_hash(manifest)
        supplied_hash = str(manifest["manifest_hash"])
        if supplied_hash != expected_hash:
            raise ValueError(
                "Dataset snapshot manifest hash mismatch: "
                f"expected {expected_hash}, received {supplied_hash}"
            )

        snapshot_id = str(uuid.uuid4())
        
        with self.get_connection() as conn:
            with conn.cursor() as cursor:
                # Check for existing
                cursor.execute(
                    """
                    SELECT snapshot_id
                    FROM ml.dataset_snapshots
                    WHERE snapshot_name = %s AND version = %s AND manifest_hash = %s
                    """,
                    (manifest["snapshot_name"], manifest["version"], supplied_hash),
                )
                row = cursor.fetchone()
                if row:
                    return str(row[0])

                cursor.execute(
                    """
                    INSERT INTO ml.dataset_snapshots
                    (snapshot_id, snapshot_name, version, status, manifest_hash)
                    VALUES (%s, %s, %s, %s, %s)
                    """,
                    (
                        snapshot_id,
                        manifest["snapshot_name"],
                        manifest["version"],
                        manifest["status"],
                        supplied_hash,
                    ),
                )
            conn.commit()
            return snapshot_id

    def register_sequence(self, raw_sequence: str, molecule_class: str = "CDS") -> str:
        """Registers a canonical sequence in common.sequences.
        
        Note: Requires an underlying common.artifacts entry as per v0.3 spec invariant.
        """
        seq_hash = hashlib.sha256(raw_sequence.encode("utf-8")).hexdigest()
        
        with self.get_connection() as conn:
            with conn.cursor(cursor_factory=DictCursor) as cursor:
                cursor.execute(
                    "SELECT sequence_id FROM common.sequences WHERE canonical_sequence_sha256 = %s",
                    (seq_hash,),
                )
                row = cursor.fetchone()
                if row:
                    return str(row["sequence_id"])

                # Insert dummy artifact to satisfy FK invariant (in real app, this points to actual fasta file in GCS/S3)
                artifact_id = str(uuid.uuid4())
                cursor.execute(
                    """
                    INSERT INTO common.artifacts (artifact_id, artifact_uri, sha256_hash, artifact_type)
                    VALUES (%s, %s, %s, 'fasta')
                    """,
                    (artifact_id, f"internal://generated/{seq_hash}.fasta", seq_hash)
                )

                seq_id = str(uuid.uuid4())
                cursor.execute(
                    """
                    INSERT INTO common.sequences (sequence_id, molecule_class, canonical_sequence_sha256, artifact_id, length)
                    VALUES (%s, %s, %s, %s, %s)
                    """,
                    (seq_id, molecule_class, seq_hash, artifact_id, len(raw_sequence)),
                )
            conn.commit()
            return seq_id

    def save_computational_provenance(
        self,
        run_metadata: Dict[str, Any],
        candidate_data: Dict[str, Any],
        check_results: List[Dict[str, Any]],
    ) -> str:
        """
        Atomically saves the execution provenance (Phase 1):
        1. design_run
        2. candidate
        3. check_results
        
        Rolls back the entire transaction if any step fails. Does NOT create a design_package.
        """
        raw_seq = candidate_data["optimized_sequence"]
        # Sequence registration could technically be part of the transaction, but it's idempotent.
        seq_id = self.register_sequence(raw_seq)
        
        run_id = str(uuid.uuid4())
        candidate_id = str(uuid.uuid4())

        # The context manager automatically commits on success, and rolls back on exception
        with self.get_connection() as conn:
            with conn.cursor() as cursor:
                # 1. Insert Design Run
                cursor.execute(
                    """
                    INSERT INTO factorforge.design_runs 
                    (design_run_id, execution_origin, actual_engine_name, actual_profile_name, 
                     generation_performed, analysis_mode, evidence_level, runner_entrypoint)
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                    """,
                    (
                        run_id,
                        run_metadata.get("execution_origin", "standalone_cli"),
                        run_metadata["actual_engine_name"],
                        run_metadata.get("actual_profile_name"),
                        run_metadata.get("generation_performed", True),
                        run_metadata.get("analysis_mode", "cli_generation"),
                        run_metadata.get("evidence_level", "prospective_factorforge"),
                        run_metadata.get("runner_entrypoint", "factorforge.cli.main"),
                    )
                )

                # 2. Insert Candidate (Metrics)
                overall_status = candidate_data.get("computational_status", "passed")
                cursor.execute(
                    """
                    INSERT INTO factorforge.candidates
                    (candidate_id, design_run_id, sequence_id, computational_status, cai_value, gc_percent)
                    VALUES (%s, %s, %s, %s, %s, %s)
                    """,
                    (
                        candidate_id,
                        run_id,
                        seq_id,
                        overall_status,
                        candidate_data.get("cai", None),
                        candidate_data.get("gc_percent", None),
                    ),
                )

                # 3. Insert Check Results (Constraints)
                for chk in check_results:
                    check_result_id = str(uuid.uuid4())
                    cursor.execute(
                        """
                        INSERT INTO factorforge.check_results
                        (check_result_id, candidate_id, check_type, result, observed_value, details_json)
                        VALUES (%s, %s, %s, %s, %s, %s)
                        """,
                        (
                            check_result_id,
                            candidate_id,
                            chk["check_domain"],  # e.g., sequence_integrity, advisory_sequence_risk
                            chk["result"],        # PASS, WARN, FAIL
                            chk.get("observed_value"),
                            json.dumps(chk.get("details_json", {})),
                        ),
                    )
            # Automatic commit here; if exception occurs, it rolls back.
            conn.commit()

        return run_id
