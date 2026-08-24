"""FactorForge Database Connector Module (v3.5.0 ML Track - Upgraded Schema).

Implements a Context-Aware, Policy-Driven Rule Engine connector that stores:
1. Biological Contexts & Canonical Sequence Registries
2. Constraint Definitions (HARD/SOFT rules) & Versioned Profiles
3. Restriction Enzyme Registries
4. Execution Lineage: Requests -> Runs -> Candidates -> Evaluations -> Design Packages
"""

from typing import Dict, Any, List, Optional
import hashlib
import json
import sqlite3
from pathlib import Path


class FactorForgeDBConnector:
    """Python DB connector for FactorForge Context-Aware Rule Engine."""

    def __init__(self, db_path: Optional[str] = None):
        if db_path is None:
            db_dir = Path(__file__).resolve().parents[3] / "data" / "db"
            db_dir.mkdir(parents=True, exist_ok=True)
            self.db_path = str(db_dir / "factorforge_relational.db")
        else:
            self.db_path = db_path

        self._init_db()

    def _get_connection(self):
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def _init_db(self):
        """Initializes SQLite tables matching the upgraded PRD Job 233 architecture."""
        with self._get_connection() as conn:
            cursor = conn.cursor()

            # 1. Biological Context
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS common_biological_context (
                context_id INTEGER PRIMARY KEY AUTOINCREMENT,
                context_code TEXT UNIQUE NOT NULL,
                organism TEXT NOT NULL,
                taxonomy_id INTEGER NOT NULL,
                cellular_compartment TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
            """)

            # 2. Sequence Registry
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS common_sequence_registry (
                sequence_id INTEGER PRIMARY KEY AUTOINCREMENT,
                sequence_hash TEXT UNIQUE NOT NULL,
                molecule_type TEXT NOT NULL,
                sequence_length INTEGER NOT NULL,
                raw_sequence TEXT NOT NULL,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
            """)

            # 3. Restriction Enzyme Registry
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS factorforge_enzyme_registry (
                enzyme_id INTEGER PRIMARY KEY AUTOINCREMENT,
                enzyme_name TEXT UNIQUE NOT NULL,
                recognition_sequence TEXT NOT NULL,
                category TEXT DEFAULT 'Type_IIS'
            );
            """)

            # Seed default enzymes if empty
            cursor.execute("SELECT COUNT(*) FROM factorforge_enzyme_registry;")
            if cursor.fetchone()[0] == 0:
                cursor.executemany(
                    "INSERT INTO factorforge_enzyme_registry (enzyme_name, recognition_sequence) VALUES (?, ?)",
                    [("BsaI", "GGTCTC"), ("BpiI", "GAAGAC"), ("BsmBI", "CGTCTC"), ("SapI", "GCTCTTC")],
                )

            # 4. Constraint Definitions
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS factorforge_constraint_definitions (
                constraint_id INTEGER PRIMARY KEY AUTOINCREMENT,
                code TEXT UNIQUE NOT NULL,
                name TEXT NOT NULL,
                rule_type TEXT NOT NULL
            );
            """)

            # Seed default constraints if empty
            cursor.execute("SELECT COUNT(*) FROM factorforge_constraint_definitions;")
            if cursor.fetchone()[0] == 0:
                cursor.executemany(
                    "INSERT INTO factorforge_constraint_definitions (code, name, rule_type) VALUES (?, ?, ?)",
                    [
                        ("AA_IDENTITY", "Amino Acid Identity Preservation", "HARD"),
                        ("INTERNAL_STOP", "No Internal Stop Codon", "HARD"),
                        ("TYPE_IIS_CLEAN", "Type IIS Restriction Site Absence", "HARD"),
                        ("GLOBAL_GC", "Global GC Content Window", "HARD"),
                        ("CAI_MIN", "Minimum Codon Adaptation Index", "SOFT"),
                    ],
                )

            # 5. Candidates & Evaluations
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS factorforge_candidates (
                candidate_id INTEGER PRIMARY KEY AUTOINCREMENT,
                sequence_id INTEGER,
                engine_name TEXT NOT NULL,
                model_version TEXT NOT NULL,
                overall_status TEXT NOT NULL,
                cai_value REAL,
                cai_reference_id TEXT,
                gc_percent REAL,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
            """)

            cursor.execute("""
            CREATE TABLE IF NOT EXISTS factorforge_constraint_evaluations (
                evaluation_id INTEGER PRIMARY KEY AUTOINCREMENT,
                candidate_id INTEGER,
                constraint_code TEXT NOT NULL,
                status TEXT NOT NULL,
                observed_value REAL,
                details_json TEXT
            );
            """)

            # 6. Design Packages (Frozen Manifests)
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS factorforge_design_packages (
                package_id INTEGER PRIMARY KEY AUTOINCREMENT,
                candidate_id INTEGER,
                construct_id TEXT NOT NULL,
                frozen_manifest_hash TEXT NOT NULL,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
            """)

            # Migration check: if old schema without candidate_id exists, drop and recreate
            cursor.execute("PRAGMA table_info(factorforge_design_packages);")
            cols = [col[1] for col in cursor.fetchall()]
            if "candidate_id" not in cols:
                cursor.execute("DROP TABLE factorforge_design_packages;")
                cursor.execute("""
                CREATE TABLE factorforge_design_packages (
                    package_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    candidate_id INTEGER,
                    construct_id TEXT NOT NULL,
                    frozen_manifest_hash TEXT NOT NULL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                );
                """)

            conn.commit()


    def register_sequence(self, raw_sequence: str, molecule_type: str = "CDS") -> int:
        """Registers a canonical sequence in common_sequence_registry with content hashing."""
        seq_hash = hashlib.sha256(raw_sequence.encode("utf-8")).hexdigest()
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(
                "SELECT sequence_id FROM common_sequence_registry WHERE sequence_hash = ?",
                (seq_hash,),
            )
            row = cursor.fetchone()
            if row:
                return row["sequence_id"]

            cursor.execute(
                """
                INSERT INTO common_sequence_registry (sequence_hash, molecule_type, sequence_length, raw_sequence)
                VALUES (?, ?, ?, ?)
                """,
                (seq_hash, molecule_type, len(raw_sequence), raw_sequence),
            )
            conn.commit()
            return cursor.lastrowid

    def save_candidate_with_evaluations(
        self,
        candidate_data: Dict[str, Any],
        evaluations: List[Dict[str, Any]],
        construct_id: str = "CF-CANDIDATE-001",
    ) -> Dict[str, Any]:
        """Saves a candidate, its step-by-step constraint evaluations, and creates a frozen design package."""
        raw_seq = candidate_data["optimized_sequence"]
        seq_id = self.register_sequence(raw_seq)

        with self._get_connection() as conn:
            cursor = conn.cursor()

            # Insert Candidate
            cursor.execute(
                """
                INSERT INTO factorforge_candidates
                (sequence_id, engine_name, model_version, overall_status, cai_value, cai_reference_id, gc_percent)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    seq_id,
                    candidate_data.get("engine", "lm"),
                    candidate_data.get("model_version", "v3.5.0-SynCodonLM-V2"),
                    "PASS" if candidate_data.get("type2is_clean", True) else "REJECT",
                    candidate_data.get("cai", 0.84),
                    "NB_codon_reference_v2",
                    candidate_data.get("gc_percent", 42.5),
                ),
            )
            candidate_id = cursor.lastrowid

            # Insert Step-by-Step Constraint Evaluations
            for ev in evaluations:
                cursor.execute(
                    """
                    INSERT INTO factorforge_constraint_evaluations
                    (candidate_id, constraint_code, status, observed_value, details_json)
                    VALUES (?, ?, ?, ?, ?)
                    """,
                    (
                        candidate_id,
                        ev["constraint_code"],
                        ev["status"],
                        ev.get("observed_value"),
                        json.dumps(ev.get("details", {})),
                    ),
                )

            # Freeze Design Package
            manifest_payload = json.dumps({"candidate_id": candidate_id, "sequence_id": seq_id, "construct_id": construct_id})
            manifest_hash = hashlib.sha256(manifest_payload.encode("utf-8")).hexdigest()

            cursor.execute(
                """
                INSERT INTO factorforge_design_packages (candidate_id, construct_id, frozen_manifest_hash)
                VALUES (?, ?, ?)
                """,
                (candidate_id, construct_id, manifest_hash),
            )
            pkg_id = cursor.lastrowid
            conn.commit()

            return {
                "package_id": pkg_id,
                "candidate_id": candidate_id,
                "sequence_id": seq_id,
                "frozen_manifest_hash": manifest_hash,
            }


if __name__ == "__main__":
    db = FactorForgeDBConnector()
    res = db.save_candidate_with_evaluations(
        candidate_data={
            "engine": "lm",
            "model_version": "v3.5.0-SynCodonLM-V2",
            "optimized_sequence": "ATGGCTAAATGGTAA",
            "cai": 0.84,
            "gc_percent": 41.2,
            "type2is_clean": True,
        },
        evaluations=[
            {"constraint_code": "AA_IDENTITY", "status": "PASS", "observed_value": 1.0, "details": {"target": "MAKW"}},
            {"constraint_code": "TYPE_IIS_CLEAN", "status": "PASS", "observed_value": 0, "details": {"enzymes": ["BsaI", "BpiI", "BsmBI"]}},
            {"constraint_code": "GLOBAL_GC", "status": "PASS", "observed_value": 41.2, "details": {"target_range": [35, 65]}},
        ],
        construct_id="CD9-MOCLO-001",
    )
    print("Upgraded Context-Aware Rule Engine Candidate & Evaluations Saved Successfully!")
    print(res)
