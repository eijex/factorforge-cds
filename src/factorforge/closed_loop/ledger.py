"""Append-only Evidence Ledger for ValidationHub / FactorForge closed loop."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, List, Optional

from factorforge.closed_loop.contracts import (
    PolicyRecord,
    ConstructRecord,
    ConstructSetRecord,
    ExperimentRecord,
    MeasurementRecord,
)
from factorforge.closed_loop.readiness import (
    ReadinessStatus,
    ReadinessReport,
    ReadinessEvaluator,
)


class EvidenceLedger:
    """Local file-based append-only registry storing closed-loop DBTL evidence."""

    def __init__(self, ledger_dir: Optional[Path] = None):
        if ledger_dir is None:
            ledger_dir = Path(__file__).resolve().parents[3] / "data" / "evidence_ledger"
        self.ledger_dir = Path(ledger_dir)
        self.ledger_dir.mkdir(parents=True, exist_ok=True)

        self.policies_file = self.ledger_dir / "policies.jsonl"
        self.constructs_file = self.ledger_dir / "constructs.jsonl"
        self.construct_sets_file = self.ledger_dir / "construct_sets.jsonl"
        self.experiments_file = self.ledger_dir / "experiments.jsonl"
        self.measurements_file = self.ledger_dir / "measurements.jsonl"

    def _append_jsonl(self, file_path: Path, record: dict):
        with open(file_path, "a", encoding="utf-8") as f:
            f.write(json.dumps(record, ensure_ascii=False) + "\n")

    def _read_jsonl(self, file_path: Path) -> List[dict]:
        if not file_path.exists():
            return []
        records = []
        with open(file_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    records.append(json.loads(line))
        return records

    def register_policy(self, record: PolicyRecord):
        self._append_jsonl(self.policies_file, record.model_dump())

    def register_construct(self, record: ConstructRecord):
        self._append_jsonl(self.constructs_file, record.model_dump())

    def register_construct_set(self, record: ConstructSetRecord):
        self._append_jsonl(self.construct_sets_file, record.model_dump())

    def register_experiment(self, record: ExperimentRecord):
        self._append_jsonl(self.experiments_file, record.model_dump())

    def record_measurement(self, record: MeasurementRecord):
        self._append_jsonl(self.measurements_file, record.model_dump())

    def list_measurements(self, include_synthetic: bool = True) -> List[MeasurementRecord]:
        raw = self._read_jsonl(self.measurements_file)
        measurements = [MeasurementRecord(**r) for r in raw]
        if not include_synthetic:
            return [m for m in measurements if not m.is_synthetic]
        return measurements

    def get_readiness(self) -> ReadinessReport:
        all_measurements = self.list_measurements(include_synthetic=True)
        return ReadinessEvaluator.evaluate(all_measurements)

    def verify_lineage(self, measurement_id: str) -> Dict[str, object]:
        """Verify unbroken provenance: Measurement -> Experiment -> ConstructSet -> Construct -> Policy."""
        measurements = {m["measurement_id"]: m for m in self._read_jsonl(self.measurements_file)}
        experiments = {e["experiment_id"]: e for e in self._read_jsonl(self.experiments_file)}
        construct_sets = {s["construct_set_id"]: s for s in self._read_jsonl(self.construct_sets_file)}
        constructs = {c["construct_id"]: c for c in self._read_jsonl(self.constructs_file)}
        policies = {p["policy_id"]: p for p in self._read_jsonl(self.policies_file)}

        if measurement_id not in measurements:
            raise KeyError(f"Measurement '{measurement_id}' not found in ledger")

        m = measurements[measurement_id]
        exp_id = m.get("experiment_id")
        set_id = m.get("construct_set_id")

        exp = experiments.get(exp_id)
        cset = construct_sets.get(set_id)

        lc_construct = constructs.get(cset.get("light_construct_id")) if cset else None
        hc_construct = constructs.get(cset.get("heavy_construct_id")) if cset else None

        policy_id = lc_construct.get("policy_id") if lc_construct else None
        pol = policies.get(policy_id) if policy_id else None

        lineage_intact = (
            exp is not None
            and cset is not None
            and lc_construct is not None
            and hc_construct is not None
            and pol is not None
        )

        return {
            "measurement_id": measurement_id,
            "lineage_intact": lineage_intact,
            "is_synthetic": m.get("is_synthetic", False),
            "experiment": exp,
            "construct_set": cset,
            "light_construct": lc_construct,
            "heavy_construct": hc_construct,
            "policy": pol,
        }

    def sync_to_gcs(
        self,
        bucket_name: str,
        destination_folder: str = "validationhub/ledger_v0",
        gcp_key_path: Optional[str] = None,
    ) -> Dict[str, str]:
        """Sync append-only ledger snapshot to Google Cloud Storage."""
        from google.cloud import storage

        if gcp_key_path:
            os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = gcp_key_path

        client = storage.Client()
        bucket = client.bucket(bucket_name)

        uploaded = {}
        for file_path in [
            self.policies_file,
            self.constructs_file,
            self.construct_sets_file,
            self.experiments_file,
            self.measurements_file,
        ]:
            if file_path.exists():
                blob_name = f"{destination_folder}/{file_path.name}"
                blob = bucket.blob(blob_name)
                blob.upload_from_filename(str(file_path))
                uploaded[file_path.name] = f"gs://{bucket_name}/{blob_name}"

        return uploaded
