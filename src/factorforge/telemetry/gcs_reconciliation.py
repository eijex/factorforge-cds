"""Durable Reconciliation Contract for Telemetry and Checkpoints (Job 448).

Local reconciliation checks (not a deployed lifecycle integration):
1. Reconciles local training telemetry (SQLite / JSONL) and model checkpoints against
   storage-client supplied object digests. No ValidationHub adapter is implemented.
2. Content verification: exact SHA-256 matching for every artifact file.
3. Fail-Closed VM Termination Gate (`assert_safe_to_terminate`):
   Raises when called with unverified or stale artifacts; VM runners must call
   this gate explicitly. Writers must be quiescent during reconciliation.
4. GCP Account & Project Policy Enforcement:
   Checks configured credentials, not authenticated human ownership.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Protocol

from factorforge.telemetry.training_logger import (
    AUTHORIZED_GCP_PROJECT_ID,
    AUTHORIZED_GCS_BUCKET,
    REQUIRED_GCP_ACCOUNT,
    verify_gcp_account_policy,
)


class DurableReconciliationError(Exception):
    """Raised when reconciliation fails or when unsafe termination is attempted."""


class StorageClientProtocol(Protocol):
    """Protocol for remote storage checks (e.g. Google Cloud Storage or mock)."""

    def object_exists(self, bucket: str, object_name: str) -> bool: ...

    def get_object_sha256(self, bucket: str, object_name: str) -> str | None: ...


@dataclass(frozen=True)
class ArtifactAuditItem:
    relative_path: str
    sha256: str
    size_bytes: int
    remote_verified: bool
    remote_sha256: str | None = None
    status: str = "PENDING"


class DurableReconciler:
    """Reconciles local SQLite/JSONL telemetry and model artifacts with durable GCS storage."""

    def __init__(
        self,
        telemetry_db_path: Path | str,
        artifacts_dir: Path | str,
        gcs_bucket: str = AUTHORIZED_GCS_BUCKET,
        gcp_project: str = AUTHORIZED_GCP_PROJECT_ID,
        account_email: str = REQUIRED_GCP_ACCOUNT,
        key_path: Path | str | None = None,
    ) -> None:
        self.telemetry_db_path = Path(telemetry_db_path)
        self.artifacts_dir = Path(artifacts_dir)
        self.gcs_bucket = gcs_bucket
        self.gcp_project = gcp_project
        self.account_email = account_email
        self.key_path = key_path
        self._verified_result_digest: str | None = None
        self._storage_client: StorageClientProtocol | None = None

        # Enforce account and project policy
        verify_gcp_account_policy(
            account_email=self.account_email,
            project_id=self.gcp_project,
            bucket_name=self.gcs_bucket,
            key_path=self.key_path,
        )

    @staticmethod
    def compute_sha256(file_path: Path) -> str:
        h = hashlib.sha256()
        with open(file_path, "rb") as f:
            while chunk := f.read(65536):
                h.update(chunk)
        return h.hexdigest()

    def audit_local_artifacts(self) -> dict[str, dict[str, Any]]:
        """Audits all local artifacts and telemetry files."""
        items: dict[str, dict[str, Any]] = {}

        if not self.telemetry_db_path.is_file() or self.telemetry_db_path.is_symlink():
            raise DurableReconciliationError("Required telemetry database is missing or linked")
        if not self.artifacts_dir.is_dir() or self.artifacts_dir.is_symlink():
            raise DurableReconciliationError("Required artifact directory is missing or linked")
        for suffix in ("-wal", "-journal"):
            sidecar = Path(str(self.telemetry_db_path) + suffix)
            if sidecar.exists() and sidecar.stat().st_size:
                raise DurableReconciliationError("Quiesce/checkpoint SQLite before reconciliation")

        # 1. Telemetry DB
        if self.telemetry_db_path.exists():
            db_sha = self.compute_sha256(self.telemetry_db_path)
            items["telemetry.db"] = {
                "relative_path": "telemetry.db",
                "sha256": db_sha,
                "size_bytes": self.telemetry_db_path.stat().st_size,
            }

        # 2. Artifacts directory
        if self.artifacts_dir.exists():
            for p in self.artifacts_dir.rglob("*"):
                if p.is_symlink():
                    raise DurableReconciliationError("Linked artifacts are not accepted")
                if p.is_file():
                    rel = p.relative_to(self.artifacts_dir).as_posix()
                    if rel in items:
                        raise DurableReconciliationError(
                            "Artifact name collides with telemetry database"
                        )
                    items[rel] = {
                        "relative_path": rel,
                        "sha256": self.compute_sha256(p),
                        "size_bytes": p.stat().st_size,
                    }

        if not items:
            raise DurableReconciliationError(
                "No local artifacts or telemetry database found to reconcile."
            )

        return items

    def reconcile_with_storage(
        self,
        storage_client: StorageClientProtocol,
    ) -> dict[str, Any]:
        """Compares local artifacts against durable remote storage."""
        self._verified_result_digest = None
        self._storage_client = None
        local_items = self.audit_local_artifacts()
        audited_items: list[dict[str, Any]] = []
        errors: list[str] = []

        all_ok = True
        for rel_path, meta in local_items.items():
            remote_name = f"runs/{rel_path}"
            exists = storage_client.object_exists(self.gcs_bucket, remote_name)
            if not exists:
                all_ok = False
                err = f"Remote object missing: gs://{self.gcs_bucket}/{remote_name}"
                errors.append(err)
                audited_items.append(
                    {
                        "relative_path": rel_path,
                        "sha256": meta["sha256"],
                        "size_bytes": meta["size_bytes"],
                        "remote_verified": False,
                        "remote_sha256": None,
                        "status": "MISSING_ON_REMOTE",
                    }
                )
                continue

            remote_sha = storage_client.get_object_sha256(self.gcs_bucket, remote_name)
            if remote_sha != meta["sha256"]:
                all_ok = False
                err = (
                    f"SHA-256 digest mismatch for gs://{self.gcs_bucket}/{remote_name}: "
                    f"local={meta['sha256']}, remote={remote_sha}"
                )
                errors.append(err)
                audited_items.append(
                    {
                        "relative_path": rel_path,
                        "sha256": meta["sha256"],
                        "size_bytes": meta["size_bytes"],
                        "remote_verified": False,
                        "remote_sha256": remote_sha,
                        "status": "DIGEST_MISMATCH",
                    }
                )
                continue

            audited_items.append(
                {
                    "relative_path": rel_path,
                    "sha256": meta["sha256"],
                    "size_bytes": meta["size_bytes"],
                    "remote_verified": True,
                    "remote_sha256": remote_sha,
                    "status": "RECONCILED",
                }
            )

        if self.audit_local_artifacts() != local_items:
            raise DurableReconciliationError("Local artifacts changed during reconciliation")
        manifest = {
            "reconciliation_time_utc": datetime.now(timezone.utc).isoformat(),
            "status": "RECONCILED" if all_ok else "RECONCILIATION_FAILED",
            "gcs_bucket": self.gcs_bucket,
            "gcp_project": self.gcp_project,
            "account_email": self.account_email,
            "total_artifacts": len(local_items),
            "reconciled_count": sum(1 for x in audited_items if x["status"] == "RECONCILED"),
            "errors": errors,
            "artifacts": audited_items,
        }
        if all_ok:
            self._verified_result_digest = self._result_digest(manifest)
            self._storage_client = storage_client
        return manifest

    @staticmethod
    def _result_digest(result: dict[str, Any]) -> str:
        return hashlib.sha256(
            json.dumps(result, sort_keys=True, allow_nan=False).encode()
        ).hexdigest()

    def assert_safe_to_terminate(self, reconciliation_result: dict[str, Any]) -> None:
        """Fail-Closed gate: raises DurableReconciliationError if unsafe to terminate instance."""
        if (
            self._verified_result_digest is None
            or self._result_digest(reconciliation_result) != self._verified_result_digest
            or reconciliation_result.get("status") != "RECONCILED"
        ):
            errors = reconciliation_result.get("errors", [])
            err_msg = "; ".join(errors) if errors else "Unspecified reconciliation failure"
            raise DurableReconciliationError(
                f"UNSAFE_TERMINATION_PREVENTED: Cannot terminate ephemeral VM. "
                f"Artifacts have not been reconciled with durable storage: {err_msg}"
            )
        # Recheck both local and remote contents; a previously exported JSON is not authority.
        prior = {
            x["relative_path"]: (x["sha256"], x["size_bytes"])
            for x in reconciliation_result["artifacts"]
        }
        current = self.audit_local_artifacts()
        if {k: (v["sha256"], v["size_bytes"]) for k, v in current.items()} != prior:
            self._verified_result_digest = None
            raise DurableReconciliationError(
                "UNSAFE_TERMINATION_PREVENTED: Local artifacts changed"
            )
        fresh = self.reconcile_with_storage(self._storage_client)
        if fresh["status"] != "RECONCILED":
            raise DurableReconciliationError(
                "UNSAFE_TERMINATION_PREVENTED: Remote artifacts changed"
            )

    def write_manifest(
        self,
        output_path: Path | str,
        reconciliation_result: dict[str, Any],
    ) -> str:
        """Writes an unsigned manifest and returns its content digest (not a signature)."""
        out = Path(output_path).resolve()
        if out == self.telemetry_db_path.resolve() or self.artifacts_dir.resolve() in out.parents:
            raise DurableReconciliationError(
                "Manifest output must not overwrite or join the input inventory"
            )
        out.parent.mkdir(parents=True, exist_ok=True)
        content = json.dumps(reconciliation_result, indent=2, sort_keys=True, allow_nan=False)
        out.write_text(content, encoding="utf-8")
        return hashlib.sha256(content.encode("utf-8")).hexdigest()
