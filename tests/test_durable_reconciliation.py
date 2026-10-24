"""Unit tests for Durable Telemetry and Checkpoint Reconciliation Contract (Job 448)."""

import json
from pathlib import Path

import pytest

from factorforge.telemetry.gcs_reconciliation import (
    DurableReconciler,
    DurableReconciliationError,
)
from factorforge.telemetry.training_logger import (
    AUTHORIZED_GCP_PROJECT_ID,
    AUTHORIZED_GCS_BUCKET,
    REQUIRED_GCP_ACCOUNT,
    GCPAccountPolicyError,
)


@pytest.fixture(autouse=True)
def synthetic_credentials(tmp_path, monkeypatch):
    key = tmp_path / "synthetic-key.json"
    key.write_text(
        json.dumps(
            {
                "type": "service_account",
                "project_id": AUTHORIZED_GCP_PROJECT_ID,
                "client_email": f"factorforge-runner@{AUTHORIZED_GCP_PROJECT_ID}.iam.gserviceaccount.com",
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setenv("EIJEX_GCP_KEY_PATH", str(key))


class MockStorageClient:
    """Mock GCS storage client for testing reconciliation."""

    def __init__(self, objects: dict[str, str] | None = None) -> None:
        # mapping of (bucket, object_name) -> sha256
        self.objects: dict[tuple[str, str], str] = {}
        if objects:
            for (bucket, obj), digest in objects.items():
                self.objects[(bucket, obj)] = digest

    def object_exists(self, bucket: str, object_name: str) -> bool:
        return (bucket, object_name) in self.objects

    def get_object_sha256(self, bucket: str, object_name: str) -> str | None:
        return self.objects.get((bucket, object_name))


def test_successful_reconciliation_allows_termination(tmp_path: Path):
    db_file = tmp_path / "telemetry.db"
    db_file.write_text("sqlite mock content", encoding="utf-8")
    art_dir = tmp_path / "artifacts"
    art_dir.mkdir()
    chk_file = art_dir / "checkpoint.pt"
    chk_file.write_text("weights mock", encoding="utf-8")

    db_sha = DurableReconciler.compute_sha256(db_file)
    chk_sha = DurableReconciler.compute_sha256(chk_file)

    storage = MockStorageClient(
        {
            (AUTHORIZED_GCS_BUCKET, "runs/telemetry.db"): db_sha,
            (AUTHORIZED_GCS_BUCKET, "runs/checkpoint.pt"): chk_sha,
        }
    )

    reconciler = DurableReconciler(
        telemetry_db_path=db_file,
        artifacts_dir=art_dir,
        gcs_bucket=AUTHORIZED_GCS_BUCKET,
        gcp_project=AUTHORIZED_GCP_PROJECT_ID,
        account_email=REQUIRED_GCP_ACCOUNT,
    )

    manifest = reconciler.reconcile_with_storage(storage)
    assert manifest["status"] == "RECONCILED"
    assert manifest["reconciled_count"] == 2
    assert not manifest["errors"]

    # Termination check passes without raising
    reconciler.assert_safe_to_terminate(manifest)

    # Manifest export
    manifest_path = tmp_path / "reconciliation_manifest.json"
    manifest_sha = reconciler.write_manifest(manifest_path, manifest)
    assert len(manifest_sha) == 64
    assert manifest_path.exists()


def test_missing_remote_artifact_blocks_termination(tmp_path: Path):
    db_file = tmp_path / "telemetry.db"
    db_file.write_text("sqlite mock", encoding="utf-8")
    art_dir = tmp_path / "artifacts"
    art_dir.mkdir()
    chk_file = art_dir / "checkpoint.pt"
    chk_file.write_text("weights mock", encoding="utf-8")

    db_sha = DurableReconciler.compute_sha256(db_file)
    # Storage is missing checkpoint.pt!
    storage = MockStorageClient(
        {
            (AUTHORIZED_GCS_BUCKET, "runs/telemetry.db"): db_sha,
        }
    )

    reconciler = DurableReconciler(
        telemetry_db_path=db_file,
        artifacts_dir=art_dir,
        gcs_bucket=AUTHORIZED_GCS_BUCKET,
        gcp_project=AUTHORIZED_GCP_PROJECT_ID,
        account_email=REQUIRED_GCP_ACCOUNT,
    )

    manifest = reconciler.reconcile_with_storage(storage)
    assert manifest["status"] == "RECONCILIATION_FAILED"
    assert len(manifest["errors"]) == 1

    with pytest.raises(DurableReconciliationError, match="UNSAFE_TERMINATION_PREVENTED"):
        reconciler.assert_safe_to_terminate(manifest)


def test_digest_mismatch_blocks_termination(tmp_path: Path):
    db_file = tmp_path / "telemetry.db"
    db_file.write_text("sqlite mock", encoding="utf-8")
    art_dir = tmp_path / "artifacts"
    art_dir.mkdir()

    storage = MockStorageClient(
        {
            (AUTHORIZED_GCS_BUCKET, "runs/telemetry.db"): "0" * 64,  # wrong digest
        }
    )

    reconciler = DurableReconciler(
        telemetry_db_path=db_file,
        artifacts_dir=art_dir,
        gcs_bucket=AUTHORIZED_GCS_BUCKET,
        gcp_project=AUTHORIZED_GCP_PROJECT_ID,
        account_email=REQUIRED_GCP_ACCOUNT,
    )

    manifest = reconciler.reconcile_with_storage(storage)
    assert manifest["status"] == "RECONCILIATION_FAILED"
    assert any("mismatch" in e.lower() for e in manifest["errors"])

    with pytest.raises(DurableReconciliationError, match="UNSAFE_TERMINATION_PREVENTED"):
        reconciler.assert_safe_to_terminate(manifest)


def test_unauthorized_account_or_project_rejected(tmp_path: Path):
    db_file = tmp_path / "telemetry.db"
    db_file.write_text("content", encoding="utf-8")
    art_dir = tmp_path / "artifacts"
    art_dir.mkdir()

    with pytest.raises(GCPAccountPolicyError, match=r"is not authorized"):
        DurableReconciler(
            telemetry_db_path=db_file,
            artifacts_dir=art_dir,
            account_email="malicious@external.com",
        )

    with pytest.raises(GCPAccountPolicyError, match=r"blocked|Use"):
        DurableReconciler(
            telemetry_db_path=db_file,
            artifacts_dir=art_dir,
            gcp_project="unauthorized-gcp-project",
        )


def test_forged_success_and_missing_db_are_rejected(tmp_path):
    art = tmp_path / "artifacts"
    art.mkdir()
    (art / "checkpoint.pt").write_bytes(b"fixture")
    r = DurableReconciler(tmp_path / "missing.db", art)
    with pytest.raises(DurableReconciliationError):
        r.assert_safe_to_terminate({"status": "RECONCILED"})
    with pytest.raises(DurableReconciliationError, match="database"):
        r.audit_local_artifacts()


@pytest.mark.parametrize("change", ["local", "remote", "manifest", "wal", "collision"])
def test_stale_or_conflicting_snapshots_are_rejected(tmp_path, change):
    db = tmp_path / "ledger.db"
    db.write_bytes(b"fixture-db")
    art = tmp_path / "artifacts"
    art.mkdir()
    f = art / "checkpoint.pt"
    f.write_bytes(b"fixture-checkpoint")
    storage = MockStorageClient(
        {
            (AUTHORIZED_GCS_BUCKET, "runs/telemetry.db"): DurableReconciler.compute_sha256(db),
            (AUTHORIZED_GCS_BUCKET, "runs/checkpoint.pt"): DurableReconciler.compute_sha256(f),
        }
    )
    r = DurableReconciler(db, art)
    result = r.reconcile_with_storage(storage)
    if change == "local":
        f.write_bytes(b"changed")
    if change == "remote":
        storage.objects.clear()
    if change == "manifest":
        result["artifacts"] = []
    if change == "wal":
        Path(str(db) + "-wal").write_bytes(b"pending-write")
    if change == "collision":
        (art / "telemetry.db").write_bytes(b"collision")
    with pytest.raises(DurableReconciliationError):
        r.assert_safe_to_terminate(result)


@pytest.mark.parametrize("destination", ["db", "artifact"])
def test_manifest_export_cannot_overwrite_inventory(tmp_path, destination):
    db = tmp_path / "ledger.db"
    db.write_bytes(b"unchanged-fixture")
    art = tmp_path / "artifacts"
    art.mkdir()
    r = DurableReconciler(db, art)
    with pytest.raises(DurableReconciliationError, match="inventory"):
        r.write_manifest(db if destination == "db" else art / "manifest.json", {})
    assert db.read_bytes() == b"unchanged-fixture"
