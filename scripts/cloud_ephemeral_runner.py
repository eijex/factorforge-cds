"""FactorForge Ephemeral Cloud Runner.

Safely provisions a temporary VM on GCP Compute Engine, executes an experiment,
verifies output on GCS, and attempts teardown/deletion of the VM.
Deletion failures and residual cloud costs require independent monitoring.
Enforces munkyu.kim@eijex.com GCP account policy (blocks legacy free-credit project).
"""

from __future__ import annotations

import os
from pathlib import Path
import time

from google.cloud import storage
from google.oauth2 import service_account
from googleapiclient import discovery
from googleapiclient.errors import HttpError

from factorforge.telemetry.training_logger import verify_gcp_account_policy

GCP_ACCOUNT = os.environ.get("EIJEX_GCP_ACCOUNT", "munkyu.kim@eijex.com")
KEY_PATH = os.environ.get("EIJEX_GCP_KEY_PATH", str(Path.home() / ".gcp" / "factorforge-runner-key.json"))
PROJECT_ID = os.environ.get("EIJEX_GCP_PROJECT_ID", "fair-conduit-510507-p5")
ZONE = os.environ.get("EIJEX_GCP_ZONE", "asia-northeast3-a")
BUCKET_NAME = os.environ.get("EIJEX_GCP_BUCKET", "factorforge-fair-conduit")
INSTANCE_NAME = "factorforge-exp-runner"
SERVICE_ACCOUNT_EMAIL = f"factorforge-runner@{PROJECT_ID}.iam.gserviceaccount.com"


def get_compute_client():
    verified = verify_gcp_account_policy(
        project_id=PROJECT_ID,
        account_email=GCP_ACCOUNT,
        key_path=KEY_PATH,
        bucket_name=BUCKET_NAME,
    )
    credentials = service_account.Credentials.from_service_account_file(verified["key_path"])
    return discovery.build("compute", "v1", credentials=credentials)


def get_storage_client():
    verified = verify_gcp_account_policy(
        project_id=PROJECT_ID,
        account_email=GCP_ACCOUNT,
        key_path=KEY_PATH,
        bucket_name=BUCKET_NAME,
    )
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = verified["key_path"]
    return storage.Client(project=verified["project_id"])


def create_ephemeral_instance(compute, instance_name: str = INSTANCE_NAME):
    # Startup script with triple safety:
    # 1. Fallback self-shutdown in 3 minutes if anything hangs
    # 2. Runs the experiment and uploads proof to GCS
    # 3. Powers off OS immediately upon completion
    startup_script = f"""#!/bin/bash
set -e

# Safety Guardrail: OS-level hard shutdown in 4 minutes
shutdown -h +4 &

echo "=== FactorForge Cloud Experiment Started ===" > /tmp/result.txt
echo "Instance: $(hostname)" >> /tmp/result.txt
echo "Date: $(date -u)" >> /tmp/result.txt
echo "Kernel: $(uname -r)" >> /tmp/result.txt
echo "Status: SUCCESS - FactorForge Cloud Node is healthy and functioning" >> /tmp/result.txt

# Upload result to GCS
gsutil cp /tmp/result.txt gs://{BUCKET_NAME}/experiments/{instance_name}_output.txt

echo "=== Experiment Finished. Powering off now ===" >> /tmp/result.txt
poweroff
"""

    config = {
        "name": instance_name,
        "machineType": f"zones/{ZONE}/machineTypes/e2-micro",
        "disks": [
            {
                "boot": True,
                "autoDelete": True,
                "initializeParams": {
                    "sourceImage": "projects/debian-cloud/global/images/family/debian-12",
                    "diskSizeGb": "10",
                },
            }
        ],
        "networkInterfaces": [
            {
                "network": "global/networks/default",
                "accessConfigs": [{"type": "ONE_TO_ONE_NAT", "name": "External NAT"}],
            }
        ],
        "serviceAccounts": [
            {
                "email": SERVICE_ACCOUNT_EMAIL,
                "scopes": ["https://www.googleapis.com/auth/cloud-platform"],
            }
        ],
        "metadata": {
            "items": [
                {
                    "key": "startup-script",
                    "value": startup_script,
                }
            ]
        },
    }

    print(f"[1/4] Launching ephemeral VM '{instance_name}' in {ZONE}...")
    operation = compute.instances().insert(project=PROJECT_ID, zone=ZONE, body=config).execute()
    return operation


def wait_for_zone_operation(compute, operation_name: str):
    while True:
        result = (
            compute.zoneOperations()
            .get(project=PROJECT_ID, zone=ZONE, operation=operation_name)
            .execute()
        )
        status = result.get("status")
        if status == "DONE":
            if "error" in result:
                raise RuntimeError(result["error"])
            return result
        time.sleep(2)


def delete_instance(compute, instance_name: str = INSTANCE_NAME):
    print(
        f"[Teardown] Deleting VM '{instance_name}' (residual resources and billing require separate checks)..."
    )
    try:
        op = (
            compute.instances()
            .delete(project=PROJECT_ID, zone=ZONE, instance=instance_name)
            .execute()
        )
        wait_for_zone_operation(compute, op["name"])
        print(f"[Teardown] VM '{instance_name}' permanently DELETED.")
    except (HttpError, OSError, TimeoutError, ValueError) as e:
        print(f"[Teardown Warning] Error deleting VM: {e}")


def run_experiment_lifecycle():
    compute = get_compute_client()
    storage_client = get_storage_client()
    bucket = storage_client.bucket(BUCKET_NAME)

    # Clean previous test artifact if any
    blob = bucket.blob(f"experiments/{INSTANCE_NAME}_output.txt")
    if blob.exists():
        blob.delete()

    start_time = time.time()
    try:
        op = create_ephemeral_instance(compute, INSTANCE_NAME)
        wait_for_zone_operation(compute, op["name"])
        print(f"[2/4] VM '{INSTANCE_NAME}' is now PROVISIONED and RUNNING!")
        print(
            ">> GCP 콘솔 화면에서 [새로고침(Refresh)]을 누르시면 녹색 체크와 함께 VM이 떠 있는 것을 직접 보실 수 있습니다!"
        )

        # Wait for experiment output in GCS or VM shutdown (max 3 minutes)
        print("[3/4] Waiting for experiment completion & auto-shutdown (max 180s)...")
        completed = False
        for _ in range(36):
            time.sleep(5)
            # Check if output uploaded to GCS
            if blob.exists():
                print("[3/4] Output successfully uploaded to GCS!")
                output_content = blob.download_as_text()
                print("-------------------- VM Output --------------------")
                print(output_content.strip())
                print("---------------------------------------------------")
                completed = True
                break

            # Check instance status
            try:
                inst = (
                    compute.instances()
                    .get(project=PROJECT_ID, zone=ZONE, instance=INSTANCE_NAME)
                    .execute()
                )
                status = inst.get("status")
                if status == "TERMINATED":
                    print(f"[3/4] VM has powered off itself (status: {status}).")
                    completed = True
                    break
            except (HttpError, OSError, TimeoutError, ValueError):
                print(
                    "[Polling Warning] Cloud status check failed; retrying within the existing limit."
                )

        if not completed:
            print("[Warning] Timeout reached. Initiating forced teardown.")

    finally:
        # Attempt teardown even when execution fails
        print("[4/4] Attempting instance teardown...")
        delete_instance(compute, INSTANCE_NAME)
        elapsed = time.time() - start_time
        print(
            f"[Complete] Total cycle time: {elapsed:.1f}s. Cost not measured; verify cloud billing and residual resources."
        )
        print(
            ">> GCP 콘솔에서 다시 [새로고침(Refresh)]을 누르시면 인스턴스가 완전히 사라진 깨끗한 화면을 확인하실 수 있습니다."
        )


if __name__ == "__main__":
    run_experiment_lifecycle()
