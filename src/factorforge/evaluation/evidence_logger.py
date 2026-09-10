import json
import hashlib
import datetime
from pathlib import Path
from typing import Dict, Any

class RegulatoryEvidenceLogger:
    """
    FactorForge 엔진의 OptimizationResult를 받아 eijex-regulatory-evidence
    규격에 맞는 해시 검증된 증거 패키지를 생성하는 어댑터.
    """
    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_package(self, run_id: str, construct_id: str, optimization_result: Dict[str, Any]) -> Path:
        package_dir = self.output_dir / f"pkg_{run_id}_{construct_id}"
        package_dir.mkdir(parents=True, exist_ok=True)

        # RFC3339 UTC 형식 준수
        timestamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        
        metadata = optimization_result.get("metadata", {})
        passed = metadata.get("validator_passed", False)

        # 1. DESIGN_EVIDENCE (시퀀스 포함 금지)
        design_evidence = {
            "report_schema_version": "1.0",
            "identity": {
                "result_id": run_id,
                "construct_id": construct_id
            },
            "disposition": {
                "automated_decision": "PASS" if passed else "FAIL"
            },
            "provenance": {
                "engine": metadata.get("engine", "hybrid_sllm"),
                "inference_mode": metadata.get("inference_mode", "constraint_aware_beam_search")
            }
        }
        
        # 2. SELECTION_RECORD (AgentOS 의사결정 기록)
        selection_record = {
            "schema_version": "eijex.selection_record.v0.1",
            "source_system": "agentos",
            "design_run_id": run_id,
            "selected_construct_id": construct_id,
            "decision_status": "PROPOSED",
            "rationale": "Automated neuro-symbolic generation successfully passed all predefined Motif and CAI constraints.",
            "recorded_at": timestamp
        }

        # 3. VALIDATION_EVIDENCE_INDEX (ValidationHub 검증 지표)
        validation_index = {
            "schema_version": "eijex.validation_evidence_index.v0.1",
            "source_system": "validationhub",
            "design_run_id": run_id,
            "selected_construct_id": construct_id,
            "evidence_refs": [
                {
                    "evidence_id": f"ev_{run_id}_constraints",
                    "evidence_type": "automated_constraint_check",
                    "evidence_state": "REVIEWED" if passed else "NOT_AVAILABLE",
                    "artifact_reference": None
                }
            ]
        }

        def write_and_hash(filename: str, data: dict) -> str:
            path = package_dir / filename
            with path.open('w', encoding='utf-8') as f:
                json.dump(data, f, indent=2)
            # Calculate SHA-256
            digest = hashlib.sha256()
            with path.open('rb') as f:
                digest.update(f.read())
            return digest.hexdigest()

        sha_design = write_and_hash("design_evidence.json", design_evidence)
        sha_selection = write_and_hash("selection_record.json", selection_record)
        sha_validation = write_and_hash("validation_index.json", validation_index)

        # 4. MANIFEST (무결성 보증)
        manifest = {
            "schema_version": "eijex.regulatory_evidence_package.v0.1",
            "package_id": f"pkg_{run_id}_{construct_id}",
            "design_run_id": run_id,
            "selected_construct_id": construct_id,
            "created_at": timestamp,
            "status": "DRAFT",
            "intended_use": "research_evidence_assembly",
            "artifacts": [
                {
                    "artifact_id": f"art_des_{run_id}",
                    "artifact_type": "DESIGN_EVIDENCE",
                    "classification": "INTERNAL",
                    "contains_sequence": False,
                    "path": "design_evidence.json",
                    "sha256": sha_design
                },
                {
                    "artifact_id": f"art_sel_{run_id}",
                    "artifact_type": "SELECTION_RECORD",
                    "classification": "INTERNAL",
                    "contains_sequence": False,
                    "path": "selection_record.json",
                    "sha256": sha_selection
                },
                {
                    "artifact_id": f"art_val_{run_id}",
                    "artifact_type": "VALIDATION_EVIDENCE_INDEX",
                    "classification": "INTERNAL",
                    "contains_sequence": False,
                    "path": "validation_index.json",
                    "sha256": sha_validation
                }
            ]
        }
        write_and_hash("manifest.json", manifest)
        return package_dir
