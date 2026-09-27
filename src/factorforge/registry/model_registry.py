# factorforge/registry/model_registry.py
"""Immutable Model Registry for Production & Candidate Model Artifacts.

Decoupled from ScientificMemory (which stores biological observations & evidence).
Stores immutable release records, cryptographic weight digests, validation bundles,
and human promotion authorization records.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import hashlib
import json
import os
from typing import Any, Dict, List, Optional


@dataclass
class PromotionRecord:
    """Cryptographically indexed human authorization record for model promotion."""

    promotion_id: str
    model_id: str
    promoted_by: str
    promotion_rationale: str
    validation_bundle_sha256: str
    is_simulated: bool = False
    promoted_at_utc: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())


@dataclass
class ModelReleaseRecord:
    """Immutable metadata record for a released or staged model artifact."""

    model_id: str
    model_version: str
    model_type: str  # e.g., "ExpressionPriorRanker_v0"
    weights_sha256: str
    training_snapshot_id: str
    git_commit_sha: str
    status: str  # STAGED_CANDIDATE, BENCHMARK_CANDIDATE, ACTIVE_PRODUCTION, DEPRECATED, REJECTED, SIMULATION_ONLY
    validation_metrics: Dict[str, float]
    promotion_record: Optional[PromotionRecord] = None
    created_at_utc: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class ModelRegistry:
    """File-backed or in-memory immutable model registry."""

    def __init__(self, registry_dir: Optional[str] = None):
        self.registry_dir = registry_dir
        self._models: Dict[str, ModelReleaseRecord] = {}
        self._active_model_id: Optional[str] = None
        if self.registry_dir and os.path.exists(self.registry_dir):
            self._load_from_disk()

    def _load_from_disk(self) -> None:
        if not self.registry_dir:
            return
        manifest_path = os.path.join(self.registry_dir, "model_registry_manifest.json")
        if os.path.exists(manifest_path):
            with open(manifest_path, "r", encoding="utf-8") as f:
                data = json.load(f)
                self._active_model_id = data.get("active_model_id")
                for m_data in data.get("models", []):
                    rec = ModelReleaseRecord(**m_data)
                    self._models[rec.model_id] = rec

    def save_to_disk(self) -> None:
        if not self.registry_dir:
            return
        os.makedirs(self.registry_dir, exist_ok=True)
        manifest_path = os.path.join(self.registry_dir, "model_registry_manifest.json")
        payload = {
            "active_model_id": self._active_model_id,
            "models": [asdict(m) for m in self._models.values()],
            "updated_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)

    def register_candidate(
        self,
        model_id: str,
        model_version: str,
        weights_sha256: str,
        training_snapshot_id: str,
        validation_metrics: Dict[str, float],
        git_commit_sha: str = "HEAD",
        model_type: str = "ExpressionPriorRanker_v0",
        is_benchmark_pilot: bool = False,
    ) -> ModelReleaseRecord:
        """Stages a newly trained candidate model in the registry."""
        initial_status = "BENCHMARK_CANDIDATE" if is_benchmark_pilot else "STAGED_CANDIDATE"
        record = ModelReleaseRecord(
            model_id=model_id,
            model_version=model_version,
            model_type=model_type,
            weights_sha256=weights_sha256,
            training_snapshot_id=training_snapshot_id,
            git_commit_sha=git_commit_sha,
            status=initial_status,
            validation_metrics=validation_metrics,
        )
        self._models[model_id] = record
        self.save_to_disk()
        return record

    def stage_for_promotion(self, model_id: str, rationale: str) -> ModelReleaseRecord:
        """Sets status to AWAITING_HUMAN_PROMOTION."""
        if model_id not in self._models:
            raise KeyError(f"Model ID '{model_id}' not found in registry.")
        rec = self._models[model_id]
        rec.status = "AWAITING_HUMAN_PROMOTION"
        self.save_to_disk()
        return rec

    def promote_model(
        self,
        model_id: str,
        promoted_by: str,
        rationale: str,
        validation_bundle_sha256: str,
        is_simulated: bool = False,
    ) -> ModelReleaseRecord:
        """Promotes a candidate model with an immutable PromotionRecord.

        If is_simulated is True, status is set to SIMULATION_ONLY rather than ACTIVE_PRODUCTION.
        """
        if model_id not in self._models:
            raise KeyError(f"Model ID '{model_id}' not found in registry.")

        rec = self._models[model_id]
        promotion_id = (
            f"PROM-{hashlib.sha256(f'{model_id}:{promoted_by}'.encode('utf-8')).hexdigest()[:8]}"
        )
        prom_record = PromotionRecord(
            promotion_id=promotion_id,
            model_id=model_id,
            promoted_by=promoted_by,
            promotion_rationale=rationale,
            validation_bundle_sha256=validation_bundle_sha256,
            is_simulated=is_simulated,
        )

        if is_simulated:
            rec.status = "SIMULATION_ONLY"
        else:
            # Deprecate previous active model
            if self._active_model_id and self._active_model_id in self._models:
                self._models[self._active_model_id].status = "DEPRECATED"
            rec.status = "ACTIVE_PRODUCTION"
            self._active_model_id = model_id

        rec.promotion_record = prom_record
        self.save_to_disk()
        return rec

    def reject_model(self, model_id: str, reason: str) -> ModelReleaseRecord:
        """Marks a candidate model as REJECTED."""
        if model_id not in self._models:
            raise KeyError(f"Model ID '{model_id}' not found in registry.")
        rec = self._models[model_id]
        rec.status = "REJECTED"
        self.save_to_disk()
        return rec

    def roll_back_active_model(self, rationale: str) -> Optional[str]:
        """Rolls back active production model if unvalidated claims were identified."""
        prev = self._active_model_id
        if prev and prev in self._models:
            self._models[prev].status = "ROLLED_BACK"
            self._active_model_id = None
            self.save_to_disk()
        return prev

    def get_active_model(self) -> Optional[ModelReleaseRecord]:
        if self._active_model_id:
            return self._models.get(self._active_model_id)
        return None

    def get_model(self, model_id: str) -> Optional[ModelReleaseRecord]:
        return self._models.get(model_id)

    def list_models(self) -> List[ModelReleaseRecord]:
        return list(self._models.values())
