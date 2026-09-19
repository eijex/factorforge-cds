# factorforge/learning/ranker.py
"""Expression Prior Ranker v0 (Linear Preference Model on 8D Biophysical Traits).

Provides interpretable, noise-resilient candidate ranking:
    f(x) = w^T z(x) + b
trained with Weighted Margin Ranking Loss over empirical within-protein pairs.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import hashlib
import json
import math
from typing import Any, Dict, List, Optional


@dataclass
class RankerConfig:
    """Hyperparameters and configuration for ExpressionPriorRanker."""

    feature_names: List[str] = field(
        default_factory=lambda: [
            "cai",
            "gc_cds",
            "gc_5p",
            "mfe_5p",
            "cpi",
            "rare_codon_fraction",
            "poly_a_run_max",
            "poly_t_run_max",
        ]
    )
    margin_gamma: float = 0.1
    learning_rate: float = 0.05
    l2_regularization: float = 0.01
    max_epochs: int = 200
    early_stopping_patience: int = 15
    seed: int = 42


@dataclass
class RankerWeights:
    """Serialisable weights, normalization stats, and provenance for a trained ranker."""

    model_id: str
    model_version: str
    weights: List[float]
    bias: float
    feature_means: List[float]
    feature_stds: List[float]
    feature_names: List[str]
    training_snapshot_id: str
    trained_at_utc: str
    train_loss_history: List[float]
    train_pairwise_acc: float
    val_pairwise_acc: float
    weights_sha256: str = ""

    def __post_init__(self) -> None:
        if not self.weights_sha256:
            self.weights_sha256 = self.compute_sha256()

    def compute_sha256(self) -> str:
        payload = {
            "model_id": self.model_id,
            "weights": [round(w, 6) for w in self.weights],
            "bias": round(self.bias, 6),
            "means": [round(m, 6) for m in self.feature_means],
            "stds": [round(s, 6) for s in self.feature_stds],
            "snapshot_id": self.training_snapshot_id,
        }
        canonical_json = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        return hashlib.sha256(canonical_json).hexdigest()

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> RankerWeights:
        return cls(**data)


class ExpressionPriorRanker:
    """Interpretable Linear Expression Prior Ranker."""

    def __init__(
        self, config: Optional[RankerConfig] = None, weights: Optional[RankerWeights] = None
    ):
        self.config = config or RankerConfig()
        self.weights = weights
        self.num_features = len(self.config.feature_names)

    def is_trained(self) -> bool:
        return self.weights is not None

    def _normalize_traits(self, traits: List[float]) -> List[float]:
        if not self.weights:
            return list(traits)
        norm = []
        for val, m, s in zip(traits, self.weights.feature_means, self.weights.feature_stds):
            std = s if s > 1e-6 else 1.0
            norm.append((val - m) / std)
        return norm

    def score(self, trait_vector: List[float]) -> float:
        """Computes preference score f(x) = w^T z(x) + b for a single construct trait vector."""
        if not self.weights:
            # Fallback heuristic: standard balanced combination (CAI + 5' MFE opening)
            cai = trait_vector[0] if len(trait_vector) > 0 else 0.8
            mfe = trait_vector[3] if len(trait_vector) > 3 else -15.0
            return float(cai * 2.0 + (mfe / 50.0))

        z = self._normalize_traits(trait_vector)
        raw_score = sum(w * val for w, val in zip(self.weights.weights, z)) + self.weights.bias
        return float(raw_score)

    def rank_candidates(self, candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Ranks a list of candidate dictionaries containing 'trait_vector' descending by score."""
        scored = []
        for c in candidates:
            traits = c.get("trait_vector", [])
            s = self.score(traits)
            c_copy = dict(c)
            c_copy["expression_prior_score"] = round(s, 5)
            scored.append(c_copy)
        scored.sort(key=lambda x: x["expression_prior_score"], reverse=True)
        return scored

    def fit(
        self,
        train_pairs: List[Dict[str, Any]],
        val_pairs: Optional[List[Dict[str, Any]]] = None,
        snapshot_id: str = "snapshot_untracked",
    ) -> RankerWeights:
        """Trains the linear model using gradient descent on Weighted Margin Ranking Loss."""
        if not train_pairs:
            raise ValueError("Cannot train ExpressionPriorRanker: train_pairs is empty.")

        # 1. Compute normalization statistics across all observed trait vectors
        all_vectors: List[List[float]] = []
        for p in train_pairs:
            all_vectors.append(p["trait_vector_a"])
            all_vectors.append(p["trait_vector_b"])

        n_samples = len(all_vectors)
        means = [0.0] * self.num_features
        for vec in all_vectors:
            for d in range(self.num_features):
                means[d] += vec[d] / n_samples

        stds = [0.0] * self.num_features
        for vec in all_vectors:
            for d in range(self.num_features):
                stds[d] += (vec[d] - means[d]) ** 2
        stds = [math.sqrt(v / max(1, n_samples - 1)) if v > 1e-8 else 1.0 for v in stds]

        # 2. Normalize pair vectors
        def norm_vec(vec: List[float]) -> List[float]:
            return [(v - m) / s for v, m, s in zip(vec, means, stds)]

        # Initialize weights deterministically
        w = [0.0] * self.num_features
        b = 0.0
        lr = self.config.learning_rate
        l2 = self.config.l2_regularization
        gamma = self.config.margin_gamma

        loss_history: List[float] = []
        best_val_loss = float("inf")
        best_weights = list(w)
        best_bias = b
        patience_counter = 0

        for epoch in range(self.config.max_epochs):
            epoch_loss = 0.0
            grad_w = [0.0] * self.num_features
            grad_b = 0.0

            for p in train_pairs:
                z_a = norm_vec(p["trait_vector_a"])
                z_b = norm_vec(p["trait_vector_b"])
                y = float(p.get("preference_label", 1))
                weight = float(p.get("tier_weight", 1.0))

                s_a = sum(wi * zi for wi, zi in zip(w, z_a)) + b
                s_b = sum(wi * zi for wi, zi in zip(w, z_b)) + b

                diff = s_a - s_b
                margin_violation = gamma - y * diff

                if margin_violation > 0:
                    epoch_loss += weight * margin_violation
                    for d in range(self.num_features):
                        grad_w[d] += weight * (-y * (z_a[d] - z_b[d]))

            # Regularization
            epoch_loss += 0.5 * l2 * sum(wi**2 for wi in w)
            for d in range(self.num_features):
                grad_w[d] += l2 * w[d]

            n_pairs = len(train_pairs)
            for d in range(self.num_features):
                w[d] -= lr * (grad_w[d] / n_pairs)
            b -= lr * (grad_b / n_pairs)

            avg_loss = epoch_loss / n_pairs
            loss_history.append(avg_loss)

            # Validation evaluation
            if val_pairs:
                val_loss = 0.0
                for vp in val_pairs:
                    z_va = norm_vec(vp["trait_vector_a"])
                    z_vb = norm_vec(vp["trait_vector_b"])
                    vy = float(vp.get("preference_label", 1))
                    v_weight = float(vp.get("tier_weight", 1.0))
                    vs_a = sum(wi * zi for wi, zi in zip(w, z_va)) + b
                    vs_b = sum(wi * zi for wi, zi in zip(w, z_vb)) + b
                    v_violation = max(0.0, gamma - vy * (vs_a - vs_b))
                    val_loss += v_weight * v_violation
                avg_val_loss = val_loss / len(val_pairs)
                if avg_val_loss < best_val_loss:
                    best_val_loss = avg_val_loss
                    best_weights = list(w)
                    best_bias = b
                    patience_counter = 0
                else:
                    patience_counter += 1
                    if patience_counter >= self.config.early_stopping_patience:
                        w = list(best_weights)
                        b = best_bias
                        break

        # Compute accuracy
        def calc_acc(pairs: List[Dict[str, Any]], weights: List[float], bias: float) -> float:
            if not pairs:
                return 1.0
            correct = 0
            for p in pairs:
                z_a = norm_vec(p["trait_vector_a"])
                z_b = norm_vec(p["trait_vector_b"])
                y = float(p.get("preference_label", 1))
                s_a = sum(wi * zi for wi, zi in zip(weights, z_a)) + bias
                s_b = sum(wi * zi for wi, zi in zip(weights, z_b)) + bias
                pred_y = 1 if s_a >= s_b else -1
                if pred_y == y:
                    correct += 1
            return correct / len(pairs)

        train_acc = calc_acc(train_pairs, w, b)
        val_acc = calc_acc(val_pairs or [], w, b)

        now_utc = datetime.now(timezone.utc).isoformat()
        model_id = f"RANKER-V0-{hashlib.sha256(now_utc.encode('utf-8')).hexdigest()[:8]}"

        self.weights = RankerWeights(
            model_id=model_id,
            model_version="0.1.0",
            weights=list(w),
            bias=float(b),
            feature_means=list(means),
            feature_stds=list(stds),
            feature_names=list(self.config.feature_names),
            training_snapshot_id=snapshot_id,
            trained_at_utc=now_utc,
            train_loss_history=loss_history,
            train_pairwise_acc=train_acc,
            val_pairwise_acc=val_acc,
        )
        return self.weights
