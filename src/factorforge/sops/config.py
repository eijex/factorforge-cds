import hashlib
import json
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Tuple

from factorforge.rules.models import EnforcementLevel, AuthorizedAction

class SopProfile:
    """Validated, hash-stable laboratory review policy."""

    SCHEMA = "factorforge-sop-v1"
    MAX_RULES = 256

    def __init__(self, source: str | Path | Mapping[str, Any], known_rule_ids: Iterable[str] | None = None):
        if isinstance(source, Mapping):
            data = dict(source)
        else:
            with Path(source).open("r", encoding="utf-8") as stream:
                data = json.load(stream)
        self._validate(data, set(known_rule_ids) if known_rule_ids is not None else None)
        self.data = data
        self.profile_id = data["profile_id"]
        self.sop_name = data["sop_name"]
        self.version = data["version"]
        self.status = data["status"]
        self.default_enforcement = EnforcementLevel(data["default_enforcement"].lower())
        self.unknown_rule_policy = data["unknown_rule_policy"].upper()
        
        self.rules: Dict[str, Tuple[EnforcementLevel, AuthorizedAction]] = {}
        for rule_id, rule_data in data["rules"].items():
            if isinstance(rule_data, str):
                # Legacy mapping
                enf = EnforcementLevel(rule_data.lower())
                act = AuthorizedAction.BLOCK if enf == EnforcementLevel.HARD_FAIL else AuthorizedAction.REPORT_ONLY
                self.rules[rule_id] = (enf, act)
            else:
                self.rules[rule_id] = (
                    EnforcementLevel(rule_data["enforcement"].lower()),
                    AuthorizedAction(rule_data["authorized_action"].lower())
                )

        canonical = json.dumps(data, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        self.digest = "sha256:" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    @classmethod
    def from_dict(cls, data: Mapping[str, Any], known_rule_ids: Iterable[str] | None = None):
        return cls(data, known_rule_ids=known_rule_ids)

    @classmethod
    def _validate(cls, data: Mapping[str, Any], known_rule_ids: set[str] | None) -> None:
        required = {
            "$schema", "profile_id", "sop_name", "version", "status",
            "default_enforcement", "unknown_rule_policy", "rules",
        }
        missing = sorted(required - set(data))
        if missing:
            raise ValueError(f"SOP profile missing required fields: {', '.join(missing)}")
        if data["$schema"] != cls.SCHEMA:
            raise ValueError(f"Unsupported SOP schema: {data['$schema']}")
        for field in ("profile_id", "sop_name", "version", "status"):
            if not isinstance(data[field], str) or not data[field].strip():
                raise ValueError(f"SOP field {field} must be a non-empty string")
        if data["unknown_rule_policy"].upper() != "ERROR":
            raise ValueError("SOP unknown_rule_policy must be ERROR")
        try:
            EnforcementLevel(str(data["default_enforcement"]).lower())
        except ValueError as exc:
            raise ValueError("Invalid SOP default_enforcement") from exc
        rules = data["rules"]
        if not isinstance(rules, dict) or len(rules) > cls.MAX_RULES:
            raise ValueError("SOP rules must be an object with at most 256 entries")
        unknown = sorted(set(rules) - known_rule_ids) if known_rule_ids is not None else []
        if unknown:
            raise ValueError(f"Unknown SOP rule IDs: {', '.join(unknown)}")
        for rule_id, rule_data in rules.items():
            if not isinstance(rule_id, str) or not rule_id.strip():
                raise ValueError("SOP rule IDs must be non-empty strings")
            if isinstance(rule_data, str):
                try:
                    EnforcementLevel(rule_data.lower())
                except ValueError as exc:
                    raise ValueError(f"Invalid legacy enforcement for {rule_id}: {rule_data}") from exc
            else:
                if not isinstance(rule_data, dict) or "enforcement" not in rule_data or "authorized_action" not in rule_data:
                    raise ValueError(f"Invalid rule object for {rule_id}: must contain 'enforcement' and 'authorized_action'")
                try:
                    EnforcementLevel(rule_data["enforcement"].lower())
                except ValueError as exc:
                    raise ValueError(f"Invalid enforcement for {rule_id}: {rule_data['enforcement']}") from exc
                try:
                    AuthorizedAction(rule_data["authorized_action"].lower())
                except ValueError as exc:
                    raise ValueError(f"Invalid authorized_action for {rule_id}: {rule_data['authorized_action']}") from exc

    def get_enforcement(self, rule_id: str) -> EnforcementLevel:
        if rule_id in self.rules:
            return self.rules[rule_id][0]
        return self.default_enforcement
        
    def get_authorized_action(self, rule_id: str) -> AuthorizedAction:
        if rule_id in self.rules:
            return self.rules[rule_id][1]
        return AuthorizedAction.BLOCK if self.default_enforcement == EnforcementLevel.HARD_FAIL else AuthorizedAction.REPORT_ONLY

    def provenance(self) -> dict[str, Any]:
        return {
            "profile_id": self.profile_id,
            "profile_name": self.sop_name,
            "profile_version": self.version,
            "profile_status": self.status,
            "profile_hash": self.digest,
            "resolved_rule_policy": {key: {"enforcement": val[0].value, "authorized_action": val[1].value} for key, val in sorted(self.rules.items())},
        }

