"""Optional PostgreSQL persistence for FactorForge design provenance.

The public FactorForge package is usable without a database. This module does
not create an engine or import the integration schema until a caller explicitly
invokes a persistence operation with a configured database URL.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import math
import os
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Mapping
from urllib.parse import urlsplit, urlunsplit


DATABASE_URL_ENV = "FACTORFORGE_DATABASE_URL"
LEGACY_DATABASE_URL_ENV = "DATABASE_URL"


class LocalArtifactStore:
    """Explicit private content store, independent of the relational schema."""

    def __init__(self, root: str | Path):
        self.root = Path(root).resolve()

    def _path(self, digest: str) -> Path:
        if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("Invalid artifact digest")
        path = self.root / f"{digest}.seq"
        if path.is_symlink():
            raise RuntimeError("Artifact symlinks are not supported")
        return path

    def put(self, content: bytes) -> str:
        digest = hashlib.sha256(content).hexdigest()
        self.root.mkdir(parents=True, exist_ok=True)
        path = self._path(digest)
        try:
            with path.open("xb") as stream:
                stream.write(content)
        except FileExistsError:
            pass
        if self.get(digest) != content:
            raise RuntimeError("Artifact content mismatch")
        return f"factorforge-artifact://sha256/{digest}"

    def get(self, digest: str) -> bytes:
        try:
            content = self._path(digest).read_bytes()
        except FileNotFoundError:
            raise RuntimeError("Sequence artifact is unavailable; restore the configured artifact store") from None
        if hashlib.sha256(content).hexdigest() != digest:
            raise RuntimeError("Sequence artifact checksum mismatch")
        return content


def configured_artifact_store() -> LocalArtifactStore:
    root = os.environ.get("FACTORFORGE_ARTIFACT_DIR", "").strip()
    if not root:
        raise PersistenceNotConfiguredError(
            "Integrated persistence requires FACTORFORGE_ARTIFACT_DIR for private sequence artifacts"
        )
    return LocalArtifactStore(root)


class PersistenceNotConfiguredError(RuntimeError):
    """Raised when integrated persistence is requested without configuration."""


class PersistenceDependencyError(RuntimeError):
    """Raised when optional persistence dependencies are unavailable."""


def configured_database_url(environ: Mapping[str, str] | None = None) -> str:
    """Return the explicitly configured PostgreSQL URL."""

    values = os.environ if environ is None else environ
    database_url = (
        values.get(DATABASE_URL_ENV) or values.get(LEGACY_DATABASE_URL_ENV) or ""
    ).strip()
    if not database_url:
        raise PersistenceNotConfiguredError(
            "FactorForge persistence is disabled. Set FACTORFORGE_DATABASE_URL "
            "to an explicitly provisioned PostgreSQL database before using it."
        )

    scheme = urlsplit(database_url).scheme.lower()
    if scheme not in {"postgres", "postgresql", "postgresql+psycopg2"}:
        raise PersistenceNotConfiguredError(
            "FactorForge integrated persistence requires an explicit PostgreSQL URL."
        )
    if scheme in {"postgres", "postgresql"}:
        database_url = "postgresql+psycopg2:" + database_url.split(":", 1)[1]
    return database_url


def redact_database_url(database_url: str) -> str:
    """Return a log-safe database URL with any password removed."""

    parts = urlsplit(database_url)
    if not parts.hostname:
        return f"{parts.scheme}://<redacted>"
    user = f"{parts.username}@" if parts.username else ""
    port = f":{parts.port}" if parts.port else ""
    netloc = f"{user}{parts.hostname}{port}"
    return urlunsplit((parts.scheme, netloc, parts.path, "", ""))


@lru_cache(maxsize=4)
def _session_factory_for_url(database_url: str):
    try:
        sqlalchemy = importlib.import_module("sqlalchemy")
        orm = importlib.import_module("sqlalchemy.orm")
        importlib.import_module("psycopg2")
        importlib.import_module("eijex_db_core.models")
    except ModuleNotFoundError as exc:
        raise PersistenceDependencyError(
            "FactorForge PostgreSQL persistence requires SQLAlchemy, psycopg2, "
            "and an installable eijex-db-core package. Standalone design does not."
        ) from exc

    engine = sqlalchemy.create_engine(database_url, pool_pre_ping=True)
    return orm.sessionmaker(autocommit=False, autoflush=False, bind=engine)


def create_session_factory(database_url: str | None = None):
    """Create a cached session factory for explicitly enabled persistence."""

    explicit_url = (
        configured_database_url({DATABASE_URL_ENV: database_url})
        if database_url is not None
        else configured_database_url()
    )
    return _session_factory_for_url(explicit_url)


def persistence_status(environ: Mapping[str, str] | None = None) -> dict[str, Any]:
    """Return sequence-free configuration metadata for diagnostics."""

    try:
        database_url = configured_database_url(environ)
    except PersistenceNotConfiguredError:
        return {"enabled": False, "backend": None, "database_url": None}
    return {
        "enabled": True,
        "backend": "postgresql",
        "database_url": redact_database_url(database_url),
    }


def _models():
    try:
        return importlib.import_module("eijex_db_core.models")
    except ModuleNotFoundError as exc:
        raise PersistenceDependencyError(
            "Install eijex-db-core before using FactorForge integrated persistence."
        ) from exc


def _check_schema(session) -> None:
    """Require the integration schema before reading or writing domain records."""
    from sqlalchemy import text

    revision = session.execute(text("SELECT version_num FROM alembic_version")).scalar()
    if revision != "b331f0a6d9c1":
        raise RuntimeError("Unsupported persistence schema; apply the pinned eijex-db-core migrations")


def _canonical_sequence(sequence: str) -> str:
    canonical = "".join(sequence.split()).upper()
    if not canonical:
        raise ValueError("Sequence must not be empty")
    return canonical


def _get_or_create_system_user(session, models) -> Any:
    system_user = session.query(models.User).filter_by(display_name="FactorForge AI").first()
    if system_user:
        return system_user.user_id

    organization = (
        session.query(models.Organization).filter_by(display_name="Eijex Systems").first()
    )
    if not organization:
        organization = models.Organization(
            display_name="Eijex Systems",
            organization_type=models.OrganizationType.internal,
        )
        session.add(organization)
        session.flush()

    system_user = models.User(
        organization_id=organization.organization_id,
        display_name="FactorForge AI",
        role=models.UserRole.researcher,
    )
    session.add(system_user)
    session.flush()
    return system_user.user_id


def _get_or_create_sequence(session, models, sequence: str, molecule_class, store) -> Any:
    canonical = _canonical_sequence(sequence)
    sequence_bytes = canonical.encode("utf-8")
    sequence_hash = hashlib.sha256(sequence_bytes).hexdigest()
    artifact_uri = store.put(sequence_bytes)

    existing_sequence = (
        session.query(models.Sequence)
        .filter_by(
            molecule_class=molecule_class,
            canonical_sequence_sha256=sequence_hash,
        )
        .first()
    )
    if existing_sequence:
        return existing_sequence

    artifact = session.query(models.Artifact).filter_by(sha256_hash=sequence_hash).first()
    if not artifact:
        artifact = models.Artifact(
            artifact_uri=artifact_uri,
            sha256_hash=sequence_hash,
            mime_type="text/plain",
            size_bytes=len(sequence_bytes),
            artifact_type=models.ArtifactType.other,
            encryption_status=models.EncryptionStatus.none,
        )
        session.add(artifact)
        session.flush()

    sequence_record = models.Sequence(
        molecule_class=molecule_class,
        canonical_sequence_sha256=sequence_hash,
        artifact_id=artifact.artifact_id,
        length=len(canonical),
    )
    session.add(sequence_record)
    session.flush()
    return sequence_record


def _campaign_sequence_map(session, models, campaign_id) -> dict[str, Any]:
    rows = (
        session.query(models.CampaignSequence, models.Sequence)
        .join(models.Sequence, models.CampaignSequence.sequence_id == models.Sequence.sequence_id)
        .filter(models.CampaignSequence.campaign_id == campaign_id)
        .all()
    )
    return {association.sequence_role: sequence for association, sequence in rows}


def _assert_idempotent_campaign(
    session,
    models,
    campaign,
    input_sequence: str,
    optimized_sequence: str,
    description: str,
) -> None:
    sequences = _campaign_sequence_map(session, models, campaign.campaign_id)
    expected = {
        "input": hashlib.sha256(_canonical_sequence(input_sequence).encode("utf-8")).hexdigest(),
        "optimized": hashlib.sha256(
            _canonical_sequence(optimized_sequence).encode("utf-8")
        ).hexdigest(),
    }
    observed = {
        role: sequence.canonical_sequence_sha256 for role, sequence in sequences.items()
    }
    if observed != expected or campaign.description != description:
        raise ValueError(
            "Study number already exists with different sequence identities or provenance"
        )


def save_optimization(
    study_number: str,
    protein_name: str,
    input_sequence: str,
    optimized_sequence: str,
    metrics: Mapping[str, Any],
    algorithm_version: str = "2.1.0",
    *,
    session_factory: Callable[[], Any] | None = None,
    artifact_store: LocalArtifactStore | None = None,
) -> str:
    """Persist one computational design in the shared PostgreSQL schema."""

    if not study_number.strip():
        raise ValueError("study_number must not be empty")

    factory = session_factory or create_session_factory()
    store = artifact_store or configured_artifact_store()
    persisted_metrics = {
        name: float(metrics[name])
        for name in ("gc_content", "cai", "tm", "execution_time")
        if metrics.get(name) is not None
    }
    if not all(math.isfinite(value) for value in persisted_metrics.values()):
        raise ValueError("Computational metrics must be finite")
    description = json.dumps(
        {"protein": protein_name, "algorithm_version": algorithm_version,
         "metrics": persisted_metrics}, sort_keys=True, allow_nan=False,
    )
    models = _models()
    with factory() as session, session.begin():
        _check_schema(session)
        from sqlalchemy import text
        # Serialize adapter writes, including system-user and sequence creation.
        # The transaction releases this lock automatically on commit/rollback.
        session.execute(text("SELECT pg_advisory_xact_lock(331, 1)"))
        existing_campaign = (
            session.query(models.Campaign).filter(models.Campaign.name == study_number).first()
        )
        if existing_campaign:
            _assert_idempotent_campaign(
                session,
                models,
                existing_campaign,
                input_sequence,
                optimized_sequence,
                description,
            )
            store.get(hashlib.sha256(_canonical_sequence(input_sequence).encode()).hexdigest())
            store.get(hashlib.sha256(_canonical_sequence(optimized_sequence).encode()).hexdigest())
            return str(existing_campaign.campaign_id)

        system_user_id = _get_or_create_system_user(session, models)
        campaign = models.Campaign(
            name=study_number,
            description=description,
        )
        session.add(campaign)
        session.flush()

        input_record = _get_or_create_sequence(
            session, models, input_sequence, models.MoleculeClass.Protein, store
        )
        optimized_record = _get_or_create_sequence(
            session, models, optimized_sequence, models.MoleculeClass.CDS, store
        )
        session.add_all(
            [
                models.CampaignSequence(
                    campaign_id=campaign.campaign_id,
                    sequence_id=input_record.sequence_id,
                    sequence_role="input",
                ),
                models.CampaignSequence(
                    campaign_id=campaign.campaign_id,
                    sequence_id=optimized_record.sequence_id,
                    sequence_role="optimized",
                ),
            ]
        )

        candidate = models.Candidate(
            campaign_id=campaign.campaign_id,
            sequence_id=optimized_record.sequence_id,
            designer_user_id=system_user_id,
            design_rationale=f"Algorithm {algorithm_version} optimization",
        )
        session.add(candidate)
        session.flush()

        for metric_name, metric_value in persisted_metrics.items():
            if metric_value is not None:
                session.add(
                    models.CandidateMetric(
                        candidate_id=candidate.candidate_id,
                        metric_name=metric_name,
                        metric_value=float(metric_value),
                    )
                )

        return str(campaign.campaign_id)


def get_batch(
    study_number: str,
    *,
    session_factory: Callable[[], Any] | None = None,
    artifact_store: LocalArtifactStore | None = None,
) -> dict[str, Any] | None:
    """Retrieve a private campaign, preserving the legacy sequence-preview contract."""

    factory = session_factory or create_session_factory()
    models = _models()
    with factory() as session:
        _check_schema(session)
        campaign = (
            session.query(models.Campaign).filter(models.Campaign.name == study_number).first()
        )
        if not campaign:
            return None

        sequences = _campaign_sequence_map(session, models, campaign.campaign_id)
        candidate = (
            session.query(models.Candidate)
            .filter(models.Candidate.campaign_id == campaign.campaign_id)
            .first()
        )
        metrics: dict[str, float] = {}
        if candidate:
            metric_rows = (
                session.query(models.CandidateMetric)
                .filter_by(candidate_id=candidate.candidate_id)
                .all()
            )
            metrics = {row.metric_name: row.metric_value for row in metric_rows}

        store = artifact_store or configured_artifact_store()
        sequence_payloads = []
        for role in ("input", "optimized"):
            sequence = sequences.get(role)
            if sequence is None:
                continue
            sequence_payloads.append(
                {
                    "type": role,
                    "data": store.get(sequence.canonical_sequence_sha256).decode("utf-8")[:50] + "...",
                    "sequence_id": str(sequence.sequence_id),
                    "sha256": sequence.canonical_sequence_sha256,
                    "length": sequence.length,
                    "gc": metrics.get("gc_content") if role == "optimized" else None,
                    "cai": metrics.get("cai") if role == "optimized" else None,
                }
            )

        protein_name = ""
        if campaign.description:
            protein_name = json.loads(campaign.description)["protein"]

        return {
            "batch_id": str(campaign.campaign_id),
            "study_number": campaign.name,
            "protein": protein_name,
            "status": "completed",
            "sequences": sequence_payloads,
        }
