# Optional persistence

FactorForge is a standalone CDS design package by default. Importing the package,
running the CLI, generating a design, and exporting its artifacts do not require a
database.

## Integrated PostgreSQL mode

Eijex deployments may explicitly enable shared campaign and candidate persistence.
That mode uses PostgreSQL and the models and migrations owned by `eijex-db-core`.
It is an integration capability, not a requirement for public FactorForge use.

Install the PostgreSQL dependencies:

```bash
python -m pip install "factorforge-cds[postgres]"
```

The deployment must separately install the pinned `eijex-db-core==0.1.1` wheel
with schema revision `b331f0a6d9c1`. The
package is not currently distributed as part of `factorforge-cds`, and a local
editable monorepo path must not be assumed by a public FactorForge installation.

After installing that wheel, migrations run from any working directory:

```bash
python -m eijex_db_core.migrate upgrade head
```

Set `EIJEX_DATABASE_URL` explicitly for migrations. Database provisioning and
migration approval belong to the deployment operator, not the public web service.

Configure the database explicitly:

```bash
export FACTORFORGE_DATABASE_URL='postgresql://USER:PASSWORD@HOST:5432/DATABASE'
export FACTORFORGE_ARTIFACT_DIR='/private/factorforge-artifacts'
```

`DATABASE_URL` remains a compatibility alias. New deployments should use the
FactorForge-specific variable. If neither is present, persistence operations fail
with an actionable configuration error; standalone design remains available.

## Data boundary

The shared schema stores:

- canonical sequence UUID, SHA-256 identity, molecule class, and length;
- role-qualified input and optimized sequence links for a campaign;
- generated candidate identity and deterministic/advisory computational metrics;
- sequence-free campaign and algorithm descriptions.

The explicit private artifact directory stores canonical sequence bytes, verified by
SHA-256 on every read. It must not be web-served or committed; operators own access
control, encryption, backups, and retention. Database transactions may leave
unreferenced content-addressed files after rollback; these are not successful designs.
Private `get_batch` retrieval preserves the legacy 50-character sequence preview
and adds canonical UUID/hash metadata. It must not be exposed as unauthenticated
public intake. FactorForge does not write physical wet-lab
measurements; those belong to ValidationHub experimental records.

## Local SQLite checkpoint

`factorforge.db.connector.FactorForgeDBConnector(db_path=...)` is an explicit local
research-checkpoint utility. It is not a replica, fallback, or system of record for
the shared PostgreSQL DBTL schemas. Configuring PostgreSQL never silently falls back
to SQLite.

## Maintainer verification

Before enabling the integration in a release or deployment:

1. Build and install the FactorForge wheel in a clean environment without database
   extras; import and one representative design must succeed.
2. Build and install a pinned `eijex-db-core` wheel plus FactorForge's `postgres`
   extra in the integration environment.
3. Set an isolated test URL, apply the `eijex-db-core` Alembic migrations, and run the
   FactorForge database round-trip test.
4. Verify errors and diagnostic output never contain the database password.

Integrated persistence remains an experimental, separately provisioned capability;
this adapter does not migrate old `factorforge.batches` records or deploy a shared
production database. Missing artifacts and unsupported schema revisions fail closed.
