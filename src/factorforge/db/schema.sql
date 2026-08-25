-- FactorForge / ValidationHub / AgentOS Tri-Tier Database Schema DDL
-- Research implementation checkpoint for the v3.5.0 ML track.
--
-- This file is a version-controlled reference schema, not a migration ledger
-- or proof that a deployed PostgreSQL instance is synchronized with Git.
-- Production/release claims require captured migrations and cross-backend tests.

CREATE SCHEMA IF NOT EXISTS common;
CREATE SCHEMA IF NOT EXISTS factorforge;
CREATE SCHEMA IF NOT EXISTS validationhub;
CREATE SCHEMA IF NOT EXISTS agentops;

-- ============================================================================
-- 1. COMMON SCHEMA: Biological Context & Central Sequence Registry
-- ============================================================================
CREATE TABLE IF NOT EXISTS common.biological_context (
    context_id SERIAL PRIMARY KEY,
    context_code VARCHAR(64) UNIQUE NOT NULL, -- e.g., 'NB_TRANSIENT_NUCLEAR'
    organism VARCHAR(128) NOT NULL,            -- e.g., 'Nicotiana benthamiana'
    taxonomy_id INT NOT NULL,                  -- e.g., 4100
    host_line VARCHAR(64),                     -- e.g., 'LAB-v1.1'
    expression_system VARCHAR(64),            -- e.g., 'transient_agrobacterium'
    cellular_compartment VARCHAR(64),         -- e.g., 'nuclear'
    molecule_type VARCHAR(32) DEFAULT 'CDS',
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS common.sequence_registry (
    sequence_id SERIAL PRIMARY KEY,
    sequence_hash CHAR(64) UNIQUE NOT NULL,   -- SHA-256 Content Address
    molecule_type VARCHAR(32) NOT NULL,        -- 'DNA', 'CDS', 'RNA', 'Protein'
    sequence_length INT NOT NULL,
    raw_sequence TEXT NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================================
-- 2. FACTORFORGE SCHEMA: Context-Aware Rule Engine & Execution Provenance
-- ============================================================================

-- 2.1 Restriction Enzyme Registry
CREATE TABLE IF NOT EXISTS factorforge.enzyme_registry (
    enzyme_id SERIAL PRIMARY KEY,
    enzyme_name VARCHAR(64) UNIQUE NOT NULL,  -- e.g., 'BsaI', 'BpiI', 'BsmBI', 'SapI'
    recognition_sequence VARCHAR(64) NOT NULL, -- e.g., 'GGTCTC'
    cut_offset INT,
    category VARCHAR(32) DEFAULT 'Type_IIS',   -- 'Type_IIS', 'Standard'
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- 2.2 Constraint Definitions (Hard vs Soft Rules)
CREATE TABLE IF NOT EXISTS factorforge.constraint_definitions (
    constraint_id SERIAL PRIMARY KEY,
    code VARCHAR(64) UNIQUE NOT NULL,          -- e.g., 'AA_IDENTITY', 'TYPE_IIS_BSAI', 'GLOBAL_GC'
    name VARCHAR(128) NOT NULL,
    rule_type VARCHAR(16) NOT NULL,            -- 'HARD', 'SOFT'
    description TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- 2.3 Constraint Profiles (Versioned Policy Bundles per Biological Context)
CREATE TABLE IF NOT EXISTS factorforge.constraint_profiles (
    profile_id SERIAL PRIMARY KEY,
    profile_code VARCHAR(64) UNIQUE NOT NULL,  -- e.g., 'FF-NB-v3.5'
    context_id INT REFERENCES common.biological_context(context_id),
    version VARCHAR(32) NOT NULL,
    description TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Profile items mapping constraints to parameters
CREATE TABLE IF NOT EXISTS factorforge.constraint_profile_items (
    item_id SERIAL PRIMARY KEY,
    profile_id INT REFERENCES factorforge.constraint_profiles(profile_id),
    constraint_id INT REFERENCES factorforge.constraint_definitions(constraint_id),
    parameters_json JSONB,                     -- e.g., {"target_min": 0.35, "target_max": 0.65}
    action_on_fail VARCHAR(32) DEFAULT 'REJECT'
);

-- 2.4 Research provenance registries
CREATE TABLE IF NOT EXISTS factorforge.dataset_snapshots (
    snapshot_id SERIAL PRIMARY KEY,
    snapshot_name VARCHAR(128) NOT NULL,
    version VARCHAR(64) NOT NULL,
    status VARCHAR(64) NOT NULL,
    host_scope VARCHAR(128),
    sequence_count INT,
    deduplication_method VARCHAR(64),
    split_ratio JSONB,
    manifest_hash CHAR(64) NOT NULL,
    manifest_uri TEXT,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (snapshot_name, version, manifest_hash)
);

CREATE TABLE IF NOT EXISTS factorforge.model_registry (
    model_id SERIAL PRIMARY KEY,
    model_name VARCHAR(128) NOT NULL,
    model_version VARCHAR(64) NOT NULL,
    architecture VARCHAR(128) NOT NULL,
    tokenizer_version VARCHAR(64),
    weights_hash CHAR(64),
    config_hash CHAR(64),
    code_commit CHAR(40),
    dataset_snapshot_id INT REFERENCES factorforge.dataset_snapshots(snapshot_id)
        ON DELETE RESTRICT,
    status VARCHAR(32) NOT NULL DEFAULT 'research_scaffold',
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (model_name, model_version)
);

CREATE TABLE IF NOT EXISTS factorforge.ml_features (
    feature_id SERIAL PRIMARY KEY,
    sequence_id INT NOT NULL REFERENCES common.sequence_registry(sequence_id)
        ON DELETE RESTRICT,
    feature_name VARCHAR(128) NOT NULL,
    feature_version VARCHAR(64) NOT NULL,
    feature_value JSONB NOT NULL,
    observed_at TIMESTAMPTZ NOT NULL,
    available_at TIMESTAMPTZ NOT NULL,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    CHECK (available_at >= observed_at)
);

-- The timestamps above make feature availability explicit. Preventing leakage
-- at a prediction cutoff remains a query/pipeline responsibility and requires
-- a separate negative integration test.

-- 2.5 Design Execution Lineage: Request -> Run -> Candidate -> Evaluation -> Package
CREATE TABLE IF NOT EXISTS factorforge.design_requests (
    request_id SERIAL PRIMARY KEY,
    context_id INT REFERENCES common.biological_context(context_id),
    profile_id INT REFERENCES factorforge.constraint_profiles(profile_id),
    target_protein_sequence TEXT NOT NULL,
    requested_by VARCHAR(128),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS factorforge.design_runs (
    run_id SERIAL PRIMARY KEY,
    request_id INT REFERENCES factorforge.design_requests(request_id),
    requested_engine VARCHAR(64) NOT NULL,     -- e.g., 'lm'
    actual_engine VARCHAR(64) NOT NULL,        -- e.g., 'lm' or 'dp_fallback'
    fallback_reason TEXT,
    model_version VARCHAR(64) NOT NULL,        -- e.g., 'v3.5.0-SynCodonLM-V2'
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS factorforge.candidates (
    candidate_id SERIAL PRIMARY KEY,
    run_id INT REFERENCES factorforge.design_runs(run_id),
    sequence_id INT REFERENCES common.sequence_registry(sequence_id),
    overall_status VARCHAR(32) NOT NULL,       -- 'PASS', 'REJECT'
    cai_value NUMERIC(5, 4),
    cai_reference_id VARCHAR(64),              -- e.g., 'NB_codon_reference_v2'
    gc_percent NUMERIC(5, 2),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Step-by-step Constraint Evaluations per Candidate
CREATE TABLE IF NOT EXISTS factorforge.constraint_evaluations (
    evaluation_id SERIAL PRIMARY KEY,
    candidate_id INT REFERENCES factorforge.candidates(candidate_id),
    constraint_id INT REFERENCES factorforge.constraint_definitions(constraint_id),
    status VARCHAR(16) NOT NULL,               -- 'PASS', 'FAIL'
    observed_value NUMERIC(10, 4),
    details_json JSONB,                        -- e.g., {"detected_sites": 0, "window": "global"}
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Frozen Release Package
CREATE TABLE IF NOT EXISTS factorforge.design_packages (
    package_id SERIAL PRIMARY KEY,
    candidate_id INT REFERENCES factorforge.candidates(candidate_id),
    construct_id VARCHAR(128) NOT NULL,
    frozen_manifest_hash CHAR(64) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================================
-- 3. VALIDATIONHUB SCHEMA: Domain-Isolated Evidence Records
-- ============================================================================
CREATE TABLE IF NOT EXISTS validationhub.computational_validation_records (
    record_id SERIAL PRIMARY KEY,
    sequence_id INT REFERENCES common.sequence_registry(sequence_id),
    validator_passed BOOLEAN NOT NULL,
    aa_identity NUMERIC(5, 4) NOT NULL,
    validation_manifest JSONB,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS validationhub.observation_records (
    observation_id SERIAL PRIMARY KEY,
    construct_id VARCHAR(128) NOT NULL,
    batch_id VARCHAR(64) NOT NULL,
    measurement_type VARCHAR(64) NOT NULL,     -- e.g., 'protein_yield_g_kg'
    numeric_value NUMERIC(12, 4),
    unit VARCHAR(32),
    observed_at TIMESTAMP NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================================
-- 4. AGENTOPS SCHEMA: Policy Enforcement & Workflow Audit
-- ============================================================================
CREATE TABLE IF NOT EXISTS agentops.workflow_runs (
    workflow_id SERIAL PRIMARY KEY,
    run_code VARCHAR(64) UNIQUE NOT NULL,
    requested_engine VARCHAR(64) NOT NULL,
    actual_engine VARCHAR(64) NOT NULL,
    approved_by VARCHAR(128),
    policy_version VARCHAR(32) DEFAULT 'v1.0',
    audit_hash CHAR(64) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
