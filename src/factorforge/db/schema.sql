-- FactorForge / ValidationHub / AgentOS Tri-Tier Database Schema DDL
-- v0.3 Canonical Relational Spec (GENERATED FROM SPEC)
-- 
-- Job 234 Master Architecture -> Job 210 DB Spec -> schema.sql v0.3

CREATE SCHEMA IF NOT EXISTS common;
CREATE SCHEMA IF NOT EXISTS factorforge;
CREATE SCHEMA IF NOT EXISTS validationhub;
CREATE SCHEMA IF NOT EXISTS agentops;
CREATE SCHEMA IF NOT EXISTS ml;
CREATE SCHEMA IF NOT EXISTS bioprocess;

-- ============================================================================
-- 1. COMMON SCHEMA: Canonical Sequence Registry & Biological Context
-- ============================================================================

CREATE TABLE IF NOT EXISTS common.artifacts (
    artifact_id UUID PRIMARY KEY,
    artifact_uri TEXT NOT NULL,
    sha256_hash CHAR(64) NOT NULL,
    mime_type VARCHAR(64),
    size_bytes BIGINT,
    encryption_status VARCHAR(32) DEFAULT 'none',
    artifact_type VARCHAR(32) NOT NULL,
    supersedes_artifact_id UUID REFERENCES common.artifacts(artifact_id),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (artifact_uri, sha256_hash)
);

CREATE TABLE IF NOT EXISTS common.reference_assemblies (
    reference_assembly_id UUID PRIMARY KEY,
    organism VARCHAR(128) NOT NULL,
    taxonomy_id INT NOT NULL,
    assembly_name VARCHAR(128) NOT NULL,
    version VARCHAR(64) NOT NULL,
    sequence_artifact_id UUID REFERENCES common.artifacts(artifact_id),
    annotation_artifact_id UUID REFERENCES common.artifacts(artifact_id),
    checksum CHAR(64),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS common.biological_contexts (
    biological_context_id UUID PRIMARY KEY,
    taxonomy_id INT NOT NULL,
    organism VARCHAR(128) NOT NULL,
    host_line VARCHAR(64),
    strain VARCHAR(64),
    cultivar VARCHAR(64),
    reference_assembly_id UUID REFERENCES common.reference_assemblies(reference_assembly_id),
    annotation_release VARCHAR(64),
    expression_system VARCHAR(64),
    cellular_compartment VARCHAR(64),
    molecule_context VARCHAR(32) DEFAULT 'CDS',
    profile_version VARCHAR(32),
    production_enabled BOOLEAN DEFAULT false,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS common.sequences (
    sequence_id UUID PRIMARY KEY,
    molecule_class VARCHAR(32) NOT NULL,
    canonical_sequence_sha256 CHAR(64) NOT NULL,
    artifact_id UUID NOT NULL REFERENCES common.artifacts(artifact_id),
    length INT NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS common.sequence_operations (
    operation_id UUID PRIMARY KEY,
    operation_type VARCHAR(64) NOT NULL,
    actor_type VARCHAR(64),
    engine_reference VARCHAR(128),
    workflow_run_id UUID,
    artifact_id UUID REFERENCES common.artifacts(artifact_id),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS common.sequence_lineage (
    lineage_id UUID PRIMARY KEY,
    child_sequence_id UUID NOT NULL REFERENCES common.sequences(sequence_id),
    parent_sequence_id UUID NOT NULL REFERENCES common.sequences(sequence_id),
    operation_id UUID REFERENCES common.sequence_operations(operation_id),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================================
-- 2. FACTORFORGE SCHEMA: Context-Aware Rule Engine & Execution Provenance
-- ============================================================================

CREATE TABLE IF NOT EXISTS factorforge.protein_targets (
    protein_target_id UUID PRIMARY KEY,
    display_alias VARCHAR(128) NOT NULL,
    sequence_id UUID REFERENCES common.sequences(sequence_id),
    source VARCHAR(64) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS factorforge.codon_references (
    codon_reference_id UUID PRIMARY KEY,
    reference_name VARCHAR(128) NOT NULL,
    tier VARCHAR(64) NOT NULL,
    checksum_sha256 CHAR(64) NOT NULL,
    biological_context_id UUID REFERENCES common.biological_contexts(biological_context_id),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS factorforge.enzyme_registry (
    enzyme_id SERIAL PRIMARY KEY,
    enzyme_name VARCHAR(64) UNIQUE NOT NULL,
    recognition_sequence VARCHAR(64) NOT NULL,
    cut_offset INT,
    category VARCHAR(32) DEFAULT 'Type_IIS',
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS factorforge.constraint_definitions (
    constraint_id SERIAL PRIMARY KEY,
    code VARCHAR(64) UNIQUE NOT NULL,
    name VARCHAR(128) NOT NULL,
    rule_type VARCHAR(16) NOT NULL,
    description TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS factorforge.constraint_profiles (
    profile_id SERIAL PRIMARY KEY,
    profile_code VARCHAR(64) UNIQUE NOT NULL,
    biological_context_id UUID REFERENCES common.biological_contexts(biological_context_id),
    version VARCHAR(32) NOT NULL,
    description TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS factorforge.constraint_profile_items (
    item_id SERIAL PRIMARY KEY,
    profile_id INT REFERENCES factorforge.constraint_profiles(profile_id),
    constraint_id INT REFERENCES factorforge.constraint_definitions(constraint_id),
    parameters_json JSONB,
    action_on_fail VARCHAR(32) DEFAULT 'REJECT'
);

CREATE TABLE IF NOT EXISTS factorforge.design_specs (
    design_spec_id UUID PRIMARY KEY,
    protein_target_id UUID REFERENCES factorforge.protein_targets(protein_target_id),
    biological_context_id UUID REFERENCES common.biological_contexts(biological_context_id),
    requested_engine_name VARCHAR(64) NOT NULL,
    requested_profile_name VARCHAR(64) NOT NULL,
    requested_engine_objective VARCHAR(64),
    constraints_json JSONB,
    approved_by_approval_action_id UUID,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS factorforge.design_runs (
    design_run_id UUID PRIMARY KEY,
    design_spec_id UUID REFERENCES factorforge.design_specs(design_spec_id),
    execution_origin VARCHAR(64) NOT NULL,
    workflow_run_id UUID,
    task_run_id UUID,
    tool_call_id UUID,
    actual_engine_name VARCHAR(64) NOT NULL,
    actual_profile_name VARCHAR(64),
    actual_engine_objective VARCHAR(64),
    generation_performed BOOLEAN NOT NULL,
    analysis_mode VARCHAR(64) NOT NULL,
    evidence_level VARCHAR(64) NOT NULL,
    factorforge_commit CHAR(40),
    agentops_commit CHAR(40),
    runner_entrypoint VARCHAR(128),
    input_sequence_id UUID REFERENCES common.sequences(sequence_id),
    record_version INT NOT NULL DEFAULT 1,
    supersedes_design_run_id UUID REFERENCES factorforge.design_runs(design_run_id),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS factorforge.candidates (
    candidate_id UUID PRIMARY KEY,
    design_run_id UUID REFERENCES factorforge.design_runs(design_run_id),
    sequence_id UUID REFERENCES common.sequences(sequence_id),
    computational_status VARCHAR(32) NOT NULL,
    cai_value NUMERIC(5, 4),
    gc_percent NUMERIC(5, 2),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS factorforge.check_results (
    check_result_id UUID PRIMARY KEY,
    candidate_id UUID REFERENCES factorforge.candidates(candidate_id),
    check_type VARCHAR(64) NOT NULL,
    result VARCHAR(16) NOT NULL,
    observed_value NUMERIC(10, 4),
    details_json JSONB,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS factorforge.design_packages (
    package_id UUID PRIMARY KEY,
    design_run_id UUID REFERENCES factorforge.design_runs(design_run_id),
    frozen_manifest_hash CHAR(64) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS factorforge.design_package_candidates (
    package_candidate_id UUID PRIMARY KEY,
    package_id UUID REFERENCES factorforge.design_packages(package_id),
    candidate_id UUID REFERENCES factorforge.candidates(candidate_id),
    UNIQUE (package_id, candidate_id)
);

-- ============================================================================
-- 3. VALIDATIONHUB SCHEMA: Domain-Isolated Evidence Records
-- ============================================================================

CREATE TABLE IF NOT EXISTS validationhub.constructs (
    construct_id UUID PRIMARY KEY,
    candidate_id UUID REFERENCES factorforge.candidates(candidate_id),
    sequence_id UUID REFERENCES common.sequences(sequence_id),
    construct_alias VARCHAR(128),
    identity_status VARCHAR(32) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS validationhub.experiments (
    experiment_id UUID PRIMARY KEY,
    batch_id VARCHAR(64) NOT NULL,
    biological_context_id UUID REFERENCES common.biological_contexts(biological_context_id),
    performed_at TIMESTAMP NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS validationhub.experimental_arms (
    experimental_arm_id UUID PRIMARY KEY,
    experiment_id UUID REFERENCES validationhub.experiments(experiment_id),
    construct_id UUID REFERENCES validationhub.constructs(construct_id),
    arm_type VARCHAR(64) NOT NULL,
    control_role VARCHAR(64),
    linkage_status VARCHAR(32) NOT NULL,
    outcome_link_allowed BOOLEAN DEFAULT false,
    expected_result_direction VARCHAR(64),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS validationhub.samples (
    sample_id UUID PRIMARY KEY,
    experimental_arm_id UUID REFERENCES validationhub.experimental_arms(experimental_arm_id),
    sample_label VARCHAR(128),
    collected_at TIMESTAMP NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS validationhub.measurements (
    measurement_id UUID PRIMARY KEY,
    sample_id UUID REFERENCES validationhub.samples(sample_id),
    assay_type VARCHAR(64) NOT NULL,
    analyte VARCHAR(64),
    value_numeric NUMERIC(12, 4),
    value_text VARCHAR(256),
    unit VARCHAR(32),
    raw_unit_text VARCHAR(64),
    measured_at TIMESTAMP NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS validationhub.control_results (
    control_result_id UUID PRIMARY KEY,
    experimental_arm_id UUID REFERENCES validationhub.experimental_arms(experimental_arm_id),
    result VARCHAR(32) NOT NULL,
    interpretation_json JSONB,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS validationhub.outcomes (
    outcome_id UUID PRIMARY KEY,
    experimental_arm_id UUID REFERENCES validationhub.experimental_arms(experimental_arm_id),
    outcome_category VARCHAR(32) NOT NULL,
    outcome_result VARCHAR(32) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS validationhub.evidence_records (
    evidence_record_id UUID PRIMARY KEY,
    record_hash CHAR(64) NOT NULL,
    record_status VARCHAR(32) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS validationhub.evidence_reviews (
    evidence_review_id UUID PRIMARY KEY,
    evidence_record_id UUID REFERENCES validationhub.evidence_records(evidence_record_id),
    reviewed_record_hash CHAR(64) NOT NULL,
    decision VARCHAR(32) NOT NULL,
    reviewed_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS validationhub.claims (
    claim_id UUID PRIMARY KEY,
    construct_id UUID REFERENCES validationhub.constructs(construct_id),
    candidate_id UUID REFERENCES factorforge.candidates(candidate_id),
    claim_type VARCHAR(64) NOT NULL,
    evidence_boundary VARCHAR(64) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS validationhub.claim_evidence_links (
    link_id UUID PRIMARY KEY,
    claim_id UUID REFERENCES validationhub.claims(claim_id),
    relationship VARCHAR(32) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS validationhub.evidence_packages (
    evidence_package_id UUID PRIMARY KEY,
    snapshot_hash CHAR(64) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================================
-- 4. AGENTOPS SCHEMA: Policy Enforcement & Workflow Audit
-- ============================================================================

CREATE TABLE IF NOT EXISTS agentops.workflow_runs (
    workflow_run_id UUID PRIMARY KEY,
    run_code VARCHAR(64) UNIQUE NOT NULL,
    requested_engine VARCHAR(64) NOT NULL,
    actual_engine VARCHAR(64) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS agentops.task_runs (
    task_run_id UUID PRIMARY KEY,
    workflow_run_id UUID REFERENCES agentops.workflow_runs(workflow_run_id),
    task_name VARCHAR(128) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS agentops.agent_invocations (
    invocation_id UUID PRIMARY KEY,
    task_run_id UUID REFERENCES agentops.task_runs(task_run_id),
    agent_name VARCHAR(64) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS agentops.tool_calls (
    tool_call_id UUID PRIMARY KEY,
    task_run_id UUID REFERENCES agentops.task_runs(task_run_id),
    tool_name VARCHAR(64) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS agentops.policy_decisions (
    decision_id UUID PRIMARY KEY,
    tool_call_id UUID REFERENCES agentops.tool_calls(tool_call_id),
    decision VARCHAR(32) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS agentops.approval_requests (
    approval_request_id UUID PRIMARY KEY,
    workflow_run_id UUID REFERENCES agentops.workflow_runs(workflow_run_id),
    target_object_type VARCHAR(64) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS agentops.approval_actions (
    approval_action_id UUID PRIMARY KEY,
    approval_request_id UUID REFERENCES agentops.approval_requests(approval_request_id),
    decision VARCHAR(32) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS agentops.audit_streams (
    stream_id UUID PRIMARY KEY,
    stream_type VARCHAR(64) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS agentops.audit_events (
    event_id UUID PRIMARY KEY,
    stream_id UUID REFERENCES agentops.audit_streams(stream_id),
    previous_event_hash CHAR(64),
    payload_hash CHAR(64) NOT NULL,
    event_hash CHAR(64) NOT NULL,
    sequence_number INT NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================================
-- 4.1 AUDIT IMMUTABILITY TRIGGERS
-- ============================================================================
-- Enforce cryptographic immutability on agentops.audit_events.
-- Updates and deletes are strictly prohibited at the DB engine layer.

CREATE OR REPLACE FUNCTION agentops.prevent_audit_tampering()
RETURNS TRIGGER AS $$
BEGIN
    RAISE EXCEPTION 'Audit events are immutable. UPDATE or DELETE operations are strictly prohibited on agentops.audit_events.';
END;
$$ LANGUAGE plpgsql;

DROP TRIGGER IF EXISTS enforce_audit_immutability_update ON agentops.audit_events;
CREATE TRIGGER enforce_audit_immutability_update
BEFORE UPDATE ON agentops.audit_events
FOR EACH ROW EXECUTE FUNCTION agentops.prevent_audit_tampering();

DROP TRIGGER IF EXISTS enforce_audit_immutability_delete ON agentops.audit_events;
CREATE TRIGGER enforce_audit_immutability_delete
BEFORE DELETE ON agentops.audit_events
FOR EACH ROW EXECUTE FUNCTION agentops.prevent_audit_tampering();

-- ============================================================================
-- 5. FUTURE EXTENSION CONTRACTS (Research Scaffold)
-- ============================================================================
-- NOTE: The ML and BioProcess schemas are currently just research scaffolds
-- and are not part of the core relational spec yet. They map to future contracts.

CREATE TABLE IF NOT EXISTS ml.dataset_snapshots (
    snapshot_id UUID PRIMARY KEY,
    snapshot_name VARCHAR(128) NOT NULL,
    version VARCHAR(64) NOT NULL,
    status VARCHAR(64) NOT NULL,
    manifest_hash CHAR(64) NOT NULL,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS ml.dataset_snapshot_items (
    item_id UUID PRIMARY KEY,
    snapshot_id UUID REFERENCES ml.dataset_snapshots(snapshot_id),
    sequence_id UUID REFERENCES common.sequences(sequence_id),
    split_role VARCHAR(32) NOT NULL, -- e.g., 'train', 'validation', 'test', 'reference'
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (snapshot_id, sequence_id)
);

CREATE TABLE IF NOT EXISTS ml.training_runs (
    training_run_id UUID PRIMARY KEY,
    dataset_snapshot_id UUID REFERENCES ml.dataset_snapshots(snapshot_id),
    training_config JSONB,
    status VARCHAR(32) NOT NULL,
    started_at TIMESTAMPTZ,
    completed_at TIMESTAMPTZ,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS ml.model_registry (
    model_id UUID PRIMARY KEY,
    model_name VARCHAR(128) NOT NULL,
    model_version VARCHAR(64) NOT NULL,
    training_run_id UUID REFERENCES ml.training_runs(training_run_id),
    status VARCHAR(32) NOT NULL DEFAULT 'research_scaffold',
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS ml.evaluation_runs (
    evaluation_run_id UUID PRIMARY KEY,
    model_id UUID REFERENCES ml.model_registry(model_id),
    evaluation_snapshot_id UUID REFERENCES ml.dataset_snapshots(snapshot_id),
    engine_name VARCHAR(64) NOT NULL,
    evaluation_protocol VARCHAR(128) NOT NULL,
    status VARCHAR(32) NOT NULL,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS ml.evaluation_checks (
    check_id UUID PRIMARY KEY,
    evaluation_run_id UUID REFERENCES ml.evaluation_runs(evaluation_run_id),
    check_type VARCHAR(64) NOT NULL, -- e.g., train_test_exact_overlap, generated_output_homology
    result VARCHAR(32) NOT NULL,     -- PASS, FAIL, WARNING, INDETERMINATE, NOT_APPLICABLE
    threshold FLOAT,
    match_count INT,
    details_json JSONB,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS ml.feature_observations (
    feature_id UUID PRIMARY KEY,
    sequence_id UUID NOT NULL REFERENCES common.sequences(sequence_id),
    feature_name VARCHAR(128) NOT NULL,
    feature_value JSONB NOT NULL,
    observed_at TIMESTAMPTZ NOT NULL,
    available_at TIMESTAMPTZ NOT NULL,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    CHECK (available_at >= observed_at)
);
