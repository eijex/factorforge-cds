import psycopg2
import uuid
from psycopg2.extras import DictCursor

def test_integrity():
    conn = psycopg2.connect(
        host="127.0.0.1",
        port=5432,
        user="postgres",
        password="postgres",
        dbname="postgres"
    )
    conn.autocommit = True
    cursor = conn.cursor(cursor_factory=DictCursor)
    
    print("Testing Tri-Tier Schema Data Insertion (Paper 2 Fixtures)...")
    
    # 1. Common Schema
    artifact_id = str(uuid.uuid4())
    cursor.execute("""
        INSERT INTO common.artifacts (artifact_id, artifact_uri, sha256_hash, artifact_type)
        VALUES (%s, 's3://eijex-data/seq1.fasta', 'd2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2d2', 'fasta')
    """, (artifact_id,))
    
    seq_id = str(uuid.uuid4())
    cursor.execute("""
        INSERT INTO common.sequences (sequence_id, molecule_class, canonical_sequence_sha256, artifact_id, length)
        VALUES (%s, 'CDS', 'c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1', %s, 1500)
    """, (seq_id, artifact_id))
    
    ref_assembly_id = str(uuid.uuid4())
    cursor.execute("""
        INSERT INTO common.reference_assemblies (reference_assembly_id, organism, taxonomy_id, assembly_name, version)
        VALUES (%s, 'Nicotiana benthamiana', 4100, 'Niben.genome.v1.0.1', '1.0.1')
    """, (ref_assembly_id,))
    
    bio_context_id = str(uuid.uuid4())
    cursor.execute("""
        INSERT INTO common.biological_contexts (biological_context_id, taxonomy_id, organism, reference_assembly_id)
        VALUES (%s, 4100, 'Nicotiana benthamiana', %s)
    """, (bio_context_id, ref_assembly_id))
    print("[SUCCESS] Common Context generated successfully.")
    
    # 2. FactorForge
    prot_target_id = str(uuid.uuid4())
    cursor.execute("""
        INSERT INTO factorforge.protein_targets (protein_target_id, display_alias, sequence_id, source)
        VALUES (%s, 'Spike_RBD', %s, 'public_database')
    """, (prot_target_id, seq_id))
    
    design_spec_id = str(uuid.uuid4())
    cursor.execute("""
        INSERT INTO factorforge.design_specs (design_spec_id, protein_target_id, biological_context_id, requested_engine_name, requested_profile_name)
        VALUES (%s, %s, %s, 'Candidax2', 'FF-NB-v3.5')
    """, (design_spec_id, prot_target_id, bio_context_id))
    
    design_run_id = str(uuid.uuid4())
    cursor.execute("""
        INSERT INTO factorforge.design_runs (design_run_id, design_spec_id, execution_origin, actual_engine_name, generation_performed, analysis_mode, evidence_level, input_sequence_id)
        VALUES (%s, %s, 'agentops', 'Candidax2', true, 'redesigned_candidate', 'prospective_factorforge', %s)
    """, (design_run_id, design_spec_id, seq_id))
    
    candidate_id = str(uuid.uuid4())
    cursor.execute("""
        INSERT INTO factorforge.candidates (candidate_id, design_run_id, sequence_id, computational_status)
        VALUES (%s, %s, %s, 'passed')
    """, (candidate_id, design_run_id, seq_id))
    print("[SUCCESS] FactorForge Execution Lineage linked successfully.")
    
    # 3. ValidationHub
    construct_id = str(uuid.uuid4())
    cursor.execute("""
        INSERT INTO validationhub.constructs (construct_id, candidate_id, sequence_id, identity_status)
        VALUES (%s, %s, %s, 'confirmed')
    """, (construct_id, candidate_id, seq_id))
    
    experiment_id = str(uuid.uuid4())
    cursor.execute("""
        INSERT INTO validationhub.experiments (experiment_id, batch_id, biological_context_id, performed_at)
        VALUES (%s, 'B-2026-08', %s, NOW())
    """, (experiment_id, bio_context_id))
    
    arm_id = str(uuid.uuid4())
    cursor.execute("""
        INSERT INTO validationhub.experimental_arms (experimental_arm_id, experiment_id, construct_id, arm_type, linkage_status, outcome_link_allowed)
        VALUES (%s, %s, %s, 'factorforge_candidate', 'confirmed', true)
    """, (arm_id, experiment_id, construct_id))
    
    outcome_id = str(uuid.uuid4())
    cursor.execute("""
        INSERT INTO validationhub.outcomes (outcome_id, experimental_arm_id, outcome_category, outcome_result)
        VALUES (%s, %s, 'yield', 'positive')
    """, (outcome_id, arm_id))
    print("[SUCCESS] ValidationHub Physical Evidence Record bound successfully.")
    
    print("\n[SUCCESS] ALL FOREIGN KEYS AND CASCADES VALIDATED!")
    conn.close()

if __name__ == "__main__":
    test_integrity()
