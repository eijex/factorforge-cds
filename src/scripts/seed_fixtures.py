import os
import uuid
import psycopg2
from factorforge.db.connector import FactorForgeDBConnector

def seed_paper_fixtures():
    print("Starting Paper 2 Synthetic Fixture Seed...")
    
    # 1. Connect to PostgreSQL
    dsn = os.environ.get(
        "FACTORFORGE_DB_DSN", 
        "host=127.0.0.1 port=5432 user=postgres password=postgres dbname=postgres"
    )
    conn = psycopg2.connect(dsn)
    conn.autocommit = True
    cursor = conn.cursor()
    
    # 2. Insert Core Biological Contexts & Enymes
    print("Seeding Common / Biological Contexts...")
    
    # Generate static UUIDs for stable fixtures
    ref_id = str(uuid.uuid5(uuid.NAMESPACE_DNS, 'niben.genome'))
    bio_id = str(uuid.uuid5(uuid.NAMESPACE_DNS, 'niben.context'))
    
    cursor.execute("""
        INSERT INTO common.reference_assemblies (reference_assembly_id, organism, taxonomy_id, assembly_name, version)
        VALUES (%s, 'Nicotiana benthamiana', 4100, 'Niben.genome.v1.0.1', '1.0.1')
        ON CONFLICT DO NOTHING
    """, (ref_id,))
    
    cursor.execute("""
        INSERT INTO common.biological_contexts (biological_context_id, taxonomy_id, organism, reference_assembly_id)
        VALUES (%s, 4100, 'Nicotiana benthamiana', %s)
        ON CONFLICT DO NOTHING
    """, (bio_id, ref_id))
    
    print("Seeding FactorForge Constraints & Engines...")
    cursor.execute("""
        INSERT INTO factorforge.enzyme_registry (enzyme_name, recognition_sequence, category)
        VALUES 
            ('BsaI', 'GGTCTC', 'Type_IIS'),
            ('BpiI', 'GAAGAC', 'Type_IIS'),
            ('BsmBI', 'CGTCTC', 'Type_IIS')
        ON CONFLICT (enzyme_name) DO NOTHING
    """)
    
    cursor.execute("""
        INSERT INTO factorforge.constraint_definitions (code, name, rule_type)
        VALUES 
            ('AA_IDENTITY', 'Amino Acid Identity Preservation', 'HARD'),
            ('TYPE_IIS_CLEAN', 'Type IIS Restriction Site Absence', 'HARD'),
            ('GLOBAL_GC', 'Global GC Content Window', 'HARD')
        ON CONFLICT (code) DO NOTHING
    """)
    
    print("Seeding Base Target & Design Spec...")
    connector = FactorForgeDBConnector(dsn)
    
    target_seq = "ATGGCTAAATGGTAA"
    target_seq_id = connector.register_sequence(target_seq, molecule_class="CDS")
    
    prot_target_id = str(uuid.uuid5(uuid.NAMESPACE_DNS, 'target.spike.rbd'))
    cursor.execute("""
        INSERT INTO factorforge.protein_targets (protein_target_id, display_alias, sequence_id, source)
        VALUES (%s, 'Spike_RBD_Mock', %s, 'public_database')
        ON CONFLICT DO NOTHING
    """, (prot_target_id, target_seq_id))
    
    design_spec_id = str(uuid.uuid5(uuid.NAMESPACE_DNS, 'spec.001'))
    cursor.execute("""
        INSERT INTO factorforge.design_specs (design_spec_id, protein_target_id, biological_context_id, requested_engine_name, requested_profile_name)
        VALUES (%s, %s, %s, 'Candidax2', 'FF-NB-v3.5')
        ON CONFLICT DO NOTHING
    """, (design_spec_id, prot_target_id, bio_id))
    
    design_run_id = str(uuid.uuid5(uuid.NAMESPACE_DNS, 'run.001'))
    cursor.execute("""
        INSERT INTO factorforge.design_runs (design_run_id, design_spec_id, execution_origin, actual_engine_name, generation_performed, analysis_mode, evidence_level, input_sequence_id)
        VALUES (%s, %s, 'agentops', 'Candidax2', true, 'redesigned_candidate', 'prospective_factorforge', %s)
        ON CONFLICT DO NOTHING
    """, (design_run_id, design_spec_id, target_seq_id))
    
    print("Saving candidates via ORM/Connector...")
    # 3. Use Connector to save 5 mock candidates
    for i in range(5):
        res = connector.save_candidate_with_evaluations(
            candidate_data={
                "engine": "Candidax2",
                "model_version": "v3.5.0",
                "optimized_sequence": f"ATGGCTAAATGGTAA{'C'*i}",
                "cai": 0.84 + (i*0.01),
                "gc_percent": 41.2 + i,
                "type2is_clean": True,
            },
            evaluations=[
                {"constraint_code": "AA_IDENTITY", "status": "PASS", "observed_value": 1.0, "details": {"target": "MAKW"}},
                {"constraint_code": "TYPE_IIS_CLEAN", "status": "PASS", "observed_value": 0, "details": {"enzymes": ["BsaI", "BpiI", "BsmBI"]}},
            ],
            design_run_id=design_run_id,
        )
        print(f"  -> Candidate {i+1} saved with package {res['package_id']}")

    print("\n[SUCCESS] Paper 2 Fixture Data successfully seeded into PostgreSQL!")
    conn.close()

if __name__ == "__main__":
    seed_paper_fixtures()
