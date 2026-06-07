from benchmarks.config import load_benchmark_config
from benchmarks.scoring import score_cds

CFG = load_benchmark_config()


def test_score_cds_full_schema():
    row = score_cds("random_synonymous", "baseline", "p1", "MKT", "ATGAAAACC", CFG, runtime_seconds=0.01)
    expected = {"method","method_type","sequence_id","aa_identity","internal_stop_count",
                "invalid_codon_count","length_multiple_of_three","cai","gc_percent",
                "gc_in_target_range","forbidden_type_iis_site_count","biological_pass",
                "assembly_pass","multi_constraint_pass","runtime_seconds"}
    assert expected.issubset(row.keys())
    assert row["biological_pass"] is True
    assert row["aa_identity"] == 1.0
