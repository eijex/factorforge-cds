"""Tests for FactorForge Plant RNA Risk Guards and Invariant Verification (Job 394)."""

import pytest

from factorforge.guards.plant_rna_guards import (
    DETECTOR_VERSION,
    EVIDENCE_CLASS,
    SENSE_SCOPE,
    compare_candidate_plant_rna_risks,
    scan_plant_rna_risks,
    validate_candidate_hard_invariants,
)


def test_plant_rna_scanner_clean_sequence() -> None:
    # Construct sequence with NO cryptic donors, polyA motifs, or ARE elements
    clean_seq = "ATGGCATCAACTCAATCTTCTACTGCATCAACAGCACGTGTGGACTCTGTGAAGTAG"
    res = scan_plant_rna_risks(clean_seq)
    
    assert res.summary_status == "CLEAN_SCAN"
    assert res.donor_count == 0
    assert res.polya_count == 0
    assert res.are_count == 0
    assert res.are_clustered is False
    assert len(res.findings) == 0
    assert res.detector_version == DETECTOR_VERSION


def test_plant_rna_scanner_detects_cryptic_donor() -> None:
    # [AC]AGGT[AG]AGT
    seq_with_donors = "ATGGCGAAGGTAAGTCCGCAGGTGAGTTAG"
    res = scan_plant_rna_risks(seq_with_donors)
    
    assert res.donor_count == 2
    assert res.summary_status == "ADVISORY_FINDINGS"
    
    donor_rules = [f for f in res.findings if f.rule_id == "rna.cryptic_splice.v1"]
    assert len(donor_rules) == 2
    assert donor_rules[0].matched_sequence == "AAGGTAAGT"
    assert donor_rules[1].matched_sequence == "CAGGTGAGT"
    assert donor_rules[0].evidence_class == EVIDENCE_CLASS
    assert donor_rules[0].scope == SENSE_SCOPE


def test_plant_rna_scanner_detects_polya() -> None:
    # AATAAA, ATTAAA, AATAAT
    seq_with_polya = "ATGGCCAATAAAGGCCATTAAAGGCAATAATTGA"
    res = scan_plant_rna_risks(seq_with_polya)
    
    assert res.polya_count == 3
    polya_rules = [f for f in res.findings if f.rule_id == "rna.polya_motifs.v1"]
    assert len(polya_rules) == 3
    motifs_found = [f.matched_sequence for f in polya_rules]
    assert motifs_found == ["AATAAA", "ATTAAA", "AATAAT"]


def test_plant_rna_scanner_detects_are_and_clustering() -> None:
    # Isolated ARE
    seq_isolated = "ATGGCCATTTACCGGTAG"
    res_iso = scan_plant_rna_risks(seq_isolated)
    assert res_iso.are_count == 1
    assert res_iso.are_clustered is False

    # Clustered ARE (< 30 nt apart)
    seq_clustered = "ATGGCCATTTACGCATTTAGCTAG"
    res_clust = scan_plant_rna_risks(seq_clustered)
    assert res_clust.are_count == 2
    assert res_clust.are_clustered is True
    assert len(res_clust.cluster_details) >= 1


def test_validate_candidate_hard_invariants() -> None:
    protein = "MAS"
    cds = "ATGGCATCA" # M=ATG, A=GCA, S=TCA
    
    # 1. Perfectly valid
    val_pass = validate_candidate_hard_invariants(
        cds=cds,
        protein=protein,
        gc_min=0.30,
        gc_max=0.60,
    )
    assert val_pass["passed"] is True
    assert val_pass["translation_identity_100pct"] is True
    assert val_pass["internal_stops"] == 0

    # 2. Translation mismatch
    val_trans_fail = validate_candidate_hard_invariants(
        cds="ATGGCACCA", # P instead of S
        protein=protein,
        gc_min=0.30,
        gc_max=0.60,
    )
    assert val_trans_fail["passed"] is False
    assert val_trans_fail["translation_identity_100pct"] is False

    # 3. Internal stop codon
    val_stop_fail = validate_candidate_hard_invariants(
        cds="ATGTAATCA", # Stop in middle
        protein=protein,
        gc_min=0.30,
        gc_max=0.60,
    )
    assert val_stop_fail["passed"] is False
    assert val_stop_fail["internal_stops"] == 1

    # 4. GC out of bounds
    val_gc_fail = validate_candidate_hard_invariants(
        cds=cds,
        protein=protein,
        gc_min=0.60, # cds has 4/9 = 44.4% GC
        gc_max=0.80,
    )
    assert val_gc_fail["passed"] is False
    assert val_gc_fail["gc_feasible"] is False

    # 5. Type IIS BsaI in construct
    val_bsai_fail = validate_candidate_hard_invariants(
        cds=cds,
        protein=protein,
        gc_min=0.30,
        gc_max=0.60,
        left_flank="GGTCTC", # BsaI
    )
    assert val_bsai_fail["passed"] is False
    assert "GGTCTC" in val_bsai_fail["forbidden_motifs_found"]


def test_compare_candidate_plant_rna_risks() -> None:
    cds_1dp = "ATGGCATCAACTAATAAATGA"  # contains 1 polyA
    cds_3dp = "ATGGCATCAACTCAATCTTGA"  # clean
    
    comp = compare_candidate_plant_rna_risks(cds_1dp, cds_3dp, protein_name="TestProtein")
    
    assert comp["1dp_scan"]["polya_count"] == 1
    assert comp["3dp_scan"]["polya_count"] == 0
    assert comp["delta_polya"] == -1
    assert "Plant RNA Risk Comparison (TestProtein):" in comp["summary_table"]
    assert "Note: A sequence scan hit indicates a motif pattern match" in comp["summary_table"]
