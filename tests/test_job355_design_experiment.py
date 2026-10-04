"""
Unit and Integration Tests for Job 355 (Phase 2.5):
DesignExperiment, Variant Lineage, Structured EvidenceRecord, CodonPolicy,
Deterministic DesignPipeline, and Full Humira LC/HC Construct Verification.

Validates the full set of user-specified criteria:
1. Variant supports an arbitrary number of candidate designs.
2. Lineage is traceable via parent_variant_id up to the root ancestor.
3. EvidenceRecord is a structured object (method, confidence, level) rather than a string log.
4. DesignPipeline is 100% deterministic (identical inputs -> identical output hashes without LLM).
5. Serialization round-trips preserve complete schema fidelity for package/UI consumption.
6. HARD_FAIL gating enforces valid_for_release=False, while Type IIS auto-repair cleanly eliminates sites.
7. Structured CodonPolicy object uses Hamilton Largest Remainder allocation and repeat avoidance.
8. Splice remediation performs motif-local synonymous mutations, preserving global codon policy balance.
9. Fake CAI placeholders removed (cai_score is None, NOT_COMPUTED) and ValidationLevel consistency.
10. Fixed timestamp support enables bit-for-bit JSON serialization determinism.
11. Full Humira LC (235 aa, 708 nt) and HC (472 aa, 1419 nt) regression verification with real seed DNA.
"""

import json
import pytest

from factorforge.core.experiment import (
    ComparisonMatrix,
    DesignExperiment,
    EvidenceRecord,
    Severity,
    ValidationLevel,
    Variant,
)
from factorforge.core.policy import (
    CodonPolicy,
    ExcludedCodon,
    PLANTFORM_BALANCED_CODON_V1,
    get_policy,
)
from factorforge.pipeline.design_pipeline import DesignPipeline, translate_cds, calculate_gc

# Real Full-Length Adalimumab constructs from chimeric_humira.fasta
FULL_HUMIRA_LC_DNA = (
    "ATGGCAAAAACAAATTTATTTTTATTTTTAATATTTTCATTATTACTGTCTCTGTCCTCTGCT"
    "GACATTCAAATGACTCAAAGCCCTTCTTCTCTTTCTGCTTCTGTTGGAGACAGAGTCACCATCACTTGCAGAGCTTCTCAAGGAATCAGAAATTACCTTGCT"
    "TGGTATCAACAGAAGCCTGGAAAAGCTCCCAAACTTTTGATTTATGCTGCTTCCACTCTTCAATCTGGAGTTCCATCCAGATTTTCTGGATCTGGATCTGGAA"
    "CTGATTTCACTTTGACCATTTCTTCTTTGCAACCTGAAGATGTTGCCACTTACTATTGTCAAAGATACAACAGAGCACCATACACTTTTGGACAGGGAACCAA"
    "AGTTGAAATCAAGAGAACTGTTGCTGCTCCTTCTGTTTTCATTTTTCCTCCTTCTGATGAACAATTGAAATCTGGAACTGCCTCTGTTGTTTGTCTTTTGAAC"
    "AATTTTTATCCCAGAGAAGCCAAAGTTCAATGGAAGGTTGACAATGCTCTTCAATCTGGAAACAGCCAGGAATCTGTCACTGAACAAGATTCCAAAGATTCCA"
    "CTTATTCTCTCTCTTCCACTTTGACACTTTCCAAGGCTGATTATGAGAAACACAAAGTTTATGCATGTGAAGTCACTCATCAAGGACTTTCTTCTCCTGTCAC"
    "CAAATCATTCAACAGAGGAGAATGT"
)
FULL_HUMIRA_LC_PROTEIN = translate_cds(FULL_HUMIRA_LC_DNA)

FULL_HUMIRA_HC_DNA = (
    "ATGGCAAAAACAAATTTATTTTTATTTTTAATATTTTCATTATTACTGTCTCTGTCCTCTGCT"
    "GAAGTGCAATTGGTTGAGTCTGGAGGAGGACTGGTTCAACCTGGAAGATCTTTGAGACTTTCTTGTGCTGCTTCTGGATTCACTTTTGATGATTATGCCATGC"
    "ATTGGGTCAGGCAGGCTCCTGGAAAGGGACTTGAATGGGTTTCTGCCATCACCTGGAATTCAGGGCACATTGATTATGCTGACTCTGTTGAAGGAAGATTCAC"
    "CATTTCCAGAGACAATGCCAAAAATTCTCTTTATCTTCAGATGAATTCTTTGAGAGCTGAGGACACTGCTGTTTATTATTGTGCCAAAGTTTCTTATCTTTCC"
    "ACTGCCAGCTCTTTGGATTATTGGGGACAAGGAACTCTTGTCACTGTTTCTTCTGCTTCCACAAAAGGGCCTTCTGTTTTTCCTCTTGCTCCTTCATCCAAAT"
    "CCACTTCTGGTGGAACTGCTGCTCTTGGATGTCTTGTCAAGGATTATTTTCCTGAACCTGTCACTGTGTCTTGGAATTCTGGAGCTTTGACTTCTGGAGTTCA"
    "CACATTTCCTGCAGTGCTTCAGTCCTCTGGCCTTTACAGTCTTTCTTCTGTTGTCACTGTTCCATCTTCTTCTTTGGGAACTCAAACTTACATTTGCAATGTC"
    "AATCACAAACCATCCAACACCAAAGTTGACAAAAAAGTTGAACCCAAGTCTTGTGACAAGACTCACACTTGTCCTCCTTGTCCTGCACCTGAACTTCTTGGAG"
    "GACCTTCTGTTTTTCTTTTCCCTCCCAAGCCAAAAGACACTTTGATGATTTCCAGAACTCCTGAAGTCACTTGTGTTGTTGTTGATGTTTCTCATGAAGATCC"
    "TGAAGTCAAATTCAACTGGTATGTTGATGGAGTGGAAGTTCACAATGCCAAAACCAAACCCAGAGAAGAACAATACAATTCCACATACAGGGTTGTTTCTGTT"
    "TTGACTGTTCTTCATCAGGATTGGTTGAATGGAAAAGAGTACAAATGCAAAGTTTCCAACAAAGCTCTTCCTGCTCCCATTGAAAAAACCATTTCCAAAGCCA"
    "AAGGACAACCCAGAGAGCCTCAAGTCTACACTCTTCCTCCTTCCAGAGATGAATTGACCAAGAATCAGGTTTCTTTGACCTGTCTTGTCAAAGGATTTTATCC"
    "TTCTGACATTGCTGTTGAGTGGGAATCCAATGGACAACCTGAAAACAATTACAAAACCACTCCCCCTGTTCTTGACAGTGATGGATCTTTTTTTCTTTATTCC"
    "AAATTGACTGTTGACAAATCCAGATGGCAACAAGGCAATGTTTTTTCTTGCAGTGTCATGCATGAGGCTCTTCACAATCATTACACCCAAAAATCTCTTTCTT"
    "TGTCTCCAGGAAAA"
)
FULL_HUMIRA_HC_PROTEIN = translate_cds(FULL_HUMIRA_HC_DNA)


class TestJob355DesignExperiment:
    """Test Suite for Job 355 Core Models, CodonPolicy, and Pipeline Verification."""

    def test_criterion_1_arbitrary_variant_capacity(self):
        """Criterion 1: DesignExperiment supports an arbitrary number of variants."""
        exp = DesignExperiment(
            experiment_id="EXP-TEST-01",
            target_name="TestTarget",
            host_organism="nbenthamiana",
            hypothesis="Arbitrary variant capacity verification",
        )

        assert len(exp.variants) == 0

        # Add 10 arbitrary variants with different intents
        for i in range(1, 11):
            var = Variant(
                variant_id=f"VAR-{i:02d}",
                parent_variant_id="VAR-01" if i > 1 else None,
                design_intent=f"Experimental condition {i}",
                interventions=[f"intervention_{i}"],
                sequence=f"ATG{'GCT' * i}TAA",
                protein_sequence=f"M{'A' * i}",
                gc_percent=50.0,
                cai_score=None,
                policy_version="plantform_balanced_codon_v1",
                valid_for_release=True,
            )
            exp.add_variant(var)

        assert len(exp.variants) == 10
        assert exp.get_variant("VAR-05").design_intent == "Experimental condition 5"

        with pytest.raises(ValueError, match="already exists"):
            exp.add_variant(
                Variant(
                    variant_id="VAR-01",
                    parent_variant_id=None,
                    design_intent="Duplicate",
                    interventions=[],
                    sequence="ATGTAA",
                    protein_sequence="M",
                    gc_percent=50.0,
                    cai_score=None,
                    policy_version="test",
                )
            )

    def test_criterion_2_lineage_traceability(self):
        """Criterion 2: Lineage is traceable via parent_variant_id back to root ancestor."""
        exp = DesignExperiment(
            experiment_id="EXP-LINEAGE-01",
            target_name="Adalimumab LC",
            host_organism="nbenthamiana",
            hypothesis="Traceable multi-step lineage",
        )

        v1 = Variant(
            variant_id="LC-VAR-01-BASELINE",
            parent_variant_id=None,
            design_intent="Unmodified baseline",
            interventions=["none"],
            sequence="ATGGCTTAA",
            protein_sequence="MA",
            gc_percent=50.0,
            cai_score=None,
            policy_version="plantform_balanced_codon_v1",
        )
        v2 = Variant(
            variant_id="LC-VAR-02-WATERMARKED",
            parent_variant_id="LC-VAR-01-BASELINE",
            design_intent="Provenance encoding",
            interventions=["watermark"],
            sequence="ATGCCCTAA",
            protein_sequence="MA",
            gc_percent=66.67,
            cai_score=None,
            policy_version="plantform_balanced_codon_v1",
        )
        v3 = Variant(
            variant_id="LC-VAR-03-SPLICE-REMEDIATED",
            parent_variant_id="LC-VAR-02-WATERMARKED",
            design_intent="Splice remediation on top of watermark",
            interventions=["watermark", "splice_remediation"],
            sequence="ATGCCGTAA",
            protein_sequence="MA",
            gc_percent=66.67,
            cai_score=None,
            policy_version="plantform_balanced_codon_v1",
        )

        exp.add_variant(v1)
        exp.add_variant(v2)
        exp.add_variant(v3)

        lineage = exp.get_lineage("LC-VAR-03-SPLICE-REMEDIATED")
        assert len(lineage) == 3
        assert [v.variant_id for v in lineage] == [
            "LC-VAR-01-BASELINE",
            "LC-VAR-02-WATERMARKED",
            "LC-VAR-03-SPLICE-REMEDIATED",
        ]

    def test_criterion_3_structured_evidence_record_schema(self):
        """Criterion 3: EvidenceRecord is a structured object with method, confidence, and validation_level."""
        rec = EvidenceRecord(
            check_id="rna.cryptic_splice",
            category="rna_stability",
            result="WARNING",
            severity=Severity.WARNING,
            method="NetGene2PredictorBridge",
            tool_version="v2.42-plant",
            confidence=0.74,
            evidence_detail="Potential donor motif flagged at nt 441 (score: 0.74)",
            validation_level=ValidationLevel.COMPUTATIONAL_PREDICTION,
        )

        assert rec.check_id == "rna.cryptic_splice"
        assert rec.severity == Severity.WARNING
        assert rec.confidence == 0.74
        assert rec.validation_level == ValidationLevel.COMPUTATIONAL_PREDICTION
        assert "441" in rec.evidence_detail

        d = rec.to_dict()
        assert isinstance(d, dict)
        assert d["severity"] == "WARNING"
        assert d["validation_level"] == "computational_prediction"

        reconstructed = EvidenceRecord.from_dict(d)
        assert reconstructed == rec

    def test_criterion_4_deterministic_pipeline_reproducibility(self):
        """Criterion 4: DesignPipeline is 100% deterministic without LLM across runs."""
        pipeline = DesignPipeline(
            host="nbenthamiana",
            policy="plantform_balanced_codon_v1",
        )

        hashes_run_1: list[str] = []

        # Run 1
        exp1 = pipeline.run(
            target_name="Adalimumab_LC",
            protein_sequence=FULL_HUMIRA_LC_PROTEIN,
            seed_dna=FULL_HUMIRA_LC_DNA,
            frozen_boundary_nt=60,
            interventions=["baseline", "watermark", "splice_remediation"],
            fixed_timestamp="2026-10-04T12:00:00Z",
        )
        for v in exp1.variants:
            hashes_run_1.append(v.provenance_hash)

        # Runs 2-5
        for _ in range(4):
            exp_n = pipeline.run(
                target_name="Adalimumab_LC",
                protein_sequence=FULL_HUMIRA_LC_PROTEIN,
                seed_dna=FULL_HUMIRA_LC_DNA,
                frozen_boundary_nt=60,
                interventions=["baseline", "watermark", "splice_remediation"],
                fixed_timestamp="2026-10-04T12:00:00Z",
            )
            hashes_n = [v.provenance_hash for v in exp_n.variants]
            assert hashes_n == hashes_run_1, "Non-deterministic output detected across runs!"

        assert len(exp1.variants) == 3
        for v in exp1.variants:
            assert v.sequence[:60] == FULL_HUMIRA_LC_DNA[:60].upper()
            assert translate_cds(v.sequence) == FULL_HUMIRA_LC_PROTEIN

    def test_criterion_5_end_to_end_serialization_fidelity(self):
        """Criterion 5: Validation and experiment schema survives complete JSON round-trip."""
        pipeline = DesignPipeline(
            host="nbenthamiana",
            policy="plantform_balanced_codon_v1",
        )
        exp = pipeline.run(
            target_name="Adalimumab_LC",
            protein_sequence=FULL_HUMIRA_LC_PROTEIN,
            seed_dna=FULL_HUMIRA_LC_DNA,
            frozen_boundary_nt=60,
            interventions=["baseline", "watermark", "splice_remediation"],
            fixed_timestamp="2026-10-04T12:00:00Z",
        )

        json_output = exp.to_json()
        assert isinstance(json_output, str)

        reloaded = DesignExperiment.from_json(json_output)
        assert reloaded.experiment_id == exp.experiment_id
        assert len(reloaded.variants) == len(exp.variants)
        assert len(reloaded.comparisons) == len(exp.comparisons)

        for orig_v, rel_v in zip(exp.variants, reloaded.variants):
            assert orig_v.variant_id == rel_v.variant_id
            assert orig_v.provenance_hash == rel_v.provenance_hash
            assert orig_v.valid_for_release == rel_v.valid_for_release
            assert orig_v.cai_score == rel_v.cai_score
            for orig_r, rel_r in zip(orig_v.validation_records, rel_v.validation_records):
                assert orig_r.check_id == rel_r.check_id
                assert orig_r.validation_level == rel_r.validation_level
                assert orig_r.result == rel_r.result

    def test_criterion_6_hard_fail_gating_and_type_iis_autorepair(self):
        """Criterion 6: HARD_FAIL blocks release (valid_for_release=False) while auto-repair clears Type IIS."""
        # A variant with a HARD_FAIL must have valid_for_release=False
        fail_record = EvidenceRecord(
            check_id="assembly.type_iis.bpii",
            category="cloning_hygiene",
            result="FAIL",
            severity=Severity.HARD_FAIL,
            method="TypeIISRegexScanner",
            tool_version="v3.5.4",
            confidence=1.0,
            evidence_detail="1 recognition site detected",
            validation_level=ValidationLevel.DETERMINISTIC_CHECK,
        )
        v_failed = Variant(
            variant_id="VAR-FAIL-TEST",
            parent_variant_id=None,
            design_intent="Test failure gating",
            interventions=["none"],
            sequence="ATGGAAGACTAA",  # Contains GAAGAC (BpiI)
            protein_sequence="MED*",
            gc_percent=40.0,
            cai_score=None,
            policy_version="plantform_balanced_codon_v1",
            valid_for_release=True,  # Should be forced to False by post_init
            validation_records=[fail_record],
        )
        assert v_failed.valid_for_release is False

        # In pipeline execution with mature sequence containing synthetic BpiI/BsaI:
        pipeline = DesignPipeline(host="nbenthamiana", policy="plantform_balanced_codon_v1")
        exp = pipeline.run(
            target_name="Adalimumab_LC",
            protein_sequence=FULL_HUMIRA_LC_PROTEIN,
            seed_dna=FULL_HUMIRA_LC_DNA,
            frozen_boundary_nt=60,
            interventions=["baseline", "watermark", "splice_remediation"],
        )

        for v in exp.variants:
            # Verify all Type IIS sites are 0 and result is PASS
            for r in v.validation_records:
                if "assembly.type_iis" in r.check_id:
                    assert r.result == "PASS"
                    assert r.severity == Severity.HARD_FAIL
                    assert "0 recognition sites detected" in r.evidence_detail
            # Must be valid for release
            assert v.valid_for_release is True

    def test_criterion_7_structured_codon_policy_object_and_hamilton_allocation(self):
        """Criterion 7: Structured CodonPolicy uses Hamilton allocation and excludes CpG codons."""
        policy = PLANTFORM_BALANCED_CODON_V1
        assert policy.name == "plantform_balanced_codon_v1"
        assert len(policy.excluded_codons) >= 8

        # Verify excluded codons have verifiable rationales and sources
        cga = next(e for e in policy.excluded_codons if e.codon == "CGA")
        assert "silencing" in cga.reason.lower()
        assert "Doug" in cga.evidence_source or "PlantForm" in cga.evidence_source

        # Test Hamilton allocation on full LC mature sequence
        frozen_aa = len(FULL_HUMIRA_LC_DNA[:60]) // 3
        mature_aa = FULL_HUMIRA_LC_PROTEIN[frozen_aa:]
        allocated = policy.allocate_sequence(mature_aa)

        assert len(allocated) == len(mature_aa)
        # Check translation
        for aa, codon in zip(mature_aa, allocated):
            assert translate_cds(codon) == aa

        # Check Arg 1:1 balance in mature chain
        arg_codons = [c for c in allocated if translate_cds(c) == "R"]
        if len(arg_codons) > 1:
            aga_count = arg_codons.count("AGA")
            agg_count = arg_codons.count("AGG")
            # Hamilton guarantees integer difference <= 1
            assert abs(aga_count - agg_count) <= 1

        # Check that NO CpG dinucleotide codons are present in allocated codons
        cpg_codons = {"CGA", "CGC", "CGG", "CGT", "TCG", "CCG", "GCG", "ACG"}
        for c in allocated:
            assert c not in cpg_codons, f"CpG codon {c} unexpectedly allocated!"

    def test_criterion_8_motif_local_splice_remediation(self):
        """Criterion 8: Splice remediation mutates ONLY the specific motif codons, preserving global codon ratios."""
        pipeline = DesignPipeline(host="nbenthamiana", policy="plantform_balanced_codon_v1")
        exp = pipeline.run(
            target_name="Adalimumab_LC",
            protein_sequence=FULL_HUMIRA_LC_PROTEIN,
            seed_dna=FULL_HUMIRA_LC_DNA,
            frozen_boundary_nt=60,
            interventions=["baseline", "splice_remediation"],
        )

        baseline = exp.get_variant("Adalimumab_LC-VAR-01-BASELINE")
        remediated = exp.get_variant("Adalimumab_LC-VAR-02-SPLICE-REMEDIATED")

        # Translation must remain identical
        assert translate_cds(baseline.sequence) == FULL_HUMIRA_LC_PROTEIN
        assert translate_cds(remediated.sequence) == FULL_HUMIRA_LC_PROTEIN

        # If baseline had splice sites, divergence should be minimal (motif-local only, not whole-gene global replace)
        comp = next(c for c in exp.comparisons if c.variant_b_id == remediated.variant_id)
        # Sequence divergence must be low (< 15 nt out of 708 nt)
        assert comp.sequence_divergence_nt < 15
        assert abs(comp.gc_shift_percent) < 2.0

    def test_criterion_9_cai_placeholder_eliminated_and_validation_level_consistency(self):
        """Criterion 9: CAI placeholder removed (None / NOT_COMPUTED) and ValidationLevel is NOT_VALIDATED."""
        pipeline = DesignPipeline(host="nbenthamiana", policy="plantform_balanced_codon_v1")
        exp = pipeline.run(
            target_name="Adalimumab_LC",
            protein_sequence=FULL_HUMIRA_LC_PROTEIN,
            seed_dna=FULL_HUMIRA_LC_DNA,
            frozen_boundary_nt=60,
        )

        for v in exp.variants:
            # Variant cai_score must be None
            assert v.cai_score is None

            # CAI record must be NOT_COMPUTED and NOT_VALIDATED
            cai_rec = next(r for r in v.validation_records if r.check_id == "genetics.cai_score")
            assert cai_rec.result == "NOT_COMPUTED"
            assert cai_rec.confidence is None
            assert cai_rec.validation_level == ValidationLevel.NOT_VALIDATED

            # Wet lab record must be NOT_TESTED and NOT_VALIDATED (resolving semantic conflict)
            wet_rec = next(r for r in v.validation_records if r.check_id == "wet_lab.expression_status")
            assert wet_rec.result == "NOT_TESTED"
            assert wet_rec.validation_level == ValidationLevel.NOT_VALIDATED

    def test_criterion_10_deterministic_fixed_timestamp_reproducibility(self):
        """Criterion 10: Fixed timestamp yields bit-for-bit identical JSON output across runs."""
        pipeline = DesignPipeline(host="nbenthamiana", policy="plantform_balanced_codon_v1")
        ts = "2026-10-04T16:00:00Z"

        exp_a = pipeline.run(
            target_name="Adalimumab_LC",
            protein_sequence=FULL_HUMIRA_LC_PROTEIN,
            seed_dna=FULL_HUMIRA_LC_DNA,
            frozen_boundary_nt=60,
            fixed_timestamp=ts,
        )
        exp_b = pipeline.run(
            target_name="Adalimumab_LC",
            protein_sequence=FULL_HUMIRA_LC_PROTEIN,
            seed_dna=FULL_HUMIRA_LC_DNA,
            frozen_boundary_nt=60,
            fixed_timestamp=ts,
        )

        assert exp_a.to_json() == exp_b.to_json()

    def test_criterion_11_full_humira_lc_and_hc_regression_fixtures(self):
        """Criterion 11: Full Humira LC (235 aa, 708 nt) and HC (472 aa, 1419 nt) regression verification."""
        pipeline = DesignPipeline(host="nbenthamiana", policy="plantform_balanced_codon_v1")

        # 1. Full Light Chain (235 aa, 708 nt CDS)
        exp_lc = pipeline.run(
            target_name="Humira_Full_LC",
            protein_sequence=FULL_HUMIRA_LC_PROTEIN,
            seed_dna=FULL_HUMIRA_LC_DNA,
            frozen_boundary_nt=60,
            interventions=["baseline", "watermark", "splice_remediation"],
            fixed_timestamp="2026-10-04T12:00:00Z",
        )
        assert len(exp_lc.variants) == 3
        for v in exp_lc.variants:
            assert len(v.sequence) == 708  # 705 nt + TAA
            assert len(v.protein_sequence) == 235
            assert v.sequence[:60] == FULL_HUMIRA_LC_DNA[:60]
            assert translate_cds(v.sequence) == FULL_HUMIRA_LC_PROTEIN
            assert v.valid_for_release is True

        # Check Watermark Confounder Warning
        comp_wm_lc = next(c for c in exp_lc.comparisons if "WATERMARKED" in c.variant_b_id)
        if abs(comp_wm_lc.gc_shift_percent) > 5.0:
            assert "confounder_warning" in comp_wm_lc.metadata

        # 2. Full Heavy Chain (472 aa, 1419 nt CDS)
        exp_hc = pipeline.run(
            target_name="Humira_Full_HC",
            protein_sequence=FULL_HUMIRA_HC_PROTEIN,
            seed_dna=FULL_HUMIRA_HC_DNA,
            frozen_boundary_nt=60,
            interventions=["baseline", "watermark", "splice_remediation"],
            fixed_timestamp="2026-10-04T12:00:00Z",
        )
        assert len(exp_hc.variants) == 3
        for v in exp_hc.variants:
            assert len(v.sequence) == 1419  # 1416 nt + TAA
            assert len(v.protein_sequence) == 472
            assert v.sequence[:60] == FULL_HUMIRA_HC_DNA[:60]
            assert translate_cds(v.sequence) == FULL_HUMIRA_HC_PROTEIN
            assert v.valid_for_release is True

        comp_wm_hc = next(c for c in exp_hc.comparisons if "WATERMARKED" in c.variant_b_id)
        if abs(comp_wm_hc.gc_shift_percent) > 5.0:
            assert "confounder_warning" in comp_wm_hc.metadata

    def test_criterion_12_evidence_package_compiler(self, tmp_path):
        """Criterion 12: EvidencePackageCompiler produces complete, checksum-verified ZIP archives."""
        from factorforge.io import EvidencePackageCompiler
        import zipfile
        import io

        pipeline = DesignPipeline(host="nbenthamiana", policy="plantform_balanced_codon_v1")
        exp = pipeline.run(
            target_name="Humira_Full_LC",
            protein_sequence=FULL_HUMIRA_LC_PROTEIN,
            seed_dna=FULL_HUMIRA_LC_DNA,
            frozen_boundary_nt=60,
            interventions=["baseline", "watermark", "splice_remediation"],
            fixed_timestamp="2026-10-04T12:00:00Z",
        )

        compiler = EvidencePackageCompiler()

        # 1. Test In-Memory compilation (for Web API download endpoints)
        zip_bytes = compiler.compile_in_memory(exp)
        assert isinstance(zip_bytes, bytes)
        assert len(zip_bytes) > 1000

        # Inspect ZIP contents
        with zipfile.ZipFile(io.BytesIO(zip_bytes), mode="r") as zf:
            namelist = zf.namelist()
            assert "experiment.json" in namelist
            assert "DESIGN_DOSSIER.md" in namelist
            assert "dashboard.html" in namelist
            assert "SHA256SUMS.txt" in namelist
            assert "sequences/all_variants.fasta" in namelist
            assert "sequences/Humira_Full_LC-VAR-01-BASELINE.fasta" in namelist

            # Verify Checksums
            sums_content = zf.read("SHA256SUMS.txt").decode("utf-8")
            for line in sums_content.strip().splitlines():
                expected_hash, rel_path = line.split("  ")
                actual_data = zf.read(rel_path)
                import hashlib
                actual_hash = hashlib.sha256(actual_data).hexdigest()
                assert actual_hash == expected_hash, f"Checksum mismatch for {rel_path}!"

        # 2. Test File-on-disk compilation
        disk_zip = tmp_path / "evidence_package.zip"
        res_path = compiler.compile_package(exp, str(disk_zip))
        assert disk_zip.exists()
        assert disk_zip.stat().st_size > 1000

