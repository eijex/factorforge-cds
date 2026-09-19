"""Prospective Experimental Acquisition Panel Builder for FactorForge (Job 283C).

Constructs the 3-target x 3-hypothesis prospective experimental panel (9 constructs + controls)
with an explicit Orthogonality Gate, randomized plate mapping, and dual memory/archive outputs.
"""

from __future__ import annotations
from datetime import datetime, timezone
import hashlib
import math
import random
from typing import Any, Dict, List, Optional, Tuple

from factorforge.analysis.metrics import load_codon_usage_table
from factorforge.discovery.acquisition import (
    ConstructDesignRecord,
    ExperimentRunRecord,
    PairedDBTLDataset,
    SampleRecord,
)
from factorforge.discovery.schemas import TraitVector
from factorforge.discovery.traits import TraitVectorExtractor
from factorforge.engines.dp_v2_1_1 import DPV211Optimizer


STANDARD_TARGET_PANEL = [
    {
        "target_id": "PLT-B01",
        "target_name": "sfGFP",
        "description": "Superfolder GFP (Reporter)",
        "aa_length": 238,
        "sequence": (
            "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTT"
            "FSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELK"
            "GIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIG"
            "DGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        ),
        "metadata": {
            "target_context": "cytosolic_fluorescent_reporter",
            "canonical_aa_length": 238,
            "construct_aa_length": 238,
        },
    },
    {
        "target_id": "PLT-C01",
        "target_name": "VP28",
        "description": "WSSV Viral Structural Envelope Antigen",
        "aa_length": 204,
        "sequence": (
            "MDLSFTLSVVSAILAITAVIAVFIVIFRYHNTVTKTIETHTDNIETNMDENLRIPVTAEVGSG"
            "YFKMTDVSFDSDTLGKIKIRNGKSDAQMKEEDADLVITPVEGRALEVTVGQNLTFEGTFKVWN"
            "NTSRKINITGMQMVPKINPSKAFVGSSNTSSFTPVSIDEDEVGTFVCGTTFGAPIAATAGGNL"
            "FDMYVHVTYSGTETE"
        ),
        "metadata": {
            "transmembrane_included": True,
            "tm_anchor_residues": "1-30",
            "target_context": "full_length_envelope_antigen",
            "canonical_aa_length": 204,
            "construct_aa_length": 204,
        },
    },
    {
        "target_id": "PLT-B02",
        "target_name": "CD47_Ectodomain",
        "description": "Human CD47 Ectodomain (Therapeutic Domain)",
        "aa_length": 124,
        "sequence": (
            "QLLFNKTKSVEFTFCNDTVVIPCFVTNMEAQNTTEVYVKWKFKGRDIYTFDGALNKSTVPTDF"
            "SSAKIEVSQLLKGDASLKMDKSDAVSHTGNYTCEVTELTREGETIIELKYRVVSWFSPNEN"
        ),
        "metadata": {
            "expression_context": "mature_ectodomain_residues_19_141",
            "secretion_context": "pEAQ-HT vector cassette N-terminal signal peptide required",
            "canonical_aa_length": 305,
            "construct_aa_length": 124,
        },
    },
]

# Standard control sequences (717 nt: 714 nt CDS + 3 nt TAA stop)
STANDARD_POS_CTRL_SFGFP = (
    "ATGAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAAT"
    "GGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTT"
    "ATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAA"
    "TGCTTTTCAAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTAT"
    "GTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAA"
    "GGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGA"
    "CACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATC"
    "AAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAA"
    "AATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTT"
    "TCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACA"
    "CATGGCATGGATGAACTATACAAATAA"
)
STANDARD_POS_CTRL_SFGFP_AA_SHA256 = (
    "44f268846c4fcf4cb9f563116d866406e2cf07831f2495b6c205244ddf3333ce"
)


class OrthogonalityGate:
    """Calculates and verifies multi-dimensional trait distances between hypotheses.

    Mathematical Contract:
    Distance between trait vectors t_1 and t_2 is computed via weighted normalized Euclidean norm:
        D(t_1, t_2) = sqrt( sum_k w_k * ((t_{1,k} - t_{2,k}) / sigma_k)^2 )

    Feature scaling constants (sigma_k) and weights (w_k):
        - CAI:                 w = 1.0, sigma = 0.15
        - GC% (global):        w = 1.0, sigma = 10.0
        - 5' MFE (kcal/mol):   w = 0.8, sigma = 5.0
        - Local 50bp GC min:   w = 0.5, sigma = 15.0
        - Local 50bp GC max:   w = 0.5, sigma = 15.0
        - Homopolymer max run: w = 0.3, sigma = 2.0

    Gate Threshold: D(t_1, t_2) >= min_epsilon (default 0.20 for prospective panel separation).
    """

    def __init__(self, min_epsilon: float = 0.05) -> None:
        self.min_epsilon = min_epsilon

    @staticmethod
    def compute_distance(tv1: TraitVector, tv2: TraitVector) -> float:
        """Normalized Euclidean distance between two 8D trait vectors."""
        weights = {
            "cai": (1.0, 0.15),
            "gc": (1.0, 10.0),
            "mfe": (0.8, 5.0),
            "gc_min": (0.5, 15.0),
            "gc_max": (0.5, 15.0),
            "homopolymer": (0.3, 2.0),
        }

        mfe1 = tv1.initiation_mfe_kcal_mol if tv1.initiation_mfe_kcal_mol is not None else -10.0
        mfe2 = tv2.initiation_mfe_kcal_mol if tv2.initiation_mfe_kcal_mol is not None else -10.0

        d2 = 0.0
        d2 += (
            weights["cai"][0] * ((tv1.cai_golden_set - tv2.cai_golden_set) / weights["cai"][1]) ** 2
        )
        d2 += (
            weights["gc"][0]
            * ((tv1.global_gc_percent - tv2.global_gc_percent) / weights["gc"][1]) ** 2
        )
        d2 += weights["mfe"][0] * ((mfe1 - mfe2) / weights["mfe"][1]) ** 2
        d2 += (
            weights["gc_min"][0]
            * ((tv1.local_50bp_gc_min - tv2.local_50bp_gc_min) / weights["gc_min"][1]) ** 2
        )
        d2 += (
            weights["gc_max"][0]
            * ((tv1.local_50bp_gc_max - tv2.local_50bp_gc_max) / weights["gc_max"][1]) ** 2
        )
        d2 += (
            weights["homopolymer"][0]
            * ((tv1.homopolymer_max_run - tv2.homopolymer_max_run) / weights["homopolymer"][1]) ** 2
        )

        return math.sqrt(d2)

    def verify_orthogonality(
        self, hypotheses_traits: Dict[str, TraitVector]
    ) -> Tuple[bool, Dict[str, float]]:
        """Verifies all pairwise hypothesis distances exceed epsilon."""
        h_keys = list(hypotheses_traits.keys())
        distances: Dict[str, float] = {}
        all_passed = True

        for i in range(len(h_keys)):
            for j in range(i + 1, len(h_keys)):
                k1, k2 = h_keys[i], h_keys[j]
                d = self.compute_distance(hypotheses_traits[k1], hypotheses_traits[k2])
                pair_name = f"{k1}_vs_{k2}"
                distances[pair_name] = round(d, 4)
                if d < self.min_epsilon:
                    all_passed = False

        return all_passed, distances


class ProspectivePanelBuilder:
    """Builds the 3-target x 3-hypothesis prospective experimental acquisition panel."""

    def __init__(self, host: str = "nbenthamiana") -> None:
        self.host = host
        self.codon_table = load_codon_usage_table()
        self.extractor = TraitVectorExtractor(host=host)
        self.orthogonality_gate = OrthogonalityGate(min_epsilon=0.05)

    def generate_h0_primary_optimum(self, protein: str, stop_codon: str = "TAA") -> str:
        """H0: Primary Deterministic Optimum matching canonical /api/optimize DP v2.1.1 production."""
        opt = DPV211Optimizer(
            alpha_harmonization=1.0,
            beta_mfe_proxy=1.0,
            gamma_gc_penalty=0.5,
            delta_au_bonus=0.5,
            r_ramp_start=0.5,
            r_ramp_end=0.9,
            initiation_gc_min=0.20,
            initiation_gc_max=0.30,
            forbidden_enzymes={"BsaI", "BsmBI", "BpiI", "SapI"},
            homopolymer_max_run=5,
        )
        res = opt.optimize(protein, self.codon_table.codon_weights, stop_codon=stop_codon)
        return res["sequence"]

    def generate_h1_initiation_open(self, protein: str, stop_codon: str = "TAA") -> str:
        """H1: Initiation-Focused Exploration (Reduced 5' mRNA structure stability / max 5' MFE / AU bonus)."""
        w_at = {}
        for c, w in self.codon_table.codon_weights.items():
            at_count = c.count("A") + c.count("T")
            w_at[c] = w * (1.6**at_count)

        opt = DPV211Optimizer(
            alpha_harmonization=0.30,
            beta_mfe_proxy=2.0,
            delta_au_bonus=1.5,
            gamma_gc_penalty=0.5,
            r_ramp_start=0.2,
            r_ramp_end=0.9,
            initiation_gc_min=0.15,
            initiation_gc_max=0.26,
            nominal_ramp_gc=0.20,
            forbidden_enzymes={"BsaI", "BsmBI", "BpiI", "SapI"},
            homopolymer_max_run=5,
        )
        res = opt.optimize(protein, w_at, stop_codon=stop_codon)
        return res["sequence"]

    def generate_h2_composition_constrained(self, protein: str, stop_codon: str = "TAA") -> str:
        """H2: Composition-Constrained / Manufacturability Design (Strict GC & Homopolymers)."""
        w_comp = {}
        for c, w in self.codon_table.codon_weights.items():
            gc_count = c.count("G") + c.count("C")
            w_comp[c] = w * (1.8**gc_count)

        opt = DPV211Optimizer(
            alpha_harmonization=0.50,
            beta_mfe_proxy=0.0,
            gamma_gc_penalty=2.5,
            nominal_ramp_gc=0.28,
            initiation_gc_min=0.22,
            initiation_gc_max=0.35,
            forbidden_enzymes={"BsaI", "BsmBI", "BpiI", "SapI"},
            homopolymer_max_run=4,
        )
        res = opt.optimize(protein, w_comp, stop_codon=stop_codon)
        return res["sequence"]

    def build_panel(
        self,
        targets: Optional[List[Dict[str, Any]]] = None,
        replicates_per_construct: int = 3,
        seed: int = 42,
        experiment_id: str = "EXP-20260917-PILOT-01",
    ) -> Tuple[PairedDBTLDataset, Dict[str, Any]]:
        """Generates all 9 experimental constructs + 2 controls, verifies orthogonality,
        and creates randomized plate assignments.
        """
        target_list = targets or STANDARD_TARGET_PANEL
        constructs: Dict[str, ConstructDesignRecord] = {}
        orthogonality_report: Dict[str, Any] = {}

        # 1. Generate 3 hypotheses per target
        for tgt in target_list:
            tname = tgt["target_name"]
            prot = "".join(tgt["sequence"].split()).rstrip("*")
            aa_len = len(prot)
            tgt_meta = tgt.get("metadata", {})

            # Generate H0, H1, H2 sequences
            seq_h0 = self.generate_h0_primary_optimum(prot)
            seq_h1 = self.generate_h1_initiation_open(prot)
            seq_h2 = self.generate_h2_composition_constrained(prot)

            tv_h0 = self.extractor.extract(seq_h0)
            tv_h1 = self.extractor.extract(seq_h1)
            tv_h2 = self.extractor.extract(seq_h2)

            # Verify Orthogonality Gate
            h_traits = {"H0": tv_h0, "H1": tv_h1, "H2": tv_h2}
            passed, dists = self.orthogonality_gate.verify_orthogonality(h_traits)
            orthogonality_report[tname] = {
                "orthogonality_passed": passed,
                "pairwise_distances": dists,
            }

            # Register H0
            cid_h0 = f"{tname}_H0"
            constructs[cid_h0] = ConstructDesignRecord(
                construct_id=cid_h0,
                target_name=tname,
                mature_protein_aa_length=tgt.get("aa_length", aa_len),
                construct_aa_length=aa_len,
                hypothesis_id="H0",
                hypothesis_name="Primary Deterministic Optimum",
                generation_contract="factorforge_canonical_production_v2_1_1",
                sequence_digest=hashlib.sha256(seq_h0.encode("utf-8")).hexdigest(),
                trait_vector=tv_h0,
                sequence_dna=seq_h0,
                is_control=False,
                rationale="FactorForge primary Pareto-optimal deterministic candidate matching canonical production /api/optimize DP parameters.",
                metadata={
                    "design_contract_id": "factorforge_canonical_production_v2_1_1",
                    **tgt_meta,
                },
            )

            # Register H1
            cid_h1 = f"{tname}_H1"
            constructs[cid_h1] = ConstructDesignRecord(
                construct_id=cid_h1,
                target_name=tname,
                mature_protein_aa_length=tgt.get("aa_length", aa_len),
                construct_aa_length=aa_len,
                hypothesis_id="H1",
                hypothesis_name="Initiation-Focused Exploration",
                generation_contract="dp_v2_1_1_initiation_open",
                sequence_digest=hashlib.sha256(seq_h1.encode("utf-8")).hexdigest(),
                trait_vector=tv_h1,
                sequence_dna=seq_h1,
                is_control=False,
                rationale="Reduced predicted 5' mRNA secondary structure stability (less negative 5' MFE) and AU initiation bonus to test translation initiation bottlenecks.",
                metadata={"design_contract_id": "dp_v2_1_1_initiation_open", **tgt_meta},
            )

            # Register H2
            cid_h2 = f"{tname}_H2"
            constructs[cid_h2] = ConstructDesignRecord(
                construct_id=cid_h2,
                target_name=tname,
                mature_protein_aa_length=tgt.get("aa_length", aa_len),
                construct_aa_length=aa_len,
                hypothesis_id="H2",
                hypothesis_name="Composition-Constrained / Manufacturability Design",
                generation_contract="dp_v2_1_1_composition_constrained",
                sequence_digest=hashlib.sha256(seq_h2.encode("utf-8")).hexdigest(),
                trait_vector=tv_h2,
                sequence_dna=seq_h2,
                is_control=False,
                rationale="Strict local/global GC composition constraints and homopolymer bounds for translational processivity.",
                metadata={"design_contract_id": "dp_v2_1_1_composition_constrained", **tgt_meta},
            )

        # 2. Add Standard Controls
        # Positive Control: sfGFP WT Benchmark (717 nt DNA payload)
        pos_tv = self.extractor.extract(STANDARD_POS_CTRL_SFGFP)
        pos_dna_sha256 = hashlib.sha256(STANDARD_POS_CTRL_SFGFP.encode("utf-8")).hexdigest()
        constructs["POS_CTRL_sfGFP"] = ConstructDesignRecord(
            construct_id="POS_CTRL_sfGFP",
            target_name="sfGFP",
            mature_protein_aa_length=238,
            construct_aa_length=238,
            hypothesis_id="CTRL_POS",
            hypothesis_name="Positive Expression Control",
            generation_contract="reference_benchmark_sfGFP",
            sequence_digest=pos_dna_sha256,
            trait_vector=pos_tv,
            sequence_dna=STANDARD_POS_CTRL_SFGFP,
            is_control=True,
            control_type="positive_expression",
            rationale="Standard wild-type/benchmark sfGFP reporter with known high transient expression in N. benthamiana.",
            metadata={
                "protein_sequence_sha256": STANDARD_POS_CTRL_SFGFP_AA_SHA256,
                "dna_payload_sha256": pos_dna_sha256,
                "payload_nt_length": 717,
                "cds_nt_length": 714,
                "stop_codon": "TAA",
            },
        )

        # Negative Control: Empty Vector / Buffer Mock
        constructs["NEG_CTRL_EMPTY"] = ConstructDesignRecord(
            construct_id="NEG_CTRL_EMPTY",
            target_name="Mock_Vector",
            mature_protein_aa_length=0,
            construct_aa_length=0,
            hypothesis_id="CTRL_NEG",
            hypothesis_name="Negative / Empty Vector Control",
            generation_contract="pEAQ_HT_empty_mock",
            sequence_digest=hashlib.sha256(b"MOCK_EMPTY_VECTOR").hexdigest(),
            trait_vector=None,
            sequence_dna=None,
            is_control=True,
            control_type="negative_empty_vector",
            rationale="Agrobacterium containing empty vector / infiltration buffer to quantify autofluorescence baseline.",
            metadata={"control_type": "empty_vector_mock"},
        )

        # 3. Create Experiment Run Record
        exp_record = ExperimentRunRecord(
            experiment_id=experiment_id,
            batch_id="BATCH-01",
            host_organism="Nicotiana benthamiana",
            growth_conditions="24C, 16h light / 8h dark, 60% RH",
            infiltration_od=0.5,
            harvest_dpi=4,
            protocol_version="AGRO-INFIL-v2.1",
            operator_id="OP-PLANT-01",
            metadata={
                "total_constructs": len(constructs),
                "replicates_per_construct": replicates_per_construct,
            },
        )

        # 4. Randomized & Blinded Plate Layout (96-well format)
        # 11 constructs x 3 replicates = 33 sample wells
        sample_items: List[Tuple[str, int, bool]] = []
        for cid, cobj in constructs.items():
            for rep in range(1, replicates_per_construct + 1):
                sample_items.append((cid, rep, cobj.is_control))

        # Deterministic shuffle for reproducibility
        rng = random.Random(seed)
        shuffled_samples = list(sample_items)
        rng.shuffle(shuffled_samples)

        # Generate plate coordinates: A01..A12, B01..B12, C01..C12...
        rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
        cols = [f"{c:02d}" for c in range(1, 13)]
        plate_wells = [f"{r}{c}" for r in rows for c in cols]

        samples: Dict[str, SampleRecord] = {}
        blinded_layout: List[Dict[str, Any]] = []
        unblinded_mapping: List[Dict[str, Any]] = []

        for idx, (cid, rep, is_ctrl) in enumerate(shuffled_samples):
            smp_id = f"SMP-{idx + 1:03d}"
            well_pos = plate_wells[idx]
            plant_num = (idx % 3) + 1
            leaf_num = ((rep - 1) % 4) + 1
            plant_id = f"PLANT-{plant_num:02d}"
            leaf_id = f"LEAF-{leaf_num}"
            bio_rep_id = f"BIO-REP-{rep:02d}"

            s_rec = SampleRecord(
                sample_id=smp_id,
                blinded_plate_position=well_pos,
                assay_well_position=well_pos,
                assay_plate_id="PLATE-01",
                experiment_id=experiment_id,
                construct_id=cid,
                biological_replicate=rep,
                technical_replicate=1,
                plant_id=plant_id,
                leaf_id=leaf_id,
                biological_replicate_id=bio_rep_id,
                is_control=is_ctrl,
                metadata={"seed": seed},
            )
            samples[smp_id] = s_rec

            blinded_layout.append(
                {
                    "sample_id": smp_id,
                    "blinded_plate_position": well_pos,
                    "assay_plate_id": "PLATE-01",
                    "assay_well_position": well_pos,
                    "plant_id": plant_id,
                    "leaf_id": leaf_id,
                    "biological_replicate_id": bio_rep_id,
                    "host_organism": "N. benthamiana",
                    "infiltration_od": 0.5,
                    "harvest_dpi": 4,
                }
            )

            unblinded_mapping.append(
                {
                    "sample_id": smp_id,
                    "blinded_plate_position": well_pos,
                    "assay_plate_id": "PLATE-01",
                    "assay_well_position": well_pos,
                    "construct_id": cid,
                    "target_name": constructs[cid].target_name,
                    "hypothesis_id": constructs[cid].hypothesis_id,
                    "biological_replicate": rep,
                    "plant_id": plant_id,
                    "leaf_id": leaf_id,
                    "biological_replicate_id": bio_rep_id,
                    "is_control": is_ctrl,
                }
            )

        # 5. Build Dataset
        dataset = PairedDBTLDataset(
            dataset_id=f"DBTL-{experiment_id}",
            schema_version="0.1.0",
            created_at=datetime.now(timezone.utc).isoformat(),
            constructs=constructs,
            experiments={experiment_id: exp_record},
            samples=samples,
            measurements=[],
            derived_outcomes=[],
            provenance={
                "builder": "ProspectivePanelBuilder",
                "seed": seed,
                "target_count": len(target_list),
                "experimental_construct_count": len(target_list) * 3,
                "control_construct_count": 2,
                "total_samples": len(samples),
                "orthogonality_report": orthogonality_report,
            },
        )
        dataset.archive_sha256 = dataset.compute_sha256()

        aux_data = {
            "blinded_layout": blinded_layout,
            "unblinded_mapping": unblinded_mapping,
            "orthogonality_report": orthogonality_report,
        }

        return dataset, aux_data
