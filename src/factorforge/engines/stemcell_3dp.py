"""Experimental human-reference CDS design with provisional post-generation checks.

Historical StemCell names are retained for import compatibility. This is not a
cell-specific tRNA model, approved laboratory SOP, PKR assay, or yield predictor.
Thermodynamics are observations after DP, not an exact RNA-constrained DP objective.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE
from factorforge.constraints.type_iis import get_canonical_forbidden_motifs
from factorforge.engines.structural_pause import StructuralPauseDPOptimizer
from factorforge.rules.registry import RuleRegistry
from factorforge.sops.config import SopProfile

POLICY_PATH = Path(__file__).resolve().parents[1] / "data/profiles/human_stemcell_policy.json"

STEM_CELL_DEFAULT_FORBIDDEN_MOTIFS: list[str] = [
    # Type IIS enzymes (Golden Gate & MoClo compatibility)
    "GGTCTC",
    "GAGACC",  # BsaI / Eco31I
    "CGTCTC",
    "GAGACG",  # BsmBI / Esp3I
    "GAAGAC",
    "GTCTTC",  # BpiI / BbsI
    "CACCTGC",
    "GCAGGTG",  # PaqCI / AarI
    "GCTCTTC",
    "GAAGAGC",  # SapI
    # Common Mammalian Cloning MCS Sites
    "GAATTC",  # EcoRI
    "GGATCC",  # BamHI
    "AAGCTT",  # HindIII
    "GCGGCCGC",  # NotI
    "CTCGAG",  # XhoI
    "GCTAGC",  # NheI
    # Premature PolyA and Transcriptional Disruptions
    "AATAAA",
    "TTTATT",
    "ATTAAA",
    "TTTAAT",
    "AGTAAA",
    "TTTACT",
    "TATAAA",
    "TTTATA",  # TATA-box promoter mimic
    # Cryptic Splice Signals
    "GTAAGT",
    "ACTTAC",
    "GTGAGA",
    "TCTCAC",
    # Provisional exact-pattern exclusions; not a miRNA target or silencing predictor
    "AGCACTT",
    "AAGTGCT",
    "GTGCAAT",
    "ATTGCAC",
    "CAGTATT",
    "AATACTG",
    "ACTGGAA",
    "TTCCAGT",
    # Homopolymer Runs (>= 7nt)
    "AAAAAAA",
    "TTTTTTT",
    "CCCCCCC",
    "GGGGGGG",
]


def load_stemcell_codon_weights(codon_table_path: Path | None = None) -> dict[str, float]:
    """Normalize the supplied generic human reference; retain prior weight formula."""
    path = (
        codon_table_path or Path(__file__).resolve().parents[1] / "data/human_stemcell_codons.json"
    )
    data = json.loads(path.read_text(encoding="utf-8"))
    codons = data.get("codons", {})
    if set(codons) != set(STANDARD_GENETIC_CODE):
        raise ValueError("Codon reference must contain the complete standard genetic code")
    by_aa = {}
    for codon, info in codons.items():
        freq = info.get("frequency")
        if (
            info.get("aa") != STANDARD_GENETIC_CODE[codon]
            or isinstance(freq, bool)
            or not isinstance(freq, (int, float))
            or not math.isfinite(freq)
            or freq < 0
        ):
            raise ValueError("Invalid codon reference entry")
        by_aa.setdefault(info["aa"], []).append((codon, freq))
    weights = {}
    for items in by_aa.values():
        maximum = max(freq for _, freq in items)
        for codon, freq in items:
            weights[codon] = max(round(freq / maximum if maximum > 0 else 0.01, 4), 0.05)
    return weights


def calculate_max_duplex_stem(structure: str) -> int:
    """Count stacked pairs (i,j),(i+1,j-1); reject malformed dot-bracket input."""
    if not isinstance(structure, str) or set(structure) - set(".()"):
        raise ValueError("Invalid dot-bracket characters")
    stack, pairs = [], {}
    for i, char in enumerate(structure):
        if char == "(":
            stack.append(i)
        elif char == ")":
            if not stack:
                raise ValueError("Unbalanced dot-bracket structure")
            pairs[stack.pop()] = i
    if stack:
        raise ValueError("Unbalanced dot-bracket structure")
    runs = {}
    for i in sorted(pairs, reverse=True):
        runs[i] = 1 + (runs.get(i + 1, 0) if pairs.get(i + 1) == pairs[i] - 1 else 0)
    return max(runs.values(), default=0)


def evaluate_37c_thermodynamics(sequence: str) -> dict[str, Any]:
    """Fold sense CDS only at 37 C with isolated model details and validated output."""
    result = dict(
        status="ERROR",
        evidence_class="MODEL_PREDICTION",
        scope="sense_cds_only",
        temperature_c=37.0,
        mfe_37c=None,
        structure_37c=None,
        initiation_mfe_37c=None,
        min_local_mfe_40nt=None,
        max_duplex_stem_bp=None,
        pkr_assessment={"status": "NOT_AVAILABLE", "value": None},
    )
    clean = "".join(sequence.upper().split()) if isinstance(sequence, str) else ""
    if not clean or set(clean) - set("ACGT"):
        result["reason"] = "Nonempty unambiguous DNA is required"
        return result
    try:
        import RNA
    except ImportError:
        result.update(status="NOT_AVAILABLE", reason="ViennaRNA is not installed")
        return result
    except Exception:
        result.update(status="ERROR", reason="ViennaRNA could not be initialized")
        return result
    try:
        rna = clean.replace("T", "U")

        def fold(text):
            md = RNA.md()
            md.temperature = 37.0
            structure, energy = RNA.fold_compound(text, md).mfe()
            if (
                not isinstance(structure, str)
                or len(structure) != len(text)
                or isinstance(energy, bool)
                or not isinstance(energy, (int, float))
                or not math.isfinite(energy)
            ):
                raise ValueError("Invalid folding output")
            calculate_max_duplex_stem(structure)
            return structure, float(energy)

        struct, energy = fold(rna)
        _, initiation = fold(rna[:60])
        starts = list(range(0, len(rna) - 40 + 1, 5))
        # Include the final full window when a non-multiple-of-five tail remains.
        if len(rna) >= 40 and (not starts or starts[-1] != len(rna) - 40):
            starts.append(len(rna) - 40)
        local = [fold(rna[i : i + 40])[1] for i in starts]
        result.update(
            status="AVAILABLE",
            engine=f"ViennaRNA_{getattr(RNA, '__version__', 'unknown')}",
            mfe_37c=energy,
            structure_37c=struct,
            initiation_mfe_37c=initiation,
            initiation_window_nt=min(60, len(rna)),
            local_window_nt=40,
            local_step_nt=5,
            min_local_mfe_40nt=min(local) if local else None,
            local_window_status="AVAILABLE" if local else "NOT_AVAILABLE",
            max_duplex_stem_bp=calculate_max_duplex_stem(struct),
        )
    except Exception:
        result.update(status="ERROR", reason="Folding failed or returned invalid observations")
    return result


class StemCell3DPOptimizer:
    """Experimental generic-human prototype; human review remains required."""

    def __init__(
        self,
        forbidden_motifs=None,
        target_gc_min=0.44,
        target_gc_max=0.50,
        pause_positions=None,
        pause_target_adaptiveness=0.45,
        pause_penalty_weight=2.0,
        structure_reference=None,
        positions_reviewed=False,
        policy_path=None,
    ):
        if (
            any(
                isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v)
                for v in (target_gc_min, target_gc_max)
            )
            or not 0 <= target_gc_min <= target_gc_max <= 1
        ):
            raise ValueError("GC bounds must be finite fractions in [0, 1]")
        self.target_gc_min, self.target_gc_max = target_gc_min, target_gc_max
        self.codon_weights = load_stemcell_codon_weights()
        self.pause_positions = list(pause_positions or [])
        if self.pause_positions and (
            positions_reviewed is not True
            or not isinstance(structure_reference, str)
            or not structure_reference.strip()
        ):
            raise ValueError("Explicit reviewed positions and a structure reference are required")
        if len(self.pause_positions) != len(set(self.pause_positions)):
            raise ValueError("Duplicate structure positions")
        self.profile = SopProfile(
            policy_path or POLICY_PATH,
            known_rule_ids=[r.rule_id for r in RuleRegistry().list_rules()],
        )
        self.observation_policy = self.profile.data.get("experimental_observation_policy")
        if not isinstance(self.observation_policy, dict):
            raise ValueError("Experimental observation policy is required")
        for key in ("initiation_min_kcal_mol", "local_min_kcal_mol", "max_stacked_pairs_exclusive"):
            value = self.observation_policy.get(key)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise ValueError("Observation thresholds must be finite numbers")
        if (
            type(self.observation_policy["max_stacked_pairs_exclusive"]) is not int
            or self.observation_policy["max_stacked_pairs_exclusive"] < 1
        ):
            raise ValueError("Stacked-pair limit must be a positive integer")
        self.reference_path = (
            Path(__file__).resolve().parents[1] / "data/human_stemcell_codons.json"
        )
        self.reference_sha256 = hashlib.sha256(self.reference_path.read_bytes()).hexdigest()
        self.forbidden_motifs = sorted(
            set(
                get_canonical_forbidden_motifs()
                + list(
                    STEM_CELL_DEFAULT_FORBIDDEN_MOTIFS
                    if forbidden_motifs is None
                    else forbidden_motifs
                )
            )
        )
        if any(
            not isinstance(m, str) or not m or set(m) - set("ACGT") for m in self.forbidden_motifs
        ):
            raise ValueError("Motifs must be uppercase unambiguous DNA")
        self.optimizer = StructuralPauseDPOptimizer(
            positions=self.pause_positions,
            target_adaptiveness=pause_target_adaptiveness,
            penalty_weight=pause_penalty_weight if self.pause_positions else 0.0,
            structure_reference=structure_reference or "sequence_only_no_structural_evidence",
            forbidden_motifs=self.forbidden_motifs,
            include_reverse_complement=False,
        )

    def optimize(self, protein_sequence, left_flank="", right_flank=""):
        protein = "".join(protein_sequence.upper().split()).rstrip("*")
        result = self.optimizer.optimize(
            protein_sequence=protein,
            codon_weights=self.codon_weights,
            target_gc_min=self.target_gc_min,
            target_gc_max=self.target_gc_max,
            left_flank=left_flank,
            right_flank=right_flank,
        )
        return self.evaluate_candidate(result["sequence"], protein, left_flank, right_flank)

    def evaluate_candidate(self, dna, expected_protein, left_flank="", right_flank=""):
        """Independently check a candidate; no post-hoc modification or relaxation."""
        dna = "".join(dna.upper().split())
        protein = "".join(expected_protein.upper().split()).rstrip("*")
        if not dna or len(dna) % 3 or set(dna) - set("ACGT"):
            raise ValueError("Candidate must be a complete unambiguous sense CDS")
        translated = "".join(STANDARD_GENETIC_CODE[dna[i : i + 3]] for i in range(0, len(dna), 3))
        if translated != protein or "*" in translated:
            raise ValueError("Candidate translation mismatch")
        gc = (dna.count("G") + dna.count("C")) / len(dna)
        found = [
            m for m in self.forbidden_motifs if m in left_flank.upper() + dna + right_flank.upper()
        ]
        thermo = evaluate_37c_thermodynamics(dna)
        checks = {
            "translation_identity": True,
            "gc_in_band": self.target_gc_min <= gc <= self.target_gc_max,
            "motif_clean": not found,
        }
        for name, field, threshold, minimum in (
            ("initiation_mfe", "initiation_mfe_37c", "initiation_min_kcal_mol", True),
            ("local_mfe", "min_local_mfe_40nt", "local_min_kcal_mol", True),
            ("stacked_pair_limit", "max_duplex_stem_bp", "max_stacked_pairs_exclusive", False),
        ):
            value = thermo.get(field)
            available = (
                thermo.get("status") == "AVAILABLE"
                and isinstance(value, (int, float))
                and not isinstance(value, bool)
                and math.isfinite(value)
            )
            checks[name] = (
                (
                    (value >= self.observation_policy[threshold])
                    if minimum
                    else (value < self.observation_policy[threshold])
                )
                if available
                else None
            )
        policy = RuleRegistry(sop_profile=self.profile).evaluate_sequence(
            dna, host="human_stemcell"
        )
        checks["registered_policy"] = policy["all_passed"]
        computational_status = (
            "BLOCKED"
            if any(v is False for v in checks.values())
            else "INCOMPLETE"
            if any(v is None for v in checks.values())
            else "PASS"
        )
        return dict(
            protein_length=len(protein),
            dna_length_nt=len(dna),
            optimized_sequence=dna,
            cai=round(
                math.exp(
                    sum(math.log(self.codon_weights[dna[i : i + 3]]) for i in range(0, len(dna), 3))
                    / len(protein)
                ),
                4,
            ),
            gc_percent=round(gc * 100, 2),
            gc_in_band=checks["gc_in_band"],
            forbidden_motifs_found=found,
            motif_clean=not found,
            thermodynamics_37c=thermo,
            translation_identity=1.0,
            status="REVIEW_REQUIRED",
            computational_status=computational_status,
            checks=checks,
            registered_policy=policy,
            profile_digest=self.profile.digest,
            codon_reference_sha256=self.reference_sha256,
            reference_claim="Generic human-frequency reference; cell-specific calibration unverified",
            algorithm_mode="reviewed_position_conditioned"
            if self.pause_positions
            else "sequence_only",
            reviewed_positions=self.pause_positions,
            observation_policy=self.observation_policy,
            protein_sha256=hashlib.sha256(protein.encode()).hexdigest(),
            dna_sha256=hashlib.sha256(dna.encode()).hexdigest(),
            biological_validation={"status": "NOT_AVAILABLE", "value": None},
            relaxation_performed=False,
        )
