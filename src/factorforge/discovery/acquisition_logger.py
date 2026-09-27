"""Acquisition Logger and Ingestion Engine for Paired DBTL Datasets (Job 283C).

Manages dual-layer export (sequence-free ScientificMemory vs. Secure Reproducibility Archive)
and automated unblinded ingestion of wet-lab measurement outcomes.
"""

from __future__ import annotations
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List

from factorforge.discovery.acquisition import (
    DerivedOutcomeRecord,
    MeasurementRecord,
    PairedDBTLDataset,
)


class AcquisitionLogger:
    """Manages artifact persistence, memory export, and wet-lab ingestion."""

    def __init__(self, base_output_dir: str | Path) -> None:
        self.base_output_dir = Path(base_output_dir)
        self.base_output_dir.mkdir(parents=True, exist_ok=True)

    def export_prospective_panel(
        self,
        dataset: PairedDBTLDataset,
        aux_data: Dict[str, Any],
        subfolder: str = "prospective_panel_v3.5",
    ) -> Path:
        """Exports synthesis order manifests, FASTA, blinded plate layout,

        and sequence-free ScientificMemory package.
        """
        pkg_dir = self.base_output_dir / subfolder
        pkg_dir.mkdir(parents=True, exist_ok=True)

        # 1. Synthesis FASTA
        fasta_path = pkg_dir / "panel_sequences.fasta"
        with fasta_path.open("w", encoding="utf-8") as f:
            for cid, cobj in dataset.constructs.items():
                if cobj.sequence_dna:
                    header = f">{cid}|{cobj.target_name}|{cobj.hypothesis_id}|{cobj.generation_contract}|sha256={cobj.sequence_digest[:12]}"
                    f.write(f"{header}\n{cobj.sequence_dna}\n")

        # 2. Synthesis CSV Manifest
        csv_path = pkg_dir / "synthesis_manifest.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "construct_id",
                    "target_name",
                    "hypothesis_id",
                    "hypothesis_name",
                    "length_bp",
                    "sha256_digest",
                    "is_control",
                    "sequence_dna",
                ]
            )
            for cid, cobj in dataset.constructs.items():
                seq = cobj.sequence_dna or ""
                writer.writerow(
                    [
                        cid,
                        cobj.target_name,
                        cobj.hypothesis_id,
                        cobj.hypothesis_name,
                        len(seq),
                        cobj.sequence_digest,
                        cobj.is_control,
                        seq,
                    ]
                )

        # 3. Blinded Plate Layout (Laboratory Operator Document)
        blinded_path = pkg_dir / "blinded_plate_layout.json"
        with blinded_path.open("w", encoding="utf-8") as f:
            json.dump(
                {
                    "experiment_id": list(dataset.experiments.keys())[0],
                    "description": "Blinded prospective infiltration plate layout (96-well format, N=3 replicates)",
                    "total_wells": len(aux_data.get("blinded_layout", [])),
                    "samples": aux_data.get("blinded_layout", []),
                },
                f,
                indent=2,
            )

        # 4. Unblinded Mapping Manifest (For secure ingestion)
        unblinded_path = pkg_dir / "unblinded_mapping.json"
        with unblinded_path.open("w", encoding="utf-8") as f:
            json.dump(
                {
                    "experiment_id": list(dataset.experiments.keys())[0],
                    "unblinded_mapping": aux_data.get("unblinded_mapping", []),
                    "orthogonality_report": aux_data.get("orthogonality_report", {}),
                },
                f,
                indent=2,
            )

        # 5. Sequence-Free ScientificMemory Feature Store
        memory_path = pkg_dir / "scientific_memory_panel.json"
        with memory_path.open("w", encoding="utf-8") as f:
            json.dump(dataset.to_dict(include_sequence=False), f, indent=2)

        # 6. Full Composite Paired Dataset
        dataset_path = pkg_dir / "prospective_dbtl_dataset.json"
        with dataset_path.open("w", encoding="utf-8") as f:
            json.dump(dataset.to_dict(include_sequence=True), f, indent=2)

        # 7. Secure Archive Manifest (SHA-256 content-addressed index)
        def _calc_sha256(p: Path) -> str:
            h = hashlib.sha256()
            with p.open("rb") as f_in:
                h.update(f_in.read())
            return h.hexdigest()

        manifest = {
            "schema_version": "eijex.secure_reproducibility_archive.v0.1",
            "package_id": f"PKG-{dataset.dataset_id}",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "archive_integrity": "SHA-256 tamper-evident",
            "files": [
                {
                    "file_name": "panel_sequences.fasta",
                    "sha256": _calc_sha256(fasta_path),
                    "contains_sequence": True,
                },
                {
                    "file_name": "synthesis_manifest.csv",
                    "sha256": _calc_sha256(csv_path),
                    "contains_sequence": True,
                },
                {
                    "file_name": "blinded_plate_layout.json",
                    "sha256": _calc_sha256(blinded_path),
                    "contains_sequence": False,
                },
                {
                    "file_name": "unblinded_mapping.json",
                    "sha256": _calc_sha256(unblinded_path),
                    "contains_sequence": False,
                },
                {
                    "file_name": "scientific_memory_panel.json",
                    "sha256": _calc_sha256(memory_path),
                    "contains_sequence": False,
                },
                {
                    "file_name": "prospective_dbtl_dataset.json",
                    "sha256": _calc_sha256(dataset_path),
                    "contains_sequence": True,
                },
            ],
        }
        manifest_path = pkg_dir / "secure_archive_manifest.json"
        with manifest_path.open("w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2)

        return pkg_dir

    @staticmethod
    def ingest_wet_lab_measurements(
        dataset: PairedDBTLDataset,
        raw_measurements: List[Dict[str, Any]],
        lod_threshold: float = 0.5,
        loq_threshold: float = 2.0,
    ) -> PairedDBTLDataset:
        """Ingests empirical assay measurements, computes normalized/relative outcomes,

        and updates dataset with new measurements and derived outcomes.
        """
        # Register measurements
        measurements: List[MeasurementRecord] = []
        for m in raw_measurements:
            m_rec = MeasurementRecord(
                measurement_id=m["measurement_id"],
                sample_id=m["sample_id"],
                assay_type=m.get("assay_type", "GFP_FLUORESCENCE"),
                raw_yield_ug_g_fw=m.get("raw_yield_ug_g_fw"),
                western_blot_intensity=m.get("western_blot_intensity"),
                solubility_ratio=m.get("solubility_ratio"),
                chlorosis_score=m.get("chlorosis_score"),
                qc_status=m.get("qc_status", "PASS"),
                flags=m.get("flags", []),
                raw_data_ref=m.get("raw_data_ref"),
                measured_at=m.get("measured_at", datetime.now(timezone.utc).isoformat()),
            )
            measurements.append(m_rec)

        dataset.measurements.extend(measurements)

        # Build lookup for positive control and H0 yields to compute relative fold-change
        # First compute average yields per construct
        construct_yields: Dict[str, List[float]] = {}
        for m in measurements:
            if m.qc_status != "PASS" or m.raw_yield_ug_g_fw is None:
                continue
            s_rec = dataset.samples.get(m.sample_id)
            if s_rec:
                construct_yields.setdefault(s_rec.construct_id, []).append(m.raw_yield_ug_g_fw)

        mean_construct_yields = {
            cid: sum(vals) / len(vals) for cid, vals in construct_yields.items() if vals
        }

        # Average positive control yield
        pos_ctrl_mean = mean_construct_yields.get("POS_CTRL_sfGFP", 1.0)
        if pos_ctrl_mean <= 0.0:
            pos_ctrl_mean = 1.0

        # Compute derived outcomes
        derived_outcomes: List[DerivedOutcomeRecord] = []
        for m in measurements:
            s_rec = dataset.samples.get(m.sample_id)
            if not s_rec:
                continue

            raw_y = m.raw_yield_ug_g_fw if m.raw_yield_ug_g_fw is not None else 0.0
            cid = s_rec.construct_id
            cobj = dataset.constructs.get(cid)
            target_name = cobj.target_name if cobj else ""

            # Determine LOD/LOQ
            if raw_y < lod_threshold:
                lod_status = "BELOW_LOD"
            elif raw_y < loq_threshold:
                lod_status = "BETWEEN_LOD_LOQ"
            else:
                lod_status = "ABOVE_LOQ"

            # Relative to H0 mean
            h0_id = f"{target_name}_H0"
            h0_mean = mean_construct_yields.get(h0_id)
            rel_h0 = (raw_y / h0_mean) if (h0_mean and h0_mean > 0) else None

            # Relative to positive control mean
            rel_pos = (raw_y / pos_ctrl_mean) if pos_ctrl_mean > 0 else None

            out_rec = DerivedOutcomeRecord(
                outcome_id=f"OUT-{m.measurement_id}",
                sample_id=m.sample_id,
                construct_id=cid,
                experiment_id=s_rec.experiment_id,
                normalized_yield_ug_g_fw=raw_y,
                relative_to_h0=round(rel_h0, 4) if rel_h0 is not None else None,
                relative_to_pos_ctrl=round(rel_pos, 4) if rel_pos is not None else None,
                lod_loq_status=lod_status,
                notes=f"Assay: {m.assay_type}; QC: {m.qc_status}",
            )
            derived_outcomes.append(out_rec)

        dataset.derived_outcomes.extend(derived_outcomes)
        dataset.archive_sha256 = dataset.compute_sha256()

        return dataset
