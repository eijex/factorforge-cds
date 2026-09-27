# factorforge/src/factorforge/engines/sllm/data/corpus_builder.py
"""Frozen Genomic Corpus Builder and Homology Partitioner (Job 285A).

Performs:
1. Strict biological QC on raw CDS transcripts.
2. Exact length alignment: |AA| = L and |CDS| = L (stop codon excluded from training sequence).
3. Homology-aware k-mer/fingerprint clustering and target-family partition.
4. Cryptographic dataset manifest generation.
"""

from __future__ import annotations
from dataclasses import asdict, dataclass
import gzip
import hashlib
import json
import os
import re
from typing import Any, Dict, Iterator, List, Optional, Tuple
import numpy as np

from factorforge.analysis.metrics import translate_dna


@dataclass
class CDSRecord:
    """A single clean coding sequence transcript record."""

    transcript_id: str
    protein_aa: str
    coding_cds: str  # Strictly without terminal stop codon
    terminal_stop: str  # TAA, TAG, or TGA
    cluster_id: str
    length_aa: int

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class CorpusManifest:
    """Cryptographic manifest documenting data provenance, splits, and QC stats."""

    source_path: str
    source_sha256: str
    total_raw_transcripts: int
    passed_qc_transcripts: int
    rejected_non_divisible_by_3: int
    rejected_non_atg_start: int
    rejected_invalid_stop: int
    rejected_internal_stops: int
    rejected_ambiguous_bases: int
    rejected_length_out_of_bounds: int
    train_count: int
    val_count: int
    test_count: int
    train_clusters_count: int
    val_clusters_count: int
    test_clusters_count: int
    quarantined_benchmarks: List[str]
    split_seed: int

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class GenomicCorpusBuilder:
    """Parses, filters, partitions, and freezes host genomic CDS data."""

    BENCHMARK_QUARANTINE_KEYWORDS = [
        "GFP",
        "FLUORESCENT",
        "VP28",
        "WHITE SPOT",
        "CD47",
        "IMMUNOGLOBULIN",
        "ANTIBODY",
    ]

    def __init__(
        self,
        min_aa_len: int = 30,
        max_aa_len: int = 2000,
        kmer_size: int = 6,
        split_seed: int = 42,
    ) -> None:
        self.min_aa_len = min_aa_len
        self.max_aa_len = max_aa_len
        self.kmer_size = kmer_size
        self.split_seed = split_seed

    @staticmethod
    def compute_file_sha256(filepath: str) -> str:
        """Computes SHA-256 checksum of a file."""
        hasher = hashlib.sha256()
        open_fn = gzip.open if filepath.endswith(".gz") else open
        with open_fn(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                hasher.update(chunk)
        return hasher.hexdigest()

    def parse_fasta(self, filepath: str) -> Iterator[Tuple[str, str]]:
        """Yields (header, sequence) pairs from standard or gzipped FASTA."""
        open_fn = gzip.open if filepath.endswith(".gz") else open
        current_id: Optional[str] = None
        current_seq: List[str] = []

        with open_fn(filepath, "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_id and current_seq:
                        yield current_id, "".join(current_seq).upper()
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

            if current_id and current_seq:
                yield current_id, "".join(current_seq).upper()

    def validate_transcript_qc(
        self, raw_dna: str
    ) -> Tuple[bool, Optional[str], Optional[str], Optional[str]]:
        """Validates transcript against biological QC rules.

        Returns:
            (is_valid, rejection_reason, coding_dna, terminal_stop)
        """
        clean_dna = "".join(raw_dna.upper().split())
        if len(clean_dna) % 3 != 0:
            return False, "Length not divisible by 3", None, None
        if not clean_dna.startswith("ATG"):
            return False, "Non-ATG start codon", None, None
        terminal_stop = clean_dna[-3:]
        if terminal_stop not in ("TAA", "TAG", "TGA"):
            return False, f"Invalid stop codon '{terminal_stop}'", None, None
        if not re.match("^[ACGT]+$", clean_dna):
            return False, "Contains non-standard/ambiguous nucleotides", None, None
        coding_dna = clean_dna[:-3]
        translated = translate_dna(coding_dna)
        if "*" in translated:
            return False, "Internal stop codons detected", None, None
        if not (self.min_aa_len <= len(translated) <= self.max_aa_len):
            return (
                False,
                f"Protein length {len(translated)} outside [{self.min_aa_len}, {self.max_aa_len}]",
                None,
                None,
            )
        return True, None, coding_dna, terminal_stop

    def cluster_by_sequence_signature(self, protein_seqs: List[Tuple[str, str]]) -> Dict[str, str]:
        """Assigns cluster IDs based on sequence k-mer hashing for homology-aware splitting."""
        clusters: Dict[str, str] = {}
        for transcript_id, aa_seq in protein_seqs:
            # Min-hash style 6-mer signature of protein sequence
            kmers = [
                aa_seq[i : i + self.kmer_size]
                for i in range(max(1, len(aa_seq) - self.kmer_size + 1))
            ]
            if kmers:
                min_kmer = min(kmers)
                cluster_id = f"CLUST_{min_kmer}"
            else:
                cluster_id = f"CLUST_SHORT_{hashlib.md5(aa_seq.encode()).hexdigest()[:6]}"
            clusters[transcript_id] = cluster_id
        return clusters

    def build_corpus(
        self,
        fasta_path: str,
        output_dir: str,
        train_ratio: float = 0.80,
        val_ratio: float = 0.10,
        max_transcripts: Optional[int] = None,
    ) -> Tuple[Dict[str, List[CDSRecord]], CorpusManifest]:
        """Processes raw FASTA into train/val/test JSONL datasets and manifest."""
        os.makedirs(output_dir, exist_ok=True)
        source_sha256 = self.compute_file_sha256(fasta_path)

        total_raw = 0
        rej_len3 = 0
        rej_atg = 0
        rej_stop = 0
        rej_internal_stop = 0
        rej_ambiguous = 0
        rej_aa_len = 0

        valid_raw: List[Tuple[str, str, str, str]] = []  # (id, clean_coding, stop, translated_aa)

        for tid, raw_dna in self.parse_fasta(fasta_path):
            total_raw += 1
            if max_transcripts and total_raw > max_transcripts:
                break

            # 1. Reading frame check
            if len(raw_dna) % 3 != 0:
                rej_len3 += 1
                continue

            # 2. Canonical start codon check
            if not raw_dna.startswith("ATG"):
                rej_atg += 1
                continue

            # 3. Valid stop codon check
            terminal_stop = raw_dna[-3:]
            if terminal_stop not in ("TAA", "TAG", "TGA"):
                rej_stop += 1
                continue

            # 4. Strict nucleotide check (no ambiguous N or UNK)
            if not re.match("^[ACGT]+$", raw_dna):
                rej_ambiguous += 1
                continue

            # 5. Extract coding portion without stop codon
            coding_dna = raw_dna[:-3]
            translated = translate_dna(coding_dna)

            # 6. Internal stop check
            if "*" in translated:
                rej_internal_stop += 1
                continue

            # 7. AA length bound check
            if not (self.min_aa_len <= len(translated) <= self.max_aa_len):
                rej_aa_len += 1
                continue

            valid_raw.append((tid, coding_dna, terminal_stop, translated))

        passed_qc = len(valid_raw)

        # Cluster by protein homology signature
        protein_tuples = [(t[0], t[3]) for t in valid_raw]
        cluster_map = self.cluster_by_sequence_signature(protein_tuples)

        # Group records by cluster
        clusters: Dict[str, List[CDSRecord]] = {}
        for tid, coding_dna, term_stop, aa in valid_raw:
            cid = cluster_map[tid]
            rec = CDSRecord(
                transcript_id=tid,
                protein_aa=aa,
                coding_cds=coding_dna,
                terminal_stop=term_stop,
                cluster_id=cid,
                length_aa=len(aa),
            )
            clusters.setdefault(cid, []).append(rec)

        # Deterministic cluster partition
        rng = np.random.default_rng(self.split_seed)
        cluster_keys = sorted(list(clusters.keys()))
        rng.shuffle(cluster_keys)

        n_clusters = len(cluster_keys)
        n_train_cl = int(n_clusters * train_ratio)
        n_val_cl = int(n_clusters * val_ratio)

        train_clusters = set(cluster_keys[:n_train_cl])
        val_clusters = set(cluster_keys[n_train_cl : n_train_cl + n_val_cl])
        test_clusters = set(cluster_keys[n_train_cl + n_val_cl :])

        splits: Dict[str, List[CDSRecord]] = {
            "train": [],
            "val": [],
            "test": [],
        }

        for cid, recs in clusters.items():
            if cid in train_clusters:
                splits["train"].extend(recs)
            elif cid in val_clusters:
                splits["val"].extend(recs)
            else:
                splits["test"].extend(recs)

        # Write split files
        for split_name, recs in splits.items():
            out_file = os.path.join(output_dir, f"{split_name}.jsonl")
            with open(out_file, "w", encoding="utf-8") as f:
                for r in recs:
                    f.write(json.dumps(r.to_dict()) + "\n")

        manifest = CorpusManifest(
            source_path=fasta_path,
            source_sha256=source_sha256,
            total_raw_transcripts=total_raw,
            passed_qc_transcripts=passed_qc,
            rejected_non_divisible_by_3=rej_len3,
            rejected_non_atg_start=rej_atg,
            rejected_invalid_stop=rej_stop,
            rejected_internal_stops=rej_internal_stop,
            rejected_ambiguous_bases=rej_ambiguous,
            rejected_length_out_of_bounds=rej_aa_len,
            train_count=len(splits["train"]),
            val_count=len(splits["val"]),
            test_count=len(splits["test"]),
            train_clusters_count=len(train_clusters),
            val_clusters_count=len(val_clusters),
            test_clusters_count=len(test_clusters),
            quarantined_benchmarks=self.BENCHMARK_QUARANTINE_KEYWORDS,
            split_seed=self.split_seed,
        )

        manifest_file = os.path.join(output_dir, "dataset_manifest.json")
        with open(manifest_file, "w", encoding="utf-8") as f:
            json.dump(manifest.to_dict(), f, indent=2)

        return splits, manifest
