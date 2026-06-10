# Dataset Provenance — nbenthamiana_reference_cds

## Source

- **Organism**: *Nicotiana benthamiana* Domin — NCBI TaxID: 4100
- **Assembly**: SGN QLD183 v1.0.3 (Sol Genomics Network)
- **Annotation version**: NbQld183.v103
- **Primary file**: `NbQld183.v103.gff3.CDS.fasta.gz`
- **Source URL**: `https://solgenomics.net/ftp/genomes/Nicotiana_benthamiana/QLD183/`
- **Note**: `NbQld183.v103.gff3.CDS.fasta.AA.fasta.gz` was 0 bytes at time of download — do not use. Protein sequences are derived by translating filtered CDS.

## Fetch procedure

Run `fetch_dataset.py --file <local_path>` (download separately; raw data is not committed).
Checksum of produced files is written to `checksums.txt`.

## Filtering rules (applied by `fetch_dataset.py`)

1. Standard DNA alphabet only (ACGT — no ambiguous bases)
2. Length is a multiple of 3
3. No internal stop codon
4. Non-empty protein after stop codon removal
5. Deduplicated by exact sequence SHA-256 hash

## Isoform / representative transcript policy

Filtering is by **sequence hash**, not by gene ID. Multiple isoforms from the same gene
are all retained if their CDS sequences differ. No representative-transcript selection
is applied. Benchmark metrics therefore reflect unique sequences, not unique genes.

## Statistics (v3.2.0 formal run, 2026-06-11)

- Input records: 61,328 (gene model count per SGN annotation)
- After filtering: **49,257 unique CDS sequences**
- CDS FASTA sha256: `714a7155c50fff4240c196389ade1860550eb5b949f3cbbd44406ddd6d6cdb53`
- Protein FASTA sha256: `35dc8ec1f2adaea2be02ad843dd18ca62053595328fe35e38c04bfa596f188e2`

## License / data use

SGN QLD183 data is publicly available. Users must verify current SGN data agreement
before using this dataset in publications or redistributing any subset.
Raw FASTA files are **not committed** to this repository (see `.gitignore`).
Only checksums and this provenance record are tracked.

## Citation

When using this dataset, cite the SGN QLD183 genome assembly and annotation.
Reference: Sol Genomics Network — https://solgenomics.net/
