# Dataset Provenance — nbenthamiana_reference_cds

- **Source:** Sol Genomics Network (SGN) / NCBI N. benthamiana CDS (public).
- **Fetched by:** `fetch_dataset.py` (not committed to git — see .gitignore).
- **Download date:** filled by fetch_dataset.py run log.
- **Filters:** length multiple of 3; no internal stop; standard ACGT only; duplicates removed by sequence hash.
- **Checksum:** SHA256 of the assembled FASTA, written to `checksums.txt` on fetch.
- **License:** verify SGN/NCBI redistribution terms before committing any subset; default = do not commit raw data.
