"""Safe synthetic benchmark proteins for sequence-design checks."""

BENCHMARK_PROTEINS: dict[str, str] = {
    "short_synthetic": "MSTNPKPQR",
    "reporter_like": (
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTT"
        "FSYGVQCFSRYPDHMK"
    ),
    "antibody_like_heavy_style": (
        "MEVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYY"
        "ADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR"
    ),
    "low_complexity": "MSSSSGGGGSSSSGGGGPPPPQQQQNNNNKKKK",
    "cysteine_rich": "MCCGCCGCCNPNCCTGCKCCTGCCGCC",
}
