from __future__ import annotations

import json

from Bio.Seq import Seq

from factorforge.analysis.metrics import load_codon_usage_table
from factorforge.engines.dp_5dp_wssv import (
    WSSV5DPOptimizer,
    assemble_wssv_vlp_construct,
    detect_n_glycosylation_motifs,
)


def test_detect_n_glycosylation_motifs_skips_proline() -> None:
    motifs = detect_n_glycosylation_motifs("MANTSKNPTM")
    assert len(motifs) == 1
    assert motifs[0]["sequon"] == "NTS"
    assert motifs[0]["position_1based"] == 3


def test_wssv_5dp_generates_topologies_and_jsonl() -> None:
    table = load_codon_usage_table()
    hbcag = "MDIDPYKEFGATVELLSFLPSDFFPSVRDLLDTASALYREALESPEHCSPHHTALRQAILCWGELMTLATWV"
    vp28 = "MDLSFTLSVVSAILAITAVIAVFIVIFRYHNTVTKTIETHTDNIETNMDENLRIPVTAEVGSGYFKMTDVSFDSD"

    assembled = assemble_wssv_vlp_construct(hbcag, vp28, topology="c_terminal")
    assert assembled["has_pr1a"] is True
    assert assembled["has_kdel"] is True
    assert len(assembled["pause_positions"]) == 3

    opt = WSSV5DPOptimizer()
    candidates = opt.generate_candidates(
        hbcag_seq=hbcag,
        vp28_seq=vp28,
        codon_weights=table.codon_weights,
        topologies=("c_terminal", "mir_loop"),
        n_variants_per_topology=1,
    )
    assert len(candidates) == 2
    for cand in candidates:
        assert cand["wet_lab_yield_mg_kg"] is None
        assert cand["design_tier"] == "5DP_WSSV_VLP"
        assert cand["biological_validation"] == "NOT_AVAILABLE"
        assert "GGTCTC" not in cand["cds_sequence"]
        assert str(Seq(cand["cds_sequence"]).translate(to_stop=True)) == cand["protein_sequence"]

    jsonl = opt.to_jsonl(candidates)
    parsed = [json.loads(line) for line in jsonl.strip().splitlines()]
    assert len(parsed) == 2
