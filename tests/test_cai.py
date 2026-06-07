from factorforge.analysis.metrics import calculate_cai


def test_cai_uniform_weights_is_one():
    # All weights 1.0 -> geometric mean of 1.0 = 1.0
    weights = {c1+c2+c3: 1.0 for c1 in "ACGT" for c2 in "ACGT" for c3 in "ACGT"}
    assert round(calculate_cai("ATGAAA", weights), 6) == 1.0


def test_cai_known_two_codon_case():
    # weights: ATG=1.0, AAA=0.25 -> geomean = sqrt(1*0.25) = 0.5
    weights = {"ATG": 1.0, "AAA": 0.25}
    assert round(calculate_cai("ATGAAA", weights), 6) == 0.5
