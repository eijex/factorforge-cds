from factorforge.analysis.metrics import calculate_gc


def test_gc_all_gc():
    assert calculate_gc("GCGC") == 100.0


def test_gc_half():
    assert calculate_gc("ATGC") == 50.0


def test_gc_none():
    assert calculate_gc("ATAT") == 0.0
