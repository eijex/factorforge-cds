import pytest

# v1 is archived/frozen — no active maintenance.
# This test requires v1 extras (pip install factorforge[v1]) and the compare script.
pytest.skip("v1 is archived/frozen — skip by default", allow_module_level=True)
