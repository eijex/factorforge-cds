"""Regression checks for unsupported objective options and diagnostic claim scope."""

import pytest

from factorforge.engines.dp_v2_4 import DPV240Optimizer, TransferMatrixFeasibilityGuard
from factorforge.engines.sllm.automaton import AutomatonCompiler


@pytest.mark.parametrize("value", [1.0, -1.0, float("nan"), float("inf"), True])
@pytest.mark.parametrize("name", ["lambda_tasep_shock", "lambda_ot_cadence"])
def test_unused_objectives_rejected(name, value):
    with pytest.raises(ValueError):
        DPV240Optimizer(**{name: value})


def test_prefix_proxy_is_explicit():
    result = TransferMatrixFeasibilityGuard.check_feasibility(
        "KEL", AutomatonCompiler.compile(["GGTCTC"], include_rc=True)
    )
    assert result["scope"] == "PREFIX_REACHABILITY_ONLY"
    assert result["entropy_is_calibrated"] is False
