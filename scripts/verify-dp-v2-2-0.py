"""Run the local synthetic DP v2.2 regression gates, propagating failures."""

from pathlib import Path
import subprocess
import sys


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    raise SystemExit(
        subprocess.call(
            [
                sys.executable,
                "-m",
                "pytest",
                "tests/test_dp_v2_2_0.py",
                "tests/test_dp_v2_1_1.py",
                "tests/test_stemcell_3dp.py",
                "tests/test_structural_pause.py",
                "tests/test_3dp_validation_suite.py",
                "-q",
            ],
            cwd=root,
        )
    )
