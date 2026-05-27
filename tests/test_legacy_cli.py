"""Tests for archived CLI compatibility behavior."""

import pytest

from factorforge.cli import legacy_cli


def test_legacy_cli_exits_with_archive_guidance(monkeypatch):
    monkeypatch.setattr("sys.argv", ["legacy_cli"])

    with pytest.raises(SystemExit) as exc_info:
        legacy_cli.main()

    assert "v1 CLI is archived" in str(exc_info.value)
    assert "factorforge optimize" in str(exc_info.value)
