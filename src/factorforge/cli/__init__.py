"""CLI package."""

from __future__ import annotations


def __getattr__(name: str):
    if name == "cli":
        from .main import cli

        return cli
    raise AttributeError(name)

__all__ = ["cli"]
