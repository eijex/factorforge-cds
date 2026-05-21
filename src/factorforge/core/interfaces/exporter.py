"""Exporter interface"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any


class Exporter(ABC):
    @abstractmethod
    def export(self, data: Any, format: str) -> str:
        """Export data to the requested format."""
        raise NotImplementedError
