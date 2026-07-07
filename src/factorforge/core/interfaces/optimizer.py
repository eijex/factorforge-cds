"""
Optimizer Engine Interface

Interface that all optimization engines must implement.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any


class OptimizationResult:
    """Optimization result"""

    def __init__(
        self,
        sequence: str,
        metrics: dict[str, Any],
        metadata: dict[str, Any] | None = None,
    ) -> None:
        self.sequence = sequence
        self.metrics = metrics
        self.metadata = metadata or {}


class OptimizerEngine(ABC):
    """Abstract optimization engine interface"""

    @property
    @abstractmethod
    def name(self) -> str:
        """Engine name"""
        ...

    @property
    @abstractmethod
    def version(self) -> str:
        """Engine version"""
        ...

    @abstractmethod
    def optimize(
        self,
        sequence: str,
        profile: str | None = None,
        **kwargs: Any,
    ) -> OptimizationResult:
        """
        Optimize a sequence

        Args:
            sequence: Input protein sequence
            profile: Optimization profile (e.g., balanced, high_gc)
            **kwargs: Additional parameters

        Returns:
            OptimizationResult
        """
        ...

    @abstractmethod
    def validate(self, sequence: str) -> bool:
        """
        Validate input

        Args:
            sequence: Sequence to validate

        Returns:
            bool: True if valid
        """
        ...

    def get_metadata(self) -> dict[str, Any]:
        """Engine metadata"""
        return {
            "name": self.name,
            "version": self.version,
            "supported_profiles": self.get_supported_profiles(),
        }

    def get_supported_profiles(self) -> list[str]:
        """List of supported profiles"""
        return []
