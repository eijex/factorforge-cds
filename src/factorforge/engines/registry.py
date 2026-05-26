"""
Engine Registry

Registry that dynamically registers and manages optimization engines.
"""

from __future__ import annotations

import logging
from typing import Callable, Type

from factorforge.core.interfaces import OptimizerEngine

logger = logging.getLogger(__name__)


class EngineRegistry:
    """Optimization engine registry"""

    _engines: dict[str, Type[OptimizerEngine]] = {}
    _instances: dict[str, OptimizerEngine] = {}
    _lazy_loaders: dict[str, Callable[[], Type[OptimizerEngine]]] = {}
    _metadata: dict[str, dict[str, object]] = {}

    @classmethod
    def register(
        cls,
        name: str,
        engine_class: Type[OptimizerEngine],
        metadata: dict[str, object] | None = None,
    ) -> None:
        """
        Register an engine

        Args:
            name: Engine identifier such as "dp" or "profile"
            engine_class: Class implementing OptimizerEngine
            metadata: Optional engine metadata.
        """
        cls._engines[name] = engine_class
        if metadata is not None:
            cls._metadata[name] = dict(metadata)
        else:
            cls._metadata.setdefault(name, {})
        logger.debug("Registered engine: %s (%s)", name, engine_class.__name__)

    @classmethod
    def register_lazy(
        cls,
        name: str,
        loader: Callable[[], Type[OptimizerEngine]],
        metadata: dict[str, object] | None = None,
    ) -> None:
        """
        Register a lazy engine loader.

        Args:
            name: Engine identifier
            loader: Callable that returns the engine class on demand
            metadata: Optional engine metadata.
        """
        cls._lazy_loaders[name] = loader
        if metadata is not None:
            cls._metadata[name] = dict(metadata)
        else:
            cls._metadata.setdefault(name, {})

    @classmethod
    def get(cls, name: str) -> OptimizerEngine:
        """
        Get engine instance (singleton)

        Args:
            name: Engine identifier

        Returns:
            OptimizerEngine instance
        """
        if name not in cls._engines and name in cls._lazy_loaders:
            engine_class = cls._lazy_loaders[name]()
            cls.register(name, engine_class, metadata=cls._metadata.get(name))

        if name not in cls._engines:
            available = ", ".join(cls._engines.keys())
            raise ValueError(f"❌ Engine '{name}' not found. Available: {available}")

        # Singleton pattern
        if name not in cls._instances:
            cls._instances[name] = cls._engines[name]()

        return cls._instances[name]

    @classmethod
    def list_engines(cls) -> dict[str, dict[str, str]]:
        """
        List available engines

        Returns:
            dict: {name: {version, description}}
        """
        result: dict[str, dict[str, str]] = {}
        for name, engine_class in cls._engines.items():
            instance = cls.get(name)
            result[name] = {
                "version": instance.version,
                "name": instance.name,
            }
        for name in cls._lazy_loaders:
            if name in result:
                continue
            result[name] = {
                "version": "lazy",
                "name": "lazy (not loaded)",
            }
        return result

    @classmethod
    def list_with_metadata(cls) -> dict[str, dict[str, object]]:
        """List registered and lazy engines with metadata."""
        result: dict[str, dict[str, object]] = {}
        all_names = set(cls._engines) | set(cls._lazy_loaders) | set(cls._metadata)
        for name in sorted(all_names):
            metadata = dict(cls._metadata.get(name, {}))
            if name in cls._engines:
                instance = cls.get(name)
                metadata.setdefault("version", instance.version)
                metadata.setdefault("name", instance.name)
            elif name in cls._lazy_loaders:
                metadata.setdefault("version", "lazy")
                metadata.setdefault("name", "lazy (not loaded)")
            result[name] = metadata
        return result

    @classmethod
    def clear(cls) -> None:
        """Reset registry (for tests)"""
        cls._engines.clear()
        cls._instances.clear()
        cls._lazy_loaders.clear()
        cls._metadata.clear()
