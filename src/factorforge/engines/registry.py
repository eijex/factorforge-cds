"""
Engine Registry

Registry that dynamically registers and manages engines
No changes needed to existing code when adding new engines (v3, v4, etc.)
"""

from __future__ import annotations

from typing import Callable, Type

from factorforge.core.interfaces import OptimizerEngine


class EngineRegistry:
    """Optimization engine registry"""

    _engines: dict[str, Type[OptimizerEngine]] = {}
    _instances: dict[str, OptimizerEngine] = {}
    _lazy_loaders: dict[str, Callable[[], Type[OptimizerEngine]]] = {}

    @classmethod
    def register(cls, name: str, engine_class: Type[OptimizerEngine]) -> None:
        """
        Register an engine

        Args:
            name: Engine identifier (e.g., "v1", "v2", "v3")
            engine_class: Class implementing OptimizerEngine
        """
        cls._engines[name] = engine_class
        print(f"Registered engine: {name} ({engine_class.__name__})")

    @classmethod
    def register_lazy(cls, name: str, loader: Callable[[], Type[OptimizerEngine]]) -> None:
        """
        Register a lazy engine loader.

        Args:
            name: Engine identifier (e.g., "v1", "v3")
            loader: Callable that returns the engine class on demand
        """
        cls._lazy_loaders[name] = loader

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
            cls.register(name, engine_class)

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
    def clear(cls) -> None:
        """Reset registry (for tests)"""
        cls._engines.clear()
        cls._instances.clear()
        cls._lazy_loaders.clear()
