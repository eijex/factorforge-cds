from factorforge import __version__
from factorforge.engines import EngineRegistry, register_builtin_engines
from factorforge.registry.versioning import engine_version, public_version_metadata


def test_product_and_engine_versions_are_independently_declared() -> None:
    metadata = public_version_metadata()

    assert __version__ == "3.5.0"
    assert metadata["product"]["release_status"] == "release_candidate"
    assert metadata["engines"]["profile"]["generation"] == 1
    assert metadata["engines"]["dp"]["version"] == engine_version("dp")
    assert metadata["engines"]["slm"]["status"] == "research_preview"
    assert metadata["engines"]["slm"]["trained_model_available"] is False


def test_registered_deterministic_engines_use_manifest_versions() -> None:
    register_builtin_engines()

    assert EngineRegistry.get("profile").version == engine_version("profile")
    assert EngineRegistry.get("dp").version == engine_version("dp")
