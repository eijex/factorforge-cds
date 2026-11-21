from factorforge import __version__
from factorforge.engines import EngineRegistry, register_builtin_engines
from factorforge.registry.versioning import (
    engine_family_version,
    engine_formulation_family,
    engine_version,
    public_version_metadata,
)


def test_product_and_engine_versions_are_independently_declared() -> None:
    metadata = public_version_metadata()

    assert __version__ == "3.6.8"
    assert metadata["product"]["release_status"] == "released"
    assert metadata["product"]["previous_stable_version"] == "3.6.7"
    assert metadata["engines"]["profile"]["generation"] == 1
    assert metadata["engines"]["dp"]["version"] == engine_version("dp")
    assert metadata["engines"]["dp_v2_2"]["version"] == "2.2.0"
    assert metadata["engines"]["dp_v2_3"]["version"] == "2.3.0"
    assert metadata["engines"]["local_rewriter"]["version"] == "1.0.0"
    assert metadata["engines"]["slm"]["status"] == "research_preview"
    assert metadata["engines"]["slm"]["trained_model_available"] is False


def test_registered_deterministic_engines_use_manifest_versions() -> None:
    register_builtin_engines()

    assert EngineRegistry.get("profile").version == engine_version("profile")
    assert EngineRegistry.get("dp").version == engine_version("dp")


def test_two_axis_formulation_family_and_independent_versions() -> None:
    assert engine_formulation_family("dp") == "1DP"
    assert engine_family_version("dp") == "2.0.1"
    assert engine_formulation_family("dp_v2_1_1") == "1DP"
    assert engine_family_version("dp_v2_1_1") == "2.1.1"
    assert engine_formulation_family("dp_v2_2") == "3DP"
    assert engine_family_version("dp_v2_2") == "1.0.0"
    assert engine_formulation_family("dp_v2_3") == "4DP"
    assert engine_family_version("dp_v2_3") == "1.0.0"
    assert engine_formulation_family("dp_v2_4") == "4DP"
    assert engine_family_version("dp_v2_4") == "1.1.0"
    assert engine_formulation_family("dp_5dp") == "5DP"
    assert engine_family_version("dp_5dp") == "0.1.0"
