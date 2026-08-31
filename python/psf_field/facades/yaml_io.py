"""YAML loader isolated so PyYAML's missing types do not leak (C1B.8.6)."""

from __future__ import annotations

from psf_field.errors import InputError


def load_mapping(text: str) -> dict[str, object]:
    try:
        import yaml  # type: ignore[import-untyped]  # PyYAML ships no type hints
    except ImportError as exc:
        raise InputError(
            "PyYAML is required to read YAML sidecars",
            module="boundary",
        ) from exc
    loaded: object = yaml.safe_load(text)
    if not isinstance(loaded, dict):
        raise InputError("sidecar YAML must be an object", module="boundary")
    out: dict[str, object] = {}
    for key_obj, value in loaded.items():
        item: object = value
        out[str(key_obj)] = item
    return out
