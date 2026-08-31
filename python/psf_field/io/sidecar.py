"""ImageMeta sidecar merge (C1.4.2)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

from psf_field.errors import InputError

JsonObject = dict[str, object]


def _object_map(raw: dict[object, object]) -> JsonObject:
    out: JsonObject = {}
    for key, value in raw.items():
        item: object = value
        out[str(key)] = item
    return out


def load_sidecar(path: Path) -> JsonObject:
    suffix = path.suffix.lower()
    text = path.read_text(encoding="utf-8")
    if suffix in {".yaml", ".yml"}:
        from psf_field.facades.yaml_io import load_mapping

        return load_mapping(text)
    loaded: object = json.loads(text)
    if not isinstance(loaded, dict):
        raise InputError("sidecar JSON must be an object", module="boundary")
    return _object_map(loaded)


def merge_sidecar(
    base: Mapping[str, object],
    sidecar: Mapping[str, object],
) -> JsonObject:
    """Overlay sidecar keys onto FITS-derived fields. No silent NOR.13 defaults."""
    merged = dict(base)
    merged.update(dict(sidecar))
    return merged
