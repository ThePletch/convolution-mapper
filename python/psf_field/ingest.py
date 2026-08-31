"""Ingest validators for Contract-1 records, ImageMeta sidecars, and pupil specs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

from pydantic import ValidationError

from psf_field.errors import InputError
from psf_field.io.fits_header import (
    header_to_image_meta_fields,
    read_primary_header,
)
from psf_field.io.sidecar import load_sidecar, merge_sidecar
from psf_field.models.catalog import Catalog
from psf_field.models.common import wrap_validation
from psf_field.models.inputs import ImageMeta, PupilSpec, StarRecord

FIELD_TOL_MM = 1e-6


def ingest_star_record(data: Mapping[str, object]) -> StarRecord:
    try:
        return StarRecord.model_validate(dict(data))
    except ValidationError as exc:
        raw_id = data.get("star_id")
        star_id = raw_id if isinstance(raw_id, str) else None
        raise wrap_validation(exc, star_id=star_id) from exc


def ingest_image_meta(data: Mapping[str, object]) -> ImageMeta:
    try:
        return ImageMeta.model_validate(dict(data))
    except ValidationError as exc:
        raise wrap_validation(exc) from exc


def ingest_pupil_spec(data: Mapping[str, object]) -> PupilSpec:
    try:
        return PupilSpec.model_validate(dict(data))
    except ValidationError as exc:
        raise wrap_validation(exc) from exc


def ingest_catalog(data: Mapping[str, object]) -> Catalog:
    try:
        return Catalog.model_validate(dict(data))
    except ValidationError as exc:
        raise wrap_validation(exc) from exc


def ingest_image_meta_from_fits(
    fits_path: Path,
    sidecar_path: Path | None = None,
) -> ImageMeta:
    """C1.4.2: merge optional sidecar over FITS keywords; reject incomplete NOR.13."""
    header = read_primary_header(fits_path)
    base = header_to_image_meta_fields(header)
    sidecar: Mapping[str, object] = {}
    if sidecar_path is not None:
        sidecar = load_sidecar(sidecar_path)
    merged = merge_sidecar(base, sidecar)
    return ingest_image_meta(merged)


def ingest_session(
    stars: Sequence[StarRecord | Mapping[str, object]],
    meta: ImageMeta | Mapping[str, object],
    pupil: PupilSpec | Mapping[str, object],
) -> tuple[list[StarRecord], ImageMeta, PupilSpec]:
    image_meta = meta if isinstance(meta, ImageMeta) else ingest_image_meta(meta)
    pupil_spec = pupil if isinstance(pupil, PupilSpec) else ingest_pupil_spec(pupil)
    if not stars:
        raise InputError("session has no stars", module="boundary")

    ingested: list[StarRecord] = []
    seen_ids: set[str] = set()
    stamp_size: int | None = None
    p_mm = image_meta.mm_per_pixel()

    for raw in stars:
        star = raw if isinstance(raw, StarRecord) else ingest_star_record(raw)
        if star.exposure_id != image_meta.exposure_id:
            raise InputError(
                "star exposure_id does not match ImageMeta",
                module="boundary",
                star_id=star.star_id,
            )
        if star.session_id != image_meta.session_id:
            raise InputError(
                "star session_id does not match ImageMeta",
                module="boundary",
                star_id=star.star_id,
            )
        if star.star_id in seen_ids:
            raise InputError("duplicate star_id", module="boundary", star_id=star.star_id)
        seen_ids.add(star.star_id)
        if stamp_size is None:
            stamp_size = star.stamp_size
        elif stamp_size != star.stamp_size:
            raise InputError("mixed stamp sizes in one exposure_id", module="boundary")
        if image_meta.n_row < star.stamp_size or image_meta.n_col < star.stamp_size:
            raise InputError("n_row and n_col must be >= S", module="boundary")
        recomputed_x = (star.source_xy_px[0] - image_meta.optical_axis_pixel[0]) * p_mm
        recomputed_y = (star.source_xy_px[1] - image_meta.optical_axis_pixel[1]) * p_mm
        dx = abs(star.field_xy_mm[0] - recomputed_x)
        dy = abs(star.field_xy_mm[1] - recomputed_y)
        if max(dx, dy) > FIELD_TOL_MM:
            raise InputError(
                f"field_xy_mm does not match NOR.8 recomputation (Δ∞={max(dx, dy)} mm)",
                module="boundary",
                star_id=star.star_id,
            )
        ingested.append(star)

    return ingested, image_meta, pupil_spec
