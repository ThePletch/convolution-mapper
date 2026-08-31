"""Contract-1 input models and C1A extraction config."""

from __future__ import annotations

from typing import Literal, Self

import numpy as np
from pydantic import BaseModel, Field, field_serializer, model_validator

from psf_field.models.common import (
    MODEL_CONFIG,
    Array2F64,
    Array2U8,
    FiniteFloat,
    FlagName,
    JsonInt,
    RecordId,
    SchemaVersion,
    Vec2,
    f64_matrix_json,
    u8_matrix_json,
)

STAMP_SIZE_MIN = 15
STAMP_SIZE_MAX = 63
CENTROID_MAX_OFFSET_PX = 0.6


class StarRecord(BaseModel):
    """One star's postage stamp, noise map, centroid, and flags as presented to stage-1 fitting."""

    model_config = MODEL_CONFIG

    schema_version: SchemaVersion
    star_id: RecordId
    exposure_id: RecordId
    session_id: RecordId
    field_xy_mm: Vec2
    source_xy_px: Vec2
    stamp: Array2F64
    variance: Array2F64
    centroid_xy_px: Vec2
    pixel_mask: Array2U8
    flags: list[FlagName]
    flux_sum_adu: FiniteFloat

    @field_serializer("stamp", "variance")
    def _ser_f64(self, value: Array2F64) -> list[list[float]]:
        return f64_matrix_json(value)

    @field_serializer("pixel_mask")
    def _ser_u8(self, value: Array2U8) -> list[list[int]]:
        return u8_matrix_json(value)

    @model_validator(mode="after")
    def _ingest_shape(self) -> Self:
        s_row, s_col = self.stamp.shape
        if s_row != s_col:
            raise ValueError("stamp must be square")
        if s_row % 2 == 0:
            raise ValueError(f"even stamp size {s_row} is rejected")
        if s_row < STAMP_SIZE_MIN or s_row > STAMP_SIZE_MAX:
            raise ValueError(
                f"stamp size {s_row} is outside the allowed odd set [{STAMP_SIZE_MIN}, {STAMP_SIZE_MAX}]"
            )
        if self.variance.shape != (s_row, s_col):
            raise ValueError("variance shape must match stamp")
        if self.pixel_mask.shape != (s_row, s_col):
            raise ValueError("pixel_mask shape must match stamp")
        c_star = (s_row - 1) / 2.0
        dx = abs(self.centroid_xy_px[0] - c_star)
        dy = abs(self.centroid_xy_px[1] - c_star)
        if max(dx, dy) > CENTROID_MAX_OFFSET_PX:
            raise ValueError("centroid_xy_px is more than 0.6 px from the stamp center")
        if len(self.flags) != len(set(self.flags)):
            raise ValueError("flags must be unique")
        valid = self.pixel_mask == 0
        if np.any(~np.isfinite(self.stamp[valid])):
            raise ValueError("stamp has a non-finite value on a valid pixel")
        var = self.variance[valid]
        if np.any(~np.isfinite(var) | (var <= 0.0)):
            raise ValueError("variance must be finite and > 0 on a valid pixel")
        return self

    @property
    def stamp_size(self) -> int:
        return int(self.stamp.shape[0])


class ImageMeta(BaseModel):
    """Per-exposure camera and geometry metadata required to interpret stamps (NOR.13)."""

    model_config = MODEL_CONFIG

    schema_version: SchemaVersion
    exposure_id: RecordId
    session_id: RecordId
    n_row: JsonInt = Field(ge=1)
    n_col: JsonInt = Field(ge=1)
    wavelength_m: FiniteFloat
    pupil_diameter_m: FiniteFloat
    focal_length_m: FiniteFloat
    pixel_scale_arcsec: FiniteFloat
    optical_axis_pixel: Vec2
    gain_e_per_adu: FiniteFloat
    read_noise_e: FiniteFloat
    saturation_adu: FiniteFloat
    exptime_s: FiniteFloat
    known_defocus_waves: FiniteFloat = 0.0
    pixel_size_m: FiniteFloat | None = None
    plate_scale_warning: bool = False

    @model_validator(mode="after")
    def _ingest_physics(self) -> Self:
        for name in (
            "wavelength_m",
            "pupil_diameter_m",
            "focal_length_m",
            "pixel_scale_arcsec",
            "gain_e_per_adu",
            "saturation_adu",
            "exptime_s",
        ):
            if getattr(self, name) <= 0.0:
                raise ValueError(f"{name} must be > 0")
        if self.read_noise_e < 0.0:
            raise ValueError("read_noise_e must be >= 0")
        p_mm = mm_per_pixel(self.focal_length_m, self.pixel_scale_arcsec)
        r_field = 0.5 * float(np.hypot(self.n_col * p_mm, self.n_row * p_mm))
        if r_field == 0.0:
            raise ValueError("R_field is 0 (NOR.9)")
        warning = False
        if self.pixel_size_m is not None:
            if self.pixel_size_m <= 0.0:
                raise ValueError("pixel_size_m must be > 0")
            s_pred = self.pixel_size_m / self.focal_length_m
            s_hdr = pixel_scale_rad(self.pixel_scale_arcsec)
            delta = abs(s_pred - s_hdr) / s_hdr
            if delta > 0.05:
                raise ValueError(f"plate-scale inconsistency δ={delta} exceeds 0.05 (C1.4.1)")
            warning = delta > 0.01
        if warning != self.plate_scale_warning:
            return self.model_copy(update={"plate_scale_warning": warning})
        return self

    def mm_per_pixel(self) -> float:
        return mm_per_pixel(self.focal_length_m, self.pixel_scale_arcsec)

    def r_field_mm(self) -> float:
        p_mm = self.mm_per_pixel()
        return 0.5 * float(np.hypot(self.n_col * p_mm, self.n_row * p_mm))


def pixel_scale_rad(pixel_scale_arcsec: float) -> float:
    return pixel_scale_arcsec * np.pi / (180.0 * 3600.0)


def mm_per_pixel(focal_length_m: float, pixel_scale_arcsec: float) -> float:
    return focal_length_m * pixel_scale_rad(pixel_scale_arcsec) * 1_000.0


_ALLOWED_N_PUPIL = frozenset({128, 256, 512})
_ALLOWED_FFT_RATIO = frozenset({2, 4, 8})


class PupilSpec(BaseModel):
    """Pupil mask and FFT sampling grid used by the forward model."""

    model_config = MODEL_CONFIG

    schema_version: SchemaVersion
    mask: Array2F64
    n_pupil: JsonInt
    n_fft: JsonInt
    amplitude: Array2F64 | None = None

    @field_serializer("mask")
    def _ser_mask(self, value: Array2F64) -> list[list[float]]:
        return f64_matrix_json(value)

    @field_serializer("amplitude")
    def _ser_amp(self, value: Array2F64 | None) -> list[list[float]] | None:
        return None if value is None else f64_matrix_json(value)

    @model_validator(mode="after")
    def _ingest_pupil(self) -> Self:
        if self.n_pupil not in _ALLOWED_N_PUPIL:
            raise ValueError(f"n_pupil {self.n_pupil} is not in {{128, 256, 512}}")
        if self.n_fft % self.n_pupil != 0:
            raise ValueError("n_fft must be an integer multiple of n_pupil")
        ratio = self.n_fft // self.n_pupil
        if ratio not in _ALLOWED_FFT_RATIO:
            raise ValueError(f"n_fft / n_pupil = {ratio} is not in {{2, 4, 8}}")
        n_row, n_col = self.mask.shape
        if n_row != n_col or n_row != self.n_pupil:
            raise ValueError("mask shape must be (n_pupil, n_pupil)")
        unique = set(np.unique(self.mask).tolist())
        if not unique <= {0.0, 1.0}:
            raise ValueError("v1 pupil mask values must be in {0, 1}")
        if self.amplitude is not None and self.amplitude.shape != self.mask.shape:
            raise ValueError("amplitude shape must match mask")
        return self


class Stage1InputConfig(BaseModel):
    """Session-level stamp-size defaults stored beside the star table."""

    model_config = MODEL_CONFIG

    schema_version: SchemaVersion = "1.0.0"
    stamp_size: JsonInt = 31

    @model_validator(mode="after")
    def _check_s(self) -> Self:
        s = self.stamp_size
        if s % 2 == 0:
            raise ValueError(f"even stamp size {s} is rejected")
        if s < STAMP_SIZE_MIN or s > STAMP_SIZE_MAX:
            raise ValueError(
                f"stamp size {s} is outside the allowed odd set [{STAMP_SIZE_MIN}, {STAMP_SIZE_MAX}]"
            )
        return self


class ExtractionConfig(BaseModel):
    """DAOPHOT-style detection and star-selection knobs for the C1A extraction front-end."""

    model_config = MODEL_CONFIG

    schema_version: SchemaVersion
    box_size: JsonInt = 64
    filter_size: Literal[1, 3, 5] = 3
    sigma_clip_sigma: FiniteFloat = 3.0
    sigma_clip_maxiters: JsonInt = 5
    fwhm: FiniteFloat
    n_sigma: FiniteFloat = 5.0
    sharpness_range: Vec2 = (0.2, 1.0)
    roundness_range: Vec2 = (-1.0, 1.0)
    min_separation_fwhm: FiniteFloat = 1.0
    snr_min: FiniteFloat = 20.0
    max_selected: JsonInt = 400
    selection_mode: Literal["highest_snr", "snr_by_cell"] = "highest_snr"
    holdout_fraction: FiniteFloat = 0.0
    holdout_seed: JsonInt | None = None

    @model_validator(mode="after")
    def _bounds(self) -> Self:
        if self.fwhm <= 1.0:
            raise ValueError("fwhm must be finite and > 1.0")
        return self
