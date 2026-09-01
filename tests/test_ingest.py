"""Ingest rejection tests: unknown MAJOR, even S, incomplete ImageMeta sidecar."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from psf_field.errors import InputError
from psf_field.ingest import ingest_image_meta_from_fits, ingest_pupil_spec, ingest_star_record
from psf_field.io.fits_header import write_primary_header


def _valid_star(size: int = 15) -> dict[str, object]:
    center = (size - 1) / 2.0
    ones = np.ones((size, size), dtype=np.float64)
    return {
        "schema_version": "1.0.0",
        "star_id": "s1",
        "exposure_id": "exp1",
        "session_id": "sess1",
        "field_xy_mm": [0.0, 0.0],
        "source_xy_px": [511.5, 511.5],
        "stamp": ones,
        "variance": ones,
        "centroid_xy_px": [center, center],
        "pixel_mask": np.zeros((size, size), dtype=np.uint8),
        "flags": ["SELECTED"],
        "flux_sum_adu": float(size * size),
    }


def test_rejects_unknown_schema_major() -> None:
    payload = _valid_star()
    payload["schema_version"] = "2.0.0"
    with pytest.raises(InputError, match="unknown schema MAJOR 2") as exc_info:
        ingest_star_record(payload)
    assert exc_info.value.code == "input"
    assert exc_info.value.module == "boundary"


def test_rejects_even_stamp_size() -> None:
    payload = _valid_star(16)
    with pytest.raises(InputError, match="even stamp size"):
        ingest_star_record(payload)


def test_sidecar_incomplete_image_meta(tmp_path: Path) -> None:
    fits_path = tmp_path / "incomplete.fits"
    sidecar_path = tmp_path / "complete.json"
    write_primary_header(
        fits_path,
        {
            "SCHEMAV": "1.0.0",
            "EXPID": "exp1",
            "SESSID": "sess1",
            "NAXIS1": 1024,
            "NAXIS2": 1024,
            "EXPTIME": 30.0,
            "GAIN": 1.5,
            "RDNOISE": 5.0,
            "SATURATE": 60000.0,
        },
    )
    with pytest.raises(InputError):
        ingest_image_meta_from_fits(fits_path)

    sidecar_path.write_text(
        json.dumps(
            {
                "wavelength_m": 5.5e-7,
                "pupil_diameter_m": 0.2,
                "focal_length_m": 1.6,
                "pixel_scale_arcsec": 0.4,
                "optical_axis_pixel": [511.5, 511.5],
            }
        ),
        encoding="utf-8",
    )
    meta = ingest_image_meta_from_fits(fits_path, sidecar_path)
    assert meta.wavelength_m == pytest.approx(5.5e-7)
    assert meta.optical_axis_pixel == (511.5, 511.5)
    np.testing.assert_allclose(meta.optical_axis_pixel, (511.5, 511.5))


def test_pupil_ingest_zeros_amplitude_outside_the_mask() -> None:
    n_pupil = 128
    yy, xx = np.ogrid[:n_pupil, :n_pupil]
    center = (n_pupil - 1) / 2.0
    rho = np.hypot(xx - center, yy - center) / (n_pupil / 2.0)
    mask = (rho <= 1.0).astype(np.float64)
    pupil = ingest_pupil_spec(
        {
            "schema_version": "1.0.0",
            "mask": mask,
            "n_pupil": n_pupil,
            "n_fft": n_pupil * 4,
            "amplitude": np.ones((n_pupil, n_pupil), dtype=np.float64),
        }
    )
    assert pupil.amplitude is not None
    np.testing.assert_array_equal(pupil.amplitude[mask == 0.0], 0.0)
    np.testing.assert_array_equal(pupil.amplitude[mask == 1.0], 1.0)
