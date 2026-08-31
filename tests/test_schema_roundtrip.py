"""JSON Schema round-trip tests for pydantic models (C1.9.2)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
from jsonschema.validators import Draft202012Validator

from psf_field.ingest import ingest_catalog, ingest_image_meta, ingest_pupil_spec, ingest_star_record
from psf_field.models.common import dump_json
from psf_field.models.inputs import ExtractionConfig
from psf_field.models.outputs import (
    CoverageReport,
    FdReport,
    PsfEval,
    ScoreReport,
    Stage1Result,
    Stage2Result,
)

REPO = Path(__file__).resolve().parents[1]
SCHEMAS = REPO / "docs" / "contracts" / "schemas"


def _schema(name: str) -> dict[str, Any]:
    with (SCHEMAS / name).open(encoding="utf-8") as handle:
        loaded: object = json.load(handle)
    assert isinstance(loaded, dict)
    return loaded


def _assert_schema(name: str, payload: dict[str, object]) -> None:
    Draft202012Validator(_schema(name)).validate(payload)


def valid_meta() -> dict[str, object]:
    return {
        "schema_version": "1.0.0",
        "exposure_id": "exp1",
        "session_id": "sess1",
        "n_row": 1024,
        "n_col": 1024,
        "wavelength_m": 5.5e-7,
        "pupil_diameter_m": 0.2,
        "focal_length_m": 1.6,
        "pixel_scale_arcsec": 0.4,
        "optical_axis_pixel": [511.5, 511.5],
        "gain_e_per_adu": 1.5,
        "read_noise_e": 5.0,
        "saturation_adu": 60000.0,
        "exptime_s": 30.0,
    }


def valid_star(size: int = 15) -> dict[str, object]:
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


def valid_pupil(n_pupil: int = 128) -> dict[str, object]:
    yy, xx = np.ogrid[:n_pupil, :n_pupil]
    center = (n_pupil - 1) / 2.0
    rho = np.hypot(xx - center, yy - center) / (n_pupil / 2.0)
    mask = (rho <= 1.0).astype(np.float64)
    return {
        "schema_version": "1.0.0",
        "mask": mask,
        "n_pupil": n_pupil,
        "n_fft": n_pupil * 4,
    }


def test_star_record_round_trip() -> None:
    rec = ingest_star_record(valid_star())
    payload = dump_json(rec)
    _assert_schema("star_record.schema.json", payload)
    again = ingest_star_record(payload)
    assert again.star_id == rec.star_id
    assert again.stamp.shape == rec.stamp.shape


def test_image_meta_round_trip() -> None:
    meta = ingest_image_meta(valid_meta())
    payload = dump_json(meta)
    _assert_schema("image_meta.schema.json", payload)
    again = ingest_image_meta(payload)
    assert again.exposure_id == meta.exposure_id


def test_pupil_spec_round_trip() -> None:
    pupil = ingest_pupil_spec(valid_pupil())
    payload = dump_json(pupil)
    _assert_schema("pupil_spec.schema.json", payload)
    again = ingest_pupil_spec(payload)
    assert again.n_pupil == pupil.n_pupil


def test_default_catalog_round_trip() -> None:
    raw = json.loads((SCHEMAS / "psf_field_v1_default.catalog.json").read_text())
    catalog = ingest_catalog(raw)
    payload = dump_json(catalog)
    _assert_schema("catalog.schema.json", payload)
    again = ingest_catalog(payload)
    assert again.catalog_id == "psf_field_v1_default"
    assert len(again.terms) == len(catalog.terms)


def test_extraction_config_round_trip() -> None:
    cfg = ExtractionConfig.model_validate({"schema_version": "1.0.0", "fwhm": 3.5})
    payload = dump_json(cfg)
    _assert_schema("extraction_config.schema.json", payload)


def test_stage1_result_round_trip() -> None:
    result = Stage1Result.model_validate(
        {
            "schema_version": "1.0.0",
            "star_id": "s1",
            "theta": [0.1],
            "theta_init": [0.0],
            "param_meta": [
                {
                    "term_id": "flux",
                    "role": "flux",
                    "scope": "per_star",
                    "frozen": False,
                    "unit": "ADU",
                }
            ],
            "free_index": [0],
            "covariance": [[1.0]],
            "covariance_chi2_scaled": [[1.0]],
            "covariance_ok": True,
            "correlation": [[1.0]],
            "degenerate_pairs": [["flux", "sky", 0.95]],
            "residual_image": [[0.0]],
            "weighted_residual": [0.0],
            "chi2": 1.0,
            "chi2_reduced": 1.0,
            "n_iter": 1,
            "n_fev": 2,
            "termination": "converged_ftol",
            "converged": True,
            "defocus_sign_ambiguous": True,
            "flag_at_bound": [],
        }
    )
    payload = dump_json(result)
    _assert_schema("stage1_result.schema.json", payload)
    again = Stage1Result.model_validate(payload)
    assert again.star_id == "s1"
    assert again.degenerate_pairs[0].term_a == "flux"


def test_stage2_result_round_trip() -> None:
    result = Stage2Result.model_validate(
        {
            "schema_version": "1.0.0",
            "maps": {
                "zernike_2_0": {
                    "coefficients": [0.0],
                    "covariance": [[1.0]],
                    "cond2": 1.0,
                    "ill_conditioned": False,
                    "independence_assumption": True,
                    "residuals": [0.0],
                }
            },
            "kernel_globals": {},
            "dropped_star_ids": [],
            "n_stars_used": 1,
            "param_meta": [],
        }
    )
    payload = dump_json(result)
    _assert_schema("stage2_result.schema.json", payload)


def test_psf_eval_round_trip() -> None:
    digest = "a" * 64
    result = PsfEval.model_validate(
        {
            "schema_version": "1.0.0",
            "psf": [[1.0]],
            "zernike_vector": {"zernike_2_0": 0.1},
            "kernel_vector": {"moffat_seeing": 1.2},
            "field_xy_mm": [0.0, 0.0],
            "u_v": [0.0, 0.0],
            "image_meta_digest": digest,
            "catalog_id": "psf_field_v1_default",
            "stage2_schema_version": "1.0.0",
            "extrapolated": False,
            "outside_unit_square": False,
            "outside_hull": False,
        }
    )
    payload = dump_json(result)
    _assert_schema("psf_eval.schema.json", payload)


def test_fd_score_coverage_round_trip() -> None:
    fd = FdReport.model_validate(
        {
            "schema_version": "1.0.0",
            "column_errors": [1e-6],
            "passed": True,
            "passed_unconstrained": True,
        }
    )
    _assert_schema("fd_report.schema.json", dump_json(fd))

    score = ScoreReport.model_validate(
        {
            "schema_version": "1.0.0",
            "per_star": [
                {
                    "star_id": "s1",
                    "chi2_reduced": 1.0,
                    "structured_residual": False,
                    "centroid_leak_suspect": False,
                    "scores": [{"term_id": "zernike_5_1", "score": 0.2, "suggest_add": True}],
                }
            ],
            "session_degeneracies": [],
            "stage2_maps": {},
            "weak_phase_all_fraction": 0.0,
        }
    )
    _assert_schema("score_report.schema.json", dump_json(score))

    coverage = CoverageReport.model_validate(
        {
            "schema_version": "1.0.0",
            "n_detected": 10,
            "n_candidate": 8,
            "n_selected": 5,
            "frac_selected_of_detected": 0.5,
            "grid_3x3": [1, 1, 1, 0, 1, 0, 0, 1, 0],
            "empty_cells": 4,
            "convex_hull_area_mm2": 1.0,
            "detector_area_mm2": 4.0,
            "hull_fill": 0.25,
            "design_cond_plane": 2.0,
            "design_cond_quad": 5.0,
        }
    )
    _assert_schema("coverage_report.schema.json", dump_json(coverage))
