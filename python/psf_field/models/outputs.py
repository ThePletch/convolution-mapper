"""Output and report models 1:1 with the remaining JSON schemas."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field, field_serializer, model_validator

from psf_field.models.common import (
    MODEL_CONFIG,
    Array1F64,
    Array2F64,
    FiniteFloat,
    JsonInt,
    RecordId,
    SchemaVersion,
    Scope,
    Vec2,
    f64_matrix_json,
    f64_vector_json,
)


class ParamMeta(BaseModel):
    """Annotation sidecar for one element of the flat θ vector (term, role, scope, freeze, unit)."""

    model_config = MODEL_CONFIG

    term_id: str
    role: str
    scope: Scope
    frozen: bool
    unit: str


class DegeneratePair(BaseModel):
    """Two free parameters whose stage-1 correlation exceeds the degeneracy threshold."""

    model_config = MODEL_CONFIG

    term_a: str
    term_b: str
    rho: FiniteFloat

    @model_validator(mode="before")
    @classmethod
    def _from_seq(cls, value: object) -> object:
        if isinstance(value, list | tuple) and len(list(value)) == 3:
            items = list(value)
            return {"term_a": items[0], "term_b": items[1], "rho": items[2]}
        return value


class Stage1ErrorPayload(BaseModel):
    """Error code and message when a per-star fit does not succeed."""

    model_config = MODEL_CONFIG

    code: str
    message: str


Termination = Literal[
    "converged_ftol",
    "converged_xtol",
    "converged_gtol",
    "converged_zero_residual",
    "max_eval",
    "numerical",
    "lost_patience",
    "user",
    "unknown",
]


class Stage1Result(BaseModel):
    """Per-star LM solution: coefficients, covariance, residuals, and termination."""

    model_config = MODEL_CONFIG

    schema_version: SchemaVersion
    star_id: RecordId
    theta: list[FiniteFloat]
    theta_init: list[FiniteFloat]
    param_meta: list[ParamMeta]
    free_index: list[JsonInt]
    covariance: Array2F64
    covariance_chi2_scaled: Array2F64
    covariance_ok: bool
    correlation: Array2F64
    degenerate_pairs: list[DegeneratePair]
    residual_image: Array2F64
    weighted_residual: Array1F64
    chi2: FiniteFloat
    chi2_reduced: FiniteFloat
    n_iter: JsonInt
    n_fev: JsonInt
    termination: Termination
    converged: bool
    defocus_sign_ambiguous: bool
    flag_at_bound: list[str]
    error: Stage1ErrorPayload | None = None

    @field_serializer("covariance", "covariance_chi2_scaled", "correlation", "residual_image")
    def _ser_mat(self, value: Array2F64) -> list[list[float]]:
        return f64_matrix_json(value)

    @field_serializer("weighted_residual")
    def _ser_vec(self, value: Array1F64) -> list[float]:
        return f64_vector_json(value)

    @field_serializer("degenerate_pairs")
    def _ser_pairs(self, value: list[DegeneratePair]) -> list[list[str | float]]:
        return [[p.term_a, p.term_b, p.rho] for p in value]


class FieldMap(BaseModel):
    """Stage-2 polynomial (or uniform) fit of one catalog term across the field."""

    model_config = MODEL_CONFIG

    coefficients: list[FiniteFloat]
    covariance: Array2F64
    cond2: FiniteFloat
    ill_conditioned: bool
    independence_assumption: Literal[True]
    residuals: list[FiniteFloat]

    @field_serializer("covariance")
    def _ser_cov(self, value: Array2F64) -> list[list[float]]:
        return f64_matrix_json(value)


class Stage2Result(BaseModel):
    """All field maps and kernel globals from the second-stage regression; no pixels are re-read."""

    model_config = MODEL_CONFIG

    schema_version: SchemaVersion
    maps: dict[str, FieldMap]
    kernel_globals: dict[str, object]
    dropped_star_ids: list[str]
    n_stars_used: JsonInt
    param_meta: list[object]
    star_ids_used: list[str] | None = None


class PsfEval(BaseModel):
    """Model PSF at an arbitrary field position, plus the local coefficients that produced it."""

    model_config = MODEL_CONFIG

    schema_version: SchemaVersion
    psf: Array2F64
    zernike_vector: dict[str, FiniteFloat]
    kernel_vector: dict[str, FiniteFloat]
    field_xy_mm: Vec2
    u_v: Vec2
    image_meta_digest: str = Field(min_length=64, max_length=64)
    catalog_id: str
    stage2_schema_version: SchemaVersion
    extrapolated: bool
    outside_unit_square: bool
    outside_hull: bool

    @field_serializer("psf")
    def _ser_psf(self, value: Array2F64) -> list[list[float]]:
        return f64_matrix_json(value)


class FdReport(BaseModel):
    """Finite-difference check of analytic Jacobian columns at a stated θ."""

    model_config = MODEL_CONFIG

    schema_version: SchemaVersion
    column_errors: list[FiniteFloat]
    column_errors_unconstrained: list[FiniteFloat] | None = None
    passed: bool
    passed_unconstrained: bool


class ScoreEntry(BaseModel):
    """Score-test result for one candidate unmodeled term on one star."""

    model_config = MODEL_CONFIG

    term_id: str
    score: FiniteFloat
    suggest_add: bool
    undefined: bool | None = None


class StarScore(BaseModel):
    """Per-star residual diagnostics and candidate-term scores."""

    model_config = MODEL_CONFIG

    star_id: str
    chi2_reduced: FiniteFloat
    structured_residual: bool
    centroid_leak_suspect: bool
    scores: list[ScoreEntry]


class ScoreReport(BaseModel):
    """Session-level missing-term ranking and documented blur-degeneracy summary."""

    model_config = MODEL_CONFIG

    schema_version: SchemaVersion
    per_star: list[StarScore]
    session_degeneracies: list[object]
    stage2_maps: dict[str, object]
    weak_phase_all_fraction: FiniteFloat


class CoverageReport(BaseModel):
    """How selected stars fill the detector, for judging stage-2 field-map conditioning."""

    model_config = MODEL_CONFIG

    schema_version: SchemaVersion
    n_detected: JsonInt
    n_candidate: JsonInt
    n_selected: JsonInt
    frac_selected_of_detected: FiniteFloat
    grid_3x3: tuple[JsonInt, JsonInt, JsonInt, JsonInt, JsonInt, JsonInt, JsonInt, JsonInt, JsonInt]
    empty_cells: JsonInt
    convex_hull_area_mm2: FiniteFloat
    detector_area_mm2: FiniteFloat
    hull_fill: FiniteFloat
    design_cond_plane: FiniteFloat
    design_cond_quad: FiniteFloat

    @model_validator(mode="before")
    @classmethod
    def _grid(cls, value: object) -> object:
        if not isinstance(value, dict):
            return value
        grid = value.get("grid_3x3")
        if isinstance(grid, list) and len(grid) == 9:
            data: dict[str, object] = {str(key): item for key, item in value.items()}
            data["grid_3x3"] = tuple(grid)
            return data
        return value
