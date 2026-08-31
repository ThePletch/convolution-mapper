"""Catalog models 1:1 with `catalog.schema.json`."""

from __future__ import annotations

from typing import Annotated, Literal, Self

from collections.abc import Callable

from pydantic import BaseModel, Field, model_serializer, model_validator

from psf_field.models.common import (
    MODEL_CONFIG,
    FiniteFloat,
    JsonInt,
    PriorMean,
    SchemaVersion,
    Scope,
    TermId,
)


class FieldBasis(BaseModel):
    model_config = MODEL_CONFIG

    family: Literal["monomial"]
    degree: JsonInt
    terms: list[list[JsonInt]]

    @model_validator(mode="after")
    def _term_pairs(self) -> Self:
        for term in self.terms:
            if len(term) != 2:
                raise ValueError("field basis terms must be [i, j] pairs")
        return self


class Stage2Prior(BaseModel):
    model_config = MODEL_CONFIG

    mean: list[FiniteFloat]
    sigma: list[FiniteFloat]


class PriorSpec(BaseModel):
    model_config = MODEL_CONFIG

    family: Literal["none", "gaussian"]
    mean: PriorMean | None = None
    sigma: FiniteFloat | None = None
    sigma_rel: FiniteFloat | None = None
    stage2: Stage2Prior | None = None

    @model_validator(mode="after")
    def _gaussian_fields(self) -> Self:
        if self.family == "gaussian":
            if self.mean == "init":
                if self.sigma_rel is None:
                    raise ValueError("gaussian prior with mean 'init' requires sigma_rel")
                if self.sigma is not None:
                    raise ValueError("gaussian prior with mean 'init' must omit sigma")
            elif isinstance(self.mean, float):
                if self.sigma is None:
                    raise ValueError("gaussian prior with numeric mean requires sigma")
                if self.sigma_rel is not None:
                    raise ValueError("gaussian prior with numeric mean must omit sigma_rel")
            else:
                raise ValueError("gaussian prior requires mean")
        return self


class InitSpec(BaseModel):
    model_config = MODEL_CONFIG

    method: Literal["zero", "flux_sum", "defocus_moment", "moffat_fwhm"]


class KernelSpec(BaseModel):
    model_config = MODEL_CONFIG

    id: Literal[
        "gaussian_iso",
        "moffat_iso",
        "linear_drift",
        "field_rotation",
        "periodic_error",
    ]


class _TermFields(BaseModel):
    model_config = MODEL_CONFIG

    term_id: TermId
    name: str
    scope: Scope
    frozen: bool
    enabled: bool
    bounds: tuple[FiniteFloat, FiniteFloat] | None = None
    init: InitSpec
    prior: PriorSpec
    units: str
    report: dict[str, object] | None = None

    @model_validator(mode="before")
    @classmethod
    def _bounds_list(cls, value: object) -> object:
        if not isinstance(value, dict):
            return value
        bounds = value.get("bounds")
        if isinstance(bounds, list) and len(bounds) == 2:
            data: dict[str, object] = {str(key): item for key, item in value.items()}
            data["bounds"] = (bounds[0], bounds[1])
            return data
        return value


class PhaseTerm(_TermFields):
    kind: Literal["phase"]
    n: JsonInt
    m: JsonInt
    field_basis: FieldBasis


class KernelTerm(_TermFields):
    kind: Literal["kernel"]
    kernel: KernelSpec
    field_basis: FieldBasis | None = None


class PhotometricTerm(_TermFields):
    kind: Literal["photometric"]

    @model_validator(mode="after")
    def _photometric_id(self) -> Self:
        if self.term_id not in {"flux", "sky"}:
            raise ValueError("photometric term_id must be flux or sky")
        return self


ErrorTerm = Annotated[
    PhaseTerm | KernelTerm | PhotometricTerm,
    Field(discriminator="kind"),
]


class Bundle(BaseModel):
    model_config = MODEL_CONFIG

    bundle_id: TermId
    name: str
    term_ids: list[str]
    matrix: None = None


class FitScheduleStep(BaseModel):
    model_config = MODEL_CONFIG

    name: str
    unfrozen_term_ids: list[str]


class Catalog(BaseModel):
    model_config = MODEL_CONFIG

    schema_version: SchemaVersion
    catalog_id: str = Field(min_length=1)
    terms: list[ErrorTerm] = Field(min_length=1)
    bundles: list[Bundle]
    fit_schedule: list[FitScheduleStep]

    @model_serializer(mode="wrap")
    def _keep_null_matrix(self, nxt: Callable[[Catalog], object]) -> dict[str, object]:
        dumped: object = nxt(self)
        if not isinstance(dumped, dict):
            msg = "catalog serializer must return an object"
            raise TypeError(msg)
        out: dict[str, object] = dict(dumped)
        bundles = out.get("bundles")
        if isinstance(bundles, list):
            patched: list[object] = []
            for bundle in bundles:
                if isinstance(bundle, dict) and "matrix" not in bundle:
                    item: dict[str, object] = dict(bundle)
                    item["matrix"] = None
                    patched.append(item)
                else:
                    patched.append(bundle)
            out["bundles"] = patched
        return out
