"""Catalog models 1:1 with `catalog.schema.json`."""

from __future__ import annotations

from typing import Annotated, Any, Literal, Self

from collections.abc import Callable

from pydantic import BaseModel, Discriminator, Field, Tag, model_serializer, model_validator

from psf_field.models.common import (
    MODEL_CONFIG,
    FiniteFloat,
    JsonInt,
    SchemaVersion,
    Scope,
    TermId,
)


class FieldBasis(BaseModel):
    """Which monomials in normalized field (u, v) a term's stage-2 map may use."""

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
    """Optional Gaussian prior on stage-2 field-map coefficients, aligned with `FieldBasis.terms`."""

    model_config = MODEL_CONFIG

    mean: list[FiniteFloat]
    sigma: list[FiniteFloat]


class PriorNone(BaseModel):
    """No stage-1 Gaussian prior on the local coefficient."""

    model_config = MODEL_CONFIG

    family: Literal["none"]
    stage2: Stage2Prior | None = None


class PriorGaussian(BaseModel):
    """Gaussian prior whose mean is a catalog number and whose σ is absolute."""

    model_config = MODEL_CONFIG

    family: Literal["gaussian"]
    mean: FiniteFloat
    sigma: FiniteFloat
    stage2: Stage2Prior | None = None


class PriorGaussianFromInit(BaseModel):
    """Gaussian prior whose mean is the initialization value a₀ and whose σ scales with |a₀|."""

    model_config = MODEL_CONFIG

    family: Literal["gaussian"]
    mean: Literal["init"]
    sigma_rel: FiniteFloat
    stage2: Stage2Prior | None = None


def _prior_tag(value: Any) -> str:
    if isinstance(value, dict):
        family = value.get("family")
        if family == "gaussian":
            if value.get("mean") == "init":
                return "gaussian_from_init"
            return "gaussian"
        return str(family)
    if isinstance(value, PriorGaussianFromInit):
        return "gaussian_from_init"
    if isinstance(value, PriorGaussian):
        return "gaussian"
    if isinstance(value, PriorNone):
        return "none"
    return type(value).__name__


PriorSpec = Annotated[
    Annotated[PriorNone, Tag("none")]
    | Annotated[PriorGaussianFromInit, Tag("gaussian_from_init")]
    | Annotated[PriorGaussian, Tag("gaussian")],
    Discriminator(_prior_tag),
]


class ZeroInit(BaseModel):
    """Initialize the local coefficient to 0, or to a numeric Gaussian prior mean if present."""

    model_config = MODEL_CONFIG

    method: Literal["zero"]


class PhaseInit(BaseModel):
    """Initialize a Zernike coefficient: 0, or the defocus second-moment formula when (n, m) = (2, 0)."""

    model_config = MODEL_CONFIG

    method: Literal["zero", "defocus_moment"]


class FluxInit(BaseModel):
    """Initialize per-star flux to 0 or to the stamp sum."""

    model_config = MODEL_CONFIG

    method: Literal["zero", "flux_sum"]


class MoffatInit(BaseModel):
    """Initialize Moffat α to 0 or from the extraction-time expected FWHM."""

    model_config = MODEL_CONFIG

    method: Literal["zero", "moffat_fwhm"]


class GaussianIsoSpec(BaseModel):
    """Isotropic Gaussian kernel on the detector."""

    model_config = MODEL_CONFIG

    id: Literal["gaussian_iso"]


class MoffatIsoSpec(BaseModel):
    """Isotropic Moffat kernel on the detector."""

    model_config = MODEL_CONFIG

    id: Literal["moffat_iso"]


class LinearDriftSpec(BaseModel):
    """Linear trail kernel (length and angle on the detector)."""

    model_config = MODEL_CONFIG

    id: Literal["linear_drift"]


class FieldRotationSpec(BaseModel):
    """Field-rotation kernel; stage-1 fits a local trail, not the three globals."""

    model_config = MODEL_CONFIG

    id: Literal["field_rotation"]


KernelSpec = GaussianIsoSpec | MoffatIsoSpec | LinearDriftSpec | FieldRotationSpec


def _prior_is_from_init(prior: PriorSpec) -> bool:
    return isinstance(prior, PriorGaussianFromInit)


class _TermFields(BaseModel):
    """Shared freeze, bounds, and units for every catalog term kind."""

    model_config = MODEL_CONFIG

    name: str
    scope: Scope
    frozen: bool
    enabled: bool
    bounds: tuple[FiniteFloat, FiniteFloat] | None = None
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


class _CatalogTerm(_TermFields):
    """Catalog term with a free-form `term_id` (phase and kernel terms)."""

    term_id: TermId


class PhaseTerm(_CatalogTerm):
    """A Zernike (n, m) pupil-phase term with a field basis."""

    kind: Literal["phase"]
    n: JsonInt
    m: JsonInt
    field_basis: FieldBasis
    init: PhaseInit
    prior: PriorSpec

    @model_validator(mode="after")
    def _init_pairing(self) -> Self:
        if self.init.method == "defocus_moment" and (self.n, self.m) != (2, 0):
            raise ValueError(f"init method defocus_moment is not valid for term {self.term_id}")
        if _prior_is_from_init(self.prior) and self.init.method == "zero":
            raise ValueError(f'term {self.term_id} pairs mean "init" with init method zero')
        return self


class GaussianIsoTerm(_CatalogTerm):
    """Isotropic Gaussian kernel (jitter or charge diffusion)."""

    kind: Literal["kernel"]
    kernel: GaussianIsoSpec
    field_basis: FieldBasis | None = None
    init: ZeroInit
    prior: PriorNone | PriorGaussian


class MoffatIsoTerm(_CatalogTerm):
    """Isotropic Moffat seeing kernel."""

    kind: Literal["kernel"]
    kernel: MoffatIsoSpec
    field_basis: FieldBasis | None = None
    init: MoffatInit
    prior: PriorSpec

    @model_validator(mode="after")
    def _init_pairing(self) -> Self:
        if _prior_is_from_init(self.prior) and self.init.method == "zero":
            raise ValueError(f'term {self.term_id} pairs mean "init" with init method zero')
        return self


class LinearDriftTerm(_CatalogTerm):
    """Linear trail kernel."""

    kind: Literal["kernel"]
    kernel: LinearDriftSpec
    field_basis: FieldBasis | None = None
    init: ZeroInit
    prior: PriorNone | PriorGaussian


class FieldRotationTerm(_CatalogTerm):
    """Field-rotation kernel as a local trail in stage 1."""

    kind: Literal["kernel"]
    kernel: FieldRotationSpec
    field_basis: FieldBasis | None = None
    init: ZeroInit
    prior: PriorNone | PriorGaussian


class FluxTerm(_TermFields):
    """Per-star flux scalar."""

    kind: Literal["photometric"]
    term_id: Literal["flux"]
    init: FluxInit
    prior: PriorSpec

    @model_validator(mode="after")
    def _init_pairing(self) -> Self:
        if _prior_is_from_init(self.prior) and self.init.method == "zero":
            raise ValueError(f'term {self.term_id} pairs mean "init" with init method zero')
        return self


class SkyTerm(_TermFields):
    """Per-star residual-sky scalar."""

    kind: Literal["photometric"]
    term_id: Literal["sky"]
    init: ZeroInit
    prior: PriorNone | PriorGaussian


def _term_tag(value: Any) -> str:
    if isinstance(value, dict):
        kind = value.get("kind")
        if kind == "kernel":
            kernel = value.get("kernel")
            kernel_id = kernel.get("id") if isinstance(kernel, dict) else None
            return f"kernel:{kernel_id}"
        if kind == "photometric":
            return f"photometric:{value.get('term_id')}"
        return str(kind)
    kind = getattr(value, "kind", None)
    if kind == "kernel":
        kernel = getattr(value, "kernel", None)
        return f"kernel:{getattr(kernel, 'id', None)}"
    if kind == "photometric":
        return f"photometric:{getattr(value, 'term_id', None)}"
    return str(kind)


KernelTerm = GaussianIsoTerm | MoffatIsoTerm | LinearDriftTerm | FieldRotationTerm
PhotometricTerm = FluxTerm | SkyTerm

ErrorTerm = Annotated[
    Annotated[PhaseTerm, Tag("phase")]
    | Annotated[GaussianIsoTerm, Tag("kernel:gaussian_iso")]
    | Annotated[MoffatIsoTerm, Tag("kernel:moffat_iso")]
    | Annotated[LinearDriftTerm, Tag("kernel:linear_drift")]
    | Annotated[FieldRotationTerm, Tag("kernel:field_rotation")]
    | Annotated[FluxTerm, Tag("photometric:flux")]
    | Annotated[SkyTerm, Tag("photometric:sky")],
    Discriminator(_term_tag),
]


class Bundle(BaseModel):
    """Named linear grouping of catalog terms; v1 ships with a null sensitivity matrix."""

    model_config = MODEL_CONFIG

    bundle_id: TermId
    name: str
    term_ids: list[str]
    matrix: None


class FitScheduleStep(BaseModel):
    """One freeze/unfreeze step in the staged stage-1 fit."""

    model_config = MODEL_CONFIG

    name: str
    unfrozen_term_ids: list[str]


class Catalog(BaseModel):
    """Named aberration and kernel terms, optional bundles, and the default fit schedule."""

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
