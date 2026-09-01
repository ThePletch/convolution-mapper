"""Pydantic boundary models 1:1 with `docs/contracts/schemas/` (C1.9.2)."""

from psf_field.models.catalog import (
    Bundle,
    Catalog,
    ErrorTerm,
    FieldBasis,
    FitScheduleStep,
    KernelSpec,
    KernelTerm,
    PhaseTerm,
    PhotometricTerm,
    PriorSpec,
)
from psf_field.models.inputs import (
    ExtractionConfig,
    ImageMeta,
    PupilSpec,
    Stage1InputConfig,
    StarRecord,
)
from psf_field.models.outputs import (
    CoverageReport,
    FdReport,
    FieldMap,
    PsfEval,
    ScoreReport,
    Stage1Result,
    Stage2Result,
)

__all__ = [
    "Bundle",
    "Catalog",
    "CoverageReport",
    "ErrorTerm",
    "ExtractionConfig",
    "FdReport",
    "FieldBasis",
    "FieldMap",
    "FitScheduleStep",
    "ImageMeta",
    "KernelSpec",
    "KernelTerm",
    "PhaseTerm",
    "PhotometricTerm",
    "PriorSpec",
    "PsfEval",
    "PupilSpec",
    "ScoreReport",
    "Stage1InputConfig",
    "Stage1Result",
    "Stage2Result",
    "StarRecord",
]
