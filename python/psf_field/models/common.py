"""Shared pydantic config, schema versioning, and array coercions (C1B.8)."""

from __future__ import annotations

import math
import re
from typing import Annotated, Literal, TypeAlias

import numpy as np
import numpy.typing as npt
from pydantic import AfterValidator, BaseModel, BeforeValidator, ConfigDict, ValidationError

from psf_field.errors import InputError

MODEL_CONFIG = ConfigDict(
    extra="forbid",
    strict=True,
    frozen=True,
    arbitrary_types_allowed=True,
)

SCHEMA_VERSION: str = "1.0.0"
SCHEMA_MAJOR: int = 1

_ID_RE = re.compile(r"^[A-Za-z0-9._:-]{1,128}$")
_TERM_ID_RE = re.compile(r"^[a-z][a-z0-9_]*$")

JsonScalar: TypeAlias = str | int | float | bool | None
JsonValue: TypeAlias = JsonScalar | list["JsonValue"] | dict[str, "JsonValue"]

Scope: TypeAlias = Literal["per_star", "per_exposure", "per_session"]
FlagName: TypeAlias = Literal[
    "SATURATED",
    "BLENDED",
    "EDGE",
    "SHAPE",
    "UNDERDETERMINED",
    "USER_EXCLUDE",
    "SELECTED",
]


def check_schema_version(value: object) -> str:
    if not isinstance(value, str):
        raise ValueError("schema_version must be a string")
    parts = value.split(".")
    if len(parts) != 3 or not all(p.isdigit() for p in parts):
        raise ValueError(f"invalid schema_version {value!r}")
    major = int(parts[0])
    if major != SCHEMA_MAJOR:
        raise ValueError(f"unknown schema MAJOR {major}")
    if value != SCHEMA_VERSION:
        raise ValueError(f"unsupported schema_version {value!r}")
    return value


SchemaVersion = Annotated[str, AfterValidator(check_schema_version)]


def check_id(value: object) -> str:
    if not isinstance(value, str) or _ID_RE.fullmatch(value) is None:
        raise ValueError(f"identifier {value!r} does not match NOR.5")
    return value


def check_term_id(value: object) -> str:
    if not isinstance(value, str) or _TERM_ID_RE.fullmatch(value) is None:
        raise ValueError(f"term_id {value!r} does not match NOR.5")
    return value


RecordId = Annotated[str, AfterValidator(check_id)]
TermId = Annotated[str, AfterValidator(check_term_id)]


def _as_float(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise TypeError("expected a number")
    out = float(value)
    if not math.isfinite(out):
        raise ValueError("expected a finite number")
    return out


def _as_int(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError("expected an integer")
    return int(value)


FiniteFloat = Annotated[float, BeforeValidator(_as_float)]
JsonInt = Annotated[int, BeforeValidator(_as_int)]


def _prior_mean(value: object) -> float | Literal["init"]:
    if value == "init":
        return "init"
    return _as_float(value)


PriorMean = Annotated[float | Literal["init"], BeforeValidator(_prior_mean)]


def _vec2(value: object) -> tuple[float, float]:
    if not isinstance(value, list | tuple):
        raise TypeError("expected a length-2 array")
    items: list[object] = list(value)
    if len(items) != 2:
        raise TypeError("expected a length-2 array")
    return (_as_float(items[0]), _as_float(items[1]))


Vec2 = Annotated[tuple[float, float], BeforeValidator(_vec2)]


def _as_f64_2d(value: object) -> npt.NDArray[np.float64]:
    if isinstance(value, np.ndarray):
        arr = np.ascontiguousarray(value, dtype=np.float64)
    elif isinstance(value, list):
        arr = np.ascontiguousarray(np.asarray(value, dtype=np.float64))
    else:
        raise TypeError("expected a 2-D float64 array")
    if arr.ndim != 2:
        raise ValueError("expected a 2-D array")
    return arr


def _as_f64_1d(value: object) -> npt.NDArray[np.float64]:
    if isinstance(value, np.ndarray):
        arr = np.ascontiguousarray(value, dtype=np.float64)
    elif isinstance(value, list):
        arr = np.ascontiguousarray(np.asarray(value, dtype=np.float64))
    else:
        raise TypeError("expected a 1-D float64 array")
    if arr.ndim != 1:
        raise ValueError("expected a 1-D array")
    return arr


def _as_u8_2d(value: object) -> npt.NDArray[np.uint8]:
    if isinstance(value, np.ndarray):
        arr = np.ascontiguousarray(value, dtype=np.uint8)
    elif isinstance(value, list):
        arr = np.ascontiguousarray(np.asarray(value, dtype=np.uint8))
    else:
        raise TypeError("expected a 2-D uint8 array")
    if arr.ndim != 2:
        raise ValueError("expected a 2-D array")
    return arr


Array2F64 = Annotated[npt.NDArray[np.float64], BeforeValidator(_as_f64_2d)]
Array1F64 = Annotated[npt.NDArray[np.float64], BeforeValidator(_as_f64_1d)]
Array2U8 = Annotated[npt.NDArray[np.uint8], BeforeValidator(_as_u8_2d)]


def f64_matrix_json(arr: npt.NDArray[np.float64]) -> list[list[float]]:
    return [[float(v) for v in row] for row in arr]


def f64_vector_json(arr: npt.NDArray[np.float64]) -> list[float]:
    return [float(v) for v in arr]


def u8_matrix_json(arr: npt.NDArray[np.uint8]) -> list[list[int]]:
    return [[int(v) for v in row] for row in arr]


def wrap_validation(exc: ValidationError, *, star_id: str | None = None) -> InputError:
    messages = []
    for err in exc.errors():
        loc = ".".join(str(part) for part in err.get("loc", ()))
        msg = str(err.get("msg", exc))
        messages.append(f"{loc}: {msg}" if loc else msg)
    joined = "; ".join(messages) if messages else str(exc)
    return InputError(joined, module="boundary", star_id=star_id)


def dump_json(model: BaseModel) -> dict[str, object]:
    dumped: dict[str, object] = model.model_dump(mode="json", exclude_none=True)
    return dumped
