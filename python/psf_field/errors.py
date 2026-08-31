"""Python exception hierarchy matching C1B.6."""

from __future__ import annotations

from typing import ClassVar, Literal

ErrorCode = Literal["input", "convergence", "numerics", "internal"]
ErrorModuleName = Literal["zernike", "pipeline", "lm", "stage2", "eval", "boundary"]

_MODULES: frozenset[str] = frozenset(
    {"zernike", "pipeline", "lm", "stage2", "eval", "boundary"}
)


class PsfFieldError(Exception):
    """Base exception for PSF field modeler failures (ingest, fit, eval, and the language boundary).

    `code` is a closed C1B.6 vocabulary (`input`, `convergence`, `numerics`, `internal`),
    not a free-form string; each subclass fixes it. `module` is one of `zernike`,
    `pipeline`, `lm`, `stage2`, `eval`, or `boundary`.
    """

    code: ClassVar[ErrorCode] = "internal"

    def __init__(
        self,
        message: str,
        *,
        module: ErrorModuleName = "boundary",
        star_id: str | None = None,
    ) -> None:
        if module not in _MODULES:
            msg = f"unknown error module {module!r}"
            raise ValueError(msg)
        super().__init__(message)
        self.message = message
        self.module = module
        self.star_id = star_id


class InputError(PsfFieldError):
    """Rejected input: schema, array shape, units, or closed-vocabulary flags."""

    code: ClassVar[ErrorCode] = "input"


class ConvergenceError(PsfFieldError):
    """Stage-1 LM finished without a successful termination (C5.8)."""

    code: ClassVar[ErrorCode] = "convergence"


class NumericsError(PsfFieldError):
    """NaN or Inf appeared in the forward model or Jacobian."""

    code: ClassVar[ErrorCode] = "numerics"


class InternalError(PsfFieldError):
    """Invariant violation or a panic that must not cross the PyO3 boundary."""

    code: ClassVar[ErrorCode] = "internal"
