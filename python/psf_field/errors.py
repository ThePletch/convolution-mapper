"""Python exception hierarchy matching C1B.6."""

from __future__ import annotations

_MODULES = frozenset({"zernike", "pipeline", "lm", "stage2", "eval", "boundary"})


class PsfFieldError(Exception):
    """Base error. `code` / `module` / `star_id` travel with the message."""

    code: str = "internal"

    def __init__(
        self,
        message: str,
        *,
        module: str = "boundary",
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
    code = "input"


class ConvergenceError(PsfFieldError):
    code = "convergence"


class NumericsError(PsfFieldError):
    code = "numerics"


class InternalError(PsfFieldError):
    code = "internal"
