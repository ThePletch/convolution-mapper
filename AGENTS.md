# Agent guidance

This repository implements a field modeling algorithm based on star point-spread functions (PSFs) in astronomical images. It exposes a Python layer for star detection/extraction and user-facing diagnostics, which delegates to a Rust core for the main algorithm to dial in the parameters of its aberration model.

The implementation contracts in `docs/contracts/` are authoritative for the algorithm's implementation. ALWAYS review them when implementing new features to ensure that your code is compliant.

## Repository layout

| Path | Purpose |
|---|---|
| `docs/contracts/` | Normative specification (C1–C11, conventions, JSON schemas). Source of truth for requirements. |
| `crates/psf-field-core/` | Rust library: types, Zernike engine, forward pipeline, fitting, PyO3 bindings. |
| `crates/psf-field-cli/` | Command-line entry point. |
| `crates/psf-field-corpus/` | Synthetic test corpus and validation helpers. |
| `python/psf_field/` | Pure Python: extraction, I/O, models, facades, plots. |
| `python/psf_field_core/` | Thin Python package wrapping the compiled Rust extension (`_core`). |
| `tests/` | Python integration tests (pytest). |
| `.github/workflows/` | CI configuration. |

The workspace root ties Rust and Python together:

- `Cargo.toml` — Rust workspace (three crates under `crates/`).
- `pyproject.toml` — Python package (`psf-field`), maturin build config, mypy/pyright/pytest settings.
- `rust-toolchain.toml` — Pinned Rust toolchain.

When in doubt about *what* to implement, start in `docs/contracts/README.md` and follow cited requirement IDs. When in doubt about *where* code belongs, keep numerics and serde types in Rust (`psf-field-core`), keep pixel-level extraction and orchestration in Python (`psf_field`), and keep contract checks in `tests/` or `psf-field-corpus`.

## Coding style

Scientific and mathematical jargon should be used when referring to specific formulae and physical phenomena. Use acronyms only for terms that are typically referred to by acronyms within their discipline (e.g. FFT for Fast Fourier Transform). Maintain `GLOSSARY.md` under the `docs/` folder with a reference definition for all acronyms used in the project. Keep names consistent across the code base.

Write function signatures so that the purpose of each argument is clear at the call site. For instance, if a flag can be passed as an argument, it should be accepted as a keyword argument to avoid an unnamed boolean being passed as a positional argument:
```python
# bad, allows calling with the unclear signature `increment_number(3, true)`
def increment_number(number: int, logging: bool) -> int: pass

# good, requires calling as `increment_number(3, logging=true)`
def increment_number(number: int, *, logging: bool) -> int: pass
```

Define constants at the top level of the context where they are relevant rather than inline.

Functionality that is not domain-specific (e.g. argument parsing, text formatting) should live in shared libraries or modules under a directory dedicated to shared functionality. Such functionality should never exist alongside domain-specific code, even if the domain-specific code is the only place it is being used.

### Naming
Name variables, methods, and so on according to their function only. Keep names and terminology consistent with other areas of the code base that refer to the same concepts.

Avoid using analogy, metaphors, or terminology not established in general use in names. Names must be as literal as possible.

Do not use terminology invented during the planning of a feature, and do not include references to implementation plans (e.g. phases, PR numbers) in comments or names.

### Commenting
Never leave comments unless details of an implementation would surprise a first-time user, or if a first-time user would not easily understand its purpose. Always attempt to explain a function's use and purpose in its name and signature, rather than with comments.

Comments should explain non-obvious motivations if code is structured unintuitively, e.g. for efficiency reasons or to work around an external bug. They should NOT explain why a prior implementation was corrected, or refer to any part of the code base that is not present in the current version. Justifications for changes or explanations of addressed bugs belong in PR descriptions, not in the code.