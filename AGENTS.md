# Agent guidance

This repository implements a field-varying point-spread function (PSF) model from stars in astronomical images. The v1 deliverable is a PSF evaluated at an arbitrary field position, intended as the kernel for a **downstream deconvolution** of an already-in-focus exposure — not as a collimation or optical-design calibration product.

It exposes a Python layer for star detection/extraction and user-facing diagnostics, which delegates to a Rust core for the aberration model (Zernike pupil phase, convolution kernels, and field polynomials).

The implementation contracts in `docs/contracts/` are authoritative for the algorithm's implementation. ALWAYS review them when implementing new features to ensure that your code is compliant. When a test inequality and a contract disagree, that is a defect in the contracts or the tests; do not silently reinterpret a closed-form check to make CI green.

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

## Pull requests

PRs should include a brief summary of implemented functionality, as well as any removed or altered features, and any potentially surprising details that require closer review.

PR descriptions should enumerate any additional necessary testing steps, but SHOULD NOT restate testing steps covered by CI (and should not list 'passes CI' or similar as a testing step). If there are no additional testing steps needed, the testing section should be omitted.

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
Write for a reviewer who has learned the underlying optics and statistics but should not have to recall every formula from memory. Names and signatures should still make the call site clear; comments remind the reader what a quantity *is* and why a constant has its value.

- Use the full word in comments and in identifiers unless the short form is a conventional acronym (FFT, PSF, RMS) or a symbol from the governing equation (θ, ρ, Φ). Do not truncate words ("param", "init", "coeff", "idx") when the full word fits.
- Document struct and enum fields whose purpose is not obvious from the name. For example `expected_fwhm_px` should say that it is the extraction-time seeing estimate used to initialize Moffat α, not a measurement taken from this stamp.
- Explain scientific constants and variables in comments: what they represent, their units, and why a frozen numeric value is that number.
- When a comment explains purpose, describe the concept in words. Contract requirement IDs (`C3.5.1`, `NOR.10`, and so on) are citations, not the explanation: write the physics or data rule first, then the ID in parentheses.
- Comment non-obvious control flow (efficiency, an external bug, a catalog-versus-prose disagreement). Do not narrate what the next line does, why a prior implementation was wrong, or refer to code that is not in this version. Justifications for changes belong in PR descriptions.