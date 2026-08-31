# C1B — Language boundary (PyO3) and Python typing

## C1B.1 Split of responsibility

| Side | Language | Owns |
|---|---|---|
| Extraction, FITS/WCS I/O, plotting, notebook driver, pydantic ingest | Python 3.12+ | C1A, C8 plots, C7 caller |
| Zernike engine, forward pipeline, LM, stage-2, score tests, evaluator | Rust | C2, C5, C6, C7, C8 (numeric), C9 |

There is **no** Python reference implementation of C2/C5/C9. Closed-form cases (C2.8), the FD harness (C5.10), and C10 are load-bearing.

Python exists solely as an ecosystem adapter. Dynamic typing is not a design affordance.

## C1B.2 Extension module

- PyO3 + maturin.
- **Frozen module name:** `psf_field_core`.
- Import: `import psf_field_core`.
- ABI: build against CPython 3.12+; no PyPy.

## C1B.3 Exported functions (complete v1 surface)

All arguments that are arrays SHALL be `numpy.ndarray` dtype `float64`, C-contiguous. Python wrappers SHALL copy to C-contiguous `float64` if needed **before** the call, never inside Rust via silent reinterpret.

```
fit_stage1(star: StarRecordDict, pupil: PupilSpecDict, catalog: CatalogDict, options: Stage1OptionsDict) -> Stage1ResultDict
fit_stage1_batch(stars: list[StarRecordDict], ...) -> list[Stage1ResultDict]
fit_stage2(results: list[Stage1ResultDict], catalog: CatalogDict, field: FieldConfigDict) -> Stage2ResultDict
evaluate_psf(stage2: Stage2ResultDict, field_xy_mm: tuple[float, float], grid: EvalGridDict) -> PsfEvalDict
forward_psf(theta_local: ndarray, pupil: PupilSpecDict, catalog: CatalogDict, image_meta: ImageMetaDict, centroid_xy_px: tuple[float, float], stamp_size: int) -> ndarray  # shape (S,S), unit flux
score_tests(residual: ndarray, variance: ndarray, pixel_mask: ndarray, theta_local: ndarray, catalog: CatalogDict, candidates: list[str], pupil: PupilSpecDict, image_meta: ImageMetaDict) -> ScoreReportDict
validate_jacobian(star: StarRecordDict, pupil: PupilSpecDict, catalog: CatalogDict, theta: ndarray) -> FdReportDict
```

`fit_stage1_batch` SHALL be embarrassingly parallel over stars: each star owns its problem, with no shared mutable arrays. Python threads SHALL NOT run the hot loop. The parallel iterator SHOULD be `rayon`.

Dict types above are the JSON-object shapes of the schemas in `schemas/`. The Python public API SHALL accept pydantic models and convert via `.model_dump(mode="python")` at the wrapper; the PyO3 signatures MAY accept bound pyclasses that mirror those models.

## C1B.4 Versioned pyclasses (preferred over raw dicts)

The extension SHALL also expose immutable pyclasses with the same field names as the JSON schemas: `StarRecord`, `ImageMeta`, `PupilSpec`, `Catalog`, `Stage1Result`, `Stage2Result`, `PsfEval`, `ScoreReport`, `FdReport`.

`schema_version` is a required attribute. Constructing a pyclass with the wrong MAJOR SHALL raise `InputError`.

## C1B.5 NumPy buffer convention

| Rule | Requirement |
|---|---|
| Dtype | `float64` for real arrays; pupil FFT internals may use `complex128` but SHALL NOT cross the boundary as complex except `forward_complex_pupil` which is **not** v1-exported |
| Contiguity | C-contiguous on entry; Rust SHALL reject (error) non-contiguous buffers rather than silently copy |
| Writeability | Inputs borrowed read-only. Outputs are new arrays owned by Python |
| Byte order | native little-endian on the platforms we ship; big-endian SHALL be rejected |
| Shape | 2-D stamps `(S, S)`; Jacobian as `(m, n)` with `m` = valid pixels, `n` = free parameters (C5) — Jacobian is **not** exported in v1 except inside `FdReport` |
| Zero-copy | Allowed when the NumPy array is already `float64` C-contiguous; not required |

## C1B.6 Error mapping

Rust panics SHALL NOT cross the boundary. A panic SHALL be mapped to `InternalError`. Entry points SHOULD use `catch_unwind` to enforce that.

| Rust condition | Python exception type (in `psf_field_core`) | `code` string |
|---|---|---|
| Schema / shape / unit / flag vocabulary | `InputError` | `input` |
| LM termination not successful (C5.8) | `ConvergenceError` | `convergence` |
| NaN/Inf in forward or J | `NumericsError` | `numerics` |
| Panic / invariant bug | `InternalError` | `internal` |

Each exception SHALL carry:

- `code: str` as above
- `module: str` — one of `zernike`, `pipeline`, `lm`, `stage2`, `eval`, `boundary`
- `message: str` — actionable, no backtrace required
- `star_id: str | None`

Python `Exception` hierarchy:

```
PsfFieldError
├── InputError
├── ConvergenceError
├── NumericsError
└── InternalError
```

`fit_stage1_batch` SHALL NOT abort the whole batch on one star’s `ConvergenceError`; that star’s result SHALL have `converged=false` and the error payload (C5.7). `InputError` on schema SHALL abort the batch.

## C1B.7 CLI without Python

The Rust crate SHALL also be a binary `psf-field` that reads C1 FITS/Parquet and writes C5/C6/C7 FITS/Parquet. This is the same schema as PyO3. The front-end is therefore swappable in practice.

Subcommands (frozen): `stage1`, `stage2`, `eval`, `score`, `check-jacobian`, `report`.

`report` writes C8.8 artifacts from existing C5/C6/C8 outputs. Python plotting MAY render those files and SHALL NOT recompute the numbers.

## C1B.8 Python typing policy (mandatory)

1. Every public and internal function, method, and class in this repository’s Python packages SHALL have complete parameter and return annotations. Untyped Python is a defect.
2. CI SHALL run both:
   - `mypy --strict` on the Python package
   - `pyright` with `"typeCheckingMode": "strict"`
   Either failure fails CI.
3. `typing.Any` SHALL NOT appear. `cast` to `Any` is forbidden. `# type: ignore` SHALL include a trailing error code and a one-line reason citing a ticket or a third-party bug; bare ignore is forbidden.
4. Boundary models SHALL be pydantic v2 `BaseModel` with `model_config = ConfigDict(extra="forbid", strict=True, frozen=True)`.
5. Arrays SHALL be annotated as `numpy.typing.NDArray[np.float64]` (or `np.uint8` for masks), never as untyped `np.ndarray`. Rank SHALL be expressed with `Annotated` + a pydantic validator on models, and with `TypeGuard`/shape asserts at the facade when calling photutils.
6. photutils / untyped astropy SHALL be imported only inside facade modules listed in `pyproject.toml` under a named package `psf_field.facades`. All other modules SHALL import from those facades or from our pydantic models.

### C1B.8.1 Frozen mypy config

```
[tool.mypy]
python_version = "3.12"
strict = true
warn_unused_ignores = true
disallow_any_generics = true
disallow_incomplete_defs = true
disallow_untyped_defs = true
disallow_untyped_calls = true
no_implicit_reexport = true
```

### C1B.8.2 Frozen pyright config

```
[tool.pyright]
pythonVersion = "3.12"
typeCheckingMode = "strict"
```

## C1B.9 Array-library split in Rust

- Pipeline images, pupils, PSFs: `ndarray` `Array2<f64>` / `Array2<Complex64>`.
- LM interface: `nalgebra` types required by `levenberg-marquardt` (C11). Conversion happens only at the LM call boundary.
- FFT: `rustfft` on `complex64` buffers.

This is not a preference; it is the unique combination implied by C9 + C11.

## C1B.10 Parallelism contract

- Per-star stage-1: no shared mutable arrays; each star owns its problem struct. The parallel iterator SHOULD be `rayon`.
- Stage-2: single-threaded dense linear algebra (problem size is tiny: tens of stars × few basis columns).
- Python SHALL NOT hold the GIL during `fit_stage1_batch`. The wrapper SHOULD release it with PyO3 `py.allow_threads`.
