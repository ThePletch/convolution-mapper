# PSF Field Modeler — Implementation Contracts

This directory is the source of truth for the PSF field modeler. It is **not** an implementation plan. Every requirement is written so that a reader can judge correctness from the mathematics, data shapes, and stated numerical tolerances **without writing code**, and so that two independent implementations of the same requirement cannot reasonably diverge by “interpretation.”

## How to read these documents

1. **`00-normative-conventions.md`** freezes units, coordinate systems, indexing, versioning, and RFC 2119 language. Read it first. Later documents do not redefine those quantities.
2. Contracts **C1–C11** are independently evaluable. A requirement in one contract that depends on another cites the requirement ID (for example `C1.12`), never an informal paraphrase.
3. **`schemas/`** holds machine-readable JSON Schema drafts that are 1:1 with the pydantic / serde types. If prose and schema disagree, that is a defect in *this* directory, not a license to improvise.

## RFC 2119 language

| Word | Meaning here |
|---|---|
| **SHALL** | Mandatory. A contrary implementation is non-conformant. |
| **SHALL NOT** | Forbidden. |
| **SHOULD** | Recommended method or robustness margin. Ignoring it does not, by itself, make an implementation non-conformant, provided every SHALL still holds. A SHOULD SHALL NOT appear in a C2.8 or C10 inequality. |
| **MAY** | Allowed; the contract specifies the default when the option is omitted. |
| **frozen** | A numeric or enumerative choice that SHALL NOT be changed by implementers. Changing it is a contract revision (semver), not an implementation detail. |

Identity (units, composition, residuals, schemas, validation inequalities) is stated with SHALL. Named algorithms and extra safety margins that do not change that identity are stated with SHOULD. Where a quantity is configurable, the **default**, **allowed domain**, and **rejection rule** are stated.

## Conceptual evaluation (do this before implementing)

A contract set is conceptually correct if all of the following hold. Check them against the documents, not against code.

1. **Unit closure.** Every formula’s inputs and outputs have units declared in `00`. No implicit `2π` or pixel/arcsec conversion.
2. **Composition closure.** Phase terms compose only by addition in the pupil; kernels compose only by convolution after the pupil-to-PSF map; sampling is last. No term writes a PSF-space additive residual.
3. **Identifiability honesty.** Every documented degeneracy (defocus sign, blur stacking, tilt-vs-centroid) has a named mitigation that is either a frozen exclusion, a prior, a scope, or a reported covariance — never “the optimizer will figure it out.”
4. **Two-stage information flow.** Stage 1 consumes `C1` records and emits per-star coefficient vectors + covariances. Stage 2 consumes only those vectors/covariances plus the catalog field bases. The evaluator uses only stage-2 maps + the forward pipeline. Nothing in stage 2 re-reads pixels.
5. **Catalog-as-data.** Adding a named aberration never requires a new code path in the engine: it is a catalog row citing `(n,m)` or a kernel id plus a field basis.
6. **Validation load-bearing.** Every frozen numerical choice that could hide a bug has a closed-form or corpus check with an inequality. Those checks are specified *before* the pipeline they validate.

## File map

| ID | File | Pins |
|---|---|---|
| NOR | `00-normative-conventions.md` | Units, coordinates, dtypes, IDs, versioning |
| C1 | `C01-stage1-input.md` | Star records, image metadata, flags, serialization |
| C1A | `C01A-extraction.md` | DAOPHOT-style front-end, selection, coverage |
| C1B | `C01B-language-boundary.md` | PyO3 types, NumPy buffers, errors, Python typing |
| C2 | `C02-zernike-engine.md` | ANSI/OSA indexing, discrete RMS, closed-form tests |
| C3 | `C03-catalog.md` | ErrorTerm schema, v1 catalog, bundles, kernels |
| C4 | `C04-parameter-scoping.md` | Flat θ, scopes, freeze schedule, phase diversity |
| C5 | `C05-stage1-fit.md` | Weighted residual, analytic J, LM I/O, covariance |
| C6 | `C06-stage2.md` | Weighted linear field maps, kernel field fit, condition |
| C7 | `C07-psf-evaluator.md` | PSF at arbitrary field position + provenance |
| C8 | `C08-diagnostics.md` | Residuals, correlations, score tests, residual maps |
| C9 | `C09-forward-pipeline.md` | Pupil sampling, FFT padding, convolution, resample |
| C10 | `C10-test-corpus.md` | Synthetic images, ground truth, scoring inequalities |
| C11 | `C11-lm-selection.md` | Library decision, rejected alternatives, acceptance tests |
| L | `consistency.md` | Fifteen lemmas for paper evaluation of internal consistency |

## What v1 delivers

Given a field position `(x, y)` in millimetres (NOR.C), the evaluator SHALL return a model PSF on a stated detector pixel grid, plus the local Zernike vector and kernel parameters that produced it. Instrument-characterization reports are a future *consumer* of C6 outputs, not a v1 deliverable.

## Explicitly out of scope (v1)

These SHALL NOT be implemented. The contracts still reserve slots so they are non-breaking later.

- Mechanical sensitivity matrices with real optical-design data; collimation bundles populated with design numbers.
- Crowded-field deblending, PSF-fitting photometry, astrometric-frame catalogs.
- Joint multi-exposure LM (the scoping mechanism in C4 exists; v1 executes one exposure at a time, with known defocus offsets applied as frozen additives).
- Non-circular pupil *models* (the mask array exists from day one; v1’s shipped mask is circular unobstructed).
- Chromatic / multi-wavelength coherent superposition (v1 is monochromatic at `wavelength_m`).
- Hand-written LM, automatic differentiation, GPU kernels, Seidel (non-orthonormal) pupil bases in the engine.
