# C11 — Levenberg–Marquardt library selection

This is a **decision record**, not an open evaluation. The synthetic corpus (C10) remains the acceptance test; the crate choice is frozen unless C10.5 fails after the forward model and Jacobians have passed C2.8 and C5.10.

## C11.1 Decision

**v1 SHALL use the crate `levenberg-marquardt` version `0.15.x`** (rust-cv, nalgebra backend), as a **library**, not a reimplementation.

## C11.2 Evaluation against the plan’s criteria

| Candidate | Analytic user Jacobian | Damped Gauss–Newton / LM | MINPACK heritage | FD checker | Verdict |
|---|---|---|---|---|---|
| `levenberg-marquardt` 0.15 | `LeastSquaresProblem::jacobian` required for our use; we always return `Some` | Yes; crate documents exact LM, \(\min \tfrac12\|r\|^2\) | Port of MINPACK; tests claim bit-level agreement on rank-deficient cases | `differentiate_numerically` | **Select** |
| `argmin` | Jacobian trait for Gauss–Newton | **No LM** (issue #591 still open as of 2026-08) | n/a | not the point | **Reject** — would violate “library-provided LM” |
| `cminpack` FFI | yes | Yes (true MINPACK) | Original | no | **Fallback only** if C11.6 |
| Hand-written LM | n/a | n/a | n/a | n/a | **Forbidden** (plan) |

`argmin`’s Gauss–Newton **without** LM damping SHALL NOT be used: far-from-minimum defocus/kernel problems need damping (conversation: even function, blur degeneracy).

## C11.3 Objective convention

The crate minimizes \(\frac12\sum r_i^2\). Our \(r\) is C5.1. **Do not** multiply residuals by \(\sqrt{2}\) or divide \(\chi^2\) by 2 when reporting `chi2` in C5.7: `chi2` is \(\sum r_i^2\), which is `2 * report.objective_function` at the solution. Implementations SHALL assert `|chi2 - 2 * report.objective_function| ≤ 1e-9 * (1+chi2)` at success.

## C11.4 Hyperparameters (frozen)

Construct `LevenbergMarquardt::new()` then set:

| Knob | Value | Rationale |
|---|---|---|
| `xtol` | `sqrt(f64::EPSILON)` = `1.4901161193847656e-8` | MINPACK LMDER |
| `ftol` | same | MINPACK |
| `gtol` | `0.0` | MINPACK default (do not stop on tiny gradient alone) |
| `stepbound` | `100.0` | MINPACK `factor` |
| `patience` | `100 * (n_free + 1)` as max iterations / evals — set the crate’s maximum function-evaluation limit to this integer | MINPACK `maxfev` |

If a setter name in 0.15.x differs, set the documented equivalent so the **numeric values** above are the contract, not the method names.

Do **not** enable the `minpack-compat` feature (it disables zero-residual termination and uses old constants). Our residuals are noisy; modern NaN handling is required. The numeric xtol/ftol above already match MINPACK.

## C11.5 Jacobian layout

`J` is \(m_{\mathrm{tot}} \times n_{\mathrm{free}}\), row = residual (C5.1 order), column = free parameter. nalgebra `Matrix<f64, Dyn, Dyn, …>` with dynamic dimensions. Conversion from `ndarray` is a copy at the LM boundary (C1B.9).

## C11.6 Fallback trigger

Switch to `cminpack` FFI **only if** all of the following hold:

1. C2.8, C9.13, C5.10 pass on C10.1–C10.2.
2. C10.5 recovery inequalities fail **and** a driver that feeds the **same** \(r,J\) into `scipy.optimize.least_squares(..., method="lm")` via a thin debug hook recovers the truth (proving the bug is the crate, not \(r,J\)).

That debug hook is not a Python core: it is a permitted diagnostic executable. It SHALL NOT be used in production fitting.

If (2) cannot be shown, the corpus failure is a model/init/catalog bug, not a license to change crates.

## C11.7 What SHALL NOT be done

- Wrap `scipy.optimize.least_squares` as the v1 core solver (language split).
- Use finite-difference Jacobians in the crate (`jacobian()` returning `None`).
- Use `levenberg-marquardt-sparse` (PCG / sparse) — our \(J\) is dense and small.
