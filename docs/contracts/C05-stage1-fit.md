# C5 — Stage-1 fit (per star)

Levenberg–Marquardt on one star. Library: C11. Loss: weighted least squares on valid pixels (C1.6).

## C5.1 Residual vector

Let \(\mathcal{V}\) be the list of valid pixel indices in **row-major** order (frozen). \(m = |\mathcal{V}|\). \(n\) = number of free (unfrozen) parameters.

For each valid pixel \(p\):

\[
r_p = \frac{d_p - \mathrm{model}_p(\theta)}{\sigma_p}
\]

where \(\sigma_p = \sqrt{\mathrm{variance}_p}\), \(d\) is `stamp`, `model` is C9.11.

Append prior rows (C3.4) for each free parameter that has a Gaussian prior. The residual vector length is \(m + n_{\mathrm{prior}}\).

If \(m < n\), do not call LM; return `UNDERDETERMINED` (C1.7). If \(m < n + 8\), the fit SHOULD still return `UNDERDETERMINED` without calling LM (robustness margin over well-posedness).

## C5.2 Jacobian

\(J_{p j} = \partial r_p / \partial \theta_j\) with \(\theta_j\) the **unconstrained** LM variable (after C4.6 mapping, \(J\) is w.r.t. \(u\), not \(\theta\)).

For a phase coefficient \(a_k\) (waves), before bound mapping:

\[
\frac{\partial r_p}{\partial a_k} = -\frac{F}{\sigma_p}\frac{\partial m_p}{\partial a_k}
\]

where \(\partial m/\partial a_k\) is C9.8 chained through C9.9–C9.10.

For flux:

\[
\frac{\partial r_p}{\partial F} = -\frac{m_p}{\sigma_p}
\]

For sky:

\[
\frac{\partial r_p}{\partial b} = -\frac{1}{\sigma_p}
\]

For kernel parameter \(\alpha\): same as phase, using C9.9 kernel derivatives.

Prior row for parameter \(j\): \(J_{\mathrm{prior}, j} = 1/\sigma_{\mathrm{prior}}\) (and 0 elsewhere).

**Finite differences SHALL NOT be used inside LM.** They are C5.10 only.

## C5.3 Initialization

Apply C3.5 per term. Then apply C4.4 schedule. \(\theta_0\) after inits, before step 1, SHALL be written to `Stage1Result.theta_init`.

Defocus sign: positive (C3.5.1). `Stage1Result.defocus_sign_ambiguous = true` always for single-exposure v1 when `|known_defocus_waves| < 0.05`. If `|known_defocus_waves| ≥ 0.05`, set `false`.

## C5.4 What is not fitted

- Stamp centroid (C1.2.3 is an input).
- Tilt modes when disabled (C3.7).
- `known_defocus_waves`.

## C5.5 LM hyperparameters (frozen)

See C11.4. The problem SHALL implement `LeastSquaresProblem` with **user Jacobian always `Some`**. Returning `None` is non-conformant.

Max function evaluations: \(100(n+1)\) with \(n\) = free parameters. Patience / other crate knobs: C11.4.

## C5.6 Covariance

Let \(J\) be the Jacobian at the **final** \(\theta\) of the last schedule step, including prior rows, w.r.t. the physical \(\theta\) (undo C4.6 by chain rule so reported covariance is in waves / ADU / pixels).

Let \(H = J^{\top} J\) (\(n\times n\)). If \(\lambda_{\min}(H) / \lambda_{\max}(H) < 10^{-12}\) or a Cholesky fails, set `covariance_ok=false`, store `covariance` as NaN matrix, and still return \(\theta\).

Otherwise \(C = H^{-1}\) (the Gauss–Newton / Fisher approximation for unit-weight residuals \(r=(d-\mu)/\sigma\)). This \(C\) is the covariance of \(\theta\) in physical units.

Correlation:

\[
R_{ij} = \frac{C_{ij}}{\sqrt{C_{ii}C_{jj}}}
\]

If \(C_{ii}\le 0\), that row/column of \(R\) is NaN.

**Degenerate pair:** `|R_ij| ≥ 0.90` for \(i\neq j\). Listed in `degenerate_pairs`.

## C5.7 `Stage1Result` fields

| Field | Type |
|---|---|
| `schema_version` | `"1.0.0"` |
| `star_id` | string |
| `theta` | `f64[n_all]` physical, including frozen (frozen = init) |
| `theta_init` | same shape |
| `param_meta` | `ParamMeta[n_all]` |
| `free_index` | `int64[n]` indices into `theta` |
| `covariance` | `f64[n, n]` for **free** parameters, C5.6 order = `free_index` |
| `covariance_ok` | bool |
| `correlation` | `f64[n, n]` |
| `degenerate_pairs` | list of `[term_id_a, term_id_b, rho]` |
| `residual_image` | `f64[S,S]` = \(d - \mathrm{model}\) (not divided by σ); invalid pixels 0 |
| `weighted_residual` | `f64[m]` \(r_p\) |
| `chi2` | \(\sum r_p^{2}\) including priors |
| `chi2_reduced` | `chi2 / max(m + n_prior - n, 1)` |
| `n_iter` | int |
| `n_fev` | int |
| `termination` | string, C5.8 |
| `converged` | bool |
| `defocus_sign_ambiguous` | bool |
| `flag_at_bound` | list of `term_id` |
| `error` | `null` or `{code, message}` |

## C5.8 Termination mapping

Map `levenberg_marquardt::TerminationReason` as follows. `converged = termination.was_successful()` from the crate, **except** `Orthogonal` / small gradient at a point with `chi2_reduced > 50` SHALL set `converged=false` (landed in a bad basin).

Serialized `termination` strings (closed set): `converged_ftol`, `converged_xtol`, `converged_gtol`, `converged_zero_residual`, `max_eval`, `numerical`, `lost_patience`, `user`, `unknown`.

`was_successful()==true` maps to the `converged_*` family. Exact crate variant names SHALL be translated by a match that fails to compile if a new variant appears (`non_exhaustive` handled via `_ => "unknown"` plus `converged=false`).

## C5.9 Invalid residual / Jacobian

If any `model` pixel among valid pixels is non-finite, or any Jacobian entry is non-finite, `residuals()` / `jacobian()` SHALL return `None`, which the crate treats as user termination. Result: `termination="user"`, `converged=false`, `error.code="numerics"`.

## C5.10 Finite-difference harness (`validate_jacobian`)

Not used in fitting.

- Method: **central** difference, step \(h_j = 1.4901161193847656\times 10^{-8} \cdot \max(1, |\theta_j|)\) (sqrt of f64 epsilon).
- Compare columns of analytic \(J\) (physical \(\theta\), valid pixels only, no prior rows) to FD.
- Metric per column: \(\|J_j^{\mathrm{an}} - J_j^{\mathrm{fd}}\|_2 / \max(\|J_j^{\mathrm{fd}}\|_2, 10^{-12})\).
- **Pass:** every **unfrozen phase and photometric** column relative error \(< 10^{-4}\). Kernel columns: \(< 5\times 10^{-4}\) (less smooth truncation).
- Defocus column at \(a_{2,0}=0\) is exempt from the relative test; it SHALL instead satisfy C2.8.4.
- `FdReport` stores per-column errors and `passed: bool`.

This harness SHALL be run in CI on C10.1 unaberrated and on C10.2 single-mode cases (defocus 0.3, coma 0.3).

## C5.11 Parallelism

`fit_stage1_batch` fits stars independently. Covariances are per-star; no cross-star terms in v1.
