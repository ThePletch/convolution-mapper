# C6 — Stage-2 field regression

Stage 2 does **not** read pixels. It fits field maps to the cloud of stage-1 local coefficients.

## C6.1 Inputs

- All `Stage1Result` with `converged=true` and whose stars are fit-eligible (C1.7).
- Catalog field bases (C3.3).
- Stars’ `field_xy_mm` → \((u,v)\) via NOR.9.

Stars with `covariance_ok=false` SHALL be excluded from stage 2 (listed in `dropped_star_ids`).

## C6.2 Per-coefficient linear problem (phase terms)

For each enabled phase `term_id` with field basis monomials \(B_j(u,v)=u^{i}v^{j}\), \(j=1\ldots q\):

Let \(a_s\) be the local coefficient of star \(s\), \(\sigma_s^{2} = C_{kk}\) the corresponding diagonal of that star’s free-parameter covariance (C5.6). If that term was frozen in stage 1, skip the term in stage 2.

Weight \(w_s = 1/\sigma_s^{2}\). If \(\sigma_s^{2}\) is non-finite or \(\le 0\), drop the star for this term.

Design matrix \(X_{s j} = B_j(u_s, v_s)\). Solve weighted least squares

\[
\min_c \sum_s w_s \bigl(a_s - X_{s\cdot} c\bigr)^{2}
\]

Normal equations \( (X^{\top} W X) c = X^{\top} W a \) with \(W=\mathrm{diag}(w)\). Solve by Cholesky of \(X^{\top}WX\). If that fails, use SVD with relative singular-value cutoff \(10^{-8}\); set `used_svd=true`.

Covariance of \(c\): \(C_c = (X^{\top} W X)^{-1}\) (Gauss–Markov, assuming \(a_s\) independent — v1 **ignores** stage-1 off-diagonal covariances between different terms; this is frozen). Cross-term correlations at stage 1 are reported in C8, not folded into stage 2 in v1.

## C6.3 Conditioning report (required)

For each term:

| Field | Definition |
|---|---|
| `cond2` | \(\kappa_2(X^{\top}WX)\) |
| `singular_values` | of \(W^{1/2}X\), descending |
| `separability_notes` | if `zernike_2_0` and `cond2 > 1e8`, set `ill_conditioned=true` |

**Frozen flag:** `ill_conditioned` iff `cond2 > 1e8` or SVD dropped a mode.

v1 SHALL still emit the evaluator using the SVD-truncated \(c\) (truncated modes set to 0). The flag tells the user the tilt/curvature/piston split may be meaningless.

## C6.4 Kernel field parameters

**Uniform kernels** (`moffat_seeing`, `gaussian_jitter`, `charge_diffusion`, `linear_drift` if enabled): stage 2 takes the **inverse-variance weighted mean** of the per-star values. No spatial polynomial. Report the weighted RMS about that mean as `scatter`.

**`field_rotation`:** fit `(center_x_mm, center_y_mm, omega)` to per-star trail vectors \((L_s, \phi_s)\) from stage 1’s local `linear_drift`-equivalent parameters.

Model: predicted length \(L(x,y)=|\omega| T R_\perp / p_{\mathrm{mm}}\), angle = tangent (C3.6.4). Residual: \((L_s\cos\phi_s - L_{\mathrm{pred}}\cos\phi_{\mathrm{pred}},\; L_s\sin\phi_s - \ldots) / \sigma_{L,s}\).

This is a **3-parameter** nonlinear least squares, same LM crate (C11), max 20 iterations, init: center = flux-weighted mean of star positions, \(\omega=0\). If `linear_drift` is the enabled local model instead of `field_rotation`, skip this fit.

## C6.5 Bundles

If `bundle.matrix` is null (v1), skip. Non-null is rejected in v1 (C3.7.1).

## C6.6 `Stage2Result`

| Field | Content |
|---|---|
| `maps` | dict `term_id` → `{ "coefficients": f64[q], "covariance": f64[q,q], "cond2": f64, "ill_conditioned": bool, "residuals": f64[n_stars] }` |
| `kernel_globals` | dict of weighted means / rotation fit |
| `dropped_star_ids` | list |
| `n_stars_used` | int |
| `param_meta` | copy from catalog assembly |

`residuals` in each map: \(a_s - \hat a(u_s,v_s)\), same order as `star_ids_used`.

## C6.7 What stage 2 SHALL NOT do

- Re-run FFTs or LM on pixels.
- Fit a joint nonlinear field model for phase terms (linear only).
- Use unweighted OLS (weights are mandatory).
