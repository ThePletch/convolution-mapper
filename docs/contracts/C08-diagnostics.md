# C8 — Diagnostics

First-class v1 deliverable. Numeric outputs are produced in Rust; Python may plot them but SHALL NOT recompute the numbers.

## C8.1 Per-star residual image

`Stage1Result.residual_image` (C5.7): \(d - \mathrm{model}\) in ADU. Invalid pixels 0.

Also write `residual_weighted_image`: \(r_p\) placed back into an \(S\times S\) array, 0 on invalid pixels.

Python plotting is out of contract except: a saved FITS HDU `RESID` / `WRESID` with the same shapes.

## C8.2 Coefficient correlation report

From C5.6 at each star, plus a **session summary**:

For each pair of `term_id`s that appear in `degenerate_pairs` in **more than 20%** of converged stars, add a session-level `systematic_degeneracy` entry `{term_a, term_b, fraction, median_rho}`.

Always include the session-median of `|R_ij|` for the pairs:

- (`gaussian_jitter`, `moffat_seeing`)
- (`gaussian_jitter`, `charge_diffusion`)
- (`moffat_seeing`, `charge_diffusion`)
- (`zernike_2_0`, `moffat_seeing`)
- (`zernike_3_1`, `zernike_1_1`)
- (`zernike_3_m1`, `zernike_1_m1`)
- (`gaussian_aniso` `sigma_a_px`, `moffat_seeing` `alpha_px`)
- (`gaussian_aniso` `sigma_a_px`, `gaussian_aniso` `sigma_b_px`)

even if they are not flagged. These are the documented blur and tilt–coma degeneracies.

Session text in `diagnostics.json` / `psf-field report` SHALL state that stage-2 map \(\sigma\) ignore stage-1 off-diagonal covariances (C6.2) and are not \(\chi^2\)-inflated unless `covariance_chi2_scaled` is plotted.

## C8.3 Stage-2 residual maps

For each phase `term_id`, a table with columns `u`, `v`, `a_stage1`, `a_hat`, `residual`, `sigma`. This is C6.6 `residuals` plus coordinates.

**Structure test (optional flag, computed always):** fit a **next-degree** monomial (if current degree is 0, fit degree 1; if 1, fit degree 2; if 2, skip) to the stage-2 residuals with the same weights. If that extra fit reduces \(\chi^2\) by a p-value < 0.01 under a χ² difference test with \(\Delta\mathrm{dof} =\) number of new monomials, set `basis_underspecified=true`.

χ² difference test: \(\Delta\chi^2 \sim \chi^2_{\Delta\mathrm{dof}}\) under the null. p-value from the regularized gamma function. Frozen α = 0.01.

## C8.4 Score tests (missing-term ranking)

### C8.4.1 Candidate set

The candidate list is the set of phase `(n,m)` with \(n \le n_{\max}\), valid C2.1, **not** already enabled in the catalog, plus kernel ids in C3.6 that are `enabled=false`.

Frozen \(n_{\max}=7\). Tilt `(1,±1)` is included as a candidate **only if** those terms are `enabled=false` in the catalog (they diagnose centroid error leaking into the model when tilt is not already a free parameter).

### C8.4.1a Which stars are scored (`ScoreOptions`)

`ScoreOptions.max_stars` default **50**, allowed `[1, max_selected]`. If more stars converged, score `max_stars` of them: highest SNR **subject to filling as many C1A.12 cells as possible** (same per-cell quota idea as C1A.11.2). Session `weak_phase_all_fraction` is over **scored** stars, not all stars. Setting `max_stars` to the number of converged stars scores everyone.

Score tests SHALL reuse the fitted \(U\) and \(Z^{(k)}\) from the last forward at that \(\theta\). They SHALL NOT rerun LM.

### C8.4.2 Sensitivity image

For a candidate phase mode \(k\), evaluate \(\partial m / \partial a_k\) at the **fitted** \(\theta\) (C9.8–C9.10), unit flux derivative, shape \(S\times S\). For a candidate kernel, \(\partial m / \partial\alpha\) at a 1-pixel (or 0.1 rad) unit of that kernel’s first parameter, with other candidate params at 0.

Restrict to valid pixels, weight:

\[
s_p = \frac{1}{\sigma_p}\frac{\partial m_p}{\partial\alpha}
\]

Residual \(r_p\) from C5.1 (data side only, no prior rows).

### C8.4.3 Score

\[
\mathrm{score}_k = \frac{\sum_p r_p s_p}{\sqrt{\sum_p r_p^2}\,\sqrt{\sum_p s_p^2}}
\]

This is the cosine similarity in \(\mathbb{R}^m\). Range \([-1,1]\).

If \(\|s\|_2 < 10^{-15}\), `score=0`, `undefined=true`.

### C8.4.4 Ranking and flags

Rank candidates by `|score|` descending.

| Flag | Condition |
|---|---|
| `suggest_add` | `|score| ≥ 0.15` |
| `weak_phase_all` | every **phase** candidate has `|score| < 0.08` **and** \(\chi^2_{\mathrm{red}} > 2\) **and** residual image has C8.5 structure |

`weak_phase_all` means: look at kernels or amplitude, not more Zernikes.

If tilt is **enabled**, set `centroid_leak_suspect=true` when that star's `|R|` for (`zernike_3_1`, `zernike_1_1`) or (`zernike_3_m1`, `zernike_1_m1`) is \(\ge 0.90\), or when \(\sqrt{a_{1,1}^{2}+a_{1,-1}^{2}} > 3\sigma\) of the tilt prior. If tilt is **disabled**, the same flag is set when the tilt candidate `zernike_1_1` or `zernike_1_m1` has \(\lvert\mathrm{score}\rvert \ge 0.15\). Selection SHALL NOT change the catalog in response.

## C8.5 Residual spatial structure

Let \(e_{ij}\) be the weighted residual image. Compute the azimuthal RMS in 4 radial bins (equal width from 0 to \(S/2\)). If the ratio of the max bin RMS to the min bin RMS is \(> 2\) and \(\chi^2_{\mathrm{red}}>1.5\), set `structured_residual=true`.

## C8.6 Gram matrix

C2.5 matrix, written once per session as `zernike_gram.json`.

## C8.7 `ScoreReport` / `DiagnosticsBundle`

One JSON object per session:

- `per_star`: list of `{star_id, chi2_reduced, structured_residual, centroid_leak_suspect, scores: [{term_id, score, suggest_add}]}`
- `session_degeneracies`: C8.2
- `stage2_maps`: C8.3 flags
- `weak_phase_all_fraction`: fraction of stars with `weak_phase_all`

## C8.8 Report files (`psf-field report`)

The CLI subcommand `report` (C1B.7) SHALL write, from already-computed C5/C6/C8 numbers (Python SHALL NOT recompute them):

- FITS HDUs `RESID` / `WRESID` (C8.1)
- `zernike_gram.json` (C8.6)
- `diagnostics.json` (`DiagnosticsBundle`, including C7.5 flags when an eval grid is supplied)
- a Markdown summary: `n_converged`, median \(\chi^2_{\mathrm{red}}\), C8.2 pairs, `empty_cells`, extrapolation flags

## C8.9 What diagnostics SHALL NOT do

- Auto-enable catalog terms (no closed-loop catalog mutation).
- Thresholds other than those above.
- Score tests in unweighted pixel space (must use \(1/\sigma\)).
