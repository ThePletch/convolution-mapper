# C4 — Parameter assembly and scoping

## C4.1 Flat vector \(\theta\)

Concatenate **enabled** terms (C3) in **catalog array order**, then within each phase term there is no sub-ordering (one local scalar). Photometric `flux` then `sky` if enabled. Then kernel parameters in catalog order, each kernel’s parameters in the order listed in C3.6.

Stage 1’s \(\theta\) is **per star** and contains:

- one local coefficient \(a_k\) per enabled unfrozen **or frozen** phase term that is `enabled=true` (frozen values are held at init and omitted from the LM free set — see C4.4)
- photometric
- local kernel parameters (for kernels with `field_basis` uniform or for the local trail of `field_rotation` / `linear_drift`)

**Frozen default length** for `psf_field_v1_default` with default enable flags: 12 unfrozen phase slots (2 tilt + defocus + 2 astig + 2 coma + 2 trefoil + spherical + 2 sec-astig) + flux + sky + moffat α + moffat β + anisotropic \((\sigma_a,\sigma_b,\phi)\) + jitter σ + diffusion σ = **21** free parameters (piston frozen). (Linear drift and field rotation disabled, not counted.)

The implementation SHALL **not** hard-code 21; it SHALL count from the catalog. The number 21 is a conformance check for the default catalog.

## C4.2 Scope meanings

| Scope | Shared among | v1 fitting |
|---|---|---|
| `per_star` | one `star_id` | stage 1 free (unless frozen). **Not** field-mapped in stage 2. v1 uses this for tilt, flux, and sky. |
| `per_exposure` | all stars with the same `exposure_id` | stage 1 still fits **local** copies; stage 2 (or a post-average) combines them. v1 does **not** run a joint multi-star LM. |
| `per_session` | all stars in `session_id` | same: local in stage 1, field map in stage 2 |

v1 SHALL NOT implement joint LM across stars. Scoping determines **stage-2 grouping** and which parameters the evaluator interpolates.

## C4.3 Photometric parameters

- `flux`: `per_star`, always enabled, never field-mapped. Evaluated PSF (C7) is unit flux unless the caller passes a flux.
- `sky`: `per_star`, enabled, unfrozen, zero-mean Gaussian prior σ=5 ADU (C3.7). Residual sky after C1A subtraction; deconvolution is sensitive to wing bias.

## C4.4 Freeze / unfreeze

`Stage1Options.freeze_mask`: optional `bool` array aligned with the **free** parameter order of C4.1 **after** dropping `enabled=false` terms. If omitted, use each term’s `frozen` flag. If present, its length SHALL equal \(n_{\mathrm{free}}\). Length mismatch SHALL raise `InputError`. Implementations SHALL NOT zip-truncate.

`fit_schedule` (catalog) is a list of steps:

```
{ "name": "coarse", "unfrozen_term_ids": ["flux", "sky", "zernike_1_1", "zernike_1_m1", "zernike_2_0", "zernike_2_2", "zernike_2_m2", "moffat_seeing", "gaussian_aniso"] }
```

v1 default schedule (frozen):

1. `coarse`: `flux`, `sky`, `zernike_1_1`, `zernike_1_m1`, `zernike_2_0`, `zernike_2_2`, `zernike_2_m2`, `moffat_seeing`, `gaussian_aniso`
2. `mid`: add `zernike_3_1`, `zernike_3_m1`, `zernike_3_3`, `zernike_3_m3`, `gaussian_jitter`
3. `full`: add `zernike_4_0`, `zernike_4_2`, `zernike_4_m2`, `charge_diffusion`

Each step **starts from the previous step’s \(\theta\)**. Frozen-at-a-step parameters are held at their current value. After the last step, the recorded covariance (C5.6) is from the **full** unfrozen set’s final Jacobian, not from an intermediate step.

If `Stage1Options.use_schedule` is `false` (frozen default **true**), a single LM is run with all `frozen=false` terms free.

## C4.5 Phase diversity

`ImageMeta.known_defocus_waves` is added to the local \(a_{2,0}\) **inside the forward model** and is **not** a free parameter:

\[
a_{2,0}^{\mathrm{used}} = a_{2,0}^{\mathrm{fitted}} + \texttt{known\_defocus\_waves}
\]

The Jacobian column is still \(\partial/\partial a_{2,0}^{\mathrm{fitted}}\). This is the v1 diversity mechanism. Joint multi-exposure LM is out of scope.

## C4.6 Bounds

If `bounds` is non-null, LM parameters SHALL be reparameterized as

\[
\theta = \mathrm{lo} + (\mathrm{hi}-\mathrm{lo})\,\sigma(u), \quad \sigma(u)=\frac{1}{1+e^{-u}}
\]

and the Jacobian SHALL include \(\partial\theta/\partial u = (\mathrm{hi}-\mathrm{lo})\,\sigma(1-\sigma)\). If a bound is hit within \(10^{-12}\) of either end after mapping back, C5 sets `flag_at_bound` for that parameter.

Kernel `sigma_px`, `sigma_a_px`, `sigma_b_px`, `alpha_px`, `length_px` SHALL have `lo=0`, `hi=20` (pixels) in the default catalog. Moffat β: `lo=1.1`, `hi=8`. Kernel `angle_rad`: `lo=-\pi`, `hi=\pi`. Phase coefficients: `lo=-5`, `hi=5` waves. Flux: `lo=0`, `hi=1e12`. Sky: `lo=-1000`, `hi=1000` ADU.

## C4.7 Annotation sidecar

Every element of \(\theta\) has a `ParamMeta`: `{term_id, role, scope, frozen, unit}` where `role` is `"local_value"` | `"flux"` | `"sky"` | kernel parameter name. This sidecar is part of `Stage1Result` so C6 can assemble columns without guessing order.
