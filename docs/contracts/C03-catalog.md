# C3 — ErrorTerm semantic catalog

The engine never switches on a string like `"coma"`. The catalog maps names to `(pupil | kernel)` slots, field bases, priors, and freeze flags.

## C3.1 Catalog file

- JSON, UTF-8. Schema: `schemas/catalog.schema.json`.
- Required top-level: `schema_version`, `catalog_id`, `terms` (array), `bundles` (array, may be empty), `fit_schedule` (C4.5).
- **v1 default catalog id:** `psf_field_v1_default`.
- Duplicate `term_id` SHALL be rejected.

## C3.2 Discriminated union: `ErrorTerm`

Every term has:

| Field | Type | Rule |
|---|---|---|
| `term_id` | string | NOR.5 |
| `name` | string | human label, e.g. `"Primary coma (cos)"` |
| `kind` | `"phase"` \| `"kernel"` \| `"photometric"` | |
| `scope` | `"per_star"` \| `"per_exposure"` \| `"per_session"` | C4 |
| `frozen` | bool | if true, held at init in LM (still present in \(\theta\)) |
| `enabled` | bool | if false, omitted from \(\theta\) and from C9 entirely |
| `bounds` | `[f64, f64] \| null` | inclusive; null = unbounded |
| `init` | `InitSpec` | C3.5 |
| `prior` | `PriorSpec` | C3.4 |
| `units` | string | `"waves"` for phase; kernel-specific (C3.6) |
| `report` | `ReportSpec` | C3.8 |

### C3.2.1 `kind = "phase"` additional fields

| Field | Type | Rule |
|---|---|---|
| `n` | int | C2.1 |
| `m` | int | C2.1 |
| `field_basis` | `FieldBasis` | C3.3 |

Exactly one `(n,m)` per term. Cosine and sine are **two terms** (`zernike_3_1` and `zernike_3_m1`), not a polar pair inside one term. Polar reporting is C3.8.

### C3.2.2 `kind = "kernel"` additional fields

| Field | Type | Rule |
|---|---|---|
| `kernel` | `KernelSpec` | C3.6 |
| `field_basis` | `FieldBasis` \| `null` | null ⇒ spatially uniform (still may be per_exposure) |

## C3.3 Field basis

```
FieldBasis = {
  "family": "monomial",
  "degree": int,          # 0, 1, or 2 in v1
  "terms": [ [i, j], ... ]  # include only these monomials u^i v^j
}
```

- `family` SHALL be `"monomial"` in v1. Other families (`"zernike_field"`, `"design_sensitivity"`) are reserved and SHALL be rejected in v1.
- `degree` SHALL equal \(\max(i+j)\).
- `terms` SHALL list pairs \((i,j)\) of non-negative integers with \(i+j \le\) `degree`, **without duplicates**, sorted lexicographically by `(i+j, i)` (frozen sort).
- The local physical coefficient at field \((u,v)\) is

\[
a(u,v) = \sum_{(i,j) \in \mathrm{terms}} c_{ij}\, u^i v^j
\]

- Stage 1 (C5) fits the **local** \(a\) (one scalar per term per star), not the \(c_{ij}\).
- Stage 2 (C6) fits the \(c_{ij}\).
- If `terms = [[0,0]]` only, the field is constant.

**Sensor tilt** is not a second pupil mode: it is `zernike_2_0` with `terms` containing `[0,0], [1,0], [0,1]` (and optionally quadratics). The name `"sensor_tilt"` is a **bundle or report alias** over the linear subset of the defocus map (C3.7).

## C3.4 Priors

```
PriorSpec = {
  "family": "none" | "gaussian",
  "mean": f64,      # required if gaussian; on the LOCAL coefficient a (waves or kernel units)
  "sigma": f64      # > 0; required if gaussian
}
```

Gaussian prior contributes extra residual rows in C5 and extra normal-equation terms in C6:

\[
r_{\mathrm{prior}} = \frac{a - \mu}{\sigma}
\]

with Jacobian \(1/\sigma\). If `family="none"`, no extra rows.

v1 SHALL apply priors to **stage 1 local parameters** and to **stage 2 field coefficients** as follows:

- Stage 1: prior on local \(a\).
- Stage 2: prior on each \(c_{ij}\) with the same \(\mu,\sigma\) **only if** `field_basis` is a single constant term. If the basis has more than one monomial, stage-2 priors are **off** unless `PriorSpec` includes optional `"stage2": { "mean": [..], "sigma": [..] }` aligned with `terms`. The v1 default catalog does not set `stage2` priors except on kernel uniform terms.

## C3.5 Init estimators

```
InitSpec = { "method": "zero" | "flux_sum" | "defocus_moment" | "moffat_fwhm" }
```

| method | Applies to | Formula |
|---|---|---|
| `zero` | anything | \(a_0 = 0\) if `prior.family=="none"`; else \(a_0 =\) `prior.mean` |
| `flux_sum` | flux only (C4.3) | \(F_0 = \max(\texttt{flux\_sum\_adu}, 0)\) |
| `defocus_moment` | `zernike_2_0` only | C3.5.1 |
| `moffat_fwhm` | `moffat_seeing` α | C3.5.2 |

Any other pairing SHALL be rejected at ingest.

### C3.5.1 Defocus from second moments

On valid stamp pixels, flux-weighted second moment about `centroid_xy_px`:

\[
\sigma_{\mathrm{meas}}^{2} = \frac{\sum_{p} w_p\, r_p^{2}\, \mathrm{stamp}_p^{+}}{\sum_p w_p\, \mathrm{stamp}_p^{+}}
\]

where \(\mathrm{stamp}^+ = \max(\mathrm{stamp},0)\), \(w_p=1\) on valid pixels, \(r_p\) in pixels from the stamp centroid.

Let \(\sigma_0^{2}\) be the same moment computed on the **forward model** with all aberrations 0, flux 1, \(b=0\), no kernels, same centroid and `S` (one call to `forward_psf`). Then

\[
|a_{2,0}|_0 = \mathrm{clip}\bigl( 0.35 \cdot \max(\sigma_{\mathrm{meas}}^{2} - \sigma_0^{2}, 0) / \max(\sigma_0^{2}, 10^{-12}),\; 0,\; 2.0 \bigr)
\]

**Sign: always \(+\)**. Units: waves. The factor `0.35` is frozen (maps extra second-moment width to waves for the C10.1 camera; C10.3 requires recovery from this init, not that 0.35 be physically universal). Clip at 2 waves.

### C3.5.2 Moffat α from FWHM

If extraction config `fwhm` is available on the Python side, pass it as `StarRecord` optional `aux` is **not** allowed into Rust. Instead `Stage1Options.expected_fwhm_px` (required) is used:

Moffat with frozen \(\beta=2.5\): \(\mathrm{HWHM} = \alpha \sqrt{2^{1/\beta}-1}\). Set \(\mathrm{FWHM} =\) `expected_fwhm_px`, \(\alpha_0 = \mathrm{FWHM} / (2 \sqrt{2^{1/\beta}-1})\).

## C3.6 Kernel catalog (closed set for v1)

v1 SHALL implement exactly these `kernel.id` values. Unknown ids SHALL be rejected.

Spatial radius \(R\) on the FFT grid is in FFT pixels. \(s\) below is in **detector pixels**, converted by \(r\) (C9.7).

### C3.6.1 `gaussian_iso` (jitter / diffusion)

Parameters: `sigma_px` > 0.

\[
K(x,y) = \exp\bigl(-(x^{2}+y^{2})/(2\sigma^{2})\bigr)
\]

with \(\sigma = \texttt{sigma\_px} \cdot r\) in FFT pixels. Truncate at \(5\sigma\) (set \(K=0\) beyond). \(\partial K/\partial\sigma = K \cdot (x^{2}+y^{2})/\sigma^{3}\).

### C3.6.2 `moffat_iso` (seeing)

Parameters: `alpha_px` > 0, `beta` (frozen **2.5** unless `frozen=false` on beta; v1 default catalog freezes beta).

\[
K(x,y) = \bigl(1 + (R/\alpha)^{2}\bigr)^{-\beta}, \quad R^{2}=x^{2}+y^{2}
\]

\(\alpha = \texttt{alpha\_px}\cdot r\). Truncate where \(K < 10^{-6} K(0,0)\) or \(R > 15 r\) (15 detector pixels), whichever is smaller.

Closed form:

\[
\frac{\partial K}{\partial\alpha} = 2\beta (R^{2}/\alpha^{3}) \bigl(1+(R/\alpha)^{2}\bigr)^{-\beta-1}
\]

\[
\frac{\partial K}{\partial\beta} = -K \ln\bigl(1+(R/\alpha)^{2}\bigr)
\]

### C3.6.3 `linear_drift`

Parameters: `length_px` ≥ 0, `angle_rad` (radians, detector \(+x\) toward \(+y\)).

\(K\) is a unit-sum Gaussian line segment of length \(L = \texttt{length\_px}\cdot r\), σ_perp = \(0.5 r\) (half a detector pixel, frozen), along direction \((\cos\phi, \sin\phi)\).

Construct: for FFT pixel offset \((x,y)\) from center, \(t = x\cos\phi + y\sin\phi\), \(n = -x\sin\phi + y\cos\phi\), \(K = \exp(-n^{2}/(2\sigma_\perp^{2})) \cdot \mathbf{1}_{|t|\le L/2}\). If \(L=0\), \(K\) is isotropic Gaussian with \(\sigma=0.5 r\). \(K\) SHALL be \(C^1\) in \(L\) (the indicator SHALL NOT be used raw in the Jacobian). Implementations SHOULD replace \(\mathbf{1}_{|t|\le L/2}\) by

\[
\frac12\bigl(\mathrm{erf}((L/2-t)/w)+\mathrm{erf}((L/2+t)/w)\bigr), \quad w=0.25 r.
\]

### C3.6.4 `field_rotation`

Global parameters (scope `per_exposure`): `center_x_mm`, `center_y_mm`, `omega` (rad/s).

At a star with field \((x,y)\) mm, trail length in mm is \(|\omega| \cdot T \cdot R_\perp\) where \(R_\perp\) is distance from `(center_x_mm, center_y_mm)` and \(T=\) `exptime_s`. Convert length to detector pixels via \(p_{\mathrm{mm}}\). Direction is tangential: \(\hat\theta = (-(y-y_c), x-x_c) / R_\perp\) (sign of \(\omega\) included). If \(R_\perp < 10^{-6}\,\mathrm{mm}\), length 0.

Then apply C3.6.3 with that length and angle. This kernel is **field-dependent**; stage 1 fits a local `(length_px, angle_rad)` pair; stage 2 fits the three globals (C6.4).

### C3.6.5 `periodic_error`

Parameters: `amp_px` ≥ 0, `period_s` > 0, `phase_rad`.

Trajectory: \(x(t) = A \sin(2\pi t / P + \varphi)\), \(y=0\) in detector pixels along \(+x\) (mount RA; frozen axis). \(A=\texttt{amp\_px}\). \(K\) SHALL be the time-average of unit Gaussians of σ \(=0.5 r\) along that path over \([-T/2,T/2]\), non-negative, then unit-normalized. Implementations SHOULD approximate the average by 64 equal time steps.

`period_s` is typically frozen from mount metadata. `phase_rad` is a nuisance (C8 reports large variance).

v1 default catalog: this term exists but `frozen=true` and is **inactive** (`enabled=false`).

### C3.6.6 Kernel `enabled` flag

Each kernel term has `enabled: bool`. Disabled kernels are omitted from the pipeline (not frozen-at-zero: omitted). Frozen-but-enabled kernels participate at their init/prior mean.

## C3.7 v1 default catalog (complete)

**Normative instance:** `schemas/psf_field_v1_default.catalog.json`. If this section’s table and that file disagree, the JSON file wins.

`catalog_id`: `psf_field_v1_default`. Phase terms all `scope: per_session` for field maps; stage 1 still fits **local** copies (C4).

| term_id | n,m | field `terms` | frozen | init | prior |
|---|---|---|---|---|---|
| `zernike_0_0` | 0,0 | `[[0,0]]` | **true** | zero | none |
| `zernike_1_1` | 1,1 | `[[0,0]]` | **true**, `enabled=false` | zero | none |
| `zernike_1_m1` | 1,-1 | `[[0,0]]` | **true**, `enabled=false` | zero | none |
| `zernike_2_0` | 2,0 | `[[0,0],[1,0],[0,1],[2,0],[1,1],[0,2]]` | false | defocus_moment | none |
| `zernike_2_2` | 2,2 | `[[0,0],[1,0],[0,1]]` | false | zero | none |
| `zernike_2_m2` | 2,-2 | `[[0,0],[1,0],[0,1]]` | false | zero | none |
| `zernike_3_1` | 3,1 | `[[0,0],[1,0],[0,1]]` | false | zero | none |
| `zernike_3_m1` | 3,-1 | `[[0,0],[1,0],[0,1]]` | false | zero | none |
| `zernike_3_3` | 3,3 | `[[0,0]]` | false | zero | none |
| `zernike_3_m3` | 3,-3 | `[[0,0]]` | false | zero | none |
| `zernike_4_0` | 4,0 | `[[0,0]]` | false | zero | none |
| `zernike_4_2` | 4,2 | `[[0,0]]` | false | zero | none |
| `zernike_4_m2` | 4,-2 | `[[0,0]]` | false | zero | none |
| `flux` | — | n/a | false | flux_sum | none |
| `sky` | — | n/a | **true** | zero | none |
| `moffat_seeing` | kernel moffat_iso | uniform | false | moffat_fwhm | gaussian μ=α₀, σ=0.5 α₀ on α; beta frozen 2.5 |
| `gaussian_jitter` | kernel gaussian_iso | uniform | false | zero (σ=0.1 px) | gaussian μ=0.1, σ=0.3 px |
| `charge_diffusion` | kernel gaussian_iso | uniform | false | σ=0.3 px | gaussian μ=0.3, σ=0.1 px |
| `linear_drift` | kernel | uniform | **true**, enabled=false | — | — |
| `field_rotation` | kernel | field law C3.6.4 | **true**, enabled=false | — | — |
| `periodic_error` | kernel | uniform | **true**, enabled=false | — | — |

`flux` and `sky` are not Zernikes; they are first-class `kind: "photometric"` in the schema (third kind). See C3.2 extension:

**C3.2.3 `kind = "photometric"`:** `term_id` in `{flux, sky}` only. `scope: per_star`. No field basis.

**Nuisance `enabled=false` vs `frozen=true`:** `enabled=false` removes the parameter from \(\theta\). `frozen=true` keeps it in \(\theta\) at its init value.

### C3.7.1 Bundles (mechanism only)

```
Bundle = {
  "bundle_id": "collimation",
  "name": "Collimation (unpopulated)",
  "term_ids": ["zernike_3_1", "zernike_3_m1", "zernike_2_2", "zernike_2_m2", "zernike_2_0"],
  "matrix": null
}
```

`matrix` is `null` in v1. When non-null, it is shape `(n_mech, n_field_coeffs)` mapping mechanical DOFs to the concatenated \(c_{ij}\) of the listed terms (C6.5). v1 SHALL accept `null` only. A non-null matrix SHALL be rejected in v1 (reserved).

Named views (not extra parameters):

- `"coma_linear"`: the linear \(c_{ij}\) of `zernike_3_1` and `zernike_3_m1`.
- `"sensor_tilt"`: the \(c_{10}, c_{01}\) of `zernike_2_0`.
- `"field_curvature"`: \(c_{20}, c_{11}, c_{02}\) of `zernike_2_0`.

These views are reporting aliases in C3.8, not extra LM parameters.

## C3.8 Reporting format

```
ReportSpec = {
  "polar_pair": string | null,   # term_id of the sine partner, cosine-only
  "scale": f64,                  # multiply a before printing; default 1
  "unit_label": string            # "waves RMS"
}
```

For a cosine/sine pair, report amplitude \(\sqrt{a_c^{2}+a_s^{2}}\) and angle \(\mathrm{atan2}(a_s, a_c)\) in degrees, in addition to the raw pair. Angle uses NOR.10’s \(\theta\).

## C3.9 Adding an aberration

A conformant implementation SHALL require **only** a new catalog row (and, for a new kernel shape, a new `kernel.id` in C3.6 — that is a contract revision, not a silent plugin). A new \((n,m)\) phase term SHALL NOT require a new Rust match arm beyond generic evaluation.
