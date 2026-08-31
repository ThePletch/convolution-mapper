# C7 — v1 PSF evaluator

## C7.1 Function

`evaluate_psf(stage2, field_xy_mm, grid) -> PsfEval`

`field_xy_mm` is \((x,y)\) in millimetres (NOR.8). It NEED NOT coincide with a measured star.

## C7.2 `EvalGrid`

| Field | Default | Rule |
|---|---|---|
| `stamp_size` | 31 | same allowed set as C1.2.1 |
| `centroid_xy_px` | `((S-1)/2, (S-1)/2)` | stamp-local; default = centered |
| `oversample` | 1 | integer in `{1, 2, 4}`. `oversample=k` means return a `kS × kS` array with pixel scale `pixel_scale/k`, using the same C9 pipeline with `pixel_scale` replaced by `pixel_scale/k` and `S` replaced by `k S` (must remain odd: if `kS` even, **forbidden** — so `S=31`, `k=2` is **rejected**). Allowed pairs: `(S,k)` such that `k*S` is odd, i.e. `k` odd. **Therefore frozen allowed oversample is `{1}` only** unless `S` is chosen so `kS` is odd (e.g. `k=3`, `S=21`). v1 default: `oversample=1`. |

To request a 3× oversampled 21-px stamp, the caller sets `stamp_size=63` and `pixel_scale_arcsec` divided by 3 in a **copy** of `ImageMeta` passed as `grid.pixel_scale_arcsec_override` (optional). If set, C9 uses it. Default: use session `ImageMeta`.

This removes an even-size oversampling footgun.

## C7.3 Local parameters from field maps

For each phase term, \(a = \sum c_{ij} u^i v^j\) with \((u,v)\) from NOR.9 using the **session** `R_field` (detector half-diagonal), even if `(x,y)` is outside the detector. No clamping of \((u,v)\).

For uniform kernels, use `kernel_globals`. For `field_rotation`, evaluate C3.6.4 at `(x,y)`.

`flux=1`, `sky=0`.

`known_defocus_waves` from the **session’s primary exposure** if a single exposure; if multiple, `EvalGrid.known_defocus_waves` (default 0) is the only diversity applied at eval time.

## C7.4 Output `PsfEval`

| Field | Type | Meaning |
|---|---|---|
| `psf` | `f64[S,S]` | model ADU with \(F=1\), \(b=0\) (so approximately unit-sum; C9.10) |
| `zernike_vector` | dict `term_id` → f64 | local \(a_k\) actually used, including frozen zeros |
| `kernel_vector` | dict | local kernel params |
| `field_xy_mm` | `[f64,2]` | echo |
| `u_v` | `[f64,2]` | normalized field |
| `image_meta_digest` | string | SHA-256 of canonical JSON of `ImageMeta` + `PupilSpec` + `catalog_id`, hex |
| `catalog_id` | string | |
| `stage2_schema_version` | string | |

Provenance: a consumer SHALL be able to reproduce `psf` by calling `forward_psf` with `zernike_vector` ∪ photometric ∪ kernels. C10.7 checks \(\max| \mathrm{psf} - \mathrm{forward\_psf}(\ldots) | < 10^{-12}\).

## C7.5 Out-of-field positions

Allowed. The polynomial will extrapolate. `PsfEval.extrapolated = true` iff \(|u|>1\) or \(|v|>1\) **or** the point lies outside the convex hull of stage-2 stars. Both conditions are recorded separately: `outside_unit_square`, `outside_hull`.

## C7.6 What the evaluator SHALL NOT do

- Average neighboring stars’ stamps.
- Return a Zernike vector without running C9 (the PSF is not optional).
- Change catalog freeze flags.
