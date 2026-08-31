# C1 — Stage-1 input (measurement / modeling boundary)

This is the only input the modeler is allowed to know about. The modeler SHALL NOT branch on whether a record came from C1A, SExtractor, photutils, or a hand-built corpus file.

## C1.1 Object: `StarRecord`

A stage-1 input file is an unordered set of `StarRecord` objects that share one `ImageMeta` (C1.4) per `exposure_id`, plus one `PupilSpec` (C1.5) per session.

### Required fields

| Field | Type | Unit | Constraint |
|---|---|---|---|
| `schema_version` | string | — | `"1.0.0"` |
| `star_id` | string | — | NOR.5 |
| `exposure_id` | string | — | NOR.5 |
| `session_id` | string | — | NOR.5 |
| `field_xy_mm` | `[f64; 2]` | mm | finite |
| `source_xy_px` | `[f64; 2]` | px | 0-based pixel coords of the extracted centroid on the full image |
| `stamp` | `f64[S, S]` | ADU | C1.2 |
| `variance` | `f64[S, S]` | ADU² | C1.3 |
| `centroid_xy_px` | `[f64; 2]` | px | stamp-local, 0-based; C1.2.3 |
| `pixel_mask` | `u8[S, S]` | — | C1.6 |
| `flags` | `FlagSet` | — | C1.7 |
| `flux_sum_adu` | `f64` | ADU | sum of `stamp` over pixels with `pixel_mask==0`; finite, used only as an init hint |

No other fields are part of the modeler input. Extra serialized columns MAY exist for human debugging; the core SHALL ignore unknown optional columns with names prefixed `aux_`.

## C1.2 Postage stamp

### C1.2.1 Size

- `S` SHALL be an odd integer.
- **Frozen default:** `S = 31`.
- **Allowed:** odd integers in `{15, 17, …, 63}`.
- Even `S` SHALL be rejected at ingest.
- `S < 15` or `S > 63` SHALL be rejected at ingest.
- All stars in one `exposure_id` SHALL share the same `S`. Mixed sizes SHALL be rejected.

### C1.2.2 Content

- `stamp[j, i]` is background-subtracted ADU at stamp pixel `(x,y) = (i, j)`.
- The stamp SHALL already have a constant (or spatially interpolated) background removed. The modeler SHALL NOT re-estimate a 2-D background. A single scalar residual-sky parameter is a scoped optional fit parameter (C4), frozen to 0 by default, not a substitute for C1A background subtraction.
- Cutout geometry (frozen): the stamp is centered on the **nearest integer pixel** to `source_xy_px`. Let

\[
i_0 = \operatorname{round}(x_{\mathrm{src}}), \qquad j_0 = \operatorname{round}(y_{\mathrm{src}})
\]

using half-away-from-zero rounding (`round` as IEEE-754 `rint` / Python `round` ties-to-even is **forbidden** here because it is platform-ambiguous). **Frozen rounding:** `i_0 = floor(x_src + 0.5)` for \(x \ge 0\) (pixel coords are non-negative). Then the stamp covers full-image pixels

\[
i \in [i_0 - \lfloor S/2\rfloor,\; i_0 + \lfloor S/2\rfloor], \quad
j \in [j_0 - \lfloor S/2\rfloor,\; j_0 + \lfloor S/2\rfloor].
\]

If that window is not fully inside the image, the star SHALL be flagged `EDGE` and SHALL NOT be passed to the modeler as a fit-eligible record (it MAY still be serialized with the flag for diagnostics).

### C1.2.3 Stamp-local centroid

\[
c_x = x_{\mathrm{src}} - (i_0 - \lfloor S/2\rfloor), \qquad
c_y = y_{\mathrm{src}} - (j_0 - \lfloor S/2\rfloor)
\]

The geometric center of the stamp is \(c_\star = (S-1)/2\). The Fourier-shift applied by C9 is \(\Delta = (c_x - c_\star,\, c_y - c_\star)\).

`centroid_xy_px` SHALL equal \((c_x, c_y)\). Ingest SHALL reject if \(\max(|c_x-c_\star|, |c_y-c_\star|) > 0.6\) (the integer-pixel centering plus a sane subpixel residual). A centroid more than 0.6 px from the stamp center indicates a cutout/centroid mismatch.

## C1.3 Variance map

- `variance[j, i]` is the estimated variance of `stamp[j, i]` in ADU².
- Every finite stamp pixel with `pixel_mask==0` SHALL have `variance > 0` and finite. Zero or negative variance on a valid pixel SHALL be rejected.
- Variance SHALL be computed from the **pre-subtraction** image via C1A.9 (CCD equation). Using the background-subtracted stamp as the Poisson source term is forbidden (it can be negative and is biased).
- The modeler treats `variance` as given. It SHALL NOT rescale it globally except via the optional frozen-off `error_scale` parameter in C4 (default 1, not fitted in v1).

## C1.4 Object: `ImageMeta` (per exposure)

| Field | Type | Unit | Required |
|---|---|---|---|
| `exposure_id` | string | — | yes |
| `session_id` | string | — | yes |
| `n_row`, `n_col` | int64 | px | ≥ `S` |
| `wavelength_m` | f64 | m | > 0, finite; v1: monochromatic effective wavelength |
| `pupil_diameter_m` | f64 | m | > 0 |
| `focal_length_m` | f64 | m | > 0 |
| `pixel_scale_arcsec` | f64 | arcsec/px | > 0 |
| `optical_axis_pixel` | `[f64; 2]` | px | finite, 0-based |
| `gain_e_per_adu` | f64 | e⁻/ADU | > 0 |
| `read_noise_e` | f64 | e⁻ | ≥ 0 |
| `saturation_adu` | f64 | ADU | > 0 |
| `exptime_s` | f64 | s | > 0 |
| `known_defocus_waves` | f64 | waves | finite; **default 0**. Frozen additive to \(a_{2}^{0}\) (C4.5) |
| `pixel_size_m` | f64 \| null | m | optional consistency check |

**Plate-scale consistency (C1.4.1).** If `pixel_size_m` is present:

\[
s_{\mathrm{pred}} = \frac{\texttt{pixel\_size\_m}}{\texttt{focal\_length\_m}} \quad\text{(rad/px)}, \qquad
s_{\mathrm{hdr}} = \texttt{pixel\_scale\_arcsec} \cdot \frac{\pi}{180\times 3600}
\]

Let \(\delta = |s_{\mathrm{pred}} - s_{\mathrm{hdr}}| / s_{\mathrm{hdr}}\). If \(\delta > 0.05\), ingest SHALL reject. If \(0.01 < \delta \le 0.05\), ingest SHALL accept and SHALL set `ImageMeta.plate_scale_warning = true`. If `pixel_size_m` is absent, skip this check.

**C1.4.2 ImageMeta sidecar.** A YAML or JSON sidecar MAY be merged over FITS keywords for NOR.13 fields (and the rest of C1.4). After merge, any still-missing required field SHALL be rejected. Implementations SHALL NOT default `optical_axis_pixel` to the array center. A unit test SHALL show: incomplete FITS headers fail; the same file plus a complete sidecar succeeds.

**Bandpass:** a FITS `FILTER` or `BANDPASS` string MAY be stored as `aux_filter`. The core SHALL use only `wavelength_m`.

## C1.5 Object: `PupilSpec`

| Field | Type | Constraint |
|---|---|---|
| `mask` | `f64[N_p, N_p]` | C9.4; values in `{0, 1}` for v1 |
| `n_pupil` | int64 | frozen default 256; allowed `{128, 256, 512}` |
| `n_fft` | int64 | frozen default `4 * n_pupil`; allowed `n_fft / n_pupil ∈ {2, 4, 8}` |
| `amplitude` | `f64[N_p, N_p]` \| omitted | if omitted, treat as 1 where `mask==1`, else 0 |

v1 shipped mask: circular unobstructed (C9.4). A non-circular mask is valid input (the engine evaluates Zernikes on the supplied mask) but is not a v1 *catalog* change. Field-dependent amplitude (vignetting) is a reserved slot: if `amplitude` varies, C9 uses it; v1 extraction SHALL write the default uniform-on-mask amplitude.

## C1.6 Pixel mask codes

`pixel_mask` is a bitfield stored as `u8`. Bits:

| Bit | Name | Meaning | Weight in C5 |
|---|---|---|---|
| 0 | `INVALID` | ignore pixel | 0 |
| 1 | `SATURATED` | at or above saturation | 0 |
| 2 | `COSMIC` | optional; v1 unused | 0 |
| 3 | `NEIGHBOR` | neighbor flux suspected in this pixel | 0 |

A pixel is **valid for the residual** iff `pixel_mask == 0` and `variance` is finite and `> 0` and `stamp` is finite.

If \(m < n_{\mathrm{free}}\), the star SHALL be flagged `UNDERDETERMINED` and excluded from fitting. If \(n_{\mathrm{free}} \le m < n_{\mathrm{free}} + 8\), the star SHOULD be flagged `UNDERDETERMINED` and excluded (robustness margin; C5.1).

## C1.7 Flags (`FlagSet`)

Flags are a set of strings from this closed vocabulary. Unknown flag names SHALL be rejected.

| Flag | Set by | Fit-eligible? |
|---|---|---|
| `SATURATED` | any stamp pixel ≥ `0.95 * saturation_adu` in the **pre-subtraction** image | no |
| `BLENDED` | neighbor within `1.5 * fwhm_px` (C1A) | no |
| `EDGE` | stamp window not fully inside image | no |
| `SHAPE` | sharpness/roundness outside C1A cuts | no |
| `UNDERDETERMINED` | too few valid pixels | no |
| `USER_EXCLUDE` | user/catalog | no |
| `SELECTED` | C1A.11 selection | yes (must also have no exclusion flags) |

Fit-eligible ⇔ `SELECTED` is present AND none of `{SATURATED, BLENDED, EDGE, SHAPE, UNDERDETERMINED, USER_EXCLUDE}` are present.

The modeler SHALL NOT apply additional rejection heuristics. Exclusion is the harness’s job (C1A or an external producer).

## C1.8 Config object: `Stage1InputConfig`

Serialized beside the stars. Frozen defaults:

| Field | Default | Allowed |
|---|---|---|
| `stamp_size` | 31 | odd in [15, 63] |
| `schema_version` | `"1.0.0"` | — |

## C1.9 Serialization (file-level contract)

The same logical tables SHALL be expressible as:

1. **FITS** (primary: `ImageMeta` in header + pupil mask image HDU; binary table HDU `STARS`).
2. **Parquet** (two files: `stars.parquet`, `exposures.parquet`, plus `pupil.fits` or `pupil.npy`).

### C1.9.1 FITS `STARS` columns

| TTYPE | TFORM | Content |
|---|---|---|
| `STAR_ID` | 128A | |
| `EXPOSURE_ID` | 128A | |
| `SESSION_ID` | 128A | |
| `FIELD_X_MM` | D | |
| `FIELD_Y_MM` | D | |
| `SRC_X_PX` | D | |
| `SRC_Y_PX` | D | |
| `CEN_X_PX` | D | stamp-local |
| `CEN_Y_PX` | D | |
| `FLUX_SUM` | D | |
| `FLAGS` | 256A | comma-separated, sorted alphabetically |
| `STAMP` | `S*S` D, or variable-length | row-major `S*S` |
| `VARIANCE` | `S*S` D | |
| `PIXMASK` | `S*S` B | |

Primary header SHALL contain `SCHEMAV = '1.0.0'` and all C1.4 fields with the names:

`EXPTIME`, `GAIN`, `RDNOISE`, `SATURATE`, `LAMBDA` (metres), `PUPILD` (metres), `FOCALLEN` (metres), `PIXSCAL` (arcsec/px), `OAX`, `OAY`, `KDEFOCUS` (waves), `STMPSIZ`.

**Golden fixture:** `tests/fixtures/c1_roundtrip.fits` (added with the implementation) is the round-trip oracle for these `TTYPE` / header names. Python (astropy) and the Rust CLI SHALL both read it.

### C1.9.2 JSON Schema

See `schemas/star_record.schema.json`, `schemas/image_meta.schema.json`, `schemas/pupil_spec.schema.json`. Remaining serialized objects: `schemas/stage1_result.schema.json`, `schemas/stage2_result.schema.json`, `schemas/psf_eval.schema.json`, `schemas/fd_report.schema.json`, `schemas/score_report.schema.json`, `schemas/coverage_report.schema.json`, `schemas/extraction_config.schema.json`. Pydantic models SHALL be generated from or tested 1:1 against those schemas. Rust `serde` structs SHALL deserialize the same JSON.

## C1.10 Invariants the modeler may assume

After ingest:

1. All fit-eligible stars of an exposure share `S`, `ImageMeta`, and `PupilSpec`.
2. `field_xy_mm` was computed by NOR.8 from `source_xy_px` (not independently invented). Ingest SHALL recompute and reject if \(\| \mathbf{x}_{\mathrm{stored}} - \mathbf{x}_{\mathrm{recomputed}} \|_\infty > 10^{-6}\,\mathrm{mm}\).
3. `stamp` is background-subtracted; `variance` is not.
4. The modeler does not receive the full image.

## C1.11 What the modeler SHALL NOT do

- SHALL NOT detect stars, estimate FWHM, or run DAOStarFinder.
- SHALL NOT fit a 2-D background image.
- SHALL NOT treat flags as scores; they are booleans.
- SHALL NOT interpret WCS beyond the already-reduced `field_xy_mm` (no `astropy.wcs` inside Rust).
