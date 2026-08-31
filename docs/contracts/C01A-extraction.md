# C1A — Extraction front-end (v1)

Python-only. Output is **exactly** C1. No extra columns are required by the modeler. Prior art is the requirement: this module SHALL wrap `photutils.detection.DAOStarFinder` (Stetson 1987 lineage), not a novel detector.

## C1A.1 Algorithm sequence (frozen order)

1. Load FITS image as `float64` 2-D array `data_raw` in ADU, 0-based (NOR.7). Do not apply any additional nonlinearity correction.
2. Read `ImageMeta` required keywords (C1.4), merging a sidecar if present (C1.4.2); reject if any NOR.13 field is missing after merge.
3. Background (C1A.2) → `bkg`, `bkg_rms`.
4. `data_sub = data_raw - bkg`.
5. Detect (C1A.3).
6. Reject on sharpness/roundness (already inside DAOStarFinder) and on saturation, blending, edge (C1A.5–C1A.7).
7. Refine centroid (C1A.8).
8. Cut out stamp + variance + pixel_mask (C1A.9–C1A.10).
9. Select modeling subset and compute coverage (C1A.11).
10. Write C1 serialization.

## C1A.2 Background

SHALL use `photutils.background.Background2D`.

| Parameter | Frozen default | Allowed |
|---|---|---|
| `box_size` | 64 | odd or even integers in [32, 128] |
| `filter_size` | 3 | odd integers in {1, 3, 5} |
| `sigma_clip.sigma` | 3.0 | [2.0, 5.0] |
| `sigma_clip.maxiters` | 5 | [3, 10] |
| `bkg_estimator` | `photutils.background.SExtractorBackground` | that class only |
| `bkgrms_estimator` | `photutils.background.StdBackgroundRMS` | that class only |

`bkg` is the 2-D interpolated background evaluated on the full image. `bkg_rms` is the 2-D RMS map from the same object (`background_rms` property).

If `Background2D` raises, extraction SHALL fail the exposure (not silently fall back to a scalar median).

## C1A.3 Detection

SHALL use `photutils.detection.DAOStarFinder`.

| Parameter | Frozen default | Allowed |
|---|---|---|
| `fwhm` | **required input**, no default | finite, > 1.0 and < `S` |
| `threshold` | `n_sigma * median(bkg_rms)` | see `n_sigma` |
| `n_sigma` | 5.0 | [3.0, 10.0] |
| `sharpness_range` | `(0.2, 1.0)` | inclusive bounds; SHALL use the photutils ≥3.0 tuple API, not deprecated `sharplo`/`sharphi` |
| `roundness_range` | `(-1.0, 1.0)` | same |
| `exclude_border` | `True` | frozen True |
| `min_separation` | `1.0 * fwhm` | [0.5, 3.0] × fwhm |
| `ratio` | 1.0 | frozen (circular kernel) |
| `theta` | 0.0 | frozen |
| `sigma_radius` | 1.5 | frozen |
| `peak_max` | `0.95 * saturation_adu` | frozen formula |
| `n_brightest` | `None` | no cap |

`fwhm` is a user/config value in pixels. It is **not** estimated inside v1. The corpus generator writes the known FWHM into the corpus manifest (C10). For real images the operator supplies it.

DAOStarFinder is called on `data_sub`. Threshold is an absolute ADU value equal to `n_sigma * median(bkg_rms)` computed over finite `bkg_rms` pixels.

## C1A.4 Finder table columns consumed

From the returned table, the following columns SHALL be read (photutils names): `xcentroid`, `ycentroid`, `sharpness`, `roundness1`, `roundness2`, `peak`, `flux`. Other columns MAY be stored under `aux_`.

`xcentroid`/`ycentroid` are 0-based pixel coordinates matching NOR.7 (photutils convention). No FITS 1-based conversion at this step.

## C1A.5 Saturation flag

For each detection, form the C1.2.2 window on **`data_raw`**. If any pixel in that window satisfies `data_raw >= 0.95 * saturation_adu`, set `SATURATED`.

`peak_max` in C1A.3 is an additional finder-level cut; C1A.5 is the stamp-level cut that travels with the record. Both apply.

## C1A.6 Blended / neighbor flag

Let \(d_{ab}\) be the Euclidean distance in pixels between detections \(a\) and \(b\). If there exists \(b \neq a\) with

\[
d_{ab} < 1.5 \times \texttt{fwhm}
\]

then **both** `a` and `b` receive `BLENDED`.

This is a pairing cut, not deblending. v1 SHALL NOT subtract neighbors.

## C1A.7 Edge flag

If the C1.2.2 window is not a subset of `[0, n_col) × [0, n_row)`, set `EDGE`. Do not crop. Do not pad.

## C1A.8 Centroid refinement

After DAOStarFinder, refine each unflagged-or-not-yet-cutout detection with `photutils.centroids.centroid_2dg` on a box of width

\[
B = \max(5,\; 2 \cdot \operatorname{floor}(\texttt{fwhm}/2) + 1)
\]

(`B` odd) centered on the nearest integer pixel to the finder centroid, on `data_sub`. If `centroid_2dg` returns NaN or raises, keep the finder centroid and set `SHAPE`.

The refined `(x, y)` becomes `source_xy_px`.

## C1A.9 Variance (CCD equation)

Let `gain` = `gain_e_per_adu`, `rn` = `read_noise_e`. For each stamp pixel corresponding to full-image `(i,j)`:

\[
\mathrm{var_{ADU}}[j,i] = \frac{\max(\texttt{data\_raw}[j,i],\, 0)}{\mathrm{gain}} + \frac{\mathrm{rn}^{2}}{\mathrm{gain}^{2}}
\]

- Poisson term uses **`data_raw`**, never `data_sub`.
- Negative raw pixels contribute 0 to the Poisson term via the `max(·,0)`.
- Result is ADU², matching C1.3.

## C1A.10 Stamp assembly

- `stamp = data_sub[window]`.
- `variance` from C1A.9 on the same window.
- `pixel_mask`: bit `SATURATED` if that pixel’s `data_raw >= 0.95 * saturation_adu`; bit `INVALID` if `data_raw` or `data_sub` is non-finite.
- `centroid_xy_px` from C1.2.3.
- `field_xy_mm` from NOR.8.
- `flux_sum_adu` = sum of `stamp` over `pixel_mask==0`.

## C1A.11 Star selection for modeling

A detection is a **candidate** if it has none of `{SATURATED, BLENDED, EDGE, SHAPE, UNDERDETERMINED}`.

SNR (frozen):

\[
\mathrm{SNR} = \frac{\sum_{p \in \mathrm{valid}} \mathrm{stamp}_p}{\sqrt{\sum_{p \in \mathrm{valid}} \mathrm{variance}_p}}
\]

Candidates with \(\mathrm{SNR} < 20\) SHALL NOT receive `SELECTED`. **Frozen SNR minimum: 20.** Allowed config interval: [10, 50].

`ExtractionConfig.selection_mode` (frozen allowed set):

| Value | Default | Behavior |
|---|---|---|
| `highest_snr` | **yes** | C1A.11.1 |
| `snr_by_cell` | no | C1A.11.2 |

**Frozen default `max_selected = 400`.** Allowed [50, 2000]. If the number of candidates with SNR ≥ 20 is \(\le\) `max_selected`, keep all of them (no brightness cap) under either mode.

### C1A.11.1 `highest_snr` (default)

If more than `max_selected` candidates survive the SNR floor, keep the `max_selected` with highest SNR. Highest-SNR-only MAY empty detector corners; that is why C1A.11.2 exists. C10 uses this mode (and \(n_{\mathrm{truth}} <\) `max_selected`), so corpus gates are unchanged.

### C1A.11.2 `snr_by_cell`

Same candidates and SNR floor. Partition the detector with the same 3×3 field-mm bins as C1A.12. Let \(q = \lfloor \texttt{max\_selected}/9 \rfloor\). In each cell, take the \(q\) highest-SNR candidates (or all if fewer). Fill leftover slots from remaining candidates by global SNR. If `max_selected` is not a multiple of 9, the remainder after the nine quotas also goes to global SNR.

v1 SHALL NOT use continuous min-distance thinning or D-optimal design. `snr_by_cell` is the only coverage-aware ranker.

### C1A.11.3 Hold-out (after selection)

`ExtractionConfig.holdout_fraction` (default **0**, allowed [0, 0.5]) MAY assign `USER_EXCLUDE` at random among stars that just received `SELECTED`, using `holdout_seed` (`uint64`, required if `holdout_fraction > 0`). This runs **after** C1A.11.1/2. Default 0 preserves the corpus. Stage-2 SHALL NOT consume `USER_EXCLUDE` stars (C6.1).

v1 SHALL NOT spatially thin beyond the blend cut except via `snr_by_cell`. Field coverage is always reported (C1A.12).

## C1A.12 Field-coverage report

Written as `coverage.json` beside the C1 files. Required fields:

| Field | Definition |
|---|---|
| `n_detected` | rows returned by DAOStarFinder before our extra flags |
| `n_candidate` | C1A.11 candidates |
| `n_selected` | with `SELECTED` |
| `frac_selected_of_detected` | `n_selected / max(n_detected, 1)` |
| `grid_3x3` | length-9 array: counts of selected stars in a 3×3 partition of the detector in field-mm, bins equally spaced in \([x_{\min}, x_{\max}]\) and \([y_{\min}, y_{\max}]\) of the **detector**, not of the star convex hull |
| `empty_cells` | number of the 9 cells with count 0 |
| `convex_hull_area_mm2` | area of the convex hull of selected `field_xy_mm`; 0 if \(n_{\mathrm{selected}} < 3\) |
| `detector_area_mm2` | \(N_{\mathrm{col}} N_{\mathrm{row}} p_{\mathrm{mm}}^{2}\) |
| `hull_fill` | `convex_hull_area_mm2 / detector_area_mm2` |
| `design_cond_plane` | \(\kappa_2\) of the design matrix whose rows are \((1, u, v)\) for selected stars (2-norm condition number). If \(n_{\mathrm{selected}} < 3\), store `+inf`. |
| `design_cond_quad` | same for monomials \(i+j \le 2\) (6 columns). If \(n_{\mathrm{selected}} < 6\), `+inf`. |

Stage-2 SHALL still run if some cells are empty; C6 reports condition numbers. C1A.12 is informational and a corpus gate (C10.6).

**Frozen corpus gate:** a synthetic exposure fails extraction-quality (not a fitter fail) if `empty_cells > 3` or `design_cond_plane > 1e6`.

## C1A.13 Typed facade

`photutils` and some `astropy` APIs are untyped. This module SHALL expose only the following typed functions (no `Any` in their signatures):

```python
def extract_exposure(
    image_path: Path,
    config: ExtractionConfig,
) -> tuple[list[StarRecord], ImageMeta, CoverageReport]: ...

def extract_to_files(
    image_path: Path,
    config: ExtractionConfig,
    out_dir: Path,
) -> None: ...
```

`ExtractionConfig` is a pydantic `BaseModel` with the fields in C1A.2, C1A.3, C1A.11 (`selection_mode`, `max_selected`, SNR floor, `holdout_fraction`, `holdout_seed`). Third-party objects stay inside the function body. Schema: `schemas/extraction_config.schema.json`.

## C1A.14 Test-corpus ingestion

The same `extract_exposure` SHALL accept synthetic FITS produced by C10. No corpus-specific code path other than reading the same headers. Ground-truth comparison happens in C10 after the core returns, not inside extraction.
