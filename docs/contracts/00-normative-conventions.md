# NOR — Normative conventions

All other contracts inherit this file. Re-stating a convention differently elsewhere is a defect.

## NOR.1 Language of requirements

Requirements use RFC 2119 **SHALL / SHALL NOT / SHOULD / MAY** as defined in `README.md`. Every configurable quantity has a **frozen default**, an **allowed set or interval**, and a **rejection** (what happens if the value is outside the allowed set). SHOULD clauses are methods or margins; they do not relax a nearby SHALL.

## NOR.2 Schema versioning

- The data-contract version is the string `schema_version` with SemVer `MAJOR.MINOR.PATCH`.
- **Frozen for v1:** `"1.0.0"`.
- A change to any field name, unit, array layout, missingness rule, or formula that would make an old file silently mis-parse SHALL bump MAJOR.
- The shipped default catalog JSON is updated in place (tilt enabled per-star, zero-mean phase priors, quadratic astigmatism, unfrozen sky and Moffat β, anisotropic Gaussian kernel). `schema_version` remains `"1.0.0"` so existing files remain valid. A strictly additive optional field with a stated default MAY bump MINOR.
- Documentation-only clarification MAY bump PATCH.
- Every serialized object (JSON, FITS header, Parquet metadata, PyO3 capsule) SHALL include `schema_version`. Readers SHALL reject unknown MAJOR.

## NOR.3 Scalar types and IEEE-754

- All floating-point computation in the core SHALL use IEEE-754 binary64 (`f64` / `numpy.float64`).
- `float32`, `float16`, and mixed-precision accumulation SHALL NOT be used in the forward model, Jacobian, LM, or stage-2 solve.
- A NaN or Inf produced in the forward model or Jacobian SHALL cause the LM residual/Jacobian call to return “invalid” (C5.9), not a silently clipped value.
- Integers that count pixels or indices SHALL be `int64` in serialized form and Python `int` at the pydantic boundary (pydantic will coerce JSON numbers; the stored dtype of tables is `int64`).

## NOR.4 Array layout

- All 2-D arrays that cross a language or file boundary SHALL be **C-contiguous**, row-major, with index `[row, col]` = `[y, x]`.
- 1-D residual and parameter vectors SHALL be C-contiguous.
- The physical meaning of axis 0 is **row / +y in pixel space**; axis 1 is **column / +x in pixel space**.

## NOR.5 Identifiers

| Kind | Pattern | Notes |
|---|---|---|
| `star_id` | `^[A-Za-z0-9._:-]{1,128}$` | Unique within a session file |
| `exposure_id` | same | Groups stars that share one kernel scope |
| `session_id` | same | Groups exposures that share aberration field maps |
| `term_id` | `^[a-z][a-z0-9_]*$` | Catalog key, e.g. `zernike_3_1` |
| `bundle_id` | same | Catalog key, e.g. `collimation` |

Duplicate IDs within a file SHALL be rejected at ingest.

## NOR.6 Units (SI internally)

All core computation SHALL use the following internal units. Conversion from FITS/header units happens once, at the Python boundary, and is not repeated in Rust.

| Quantity | Internal unit | Symbol |
|---|---|---|
| Length (pupil diameter, focal length) | metre | m |
| Wavelength | metre | m |
| Wavefront / Zernike coefficient | **waves RMS** (dimensionless cycles) | waves |
| Pupil phase | radian | rad |
| Angle on sky / pixel scale | radian | rad |
| Field coordinates | millimetre | mm |
| Detector pixel coordinates | pixel (see NOR.7) | px |
| Time | second | s |
| Rotation rate | radian / second | rad/s |
| Flux | ADU (image native) | ADU |
| Variance | ADU² | ADU² |
| Gain | electron / ADU | e⁻/ADU |
| Read noise | electron RMS | e⁻ |

**Phase conversion (frozen):**

\[
\Phi(\rho,\theta) = 2\pi \sum_k a_k\, Z_k(\rho,\theta)
\]

where \(a_k\) is in waves RMS and \(Z_k\) is discrete-RMS-normalized (C2), so \(a_k = 1\) means 1 wave RMS of that mode over the mask.

**Sky-angle conversion (frozen):**

\[
\alpha_{\mathrm{rad}} = \alpha_{\mathrm{arcsec}} \cdot \frac{\pi}{180 \times 3600}
\]

**FITS wavelength:** if a header supplies Ångström, convert \(\lambda_m = \lambda_Å \times 10^{-10}\). If it supplies nanometres, \(\lambda_m = \lambda_{nm} \times 10^{-9}\). The pydantic model stores only metres. v1 is monochromatic: `wavelength_m` SHALL be the flux-weighted effective wavelength of the bandpass. A multi-sample incoherent sum over a filter is a reserved later slot, not a v1 requirement.

## NOR.7 Pixel coordinate system

- Pixel coordinates are **0-based**. The center of the first pixel is \((x,y) = (0, 0)\).
- The sample at array index `[j, i]` (row `j`, column `i`) is the pixel whose center is \((x,y) = (i, j)\).
- This matches NumPy indexing and **differs** from FITS 1-based convention. FITS I/O SHALL convert: \(x_{0} = \mathrm{FITS}_x - 1\), \(y_{0} = \mathrm{FITS}_y - 1\).
- Sub-pixel values are allowed and required for centroids. They refer to the same origin (pixel centers).

## NOR.8 Field coordinate system

- Field coordinates \((x, y)\) are millimetres in the focal plane.
- Origin: the **optical axis** intersection with the focal plane, **not** the array corner and **not** automatically the array center.
- \(+x\) increases with pixel \(+x\) (column). \(+y\) increases with pixel \(+y\) (row).
- Mapping (frozen):

\[
x_{\mathrm{mm}} = (x_{\mathrm{px}} - x_{\mathrm{oa}})\, p_{\mathrm{mm}}, \qquad
y_{\mathrm{mm}} = (y_{\mathrm{px}} - y_{\mathrm{oa}})\, p_{\mathrm{mm}}
\]

where \((x_{\mathrm{oa}}, y_{\mathrm{oa}})\) is `optical_axis_pixel` and \(p_{\mathrm{mm}}\) is millimetres per pixel.

- Millimetres per pixel (frozen):

\[
p_{\mathrm{mm}} = \frac{f_{\mathrm{m}} \cdot \mathrm{pixel\_scale}_{\mathrm{rad}}}{1\,\mathrm{rad}} \times 10^{3}
\]

which is the small-angle map \(s = f \theta\). Equivalently \(p_{\mathrm{mm}} = f_{\mathrm{m}} \cdot \mathrm{pixel\_scale}_{\mathrm{rad}} \times 10^{3}\).

- `optical_axis_pixel`, `focal_length_m`, `pixel_scale_arcsec` are **required** image metadata (C1.4). There is no default that silently assumes the axis is at the array center.

## NOR.9 Normalized field coordinates (for polynomial bases)

Let \(R_{\mathrm{field}}\) be the detector half-diagonal in millimetres:

\[
R_{\mathrm{field}} = \frac{1}{2}\sqrt{(N_{\mathrm{col}} p_{\mathrm{mm}})^2 + (N_{\mathrm{row}} p_{\mathrm{mm}})^2}
\]

Normalized coordinates (dimensionless):

\[
u = \frac{x_{\mathrm{mm}}}{R_{\mathrm{field}}}, \qquad v = \frac{y_{\mathrm{mm}}}{R_{\mathrm{field}}}
\]

Polynomial field bases in C3 are monomials in \((u,v)\). This keeps field-coefficient magnitudes O(waves), not O(waves·mm⁻ⁿ).

If \(R_{\mathrm{field}} = 0\) (zero-size image), ingest SHALL reject the image.

## NOR.10 Pupil polar coordinates

On the pupil grid defined in C9:

\[
\rho = \sqrt{\xi^2 + \eta^2}, \qquad \theta = \operatorname{atan2}(\eta, \xi)
\]

where \(\xi$ is the pupil \(+x\) coordinate (same handedness as field \(+x\)) and \(\eta\) is pupil \(+y\). \(\theta = 0\) on \(+\xi\), increasing toward \(+\eta\) (right-handed, standard `atan2`).

The unit disk is \(\rho \le 1\). The v1 circular mask is exactly that set (C9.4).

## NOR.11 Time origin

- Exposure midpoint is the time origin for kernel trajectories unless `date_obs_mjd` is supplied, in which case kernel phase parameters are relative to that MJD converted to seconds in the Unix epoch only for serialization; internally, \(t=0\) is the exposure midpoint and \(t \in [-T/2, T/2]\) with \(T =\) `exptime_s`.

## NOR.12 Randomness

The core forward model, fit, and diagnostics SHALL be deterministic given their inputs. The synthetic corpus generator (C10) SHALL take an explicit `seed: uint64`. There is no global RNG.

## NOR.13 Forbidden silent defaults

The following SHALL NOT be defaulted by implementations; they are required fields:

- `wavelength_m`, `pupil_diameter_m`, `focal_length_m`, `pixel_scale_arcsec`
- `gain_e_per_adu`, `read_noise_e`, `saturation_adu`
- `optical_axis_pixel`
- `stamp_size` (odd integer; default is allowed only as the catalog/config default in C1.8, not as a hidden constant in engine code)

## NOR.14 Requirement ID grammar

Requirement IDs are `C{contract}.{n}` or `NOR.{n}` or `C1A.{n}` etc. Sub-clauses are `C{contract}.{n}.{letter}`. Citations SHALL use those IDs.
