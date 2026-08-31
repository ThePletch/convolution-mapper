# Conceptual consistency lemmas

These statements can be checked against C1–C11 without writing code. If any lemma fails, the contracts are internally inconsistent and MUST be revised before implementation.

## L1. Composition

Aberration modules expose \(Z_n^m\) or \(K(\alpha)\), never a PSF residual. C9.1 is the only path from coefficients to pixels. **Check:** C3 has no “psf_add” kind.

## L2. Units of the Jacobian

\(\Phi = 2\pi \sum a_k Z_k\) with \(a_k\) in waves (NOR.6, C2.6). C9.8 includes \(2\pi\) in \(Q = i 2\pi Z P\). C5.2 weights by \(F/\sigma\). **Check:** no document uses “phase coefficients in radians” as the stored unit.

## L3. Centroid vs tilt vs coma

C1 supplies centroid; C5 SHALL NOT fit it. C3.7 disables \(Z_1^{\pm1}\). C2.8.5 requires Zernike (balanced) coma. **Check:** a Seidel \(\rho^3\cos\theta\) term is not in the catalog.

## L4. Sensor tilt is not a second defocus mode

Sensor tilt is the linear monomials of `zernike_2_0` (C3.3, C3.7.1). There is no extra \((n,m)\). Stage-2 ill-conditioning of the defocus map is how tilt vs curvature is diagnosed (C6.3), not a rivalry of modules.

## L5. Defocus sign

C2.8.3 evenness + C3.5.1 positive init + C5.3 `defocus_sign_ambiguous` + C4.5 known offset + C10.3.3 absolute-value scoring. **Check:** no requirement says LM will recover the algebraic sign from a single in-focus exposure.

## L6. Two-stage information barrier

C6.7 forbids rereading pixels. C7 uses maps + C9 only. **Check:** `Stage2Result` has no stamp arrays.

## L7. Catalog-as-data

New \((n,m)\) = new row. New kernel *shape* = contract revision (C3.6 closed set). **Check:** engine API in C2 takes `(n,m)` only.

## L8. Blur stacking

Jitter, Moffat, diffusion are three `gaussian_iso`/`moffat_iso` terms with priors (C3.7) and mandatory session correlation reporting (C8.2). **Check:** they are not merged into one “effective FWHM” parameter in v1 (they remain distinct, possibly degenerate).

## L9. Discrete vs analytic Zernikes

Runtime uses discrete RMS on the mask (C2.4). Analytic \(N_n^m\) is only the pre-normalization. Orthogonality is reported (C2.5), not enforced by Gram–Schmidt.

## L10. LM crate vs residual definition

C11.3: reported `chi2 = 2 * crate_objective`. C5.1 \(r=(d-\mathrm{model})/\sigma\). **Check:** nobody stores \(\tfrac12\chi^2\) as `chi2`.

## L11. Extraction cannot leak into the core

C1.11: no DAOStarFinder in Rust. C1A outputs C1 only. Corpus still runs C1A (C10) so centroid error is a measured bias, not an ignored one.

## L12. Amplitude slot without vignetting model

C1.5 amplitude exists; v1 catalog has no field-dependent amplitude term. C8.4 `weak_phase_all` is the diagnostic that would demand one.

## L13. Parameter count

Default catalog free count is 15 (C4.1). Implementations count from the catalog; 15 is a checksum for `psf_field_v1_default` with default flags.

## L14. Validation-before-pipeline

C2.8.8 forbids merging FFT code before closed-form tests exist. C10.8 Jacobian CI is on corpus cameras, not toy polynomials only.

## L15. Oversampling

C7.2 forbids even `kS`. Default oversample is 1. No implicit `fft` crop.
