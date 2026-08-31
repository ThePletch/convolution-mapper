# C9 — Forward pipeline numerics

The shared nonlinear map from (coefficients, kernels, centroid) to a stamp-sized model. Modules SHALL NOT produce PSF-space additive errors.

## C9.1 Stage order (frozen)

1. Phase + amplitude on the pupil (C2, C9.5)
2. Zero-pad
3. FFT → coherent field
4. Intensity \(|U|^2\)
5. Normalize to unit sum on the FFT grid
6. Convolve kernels (C9.9), each kernel already unit-sum
7. Fourier-shift by stamp centroid offset (C1.2.3)
8. Resample to detector pixel grid of size \(S\times S\) (C9.10)
9. Multiply by flux \(F\) and add residual sky \(b\) (C4)

## C9.2 Pupil grid size

- \(N_p\) = `PupilSpec.n_pupil`. Frozen default **256**. Allowed `{128, 256, 512}`.
- Coordinates of pixel `[p, q]` (row \(p\), column \(q\)):

\[
\xi_q = \frac{q - (N_p-1)/2}{N_p/2}, \qquad
\eta_p = \frac{p - (N_p-1)/2}{N_p/2}
\]

so the nominal pupil radius 1 is at a distance \(N_p/2\) samples from the array center \(((N_p-1)/2, (N_p-1)/2)\).

## C9.3 FFT size and padding

- \(N_f\) = `PupilSpec.n_fft`. Frozen default **1024** when \(N_p=256\). Constraint: \(N_f / N_p \in \{2,4,8\}\), \(N_f\) even.
- Embed the \(N_p\times N_p\) complex pupil into an \(N_f\times N_f\) array of zeros, **centered**: pupil pixel `[p,q]` maps to FFT index

\[
p' = p + (N_f - N_p)/2, \qquad q' = q + (N_f - N_p)/2
\]

(integer: \(N_f-N_p\) is even).

## C9.4 v1 circular mask

\[
M_{pq} = \begin{cases} 1 & \rho_{pq} \le 1 \\ 0 & \text{otherwise} \end{cases}
\]

\(\rho\) from NOR.10. Pixels with \(\rho=1\) (on the boundary) are **included**.

## C9.5 Complex pupil

Amplitude \(A_{pq}\): `PupilSpec.amplitude` if provided, else \(A = M\).

\[
P_{pq} = A_{pq}\, M_{pq}\, \exp\bigl(i \Phi_{pq}\bigr)
\]

with \(\Phi$ from C2.6. \(A\) SHALL be 0 where \(M=0\) (ingest shall zero it).

## C9.6 FFT definition

- `rustfft` forward DFT, **unnormalized**:

\[
U_{kl} = \sum_{p'=0}^{N_f-1}\sum_{q'=0}^{N_f-1} P^{\mathrm{pad}}_{p'q'}
\exp\bigl(-2\pi i (k p' + l q') / N_f\bigr)
\]

- Then **fftshift** so that DC (zero frequency) is at index \((N_f/2, N_f/2)\).
- Intensity \(I^{\mathrm{fft}}_{kl} = |U_{kl}|^2\).
- Let \(E = \sum_{k,l} I^{\mathrm{fft}}_{kl}\). If \(E \le 0\) or non-finite → `NumericsError`.
- Unit-sum PSF on the FFT grid: \(I^{\mathrm{fft}} \leftarrow I^{\mathrm{fft}} / E\).

The overall DFT scale is absorbed by this normalization. Implementations SHALL NOT apply extra `1/N_f` factors after this step.

## C9.7 Angular scale of FFT pixels (frozen)

Pupil sample spacing \(\Delta\xi_{\mathrm{m}} = D / N_p\) where \(D =\) `pupil_diameter_m`.

Angular size of one FFT pixel (radians):

\[
\Delta\alpha = \frac{\lambda}{N_f \cdot \Delta\xi_{\mathrm{m}}} = \frac{\lambda}{D}\cdot \frac{N_p}{N_f}
\]

Detector pixel scale \(\Delta\beta =\) `pixel_scale_arcsec` converted to radians (NOR.6).

Oversampling factor (FFT pixels per detector pixel):

\[
r = \frac{\Delta\beta}{\Delta\alpha} = \mathrm{pixel\_scale\_rad}\cdot \frac{D}{\lambda}\cdot \frac{N_f}{N_p}
\]

If \(r < 0.5\) or \(r > 64\), ingest SHALL reject (`InputError`: sampling insane). For the C10.1 camera, \(r\) is stored in the corpus header for debugging; the pipeline uses the formula, not a frozen \(r\).

## C9.8 Analytic PSF Jacobian w.r.t. a phase coefficient (on the FFT grid, before kernels)

Let \(U = \mathrm{fftshift}(\mathrm{DFT}(P^{\mathrm{pad}}))\). For mode \(k\):

\[
Q^{(k)}_{pq} = i \cdot 2\pi Z^{(k)}_{pq} P_{pq}
\]

Pad \(Q^{(k)}\) identically to \(P\). Let \(V^{(k)} = \mathrm{fftshift}(\mathrm{DFT}(Q^{(k),\mathrm{pad}}))\).

\[
\frac{\partial I^{\mathrm{fft}}}{\partial a_k} = 2\,\mathrm{Re}\bigl(\overline{U}\circ V^{(k)}\bigr)
\]

(elementwise). Then apply the same unit-sum constraint’s derivative: because \(I^{\mathrm{fft}} = J / E\) with \(E=\sum J\), \(J=|U|^2\),

\[
\frac{\partial I^{\mathrm{fft}}}{\partial a_k}
\leftarrow
\frac{1}{E}\frac{\partial J}{\partial a_k}
-
\frac{J}{E^2}\sum_{k'l'} \frac{\partial J_{k'l'}}{\partial a_k}
\]

with \(\partial J/\partial a_k = 2\mathrm{Re}(\overline U \circ V^{(k)})\) **before** the normalization in C9.6. Implementations SHALL use this exact two-term expression (the second term is the derivative of the unit-sum projection). Omitting the second term is non-conformant.

Kernels and resampling are linear; they apply to this derivative by the chain rule (C9.9–C9.10).

## C9.9 Convolution kernels

After C9.6, for each **active** kernel in catalog order (C3.6):

\[
I \leftarrow I \otimes K
\]

implemented as multiplication in Fourier space of the **already fftshifted** arrays, using the same \(N_f\) grid:

- Build \(K\) on the FFT grid in **centered** spatial coordinates (origin at \(N_f/2, N_f/2\)), in **detector-pixel units converted to FFT pixels** via \(r\) (C9.7): a kernel length of \(s\) detector pixels is \(s \cdot r\) FFT pixels.
- \(K \ge 0\), then \(K \leftarrow K / \sum K\). If \(\sum K = 0\), `NumericsError`.
- \(\hat K = \mathrm{DFT}(\mathrm{ifftshift}(K))\) (so the kernel’s spatial origin is at array index (0,0) for the DFT convolution theorem). Then \(U_I = \mathrm{DFT}(\mathrm{ifftshift}(I))\), \(I \leftarrow \mathrm{fftshift}(\mathrm{IDFT}(U_I \circ \hat K))\) taking the real part, clip imaginary residual: if \(\max|\mathrm{Im}| > 10^{-10} \max|I|\), `NumericsError`.
- Re-apply unit-sum: \(I \leftarrow I / \sum I\).

**Analytic kernel derivatives:** \(\partial I / \partial \alpha = I_0 \otimes \partial K/\partial\alpha\) with the same FFT convolution, then the unit-sum projection analogously. Closed forms of \(\partial K/\partial\alpha\) are in C3.6. After convolution, the phase Jacobian columns SHALL be convolved with \(K\) (not with \(\partial K\)).

Kernel application order is **catalog list order**. Convolution commutes for scalar kernels; sequential apply + re-unit-sum is the **reference**. Fused \(\hat K=\prod_i \hat K_i\) then one multiply is conformant iff \(\max|I_{\mathrm{fused}}-I_{\mathrm{ref}}|<10^{-12}\) on the C10.1 unaberrated FFT grid. Jacobian: phase columns still convolve with the **product** \(K\); kernel column \(i\) uses \(\partial K_i\) and \(\prod_{j\neq i} K_j\).

## C9.10 Detector resampling

Goal: produce an \(S\times S\) stamp whose pixel \((i,j)\) matches C1.2 geometry, with the optical PSF shifted so that the model centroid is at `centroid_xy_px`.

### C9.10.1 Fourier shift (subpixel)

Let \(\Delta_x = c_x - c_\star\), \(\Delta_y = c_y - c_\star\) in **detector pixels** (C1.2.3). In FFT pixels: \(\delta_x = \Delta_x \cdot r\), \(\delta_y = \Delta_y \cdot r\).

Multiply the unshifted (ifftshifted) Fourier transform of \(I\) by

\[
\exp\bigl(-2\pi i (k_x \delta_x + k_y \delta_y) / N_f\bigr)
\]

with frequencies \(k_x, k_y\) in `fftfreq` convention (0, …, \(N_f/2-1\), \(-N_f/2\), …, \(-1\)). Inverse DFT, fftshift. This is the unique band-limited shift on the FFT grid.

**Placement vs shift (SHALL):** C9.10.1 uses \(\Delta = c - c_\star\). C9.10.2 places detector pixel \(i\) using \(c_\star\) only. Substituting \(c\) into the placement formula is non-conformant (double shift).

### C9.10.2 Bin / interpolate onto detector pixels

Each detector pixel is a square of side \(r\) FFT pixels (C9.7), not necessarily integer.

**Identity:** area-weighted average (box integration) of \(I^{\mathrm{fft}}\) over the square of side \(r\) centered on the detector pixel’s conjugate location. Nodes that fall outside `[0, N_f)` contribute 0 (no wrap).

Place detector pixel \((i,j)\) (stamp-local, 0-based) at FFT coordinates:

\[
q_{\mathrm{ctr}} = N_f/2 + (i - c_\star)\, r, \qquad
p_{\mathrm{ctr}} = N_f/2 + (j - c_\star)\, r
\]

The integration domain is the axis-aligned square \([q_{\mathrm{ctr}}-r/2, q_{\mathrm{ctr}}+r/2] \times [p_{\mathrm{ctr}}-r/2, p_{\mathrm{ctr}}+r/2]\). v1 SHALL evaluate it by bilinear sampling of \(I^{\mathrm{fft}}\) at **4×4 Gauss–Legendre** nodes per detector pixel.

If the integrated weight is 0, that detector pixel is 0.

After filling the stamp \(m_{ij}\), **do not** re-normalize to unit sum (flux \(F\) is a free parameter). The unit-sum constraint already held on the FFT grid; box integration approximately conserves flux. C10.4 checks \(\sum m\) is within 2% of 1 for the C10.1 camera at zero aberration before flux scaling. That 2% check is a **pipeline** sanity bound, not the resample definition (C9.10.3).

### C9.10.3 Resample harness (load-bearing)

- **Constant image.** \(I^{\mathrm{fft}}=1/N_f^2\) everywhere. Each in-bounds stamp pixel SHALL equal the overlap area of its \(r\times r\) box with \([0,N_f)^2\), relative error \(<10^{-12}\).
- **Rectangle.** An axis-aligned rectangle of ones on the FFT grid, zeros elsewhere. Each stamp pixel SHALL equal (analytic box overlap) / \(r^2\) within \(10^{-8}\).

## C9.11 Flux and sky

\[
\mathrm{model}_{ij} = F \cdot m_{ij} + b
\]

\(F\) in ADU, \(b\) in ADU. Defaults: \(F\) initialized from `flux_sum_adu` (C1), \(b=0\) and frozen unless unfrozen (C4).

## C9.12 What is forbidden

- Additive PSF-space “error images” from aberration modules.
- Wrap-around convolution without the ifftshift convention above.
- Cropping the FFT PSF without the box integral (nearest-neighbor subsample).
- `float32` FFTs.
- Applying kernels in detector space after a lossy resample (kernels SHALL be applied on the FFT grid, C9.9, **before** C9.10.2). Linear drift/rotation kernels that are defined in detector pixels SHALL be converted to FFT pixels via \(r\) before C9.9.

## C9.13 Energy / Parseval sanity (harness)

For \(A=M\), \(\Phi=0\): \(\sum |P|^2\) over the pupil vs \(\sum |U|^2 / N_f^{2}\) SHALL match to relative \(10^{-10}\) (unnormalized DFT Parseval). This test uses **pre-intensity-normalization** \(U\). It is a pipeline unit test, not a fitter test.
