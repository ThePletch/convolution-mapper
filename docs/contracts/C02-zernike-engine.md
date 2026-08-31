# C2 — Zernike engine

A single generator parameterized by integers \((n, m)\). Named aberrations are catalog data (C3), not code.

## C2.1 Indexing convention (frozen)

**ANSI/OSA two-index** (Thibos et al. 2002 / ANSI Z80.28):

- \(n \in \mathbb{N}_0\) (radial order), \(m \in \mathbb{Z}\), \(|m| \le n\), \(n - |m|\) even.
- Invalid \((n, m)\) SHALL be rejected at catalog ingest, not at evaluation.
- Sequential OSA index (for ordering and reports only):

\[
j = \frac{n(n+2) + m}{2}
\]

- The engine API SHALL take `(n, m)`, never Noll `Z_7` integers as primary keys.
- Noll, Fringe, and Zemax “standard” numbers SHALL NOT appear in catalog files. A conversion table MAY exist in documentation only.

Sign convention for azimuthal parts (frozen, right-handed, NOR.10):

\[
\begin{aligned}
m > 0 &\Rightarrow \text{cosine mode } \cos(m\theta), \\
m < 0 &\Rightarrow \text{sine mode } \sin(|m|\theta), \\
m = 0 &\Rightarrow \text{no azimuthal factor}.
\end{aligned}
\]

## C2.2 Radial polynomial (frozen)

\[
R_n^{|m|}(\rho) = \sum_{k=0}^{(n-|m|)/2}
\frac{(-1)^k (n-k)!}{k!\,\bigl(\tfrac{n+|m|}{2}-k\bigr)!\,\bigl(\tfrac{n-|m|}{2}-k\bigr)!}
\,\rho^{n-2k}
\]

The sum is the identity; any evaluation that matches it in `f64` for \(n \le 15\) is conformant. Implementations SHOULD use a stable product/Horner form rather than naive factorial ratios. **Frozen max \(n\):** 15 for v1 evaluation; catalog ingest SHALL reject \(n > 15\).

## C2.3 Analytic orthonormalization on the unit disk

OSA/ANSI RMS factor (frozen):

\[
N_n^m = \sqrt{\frac{2(n+1)}{1+\delta_{m,0}}}
\]

Analytic function on the unit disk:

\[
\tilde Z_n^m(\rho,\theta) =
\begin{cases}
N_n^m R_n^{|m|}(\rho) \cos(m\theta) & m > 0 \\
N_n^m R_n^{|m|}(\rho) \sin(|m|\theta) & m < 0 \\
N_n^0 R_n^{0}(\rho) & m = 0
\end{cases}
\]

This \(\tilde Z\) satisfies \(\frac{1}{\pi}\int_{\rho\le 1} (\tilde Z)^2\,dA = 1\) (RMS = 1 on the continuous disk).

## C2.4 Discrete RMS normalization on the mask (frozen; overrides C2.3 at runtime)

Let \(M_{pq} \in \{0,1\}\) be the pupil mask (C9.4) on the \(N_p\times N_p\) grid. Let \(S = \sum_{p,q} M_{pq}\) (number of unmasked pixels). \(S=0\) is a fatal `InputError`.

Evaluate \(\tilde Z\) at each pixel with \(\rho\le 1\); pixels with \(M=0\) are treated as 0.

Discrete RMS:

\[
\mathrm{rms} = \sqrt{\frac{1}{S}\sum_{p,q} M_{pq}\, \tilde Z_{pq}^{2}}
\]

If \(\mathrm{rms} < 10^{-15}\), reject (would occur for a mode identically zero on the mask).

\[
Z_{pq} = \frac{\tilde Z_{pq}}{\mathrm{rms}}
\]

so \(\frac{1}{S}\sum M Z^{2} = 1\). Coefficients \(a_k\) in waves RMS are then RMS over the **supplied mask**, not the continuous disk.

For the v1 circular mask and \(N_p \ge 256\), \(|\mathrm{rms} - 1| < 0.02\) for all modes with \(n \le 8\) is a C2.8 check; the runtime path still divides by `rms`.

## C2.5 Gram matrix (reported, not used to orthogonalize further)

\[
G_{ij} = \frac{1}{S}\sum_{p,q} M_{pq} Z^{(i)}_{pq} Z^{(j)}_{pq}
\]

v1 SHALL compute \(G\) once per `PupilSpec` + catalog phase-term set. C8 reports it. v1 SHALL NOT apply Gram–Schmidt. Off-diagonal \(|G_{ij}|\) for \(i\neq j\) on the circular v1 mask with \(n\le 6\) SHALL be \(< 0.05\) (C2.8). Larger values on a non-circular mask are expected and SHALL be reported, not treated as a bug.

## C2.6 Phase screen

\[
\Phi_{pq} = 2\pi \sum_k a_k Z^{(k)}_{pq} \quad [\mathrm{rad}]
\]

(NOR.6). Piston \(Z_0^0\) does not affect \(|\mathrm{FT}\{Ae^{i\Phi}\}|^2\). Catalog v1 freezes piston (C3). The engine SHALL still evaluate it if unfrozen (for tests).

Complex pupil (C9.5) uses this \(\Phi\).

## C2.7 Analytic derivative of the phase

\[
\frac{\partial \Phi}{\partial a_k} = 2\pi Z^{(k)}
\]

The PSF Jacobian formula is C9.8; it consumes \(Z^{(k)}\) and the extra factor \(2\pi\).

## C2.8 Closed-form validation cases (load-bearing; write before the pipeline)

All checks use the v1 circular mask, \(N_p=256\), \(N_{\mathrm{fft}}=1024\), unit flux, centroid at stamp center, no kernels, `S=31`, pixel scale and \(\lambda, D, f\) from C10.1 standard camera.

Let \(I(a)\) be the stamp intensity. Let first-moment centroid be

\[
\bar x = \frac{\sum i\, I}{\sum I}, \quad \bar y = \frac{\sum j\, I}{\sum I}
\]

in stamp-local pixels. Let polar radius from stamp center \(c_\star=(S-1)/2\).

### C2.8.1 Piston independence

\(I(a_{0,0}=1)\) and \(I(a_{0,0}=0)\) SHALL satisfy \(\max |I_1-I_0| / \max I_0 < 10^{-10}\).

### C2.8.2 Zero coefficients = Airy, centered, azimuthal

With all \(a_k=0\):

- Peak pixel SHALL be the unique maximum and SHALL lie at the stamp-center pixel (index \((c_\star, c_\star)\)).
- **Azimuthal anisotropy (not radial slope).** Let \(r_p\) be the distance in stamp pixels from \(c_\star\). Exclude the center (\(r_p < 0.5\,\mathrm{px}\)). Bin the remaining pixels into annuli of width **0.25 px**. An annulus is tested only if it contains \(\ge 8\) pixels. Let \(\mu\) be the mean intensity in the annulus. Relative azimuthal RMS \(\mathrm{rms}(I-\mu)/\max(\mu, 10^{-15})\) SHALL be \(< 10^{-4}\). This catches grid/fftshift anisotropy and Seidel-like \(m=1\) leakage; it SHALL NOT be read as a bound on the Airy radial gradient.

### C2.8.3 Defocus evenness

For \(a_{2,0} = \alpha\) and \(a_{2,0}=-\alpha\) with \(\alpha=0.3\) waves, all other \(a=0\):

\[
\max |I(\alpha)-I(-\alpha)| / \max I(0) < 10^{-8}
\]

and each \(I(\pm\alpha)\) SHALL pass the same azimuthal-anisotropy test as C2.8.2 (annuli, \(10^{-4}\)).

### C2.8.4 Defocus at zero has vanishing Jacobian column

At all \(a=0\), the Jacobian column for \(a_{2,0}\) (C9.8, then resampled and weighted as C5) SHALL have \(\ell_2\) norm \(< 10^{-8}\) times that column’s norm at \(a_{2,0}=0.2\). (Even function ⇒ zero derivative at 0.)

### C2.8.5 Zernike coma is centroid-preserving

For \(a_{3,1}=0.4\) waves, all other \(a=0\):

\[
|\bar x - c_\star| < 5\times 10^{-3}\ \mathrm{px}, \qquad |\bar y - c_\star| < 5\times 10^{-3}\ \mathrm{px}
\]

Repeat for \(a_{3,-1}=0.4\). A Seidel (unbalanced \(\rho^3\cos\theta\)) implementation SHALL fail this test; that is the point.

### C2.8.6 Coma sign flips the flare

The first-moment of \(I^2\) restricted to pixels with \(I > 0.05 \max I\) SHALL reverse in \(x\) when \(a_{3,1}\) flips sign (the two centroids SHALL be on opposite sides of \(c_\star\) in \(x\), each at least 0.02 px from \(c_\star\)).

### C2.8.7 Known analytic \(\tilde Z\) samples

Split the continuous identity from the sampled grid. Do **not** stretch C9.2 \(\xi,\eta\) to force a \(\rho=1\) sample.

- **Continuous formula (no grid):** \(\tilde Z_2^0(\rho=1,\theta=0)=\sqrt{3}\) and \(\tilde Z_2^0(\rho=0)=-\sqrt{3}\) SHALL match the closed form to 1 ulp in `f64`.
- **Sampled, pre-normalization \(\tilde Z\):** at the pixel nearest \((\rho=1,\theta=0)\) on the C9.2 grid, \(|\tilde Z-\sqrt{3}|<0.03\). At the pixel(s) nearest \(\rho=0\), \(|\tilde Z+\sqrt{3}|<0.02\). (On \(N_p=256\), max on-axis \(\rho=(N_p-1)/N_p\approx 0.996\), so the rim sample is not \(\sqrt{3}\).) After discrete RMS the stored \(Z\) MAY differ.

### C2.8.8 Order of writing

These tests SHALL exist in the repository before the FFT pipeline is merged. They MAY call a temporary stub that returns zeros only if the test is marked `#[ignore]` — they SHALL NOT be deleted to make CI green.

## C2.9 Engine outputs

Given `(n,m)` and a `PupilSpec`:

- `basis_screen: f64[N_p, N_p]` — \(Z\) from C2.4
- `phase_screen(coefficients): f64[N_p, N_p]` — \(\Phi\) from C2.6

No PSF is computed here. PSF is C9.

## C2.10 Tilt modes

\(Z_1^{\pm 1}\) SHALL be implemented (the generator is generic). The **v1 default catalog SHALL freeze and exclude them** (C3.7) because they are degenerate with the extraction centroid. Tests MAY unfreeze them.
