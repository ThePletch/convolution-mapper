# Glossary

| Acronym | Definition |
|---|---|
| ADU | Analog-to-digital unit. Native detector image intensity. |
| ANSI | American National Standards Institute. Pupil-phase modes use the ANSI Z80.28 two-index scheme \((n, m)\). |
| DFT | Discrete Fourier Transform. The pupil-to-field step is rustfft's unnormalized forward DFT; overall scale is absorbed by unit-sum intensity. |
| FFT | Fast Fourier Transform. |
| FWHM | Full width at half maximum. |
| LM | Levenberg–Marquardt nonlinear least squares. |
| OSA | Optical Society of America. Sequential OSA index \(j = (n(n+2)+m)/2\) is used for ordering and reports only; the engine API is keyed by \((n, m)\). |
| PSF | Point-spread function. This crate's Zernike engine evaluates orthonormal phase on the pupil; the detector PSF is produced later by the forward pipeline. |
| RMS | Root mean square. Zernike coefficients are in waves RMS over the supplied pupil mask after discrete normalization. |
