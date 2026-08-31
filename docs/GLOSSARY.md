# Glossary

**ANSI** — American National Standards Institute. Pupil-phase modes use the ANSI Z80.28 two-index scheme \((n, m)\).

**OSA** — Optical Society of America. Sequential OSA index \(j = (n(n+2)+m)/2\) is used for ordering and reports only; the engine API is keyed by \((n, m)\).

**PSF** — Point-Spread Function. This crate's Zernike engine evaluates orthonormal phase on the pupil; the detector PSF is produced later by the forward pipeline.

**RMS** — Root Mean Square. Zernike coefficients are in waves RMS over the supplied pupil mask after discrete normalization.
