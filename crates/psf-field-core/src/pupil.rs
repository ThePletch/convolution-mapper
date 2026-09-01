//! Pupil-plane coordinates, the v1 circular aperture mask, and the complex pupil.
//! ξ and η are dimensionless pupil Cartesian coordinates with the unit disk at
//! radius 1. (C9.2, C9.4, C9.5)

use ndarray::Array2;
use rustfft::num_complex::Complex;

use crate::error::{ErrorModule, PsfFieldError};
use crate::types::{PupilSpec, SCHEMA_VERSION};

/// Horizontal pupil coordinate ξ of column `column`.
/// ξ = (column − (N_p − 1)/2) / (N_p / 2), so the array center is ξ = 0 and the
/// nominal rim is |ξ| = 1 at a distance N_p/2 samples from center. (C9.2)
#[must_use]
pub fn xi_at_column(column: usize, n_pupil: usize) -> f64 {
    let n = n_pupil as f64;
    (column as f64 - (n - 1.0) / 2.0) / (n / 2.0)
}

/// Vertical pupil coordinate η of row `row`, same centering as [`xi_at_column`].
/// Row is detector y; the polar angle is atan2(η, ξ). (C9.2, NOR.10)
#[must_use]
pub fn eta_at_row(row: usize, n_pupil: usize) -> f64 {
    let n = n_pupil as f64;
    (row as f64 - (n - 1.0) / 2.0) / (n / 2.0)
}

/// Polar (ρ, θ) at pixel `[row, column]`: ρ = √(ξ² + η²), θ = atan2(η, ξ).
/// θ = 0 on +ξ and increases toward +η (right-handed). (NOR.10)
#[must_use]
pub fn rho_theta(row: usize, column: usize, n_pupil: usize) -> (f64, f64) {
    let xi = xi_at_column(column, n_pupil);
    let eta = eta_at_row(row, n_pupil);
    (xi.hypot(eta), eta.atan2(xi))
}

/// Square pupil-grid length, requiring `mask` to match `n_pupil`.
pub fn grid_size(pupil: &PupilSpec) -> Result<usize, PsfFieldError> {
    let n_pupil = pupil.mask.len();
    if n_pupil == 0 {
        return Err(PsfFieldError::input(
            ErrorModule::Pipeline,
            "pupil mask is empty",
        ));
    }
    for row in &pupil.mask {
        if row.len() != n_pupil {
            return Err(PsfFieldError::input(
                ErrorModule::Pipeline,
                "pupil mask must be square",
            ));
        }
    }
    if pupil.n_pupil as usize != n_pupil {
        return Err(PsfFieldError::input(
            ErrorModule::Pipeline,
            "PupilSpec.n_pupil does not match mask shape",
        ));
    }
    Ok(n_pupil)
}

/// FFT length stored on the pupil spec. Must be even so centered padding and
/// fftshift land DC at (N_f/2, N_f/2). (C9.3, C9.6)
pub fn fft_size(pupil: &PupilSpec) -> Result<usize, PsfFieldError> {
    let n_fft = usize::try_from(pupil.n_fft)
        .map_err(|_| PsfFieldError::input(ErrorModule::Pipeline, "n_fft does not fit in usize"))?;
    if n_fft < 2 || n_fft % 2 != 0 {
        return Err(PsfFieldError::input(
            ErrorModule::Pipeline,
            format!("n_fft {n_fft} must be even and at least 2"),
        ));
    }
    Ok(n_fft)
}

/// Circular unobstructed mask: 1 inside and on the unit disk (ρ ≤ 1), else 0.
/// Boundary pixels with ρ exactly 1 are included. (C9.4)
#[must_use]
pub fn circular_mask(n_pupil: usize) -> Vec<Vec<f64>> {
    let mut mask = vec![vec![0.0; n_pupil]; n_pupil];
    for p in 0..n_pupil {
        let eta = eta_at_row(p, n_pupil);
        for q in 0..n_pupil {
            let xi = xi_at_column(q, n_pupil);
            if xi.hypot(eta) <= 1.0 {
                mask[p][q] = 1.0;
            }
        }
    }
    mask
}

/// `PupilSpec` whose mask is [`circular_mask`]. Does not run ingest, so unit
/// tests may use an `n_pupil` outside the production set {128, 256, 512}.
#[must_use]
pub fn circular_pupil_spec(n_pupil: i64, n_fft: i64) -> PupilSpec {
    PupilSpec {
        schema_version: SCHEMA_VERSION.to_string(),
        mask: circular_mask(n_pupil as usize),
        n_pupil,
        n_fft,
        amplitude: None,
    }
}

/// Complex pupil P = A M exp(i Φ). If `amplitude` is omitted, A = M; A is
/// multiplied by M so masked pixels are identically zero even if a caller
/// skipped ingest. Φ is in radians (C2.6). (C9.5)
pub fn complex_pupil(
    pupil: &PupilSpec,
    phase: &Array2<f64>,
) -> Result<Array2<Complex<f64>>, PsfFieldError> {
    let n_pupil = grid_size(pupil)?;
    if phase.nrows() != n_pupil || phase.ncols() != n_pupil {
        return Err(PsfFieldError::input(
            ErrorModule::Pipeline,
            "phase screen shape must match the pupil mask",
        ));
    }
    let mut field = Array2::from_elem((n_pupil, n_pupil), Complex::new(0.0, 0.0));
    for row in 0..n_pupil {
        for column in 0..n_pupil {
            let mask = pupil.mask[row][column];
            if mask == 0.0 {
                continue;
            }
            let amplitude = match &pupil.amplitude {
                Some(values) => values[row][column] * mask,
                None => mask,
            };
            let phi = phase[[row, column]];
            field[[row, column]] = Complex::from_polar(amplitude, phi);
        }
    }
    Ok(field)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn c9_2_grid_at_array_center() {
        let n = 256_usize;
        let center = (n as f64 - 1.0) / 2.0;
        assert_eq!(xi_at_column(0, n), (0.0 - center) / (n as f64 / 2.0));
        assert_eq!(eta_at_row(0, n), (0.0 - center) / (n as f64 / 2.0));
        let (rho_origin_corner, _) = rho_theta(0, 0, n);
        assert!(rho_origin_corner > 1.0);
    }

    #[test]
    fn c9_4_includes_unit_disk_boundary() {
        let mask = circular_mask(32);
        for p in 0..32 {
            for q in 0..32 {
                let (rho, _) = rho_theta(p, q, 32);
                let expected = if rho <= 1.0 { 1.0 } else { 0.0 };
                assert_eq!(mask[p][q], expected);
            }
        }
    }

    #[test]
    fn complex_pupil_is_zero_off_mask_even_if_amplitude_leaks() {
        let mut pupil = circular_pupil_spec(32, 64);
        pupil.amplitude = Some(vec![vec![0.5; 32]; 32]);
        let phase = Array2::<f64>::zeros((32, 32));
        let field = complex_pupil(&pupil, &phase).unwrap();
        for row in 0..32 {
            for column in 0..32 {
                if pupil.mask[row][column] == 0.0 {
                    assert_eq!(field[[row, column]], Complex::new(0.0, 0.0));
                } else {
                    assert!((field[[row, column]] - Complex::new(0.5, 0.0)).norm() < 1e-15);
                }
            }
        }
    }

    #[test]
    fn omitted_amplitude_equals_the_mask_at_zero_phase() {
        let pupil = circular_pupil_spec(32, 64);
        let phase = Array2::<f64>::zeros((32, 32));
        let field = complex_pupil(&pupil, &phase).unwrap();
        for row in 0..32 {
            for column in 0..32 {
                let expected = Complex::new(pupil.mask[row][column], 0.0);
                assert_eq!(field[[row, column]], expected);
            }
        }
    }
}
