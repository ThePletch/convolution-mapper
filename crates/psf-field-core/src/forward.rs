//! Forward pipeline orchestration. This module currently implements C9.1 steps
//! 1–5 (complex pupil through unit-sum intensity on the FFT grid). Kernels,
//! Fourier shift, and detector resampling are later stages.

use ndarray::Array2;
use rustfft::num_complex::Complex;

use crate::error::{ErrorModule, PsfFieldError};
use crate::fftutil::{fftshift, intensity_and_sum, pad_centered, unit_sum_intensity, Fft2D};
use crate::pupil::{complex_pupil, fft_size, grid_size};
use crate::types::PupilSpec;
use crate::zernike::{phase_screen, PhaseCoefficient};

/// Optical PSF on the FFT grid after C9.1 steps 1–5, before kernels.
#[derive(Debug, Clone)]
pub struct FftGridIntensity {
    /// fftshifted coherent field U. Intensity normalization has not been
    /// applied to this array; Parseval and the two-term Jacobian use it. (C9.6, C9.8, C9.13)
    pub coherent_field: Array2<Complex<f64>>,
    /// I^{fft} = |U|² / E, unit sum over the N_f × N_f grid. (C9.6)
    pub intensity: Array2<f64>,
    /// E = Σ |U|² before dividing. The unit-sum Jacobian's second term needs this. (C9.8)
    pub intensity_sum: f64,
}

/// C9.1 steps 1–5: Φ from the Zernike engine, P = A M exp(iΦ), centered pad,
/// unnormalized DFT, fftshift, |U|², divide by E. No extra 1/N_f factor. (C9.6)
pub fn fft_grid_intensity(
    pupil: &PupilSpec,
    terms: &[PhaseCoefficient],
    fft: &mut Fft2D,
) -> Result<FftGridIntensity, PsfFieldError> {
    let phase = phase_screen(terms, pupil)?;
    fft_grid_intensity_from_phase(pupil, &phase, fft)
}

/// Same as [`fft_grid_intensity`] with Φ already sampled on the pupil grid.
pub fn fft_grid_intensity_from_phase(
    pupil: &PupilSpec,
    phase: &Array2<f64>,
    fft: &mut Fft2D,
) -> Result<FftGridIntensity, PsfFieldError> {
    let n_pupil = grid_size(pupil)?;
    let n_fft = fft_size(pupil)?;
    if fft.n_fft() != n_fft {
        return Err(PsfFieldError::input(
            ErrorModule::Pipeline,
            "Fft2D length does not match PupilSpec.n_fft",
        ));
    }
    if n_fft < n_pupil {
        return Err(PsfFieldError::input(
            ErrorModule::Pipeline,
            "n_fft must be at least n_pupil",
        ));
    }

    let pupil_field = complex_pupil(pupil, phase)?;
    let mut padded = pad_centered(&pupil_field, n_fft, Complex::new(0.0, 0.0))?;
    fft.forward_dft(&mut padded)?;
    let coherent_field = fftshift(&padded);
    let (unnormalized, intensity_sum) = intensity_and_sum(&coherent_field);
    let intensity = unit_sum_intensity(unnormalized, intensity_sum)?;
    Ok(FftGridIntensity {
        coherent_field,
        intensity,
        intensity_sum,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pupil::circular_pupil_spec;

    /// Fast unit-test grids. Production defaults are 256 / 1024. (C9.2, C9.3)
    const UNIT_N_PUPIL: i64 = 32;
    const UNIT_N_FFT: i64 = 64;
    /// Frozen contract grids used by C2.8 / C9.13. (C9.2, C9.3)
    const CONTRACT_N_PUPIL: i64 = 256;
    const CONTRACT_N_FFT: i64 = 1024;
    /// Unnormalized 2-D DFT Parseval relative tolerance. (C9.13)
    const PARSEVAL_RELATIVE_TOLERANCE: f64 = 1e-10;

    fn unit_pupil() -> PupilSpec {
        circular_pupil_spec(UNIT_N_PUPIL, UNIT_N_FFT)
    }

    fn contract_pupil() -> PupilSpec {
        circular_pupil_spec(CONTRACT_N_PUPIL, CONTRACT_N_FFT)
    }

    fn relative_difference(left: f64, right: f64) -> f64 {
        (left - right).abs() / left.abs().max(right.abs()).max(1.0)
    }

    /// Σ |P|² over the pupil vs Σ |U|² / N_f² for A = M and Φ = 0, using
    /// pre-normalization U. (C9.13)
    fn assert_parseval(pupil: &PupilSpec) {
        let n_fft = fft_size(pupil).unwrap();
        let mut fft = Fft2D::new(n_fft).unwrap();
        let phase = Array2::<f64>::zeros((grid_size(pupil).unwrap(), grid_size(pupil).unwrap()));
        let result = fft_grid_intensity_from_phase(pupil, &phase, &mut fft).unwrap();
        let pupil_field = complex_pupil(pupil, &phase).unwrap();
        let mut pupil_energy = 0.0;
        for value in pupil_field.iter() {
            pupil_energy += value.norm_sqr();
        }
        let mut field_energy = 0.0;
        for value in result.coherent_field.iter() {
            field_energy += value.norm_sqr();
        }
        let parseval = field_energy / (n_fft as f64).powi(2);
        assert!(
            relative_difference(pupil_energy, parseval) < PARSEVAL_RELATIVE_TOLERANCE,
            "Parseval relative difference {} (pupil energy {pupil_energy}, Σ|U|²/N_f² {parseval})",
            relative_difference(pupil_energy, parseval)
        );
    }

    #[test]
    fn c9_13_parseval_tiny_grid() {
        assert_parseval(&unit_pupil());
    }

    #[test]
    fn c9_13_parseval_contract_grid() {
        assert_parseval(&contract_pupil());
    }

    #[test]
    fn unit_sum_intensity_adds_to_one() {
        let pupil = unit_pupil();
        let mut fft = Fft2D::new(fft_size(&pupil).unwrap()).unwrap();
        let result = fft_grid_intensity(&pupil, &[], &mut fft).unwrap();
        let sum: f64 = result.intensity.iter().sum();
        assert!((sum - 1.0).abs() < 1e-14);
        assert!(result.intensity.iter().all(|&value| value >= 0.0));
        assert!(result.intensity_sum > 0.0);
        assert!(result.intensity_sum.is_finite());
    }

    #[test]
    fn unaberrated_peak_is_at_fftshift_dc() {
        let pupil = unit_pupil();
        let n_fft = fft_size(&pupil).unwrap();
        let mut fft = Fft2D::new(n_fft).unwrap();
        let result = fft_grid_intensity(&pupil, &[], &mut fft).unwrap();
        let mut peak = f64::NEG_INFINITY;
        let mut peak_at = (0_usize, 0_usize);
        let mut ties = 0_usize;
        for row in 0..n_fft {
            for column in 0..n_fft {
                let value = result.intensity[[row, column]];
                if value > peak {
                    peak = value;
                    peak_at = (row, column);
                    ties = 1;
                } else if value == peak {
                    ties += 1;
                }
            }
        }
        assert_eq!(ties, 1);
        assert_eq!(peak_at, (n_fft / 2, n_fft / 2));
    }

    #[test]
    fn zero_amplitude_is_numerics_error() {
        let mut pupil = unit_pupil();
        let n_pupil = grid_size(&pupil).unwrap();
        pupil.amplitude = Some(vec![vec![0.0; n_pupil]; n_pupil]);
        let mut fft = Fft2D::new(fft_size(&pupil).unwrap()).unwrap();
        let err = fft_grid_intensity(&pupil, &[], &mut fft).unwrap_err();
        assert_eq!(err.code, crate::error::ErrorCode::Numerics);
        assert_eq!(err.module, ErrorModule::Pipeline);
        assert!(err.message.contains("intensity sum"));
    }

    #[test]
    fn non_finite_amplitude_is_numerics_error() {
        let mut pupil = unit_pupil();
        let n_pupil = grid_size(&pupil).unwrap();
        let mut amplitude = vec![vec![1.0; n_pupil]; n_pupil];
        amplitude[n_pupil / 2][n_pupil / 2] = f64::NAN;
        pupil.amplitude = Some(amplitude);
        let mut fft = Fft2D::new(fft_size(&pupil).unwrap()).unwrap();
        let err = fft_grid_intensity(&pupil, &[], &mut fft).unwrap_err();
        assert_eq!(err.code, crate::error::ErrorCode::Numerics);
    }

    #[test]
    fn piston_does_not_change_fft_grid_intensity() {
        let pupil = unit_pupil();
        let mut fft = Fft2D::new(fft_size(&pupil).unwrap()).unwrap();
        let zero = fft_grid_intensity(&pupil, &[], &mut fft).unwrap();
        let piston = fft_grid_intensity(
            &pupil,
            &[PhaseCoefficient {
                n: 0,
                m: 0,
                waves_rms: 1.0,
            }],
            &mut fft,
        )
        .unwrap();
        let mut max_abs_diff = 0.0_f64;
        for (left, right) in zero.intensity.iter().zip(piston.intensity.iter()) {
            max_abs_diff = max_abs_diff.max((left - right).abs());
        }
        assert!(max_abs_diff < 1e-12);
    }
}
