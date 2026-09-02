//! Forward pipeline orchestration. C9.1 steps 1–5 produce unit-sum intensity on
//! the FFT grid. This module then Fourier-shifts and box-resamples onto the
//! detector stamp and applies flux and sky (steps 7–9). Convolution kernels
//! (step 6) are not applied yet.

use ndarray::Array2;
use rustfft::num_complex::Complex;

use crate::error::{ErrorModule, PsfFieldError};
use crate::fftutil::{fftshift, intensity_and_sum, pad_centered, unit_sum_intensity, Fft2D};
use crate::pupil::{complex_pupil, fft_size, grid_size};
use crate::resample::{
    apply_flux_and_sky, box_resample, fourier_shift, oversampling_factor, stamp_center,
};
use crate::types::{check_stamp_size, ImageMeta, PupilSpec};
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

/// Inputs to [`forward_psf`]. Named fields keep flux, sky, and centroid distinct
/// at the call site.
pub struct ForwardPsfSpec<'a> {
    pub pupil: &'a PupilSpec,
    pub image_meta: &'a ImageMeta,
    pub phase_terms: &'a [PhaseCoefficient],
    /// Stamp-local centroid (C1.2.3). Placement uses only c_★; this enters solely
    /// through the Fourier shift as c − c_★. (C9.10.1)
    pub centroid_xy_px: [f64; 2],
    pub stamp_size: usize,
    pub flux_adu: f64,
    pub sky_adu: f64,
}

/// Detector-stamp model: C9.1 steps 1–5, skip kernels, Fourier-shift by
/// c − c_★, Gauss–Legendre box resample, then F · m + b. (C9.1, C9.10, C9.11)
pub fn forward_psf(
    spec: &ForwardPsfSpec<'_>,
    fft: &mut Fft2D,
) -> Result<Array2<f64>, PsfFieldError> {
    check_stamp_size(spec.stamp_size)?;
    if !spec.centroid_xy_px[0].is_finite() || !spec.centroid_xy_px[1].is_finite() {
        return Err(PsfFieldError::input(
            ErrorModule::Pipeline,
            "centroid_xy_px must be finite",
        ));
    }
    let oversampling = oversampling_factor(spec.image_meta, spec.pupil)?;
    let grid = fft_grid_intensity(spec.pupil, spec.phase_terms, fft)?;
    let c_star = stamp_center(spec.stamp_size);
    // Δ = c − c_★ in detector pixels, converted to FFT pixels by r. (C9.10.1)
    let delta_fft_x = (spec.centroid_xy_px[0] - c_star) * oversampling;
    let delta_fft_y = (spec.centroid_xy_px[1] - c_star) * oversampling;
    let shifted = fourier_shift(&grid.intensity, delta_fft_x, delta_fft_y, fft)?;
    let unit_flux = box_resample(&shifted, spec.stamp_size, oversampling)?;
    apply_flux_and_sky(&unit_flux, spec.flux_adu, spec.sky_adu)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pupil::circular_pupil_spec;
    use crate::types::ImageMeta;

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

    fn unique_peak(intensity: &Array2<f64>) -> (usize, usize) {
        let mut peak = f64::NEG_INFINITY;
        let mut peak_at = (0_usize, 0_usize);
        let mut ties = 0_usize;
        for row in 0..intensity.nrows() {
            for column in 0..intensity.ncols() {
                let value = intensity[[row, column]];
                if value > peak {
                    peak = value;
                    peak_at = (row, column);
                    ties = 1;
                } else if value == peak {
                    ties += 1;
                }
            }
        }
        assert_eq!(ties, 1, "peak pixel must be unique");
        peak_at
    }

    fn unaberrated_spec<'a>(
        pupil: &'a PupilSpec,
        meta: &'a ImageMeta,
        centroid_xy_px: [f64; 2],
        stamp_size: usize,
        flux_adu: f64,
        sky_adu: f64,
    ) -> ForwardPsfSpec<'a> {
        ForwardPsfSpec {
            pupil,
            image_meta: meta,
            phase_terms: &[],
            centroid_xy_px,
            stamp_size,
            flux_adu,
            sky_adu,
        }
    }

    #[test]
    fn unaberrated_airy_peak_is_at_stamp_center() {
        let pupil = unit_pupil();
        let meta = ImageMeta::c10_1_standard_camera();
        let mut fft = Fft2D::new(fft_size(&pupil).unwrap()).unwrap();
        let stamp_size = 15_usize;
        let c_star = stamp_center(stamp_size);
        let model = forward_psf(
            &unaberrated_spec(&pupil, &meta, [c_star, c_star], stamp_size, 1.0, 0.0),
            &mut fft,
        )
        .unwrap();
        let peak = unique_peak(&model);
        assert_eq!(peak, (c_star as usize, c_star as usize));
        assert!(model.iter().all(|&value| value >= -1e-18));
    }

    #[test]
    fn flux_and_sky_scale_the_resampled_stamp() {
        let pupil = unit_pupil();
        let meta = ImageMeta::c10_1_standard_camera();
        let mut fft = Fft2D::new(fft_size(&pupil).unwrap()).unwrap();
        let stamp_size = 15_usize;
        let c_star = stamp_center(stamp_size);
        let unit = forward_psf(
            &unaberrated_spec(&pupil, &meta, [c_star, c_star], stamp_size, 1.0, 0.0),
            &mut fft,
        )
        .unwrap();
        let scaled = forward_psf(
            &unaberrated_spec(&pupil, &meta, [c_star, c_star], stamp_size, 3.0, 0.5),
            &mut fft,
        )
        .unwrap();
        for (unit_value, scaled_value) in unit.iter().zip(scaled.iter()) {
            assert!((scaled_value - (3.0 * unit_value + 0.5)).abs() < 1e-12);
        }
    }

    #[test]
    fn subpixel_centroid_moves_the_first_moment() {
        let pupil = unit_pupil();
        let meta = ImageMeta::c10_1_standard_camera();
        let mut fft = Fft2D::new(fft_size(&pupil).unwrap()).unwrap();
        let stamp_size = 15_usize;
        let c_star = stamp_center(stamp_size);
        let offset = 0.3;
        let centered = forward_psf(
            &unaberrated_spec(&pupil, &meta, [c_star, c_star], stamp_size, 1.0, 0.0),
            &mut fft,
        )
        .unwrap();
        let shifted = forward_psf(
            &unaberrated_spec(
                &pupil,
                &meta,
                [c_star + offset, c_star],
                stamp_size,
                1.0,
                0.0,
            ),
            &mut fft,
        )
        .unwrap();
        let moment = |stamp: &Array2<f64>| {
            let mut wx = 0.0;
            let mut w = 0.0;
            for row in 0..stamp.nrows() {
                for column in 0..stamp.ncols() {
                    let v = stamp[[row, column]].max(0.0);
                    wx += column as f64 * v;
                    w += v;
                }
            }
            wx / w
        };
        let dx = moment(&shifted) - moment(&centered);
        assert!(
            (dx - offset).abs() < 0.05,
            "first-moment shift {dx}, requested {offset}"
        );
    }
}
