//! Fourier shift, Gauss–Legendre box resample onto the detector stamp, and
//! flux/sky scaling. Kernels are applied on the FFT grid before this module
//! runs; v1 still skips them until that stage is built. (C9.10, C9.11)

use ndarray::Array2;
use rustfft::num_complex::Complex;

use crate::error::{ErrorModule, PsfFieldError};
use crate::fftutil::{fftfreq_bin, fftshift, ifftshift, Fft2D};
use crate::types::{check_stamp_size, ImageMeta, PupilSpec};

/// Allowed oversampling r = detector pixel / FFT pixel. Outside this interval
/// the box integral is either undersampled or pathologically fine. (C9.7)
const OVERSAMPLING_MIN: f64 = 0.5;
const OVERSAMPLING_MAX: f64 = 64.0;

/// 4×4 Gauss–Legendre nodes per detector pixel. Degree-7 exact, so bilinear
/// (degree 1 in each axis) is integrated exactly. (C9.10.2)
const GAUSS_LEGENDRE_ORDER: usize = 4;

/// Same imaginary-residual clip as kernel convolution: a real intensity that
/// acquires a large imaginary part after the shift DFT is a numerics failure.
/// (C9.9, applied to C9.10.1)
const SHIFT_IMAGINARY_RELATIVE: f64 = 1e-10;

/// FFT pixels per detector pixel:
/// r = pixel_scale_rad · D / λ · N_f / N_p.
/// Rejected outside [0.5, 64] as sampling insane. (C9.7)
pub fn oversampling_factor(
    image_meta: &ImageMeta,
    pupil: &PupilSpec,
) -> Result<f64, PsfFieldError> {
    let n_pupil = pupil.n_pupil as f64;
    let n_fft = pupil.n_fft as f64;
    if n_pupil <= 0.0 || n_fft <= 0.0 {
        return Err(PsfFieldError::input(
            ErrorModule::Pipeline,
            "n_pupil and n_fft must be positive to form the oversampling factor",
        ));
    }
    let r = image_meta.pixel_scale_rad() * image_meta.pupil_diameter_m / image_meta.wavelength_m
        * n_fft
        / n_pupil;
    if !r.is_finite() || r < OVERSAMPLING_MIN || r > OVERSAMPLING_MAX {
        return Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!(
                "sampling insane: oversampling factor r={r} is outside [{OVERSAMPLING_MIN}, {OVERSAMPLING_MAX}]"
            ),
        ));
    }
    Ok(r)
}

/// Band-limited shift of an already-fftshifted intensity by `(delta_fft_x, delta_fft_y)`
/// FFT pixels. Positive `delta_fft_x` moves content toward increasing column (detector +x).
///
/// The unshifted DFT of `I` is multiplied by
/// exp(−2πi (k_x δ_x + k_y δ_y) / N_f) with k in the fftfreq integer convention,
/// then inverse-DFT and fftshift. A zero offset is the identity up to roundoff. (C9.10.1)
pub fn fourier_shift(
    intensity: &Array2<f64>,
    delta_fft_x: f64,
    delta_fft_y: f64,
    fft: &mut Fft2D,
) -> Result<Array2<f64>, PsfFieldError> {
    let n_fft = fft.n_fft();
    if intensity.nrows() != n_fft || intensity.ncols() != n_fft {
        return Err(PsfFieldError::input(
            ErrorModule::Pipeline,
            "Fourier-shift intensity must be (n_fft, n_fft)",
        ));
    }
    if !delta_fft_x.is_finite() || !delta_fft_y.is_finite() {
        return Err(PsfFieldError::numerics(
            ErrorModule::Pipeline,
            "Fourier-shift offsets must be finite",
        ));
    }

    let unshifted = ifftshift(intensity);
    let mut spectrum = Array2::from_elem((n_fft, n_fft), Complex::new(0.0, 0.0));
    for row in 0..n_fft {
        for column in 0..n_fft {
            spectrum[[row, column]] = Complex::new(unshifted[[row, column]], 0.0);
        }
    }
    fft.forward_dft(&mut spectrum)?;

    let two_pi = 2.0 * std::f64::consts::PI;
    let n = n_fft as f64;
    for row in 0..n_fft {
        let k_y = fftfreq_bin(row, n_fft);
        for column in 0..n_fft {
            let k_x = fftfreq_bin(column, n_fft);
            let phase = -two_pi * (k_x * delta_fft_x + k_y * delta_fft_y) / n;
            spectrum[[row, column]] *= Complex::from_polar(1.0, phase);
        }
    }

    fft.inverse_dft(&mut spectrum)?;
    let scale = n * n;
    let mut spatial = Array2::<f64>::zeros((n_fft, n_fft));
    let mut max_abs_real = 0.0_f64;
    let mut max_abs_imag = 0.0_f64;
    for row in 0..n_fft {
        for column in 0..n_fft {
            let value = spectrum[[row, column]] / scale;
            max_abs_real = max_abs_real.max(value.re.abs());
            max_abs_imag = max_abs_imag.max(value.im.abs());
            spatial[[row, column]] = value.re;
        }
    }
    if max_abs_imag > SHIFT_IMAGINARY_RELATIVE * max_abs_real.max(1e-30) {
        return Err(PsfFieldError::numerics(
            ErrorModule::Pipeline,
            format!(
                "Fourier-shift imaginary residual {max_abs_imag} exceeds {SHIFT_IMAGINARY_RELATIVE} of the real peak"
            ),
        ));
    }
    Ok(fftshift(&spatial))
}

/// Stamp-local geometric center c_★ = (S − 1)/2. Detector pixel (i, j) is placed
/// using c_★ only; the Fourier shift carries c − c_★. Substituting the measured
/// centroid into this placement would double-count the offset. (C1.2.3, C9.10.1–2)
#[must_use]
pub fn stamp_center(stamp_size: usize) -> f64 {
    (stamp_size as f64 - 1.0) / 2.0
}

/// FFT coordinates of the center of stamp pixel `[row, column]`.
/// q_ctr = N_f/2 + (column − c_★) r, p_ctr = N_f/2 + (row − c_★) r. (C9.10.2)
#[must_use]
pub fn detector_pixel_fft_center(
    row: usize,
    column: usize,
    n_fft: usize,
    stamp_center: f64,
    oversampling: f64,
) -> (f64, f64) {
    let origin = n_fft as f64 / 2.0;
    let p_ctr = origin + (row as f64 - stamp_center) * oversampling;
    let q_ctr = origin + (column as f64 - stamp_center) * oversampling;
    (p_ctr, q_ctr)
}

/// Area-weighted box integral of `intensity` onto an S×S detector stamp.
///
/// Each detector pixel is the square of side `oversampling` (FFT pixels) centered
/// at [`detector_pixel_fft_center`]. The integral is 4×4 Gauss–Legendre quadrature
/// of the bilinear interpolant. Nodes outside [0, N_f) contribute 0 (no wrap).
/// If every node misses the array, the pixel is 0.
///
/// The stored value is the integral, not the mean, so a unit-sum FFT image remains
/// approximately unit-sum on the stamp (C10.4). Do not re-normalize. (C9.10.2)
pub fn box_resample(
    intensity: &Array2<f64>,
    stamp_size: usize,
    oversampling: f64,
) -> Result<Array2<f64>, PsfFieldError> {
    check_stamp_size(stamp_size)?;
    let n_fft = intensity.nrows();
    if intensity.ncols() != n_fft {
        return Err(PsfFieldError::input(
            ErrorModule::Pipeline,
            "resample intensity must be square",
        ));
    }
    if !oversampling.is_finite() || oversampling <= 0.0 {
        return Err(PsfFieldError::input(
            ErrorModule::Pipeline,
            format!("oversampling factor {oversampling} must be finite and positive"),
        ));
    }

    let nodes = gauss_legendre_4();
    let c_star = stamp_center(stamp_size);
    let half = oversampling / 2.0;
    let mut stamp = Array2::<f64>::zeros((stamp_size, stamp_size));
    for row in 0..stamp_size {
        for column in 0..stamp_size {
            let (p_ctr, q_ctr) =
                detector_pixel_fft_center(row, column, n_fft, c_star, oversampling);
            // ∫∫ f dp dq over a square of side r = (r/2)² Σ_i Σ_j w_i w_j f(p_i, q_j).
            // Nodes outside [0, N_f) sample as 0, so a box that misses the array is 0. (C9.10.2)
            let mut integral = 0.0;
            for &(xi_p, w_p) in &nodes {
                let p = p_ctr + half * xi_p;
                for &(xi_q, w_q) in &nodes {
                    let q = q_ctr + half * xi_q;
                    integral += w_p * w_q * bilinear_sample(intensity, p, q);
                }
            }
            stamp[[row, column]] = integral * half * half;
        }
    }
    Ok(stamp)
}

/// model = F · m + b, F and b in ADU. Non-finite F or b is a numerics failure,
/// not a silent clip. (C9.11, NOR.3)
pub fn apply_flux_and_sky(
    unit_flux_stamp: &Array2<f64>,
    flux_adu: f64,
    sky_adu: f64,
) -> Result<Array2<f64>, PsfFieldError> {
    if !flux_adu.is_finite() || !sky_adu.is_finite() {
        return Err(PsfFieldError::numerics(
            ErrorModule::Pipeline,
            "flux and sky must be finite ADU",
        ));
    }
    Ok(unit_flux_stamp * flux_adu + sky_adu)
}

/// 4-point Gauss–Legendre nodes and weights on [−1, 1], ordered −outer … +outer.
fn gauss_legendre_4() -> [(f64, f64); GAUSS_LEGENDRE_ORDER] {
    let sqrt_six_fifths = (6.0_f64 / 5.0).sqrt();
    let inner = ((3.0 - 2.0 * sqrt_six_fifths) / 7.0).sqrt();
    let outer = ((3.0 + 2.0 * sqrt_six_fifths) / 7.0).sqrt();
    let w_inner = (18.0 + 30.0_f64.sqrt()) / 36.0;
    let w_outer = (18.0 - 30.0_f64.sqrt()) / 36.0;
    [
        (-outer, w_outer),
        (-inner, w_inner),
        (inner, w_inner),
        (outer, w_outer),
    ]
}

/// Bilinear sample of `image` at continuous FFT coordinates `(row, column)`.
/// A node outside [0, N_f) is 0; a neighbor that would wrap is also 0. (C9.10.2)
fn bilinear_sample(image: &Array2<f64>, row: f64, column: f64) -> f64 {
    let n = image.nrows() as f64;
    if row < 0.0 || row >= n || column < 0.0 || column >= n {
        return 0.0;
    }
    let y0 = row.floor();
    let x0 = column.floor();
    let ty = row - y0;
    let tx = column - x0;
    let y0i = y0 as usize;
    let x0i = x0 as usize;
    let y1i = y0i + 1;
    let x1i = x0i + 1;
    let n_usize = image.nrows();
    let v00 = image[[y0i, x0i]];
    let v01 = if x1i < n_usize {
        image[[y0i, x1i]]
    } else {
        0.0
    };
    let v10 = if y1i < n_usize {
        image[[y1i, x0i]]
    } else {
        0.0
    };
    let v11 = if y1i < n_usize && x1i < n_usize {
        image[[y1i, x1i]]
    } else {
        0.0
    };
    (1.0 - tx) * (1.0 - ty) * v00 + tx * (1.0 - ty) * v01 + (1.0 - tx) * ty * v10 + tx * ty * v11
}

/// Overlap length of two axis-aligned intervals.
#[cfg(test)]
#[must_use]
fn interval_overlap(a0: f64, a1: f64, b0: f64, b1: f64) -> f64 {
    (a1.min(b1) - a0.max(b0)).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pupil::circular_pupil_spec;

    const UNIT_N_PUPIL: i64 = 32;
    const UNIT_N_FFT: i64 = 64;
    const UNIT_STAMP: usize = 15;
    const CONTRACT_N_PUPIL: i64 = 256;
    const CONTRACT_N_FFT: i64 = 1024;

    fn unit_pupil() -> PupilSpec {
        circular_pupil_spec(UNIT_N_PUPIL, UNIT_N_FFT)
    }

    fn contract_pupil() -> PupilSpec {
        circular_pupil_spec(CONTRACT_N_PUPIL, CONTRACT_N_FFT)
    }

    fn relative_error(got: f64, expected: f64) -> f64 {
        (got - expected).abs() / expected.abs().max(1e-30)
    }

    #[test]
    fn c10_1_oversampling_is_finite_and_in_range() {
        let r =
            oversampling_factor(&ImageMeta::c10_1_standard_camera(), &contract_pupil()).unwrap();
        assert!(r.is_finite());
        assert!(r > OVERSAMPLING_MIN && r < OVERSAMPLING_MAX);
        // λ = 550 nm, D = 0.20 m, 0.50 arcsec/px, N_f/N_p = 4.
        let expected = ImageMeta::c10_1_standard_camera().pixel_scale_rad() * 0.20 / 550e-9 * 4.0;
        assert!((r - expected).abs() < 1e-15);
    }

    #[test]
    fn insane_undersampling_is_rejected() {
        let mut meta = ImageMeta::c10_1_standard_camera();
        meta.pixel_scale_arcsec = 0.01;
        let err = oversampling_factor(&meta, &unit_pupil()).unwrap_err();
        assert_eq!(err.code, crate::error::ErrorCode::Input);
        assert!(err.message.contains("sampling insane"));
    }

    #[test]
    fn insane_oversampling_is_rejected() {
        let mut meta = ImageMeta::c10_1_standard_camera();
        meta.pixel_scale_arcsec = 20.0;
        let err = oversampling_factor(&meta, &contract_pupil()).unwrap_err();
        assert_eq!(err.code, crate::error::ErrorCode::Input);
        assert!(err.message.contains("sampling insane"));
    }

    #[test]
    fn zero_fourier_shift_is_identity() {
        let n_fft = 16_usize;
        let mut fft = Fft2D::new(n_fft).unwrap();
        let mut intensity = Array2::<f64>::zeros((n_fft, n_fft));
        intensity[[n_fft / 2, n_fft / 2]] = 0.7;
        intensity[[n_fft / 2, n_fft / 2 + 1]] = 0.2;
        intensity[[n_fft / 2 + 1, n_fft / 2]] = 0.1;
        let shifted = fourier_shift(&intensity, 0.0, 0.0, &mut fft).unwrap();
        let mut max_diff = 0.0_f64;
        for (a, b) in intensity.iter().zip(shifted.iter()) {
            max_diff = max_diff.max((a - b).abs());
        }
        assert!(max_diff < 1e-14, "zero-shift residual {max_diff}");
    }

    #[test]
    fn fourier_shift_moves_a_delta_by_integer_fft_pixels() {
        let n_fft = 16_usize;
        let mut fft = Fft2D::new(n_fft).unwrap();
        let mut intensity = Array2::<f64>::zeros((n_fft, n_fft));
        intensity[[n_fft / 2, n_fft / 2]] = 1.0;
        let shifted = fourier_shift(&intensity, 2.0, -1.0, &mut fft).unwrap();
        let mut peak = f64::NEG_INFINITY;
        let mut at = (0, 0);
        for row in 0..n_fft {
            for column in 0..n_fft {
                if shifted[[row, column]] > peak {
                    peak = shifted[[row, column]];
                    at = (row, column);
                }
            }
        }
        assert_eq!(at, (n_fft / 2 - 1, n_fft / 2 + 2));
        assert!((peak - 1.0).abs() < 1e-12);
    }

    #[test]
    fn c9_10_3_constant_image_matches_overlap_area() {
        let n_fft = UNIT_N_FFT as usize;
        let r = oversampling_factor(&ImageMeta::c10_1_standard_camera(), &unit_pupil()).unwrap();
        let value = 1.0 / (n_fft as f64).powi(2);
        let intensity = Array2::from_elem((n_fft, n_fft), value);
        let stamp = box_resample(&intensity, UNIT_STAMP, r).unwrap();
        let c_star = stamp_center(UNIT_STAMP);
        let n = n_fft as f64;
        for row in 0..UNIT_STAMP {
            for column in 0..UNIT_STAMP {
                let (p_ctr, q_ctr) = detector_pixel_fft_center(row, column, n_fft, c_star, r);
                let p0 = p_ctr - r / 2.0;
                let p1 = p_ctr + r / 2.0;
                let q0 = q_ctr - r / 2.0;
                let q1 = q_ctr + r / 2.0;
                let overlap = interval_overlap(p0, p1, 0.0, n) * interval_overlap(q0, q1, 0.0, n);
                // I = 1/N_f², so N_f² times the box integral is the overlap area. (C9.10.3)
                let recovered_area = stamp[[row, column]] / value;
                if overlap == 0.0 {
                    assert_eq!(stamp[[row, column]], 0.0);
                } else {
                    assert!(
                        relative_error(recovered_area, overlap) < 1e-12,
                        "pixel [{row},{column}] recovered overlap {recovered_area}, analytic {overlap}"
                    );
                }
            }
        }
    }

    #[test]
    fn c9_10_3_constant_out_of_bounds_pixels_are_zero() {
        // r large enough that corner stamp pixels miss [0, N_f)² entirely.
        let n_fft = 16_usize;
        let r = 8.0;
        let stamp_size = 15_usize;
        let value = 1.0 / (n_fft as f64).powi(2);
        let intensity = Array2::from_elem((n_fft, n_fft), value);
        let stamp = box_resample(&intensity, stamp_size, r).unwrap();
        assert_eq!(stamp[[0, 0]], 0.0);
        assert_eq!(stamp[[0, stamp_size - 1]], 0.0);
        assert_eq!(stamp[[stamp_size - 1, 0]], 0.0);
        assert_eq!(stamp[[stamp_size - 1, stamp_size - 1]], 0.0);
        let center = stamp_center(stamp_size) as usize;
        assert!(stamp[[center, center]] > 0.0);
    }

    #[test]
    fn c9_10_3_rectangle_matches_analytic_overlap() {
        let n_fft = UNIT_N_FFT as usize;
        let r = oversampling_factor(&ImageMeta::c10_1_standard_camera(), &unit_pupil()).unwrap();
        // Ones on a block of FFT pixels. The bilinear interpolant is identically 1
        // on the rectangle connecting those pixel centers. (C9.10.3)
        let p_lo = 20_usize;
        let p_hi = 44_usize;
        let q_lo = 20_usize;
        let q_hi = 44_usize;
        let mut intensity = Array2::<f64>::zeros((n_fft, n_fft));
        for row in p_lo..=p_hi {
            for column in q_lo..=q_hi {
                intensity[[row, column]] = 1.0;
            }
        }
        let stamp = box_resample(&intensity, UNIT_STAMP, r).unwrap();
        let c_star = stamp_center(UNIT_STAMP);
        // Flat-1 region is the pixel-center rectangle; stay 1 FFT pixel inside so
        // the 4×4 nodes do not sample the bilinear ramp at the block edge.
        let flat_p0 = p_lo as f64 + 1.0;
        let flat_p1 = p_hi as f64 - 1.0;
        let flat_q0 = q_lo as f64 + 1.0;
        let flat_q1 = q_hi as f64 - 1.0;
        let mut n_interior = 0_usize;
        for row in 0..UNIT_STAMP {
            for column in 0..UNIT_STAMP {
                let (p_ctr, q_ctr) = detector_pixel_fft_center(row, column, n_fft, c_star, r);
                let p0 = p_ctr - r / 2.0;
                let p1 = p_ctr + r / 2.0;
                let q0 = q_ctr - r / 2.0;
                let q1 = q_ctr + r / 2.0;
                let fully_inside = p0 >= flat_p0 && p1 <= flat_p1 && q0 >= flat_q0 && q1 <= flat_q1;
                let fully_outside = p1 <= p_lo as f64 - 1.0
                    || p0 >= p_hi as f64 + 1.0
                    || q1 <= q_lo as f64 - 1.0
                    || q0 >= q_hi as f64 + 1.0;
                if fully_inside {
                    n_interior += 1;
                    let overlap = r * r;
                    // Contract quotes the mean (overlap / r²); the stamp stores the
                    // integral, so the mean is stamp / r². (C9.10.2–3)
                    let mean = stamp[[row, column]] / (r * r);
                    assert!(
                        (mean - overlap / (r * r)).abs() < 1e-8,
                        "interior [{row},{column}] mean {mean}, expected 1"
                    );
                } else if fully_outside {
                    assert!(
                        stamp[[row, column]].abs() < 1e-8,
                        "exterior [{row},{column}] = {}",
                        stamp[[row, column]]
                    );
                }
            }
        }
        assert!(
            n_interior > 0,
            "rectangle must cover at least one stamp pixel"
        );
    }

    #[test]
    fn gauss_legendre_weights_sum_to_two() {
        let nodes = gauss_legendre_4();
        let sum: f64 = nodes.iter().map(|(_, w)| w).sum();
        assert!((sum - 2.0).abs() < 1e-15);
    }

    #[test]
    fn bilinear_at_integer_lattice_returns_the_pixel() {
        let mut image = Array2::<f64>::zeros((8, 8));
        image[[3, 5]] = 1.25;
        assert!((bilinear_sample(&image, 3.0, 5.0) - 1.25).abs() < 1e-15);
        assert_eq!(bilinear_sample(&image, -0.1, 4.0), 0.0);
        assert_eq!(bilinear_sample(&image, 4.0, 8.0), 0.0);
    }

    #[test]
    fn flux_and_sky_are_affine() {
        let m = Array2::from_elem((3, 3), 0.5);
        let model = apply_flux_and_sky(&m, 4.0, 0.25).unwrap();
        for value in model.iter() {
            assert!((value - 2.25).abs() < 1e-15);
        }
        let err = apply_flux_and_sky(&m, f64::NAN, 0.0).unwrap_err();
        assert_eq!(err.code, crate::error::ErrorCode::Numerics);
    }

    #[test]
    fn placement_uses_stamp_center_not_measured_centroid() {
        let n_fft = 64_usize;
        let r = 2.0;
        let c_star = stamp_center(15);
        let (p, q) = detector_pixel_fft_center(7, 7, n_fft, c_star, r);
        assert!((p - 32.0).abs() < 1e-15);
        assert!((q - 32.0).abs() < 1e-15);
    }
}
