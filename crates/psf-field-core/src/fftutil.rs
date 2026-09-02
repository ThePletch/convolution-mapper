//! Centered pad, fftshift, and unnormalized 2-D DFT via rustfft. (C9.3, C9.6)

use std::sync::Arc;

use ndarray::Array2;
use rustfft::num_complex::Complex;
use rustfft::{Fft, FftPlanner};

use crate::error::{ErrorModule, PsfFieldError};

/// Planned 2-D DFT of size `n_fft × n_fft` with reusable scratch.
/// rustfft's 1-D unnormalized transform is applied along rows, then columns,
/// which is the separable 2-D DFT in C9.6. Inverse is the matching unnormalized
/// IDFT; a forward–inverse pair scales the array by N_f².
pub struct Fft2D {
    n_fft: usize,
    forward: Arc<dyn Fft<f64>>,
    inverse: Arc<dyn Fft<f64>>,
    axis_samples: Vec<Complex<f64>>,
    scratch: Vec<Complex<f64>>,
}

impl Fft2D {
    /// Plan even-length forward and inverse transforms. `n_fft` is N_f, the
    /// padded array side. (C9.3)
    pub fn new(n_fft: usize) -> Result<Self, PsfFieldError> {
        if n_fft < 2 || n_fft % 2 != 0 {
            return Err(PsfFieldError::input(
                ErrorModule::Pipeline,
                format!("n_fft {n_fft} must be even and at least 2"),
            ));
        }
        let mut planner = FftPlanner::<f64>::new();
        let forward = planner.plan_fft_forward(n_fft);
        let inverse = planner.plan_fft_inverse(n_fft);
        let scratch_len = forward
            .get_inplace_scratch_len()
            .max(inverse.get_inplace_scratch_len());
        Ok(Self {
            n_fft,
            forward,
            inverse,
            axis_samples: vec![Complex::new(0.0, 0.0); n_fft],
            scratch: vec![Complex::new(0.0, 0.0); scratch_len],
        })
    }

    #[must_use]
    pub fn n_fft(&self) -> usize {
        self.n_fft
    }

    /// Unnormalized forward DFT
    /// U_kl = Σ_{p',q'} P_{p'q'} exp(−2πi (k p' + l q') / N_f).
    /// After this call, DC (zero frequency) is at index (0, 0). (C9.6)
    pub fn forward_dft(&mut self, buffer: &mut Array2<Complex<f64>>) -> Result<(), PsfFieldError> {
        let plan = Arc::clone(&self.forward);
        self.transform_axes(buffer, &plan)
    }

    /// Unnormalized inverse DFT, the adjoint of [`Self::forward_dft`].
    /// Forward then inverse multiplies every entry by N_f²; the Fourier shift
    /// divides that factor back out so a zero shift is the identity. (C9.10.1)
    pub fn inverse_dft(&mut self, buffer: &mut Array2<Complex<f64>>) -> Result<(), PsfFieldError> {
        let plan = Arc::clone(&self.inverse);
        self.transform_axes(buffer, &plan)
    }

    fn transform_axes(
        &mut self,
        buffer: &mut Array2<Complex<f64>>,
        plan: &Arc<dyn Fft<f64>>,
    ) -> Result<(), PsfFieldError> {
        let n_fft = self.n_fft;
        if buffer.nrows() != n_fft || buffer.ncols() != n_fft {
            return Err(PsfFieldError::input(
                ErrorModule::Pipeline,
                "DFT buffer shape must be (n_fft, n_fft)",
            ));
        }
        // Rows: for each p', transform over q' (column index, pupil +x).
        for row in 0..n_fft {
            for column in 0..n_fft {
                self.axis_samples[column] = buffer[[row, column]];
            }
            plan.process_with_scratch(&mut self.axis_samples, &mut self.scratch);
            for column in 0..n_fft {
                buffer[[row, column]] = self.axis_samples[column];
            }
        }
        // Columns: for each l, transform over p' (row index, pupil +y).
        for column in 0..n_fft {
            for row in 0..n_fft {
                self.axis_samples[row] = buffer[[row, column]];
            }
            plan.process_with_scratch(&mut self.axis_samples, &mut self.scratch);
            for row in 0..n_fft {
                buffer[[row, column]] = self.axis_samples[row];
            }
        }
        Ok(())
    }
}

/// DFT frequency bin at index `k` for even length `n`: 0, 1, …, n/2−1, −n/2, …, −1.
/// This is numpy `fftfreq(n) * n`, the integer convention used by the Fourier-shift
/// phase. (C9.10.1)
#[must_use]
pub fn fftfreq_bin(k: usize, n: usize) -> f64 {
    let half = n / 2;
    if k < half {
        k as f64
    } else {
        k as f64 - n as f64
    }
}

/// Embed an N_p × N_p array in the center of an N_f × N_f zero array.
/// Pupil pixel `[p, q]` maps to `[p + (N_f − N_p)/2, q + (N_f − N_p)/2]`.
/// N_f − N_p is even. (C9.3)
pub fn pad_centered<T: Copy>(
    source: &Array2<T>,
    n_fft: usize,
    zero: T,
) -> Result<Array2<T>, PsfFieldError> {
    let n_pupil = source.nrows();
    if source.ncols() != n_pupil {
        return Err(PsfFieldError::input(
            ErrorModule::Pipeline,
            "padded source must be square",
        ));
    }
    if n_fft < n_pupil {
        return Err(PsfFieldError::input(
            ErrorModule::Pipeline,
            "n_fft must be at least n_pupil",
        ));
    }
    if (n_fft - n_pupil) % 2 != 0 {
        return Err(PsfFieldError::input(
            ErrorModule::Pipeline,
            "n_fft − n_pupil must be even so the pupil is centered",
        ));
    }
    let offset = (n_fft - n_pupil) / 2;
    let mut padded = Array2::from_elem((n_fft, n_fft), zero);
    for row in 0..n_pupil {
        for column in 0..n_pupil {
            padded[[row + offset, column + offset]] = source[[row, column]];
        }
    }
    Ok(padded)
}

/// Origin offset of the centered pupil inside the FFT array: (N_f − N_p)/2.
/// Requires N_f ≥ N_p and N_f − N_p even. (C9.3)
#[must_use]
pub fn pad_origin(n_pupil: usize, n_fft: usize) -> usize {
    debug_assert!(n_fft >= n_pupil && (n_fft - n_pupil) % 2 == 0);
    (n_fft - n_pupil) / 2
}

/// Swap halves so that index 0 (DC after the DFT) moves to (N/2, N/2).
/// `new[row, col] = old[(row + N/2) mod N, (col + N/2) mod N]`. (C9.6)
#[must_use]
pub fn fftshift<T: Copy>(array: &Array2<T>) -> Array2<T> {
    rotate_halves(array, array.nrows() / 2, array.ncols() / 2)
}

/// Inverse of [`fftshift`]. For even N the two coincide; kernel convolution
/// uses this half-swap so the spatial origin sits at array index (0, 0). (C9.9)
#[must_use]
pub fn ifftshift<T: Copy>(array: &Array2<T>) -> Array2<T> {
    rotate_halves(array, array.nrows().div_ceil(2), array.ncols().div_ceil(2))
}

fn rotate_halves<T: Copy>(array: &Array2<T>, row_shift: usize, column_shift: usize) -> Array2<T> {
    let n_row = array.nrows();
    let n_col = array.ncols();
    let mut shifted = array.clone();
    for row in 0..n_row {
        for column in 0..n_col {
            let source_row = (row + row_shift) % n_row;
            let source_column = (column + column_shift) % n_col;
            shifted[[row, column]] = array[[source_row, source_column]];
        }
    }
    shifted
}

/// I_kl = |U_kl|² and E = Σ I. E is the unnormalized intensity sum; unit-sum
/// divides by it. (C9.6)
#[must_use]
pub fn intensity_and_sum(field: &Array2<Complex<f64>>) -> (Array2<f64>, f64) {
    let mut intensity = Array2::<f64>::zeros(field.raw_dim());
    let mut intensity_sum = 0.0;
    for row in 0..field.nrows() {
        for column in 0..field.ncols() {
            let value = field[[row, column]].norm_sqr();
            intensity[[row, column]] = value;
            intensity_sum += value;
        }
    }
    (intensity, intensity_sum)
}

/// I ← I / E. E ≤ 0 or non-finite is a numerics failure, not a silent clip. (C9.6)
pub fn unit_sum_intensity(
    intensity: Array2<f64>,
    intensity_sum: f64,
) -> Result<Array2<f64>, PsfFieldError> {
    if !intensity_sum.is_finite() || intensity_sum <= 0.0 {
        return Err(PsfFieldError::numerics(
            ErrorModule::Pipeline,
            format!("FFT intensity sum {intensity_sum} is not positive and finite"),
        ));
    }
    Ok(intensity / intensity_sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pad_maps_pupil_pixel_by_half_the_size_difference() {
        let n_pupil = 32_usize;
        let n_fft = 64_usize;
        let mut source = Array2::from_elem((n_pupil, n_pupil), Complex::new(0.0, 0.0));
        source[[3, 5]] = Complex::new(1.25, -0.5);
        source[[0, 0]] = Complex::new(2.0, 0.0);
        source[[n_pupil - 1, n_pupil - 1]] = Complex::new(0.0, 3.0);
        let padded = pad_centered(&source, n_fft, Complex::new(0.0, 0.0)).unwrap();
        let origin = pad_origin(n_pupil, n_fft);
        assert_eq!(origin, 16);
        assert_eq!(padded[[3 + origin, 5 + origin]], Complex::new(1.25, -0.5));
        assert_eq!(padded[[origin, origin]], Complex::new(2.0, 0.0));
        assert_eq!(
            padded[[origin + n_pupil - 1, origin + n_pupil - 1]],
            Complex::new(0.0, 3.0)
        );
        assert_eq!(padded[[0, 0]], Complex::new(0.0, 0.0));
        assert_eq!(padded[[n_fft - 1, n_fft - 1]], Complex::new(0.0, 0.0));
        let mut nonzero = 0_usize;
        for row in 0..n_fft {
            for column in 0..n_fft {
                if padded[[row, column]] != Complex::new(0.0, 0.0) {
                    nonzero += 1;
                }
            }
        }
        assert_eq!(nonzero, 3);
    }

    #[test]
    fn fftshift_moves_dc_to_array_center() {
        let n_fft = 8_usize;
        let mut array = Array2::from_elem((n_fft, n_fft), 0.0_f64);
        array[[0, 0]] = 1.0;
        let shifted = fftshift(&array);
        assert_eq!(shifted[[n_fft / 2, n_fft / 2]], 1.0);
        assert_eq!(shifted[[0, 0]], 0.0);
        let twice = fftshift(&shifted);
        assert_eq!(twice, array);
        assert_eq!(fftshift(&array), ifftshift(&array));
    }

    #[test]
    fn fftfreq_bin_matches_numpy_integer_convention() {
        let n = 8_usize;
        let expected = [0.0, 1.0, 2.0, 3.0, -4.0, -3.0, -2.0, -1.0];
        for (k, &bin) in expected.iter().enumerate() {
            assert_eq!(fftfreq_bin(k, n), bin);
        }
    }

    #[test]
    fn inverse_undoes_forward_up_to_n_fft_squared() {
        let n_fft = 8_usize;
        let mut fft = Fft2D::new(n_fft).unwrap();
        let mut original = Array2::from_elem((n_fft, n_fft), Complex::new(0.0, 0.0));
        original[[1, 2]] = Complex::new(0.5, -0.25);
        original[[4, 4]] = Complex::new(1.0, 0.0);
        original[[7, 0]] = Complex::new(0.0, 0.75);
        let mut buffer = original.clone();
        fft.forward_dft(&mut buffer).unwrap();
        fft.inverse_dft(&mut buffer).unwrap();
        let scale = (n_fft * n_fft) as f64;
        let mut max_diff = 0.0_f64;
        for row in 0..n_fft {
            for column in 0..n_fft {
                let recovered = buffer[[row, column]] / scale;
                max_diff = max_diff.max((recovered - original[[row, column]]).norm());
            }
        }
        assert!(max_diff < 1e-14, "round-trip residual {max_diff}");
    }

    #[test]
    fn pad_rejects_odd_size_difference() {
        let source = Array2::from_elem((4, 4), Complex::new(1.0, 0.0));
        let err = pad_centered(&source, 7, Complex::new(0.0, 0.0)).unwrap_err();
        assert_eq!(err.code, crate::error::ErrorCode::Input);
        assert!(err.message.contains("even"));
    }
}
