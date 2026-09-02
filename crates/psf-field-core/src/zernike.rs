//! Generic `(n, m)` Zernike engine: orthonormal pupil-phase modes parameterized by
//! ANSI/OSA indices, sampled on a discrete mask, without forming a detector PSF. (C2)

use ndarray::Array2;

use crate::error::{ErrorModule, PsfFieldError};
use crate::pupil::{eta_at_row, xi_at_column};
use crate::types::PupilSpec;

/// Highest radial order `n` the engine will evaluate in v1. Catalog ingest also
/// rejects modes above this cap. (C2.2)
pub const MAX_RADIAL_ORDER: u32 = 15;

/// Discrete RMS below this threshold means the mode is identically zero on the
/// supplied mask (would divide by zero when normalizing). (C2.4)
const RMS_ZERO_THRESHOLD: f64 = 1e-15;

/// One term in the pupil phase Φ = 2π Σ a_k Z^{(k)}. Coefficients `a_k` are in
/// waves RMS over the supplied mask. (C2.6)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhaseCoefficient {
    /// ANSI/OSA radial order n.
    pub n: u32,
    /// ANSI/OSA azimuthal index m: cosine for m > 0, sine for m < 0, radial-only for m = 0.
    pub m: i32,
    /// Coefficient a_k in waves RMS (root-mean-square over the pupil mask).
    pub waves_rms: f64,
}

/// Reject indices that are not a Zernike mode under ANSI/OSA: `|m| ≤ n`,
/// `n − |m|` even, and `n` at most the v1 cap. (C2.1)
pub fn validate_n_m(n: u32, m: i32) -> Result<(), PsfFieldError> {
    if n > MAX_RADIAL_ORDER {
        return Err(PsfFieldError::input(
            ErrorModule::Zernike,
            format!("radial order n={n} exceeds frozen maximum {MAX_RADIAL_ORDER}"),
        ));
    }
    let m_abs = m.unsigned_abs();
    if m_abs > n || (n - m_abs) % 2 != 0 {
        return Err(PsfFieldError::input(
            ErrorModule::Zernike,
            format!("({n}, {m}) is not a valid ANSI/OSA index"),
        ));
    }
    Ok(())
}

/// Sequential OSA index j = (n(n+2) + m) / 2, used for ordering and reports only.
/// The engine API is keyed by `(n, m)`, not by this integer. (C2.1)
pub fn osa_index(n: u32, m: i32) -> Result<i32, PsfFieldError> {
    validate_n_m(n, m)?;
    Ok((n as i32 * (n as i32 + 2) + m) / 2)
}

/// Radial polynomial R_n^{|m|}(ρ) on the unit disk. ρ is the dimensionless pupil
/// radius (0 at center, 1 at the rim). (C2.2)
pub fn radial_polynomial(n: u32, m_abs: u32, rho: f64) -> Result<f64, PsfFieldError> {
    if n > MAX_RADIAL_ORDER {
        return Err(PsfFieldError::input(
            ErrorModule::Zernike,
            format!("radial order n={n} exceeds frozen maximum {MAX_RADIAL_ORDER}"),
        ));
    }
    if m_abs > n || (n - m_abs) % 2 != 0 {
        return Err(PsfFieldError::input(
            ErrorModule::Zernike,
            format!("R_{n}^{m_abs} is not a valid radial polynomial"),
        ));
    }
    Ok(radial_polynomial_horner(n, m_abs, rho))
}

/// Horner product form of the radial polynomial: successive coefficient ratios avoid
/// factorial overflow and stay stable near ρ = 0 (unlike a high-to-low ρ^n Horner).
/// `p = (n − |m|)/2` is the number of radial nodes; `q = (n + |m|)/2`.
fn radial_polynomial_horner(n: u32, m_abs: u32, rho: f64) -> f64 {
    let p = (n - m_abs) / 2;
    let q = (n + m_abs) / 2;
    let mut coefficient = 1.0;
    for i in 1..=p {
        coefficient *= f64::from(q + i) / f64::from(i);
    }
    let x = rho * rho;
    let mut acc = coefficient;
    for k in 0..p {
        coefficient *=
            -(f64::from(p - k) * f64::from(q - k)) / (f64::from(k + 1) * f64::from(n - k));
        acc = acc * x + coefficient;
    }
    if m_abs == 0 {
        acc
    } else {
        acc * rho.powi(m_abs as i32)
    }
}

/// Analytic RMS factor N_n^m = √(2(n+1) / (1 + δ_{m,0})) so that Z̃ has unit RMS
/// on the continuous unit disk. The Kronecker δ_{m,0} is 1 for rotationally
/// symmetric modes and 0 otherwise. (C2.3)
pub fn analytic_normalization(n: u32, m: i32) -> Result<f64, PsfFieldError> {
    validate_n_m(n, m)?;
    Ok(analytic_normalization_value(n, m))
}

fn analytic_normalization_value(n: u32, m: i32) -> f64 {
    // δ_{m,0}: rotationally symmetric modes (m = 0) get a smaller N so piston, defocus,
    // spherical, … still have unit RMS on the disk.
    let delta_m0 = if m == 0 { 1.0 } else { 0.0 };
    (2.0 * (f64::from(n) + 1.0) / (1.0 + delta_m0)).sqrt()
}

/// Analytic Z̃_n^m(ρ, θ) on the continuous unit disk, before discrete-mask RMS.
/// θ is the pupil polar angle (atan2(η, ξ)); m > 0 multiplies by cos(mθ),
/// m < 0 by sin(|m|θ). (C2.3, NOR.10)
pub fn analytic_zernike(n: u32, m: i32, rho: f64, theta: f64) -> Result<f64, PsfFieldError> {
    validate_n_m(n, m)?;
    Ok(analytic_zernike_at(n, m, rho, theta))
}

fn analytic_zernike_at(n: u32, m: i32, rho: f64, theta: f64) -> f64 {
    let radial = radial_polynomial_horner(n, m.unsigned_abs(), rho);
    let n_factor = analytic_normalization_value(n, m);
    if m > 0 {
        n_factor * radial * (f64::from(m) * theta).cos()
    } else if m < 0 {
        n_factor * radial * (f64::from(m.unsigned_abs()) * theta).sin()
    } else {
        n_factor * radial
    }
}

/// Discrete RMS of a screen over the pupil mask: √((1/S) Σ M_pq Z_pq²), where S is
/// the number of unmasked pixels. Coefficients in waves RMS are defined with this
/// discrete inner product, not the continuous disk. (C2.4)
pub fn discrete_rms(values: &Array2<f64>, mask: &[Vec<f64>]) -> Result<f64, PsfFieldError> {
    let (n_row, n_col) = values.dim();
    if mask.len() != n_row || mask.iter().any(|row| row.len() != n_col) {
        return Err(PsfFieldError::input(
            ErrorModule::Zernike,
            "values shape must match mask shape",
        ));
    }
    let mut sum_m = 0.0;
    let mut sum_m_z2 = 0.0;
    for p in 0..n_row {
        for q in 0..n_col {
            let m_pq = mask[p][q];
            sum_m += m_pq;
            let z = values[[p, q]];
            sum_m_z2 += m_pq * z * z;
        }
    }
    if sum_m == 0.0 {
        return Err(PsfFieldError::input(
            ErrorModule::Zernike,
            "pupil mask has no unmasked pixels",
        ));
    }
    Ok((sum_m_z2 / sum_m).sqrt())
}

/// Pre-normalization Z̃ sampled on the pupil grid. Pixels with mask 0 or ρ > 1 are
/// zero. This is the quantity compared to closed-form √3 samples before dividing
/// by discrete RMS. (C2.3, C2.4, C2.8.7)
pub fn analytic_basis_screen(
    n: u32,
    m: i32,
    pupil: &PupilSpec,
) -> Result<Array2<f64>, PsfFieldError> {
    let (screen, _, _) = analytic_screen_and_rms(n, m, pupil)?;
    Ok(screen)
}

/// Discrete-RMS-normalized Z so that (1/S) Σ M Z² = 1 on the supplied mask.
/// Runtime always divides by the sampled RMS even when it is close to 1. (C2.4)
pub fn basis_screen(n: u32, m: i32, pupil: &PupilSpec) -> Result<Array2<f64>, PsfFieldError> {
    let (mut screen, rms, _) = analytic_screen_and_rms(n, m, pupil)?;
    screen /= rms;
    Ok(screen)
}

/// Pupil phase Φ = 2π Σ a_k Z^{(k)} in radians. The 2π converts waves RMS into
/// radians of optical path. Piston (n=0, m=0) is evaluated if requested even though
/// it does not change |FT{A e^{iΦ}}|². (C2.6)
pub fn phase_screen(
    terms: &[PhaseCoefficient],
    pupil: &PupilSpec,
) -> Result<Array2<f64>, PsfFieldError> {
    let n_pupil = pupil_grid_size(pupil)?;
    let mut phi = Array2::<f64>::zeros((n_pupil, n_pupil));
    for term in terms {
        let z = basis_screen(term.n, term.m, pupil)?;
        phi += &(z * (2.0 * std::f64::consts::PI * term.waves_rms));
    }
    Ok(phi)
}

/// ∂Φ/∂a_k = 2π Z^{(k)}. The PSF Jacobian consumes this factor; this engine does
/// not form a detector image. (C2.7)
pub fn phase_derivative(n: u32, m: i32, pupil: &PupilSpec) -> Result<Array2<f64>, PsfFieldError> {
    let z = basis_screen(n, m, pupil)?;
    Ok(z * (2.0 * std::f64::consts::PI))
}

/// Gram matrix G_ij = (1/S) Σ M Z^{(i)} Z^{(j)} of discrete-normalized modes.
/// Reported as a diagnostic of mask-induced non-orthogonality; v1 does not
/// apply Gram–Schmidt. (C2.5)
pub fn gram_matrix(modes: &[(u32, i32)], pupil: &PupilSpec) -> Result<Array2<f64>, PsfFieldError> {
    let n_modes = modes.len();
    let mut screens = Vec::with_capacity(n_modes);
    let mut sum_m = 0.0;
    for &(n, m) in modes {
        let screen = basis_screen(n, m, pupil)?;
        if sum_m == 0.0 {
            sum_m = pupil.mask.iter().flatten().sum();
        }
        screens.push(screen);
    }
    let mut gram = Array2::<f64>::zeros((n_modes, n_modes));
    let n_pupil = pupil_grid_size(pupil)?;
    for i in 0..n_modes {
        for j in i..n_modes {
            let mut acc = 0.0;
            for p in 0..n_pupil {
                for q in 0..n_pupil {
                    acc += pupil.mask[p][q] * screens[i][[p, q]] * screens[j][[p, q]];
                }
            }
            let g = acc / sum_m;
            gram[[i, j]] = g;
            gram[[j, i]] = g;
        }
    }
    Ok(gram)
}

fn pupil_grid_size(pupil: &PupilSpec) -> Result<usize, PsfFieldError> {
    let n = pupil.mask.len();
    if n == 0 {
        return Err(PsfFieldError::input(
            ErrorModule::Zernike,
            "pupil mask is empty",
        ));
    }
    for row in &pupil.mask {
        if row.len() != n {
            return Err(PsfFieldError::input(
                ErrorModule::Zernike,
                "pupil mask must be square",
            ));
        }
    }
    if pupil.n_pupil as usize != n {
        return Err(PsfFieldError::input(
            ErrorModule::Zernike,
            "PupilSpec.n_pupil does not match mask shape",
        ));
    }
    Ok(n)
}

/// Sample analytic Z̃ on the pupil grid and accumulate discrete RMS in one pass.
/// `sum_m` is S, the unmasked pixel count.
fn analytic_screen_and_rms(
    n: u32,
    m: i32,
    pupil: &PupilSpec,
) -> Result<(Array2<f64>, f64, f64), PsfFieldError> {
    validate_n_m(n, m)?;
    let n_pupil = pupil_grid_size(pupil)?;
    let mut screen = Array2::<f64>::zeros((n_pupil, n_pupil));
    let mut sum_m = 0.0;
    let mut sum_m_z2 = 0.0;
    for p in 0..n_pupil {
        let eta = eta_at_row(p, n_pupil);
        for q in 0..n_pupil {
            let m_pq = pupil.mask[p][q];
            sum_m += m_pq;
            if m_pq == 0.0 {
                continue;
            }
            let xi = xi_at_column(q, n_pupil);
            let rho = xi.hypot(eta);
            // Analytic Z̃ is defined on the unit disk; outside it the mode is zero
            // even if a non-circular mask were 1 there.
            if rho > 1.0 {
                continue;
            }
            let theta = eta.atan2(xi);
            let z = analytic_zernike_at(n, m, rho, theta);
            screen[[p, q]] = z;
            sum_m_z2 += m_pq * z * z;
        }
    }
    if sum_m == 0.0 {
        return Err(PsfFieldError::input(
            ErrorModule::Zernike,
            "pupil mask has no unmasked pixels",
        ));
    }
    let rms = (sum_m_z2 / sum_m).sqrt();
    if rms < RMS_ZERO_THRESHOLD {
        return Err(PsfFieldError::input(
            ErrorModule::Zernike,
            format!("mode ({n}, {m}) is identically zero on the mask"),
        ));
    }
    Ok((screen, rms, sum_m))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pupil::{circular_pupil_spec, eta_at_row, rho_theta, xi_at_column};
    use ndarray::Array2;

    /// v1 circular-mask grid used by closed-form Zernike checks (not a fast unit grid).
    const CONTRACT_N_PUPIL: i64 = 256;
    /// FFT length paired with `CONTRACT_N_PUPIL` in the forward pipeline; unused here.
    const CONTRACT_N_FFT: i64 = 1024;
    /// Tiny grid for engine unit tests that do not need the contract-sized pupil.
    const UNIT_N_PUPIL: i64 = 32;
    /// Odd stamp side length used by the ignored PSF-level checks (filled in later).
    const STAMP_SIZE: usize = 31;

    fn contract_pupil() -> PupilSpec {
        circular_pupil_spec(CONTRACT_N_PUPIL, CONTRACT_N_FFT)
    }

    fn unit_pupil() -> PupilSpec {
        circular_pupil_spec(UNIT_N_PUPIL, UNIT_N_PUPIL * 2)
    }

    fn valid_modes_through(max_n: u32) -> Vec<(u32, i32)> {
        let mut modes = Vec::new();
        for n in 0..=max_n {
            let n_i = n as i32;
            for m in -n_i..=n_i {
                if (n_i - m.abs()) % 2 == 0 {
                    modes.push((n, m));
                }
            }
        }
        modes
    }

    fn factorial_f64(k: u32) -> f64 {
        (1..=k).map(f64::from).product()
    }

    fn radial_identity(n: u32, m_abs: u32, rho: f64) -> f64 {
        let k_max = (n - m_abs) / 2;
        let mut sum = 0.0;
        for k in 0..=k_max {
            let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
            let num = factorial_f64(n - k);
            let den = factorial_f64(k)
                * factorial_f64((n + m_abs) / 2 - k)
                * factorial_f64((n - m_abs) / 2 - k);
            sum += sign * num / den * rho.powi((n - 2 * k) as i32);
        }
        sum
    }

    fn within_one_ulp(got: f64, expected: f64) -> bool {
        if got == expected {
            return true;
        }
        let map = |x: f64| -> i64 {
            let bits = x.to_bits() as i64;
            if bits < 0 {
                i64::MIN - bits
            } else {
                bits
            }
        };
        map(got).abs_diff(map(expected)) <= 1
    }

    fn temporary_stamp_stub() -> Vec<Vec<f64>> {
        vec![vec![0.0; STAMP_SIZE]; STAMP_SIZE]
    }

    fn temporary_jacobian_column_stub() -> Vec<f64> {
        vec![0.0; STAMP_SIZE * STAMP_SIZE]
    }

    fn stamp_max(intensity: &[Vec<f64>]) -> f64 {
        intensity.iter().flatten().copied().fold(0.0_f64, f64::max)
    }

    fn first_moment(intensity: &[Vec<f64>]) -> (f64, f64) {
        let mut wx = 0.0;
        let mut wy = 0.0;
        let mut w = 0.0;
        for (j, row) in intensity.iter().enumerate() {
            for (i, &v) in row.iter().enumerate() {
                wx += i as f64 * v;
                wy += j as f64 * v;
                w += v;
            }
        }
        (wx / w, wy / w)
    }

    fn unique_peak(intensity: &[Vec<f64>]) -> Option<(usize, usize)> {
        let mut best = f64::NEG_INFINITY;
        let mut loc = (0, 0);
        let mut count = 0_usize;
        for (j, row) in intensity.iter().enumerate() {
            for (i, &v) in row.iter().enumerate() {
                if v > best {
                    best = v;
                    loc = (j, i);
                    count = 1;
                } else if v == best {
                    count += 1;
                }
            }
        }
        if count == 1 {
            Some(loc)
        } else {
            None
        }
    }

    fn max_azimuthal_relative_rms(intensity: &[Vec<f64>], c_star: f64) -> f64 {
        // 0.25 px hypot annuli mix axis (r) and diagonal (r√2) pixels, so RMS vs the
        // annulus mean is the Airy radial slope. Compare each pixel to its 90° rotate
        // instead: a circular pattern agrees; m=1 / fftshift leakage does not. (C2.8.2)
        const ANNULUS_WIDTH_PX: f64 = 0.25;
        let n = intensity.len();
        let sample = |column: f64, row: f64| -> Option<f64> {
            if column < -0.5 || row < -0.5 {
                return None;
            }
            let i = column.round() as isize;
            let j = row.round() as isize;
            if i < 0 || j < 0 {
                return None;
            }
            let (i, j) = (i as usize, j as usize);
            if j < n && i < n {
                Some(intensity[j][i])
            } else {
                None
            }
        };
        let mut max_rel = 0.0;
        let mut r = 0.5;
        while r < n as f64 {
            let r_hi = r + ANNULUS_WIDTH_PX;
            let mut n_pix = 0_usize;
            let mut mean = 0.0;
            let mut diffs = Vec::new();
            for (j, row) in intensity.iter().enumerate() {
                for (i, &v) in row.iter().enumerate() {
                    let rp = (i as f64 - c_star).hypot(j as f64 - c_star);
                    if rp < r || rp >= r_hi {
                        continue;
                    }
                    n_pix += 1;
                    mean += v;
                    let di = i as f64 - c_star;
                    let dj = j as f64 - c_star;
                    // 90° rotation in (x, y) = (column, row): (di, dj) → (−dj, di).
                    if let Some(rotated) = sample(c_star - dj, c_star + di) {
                        diffs.push(v - rotated);
                    }
                }
            }
            if n_pix >= 8 && !diffs.is_empty() {
                mean /= n_pix as f64;
                let rms = (diffs.iter().map(|d| d * d).sum::<f64>() / diffs.len() as f64).sqrt();
                let rel = rms / mean.max(1e-15);
                if rel > max_rel {
                    max_rel = rel;
                }
            }
            r = r_hi;
        }
        max_rel
    }

    fn l2_norm(column: &[f64]) -> f64 {
        column.iter().map(|v| v * v).sum::<f64>().sqrt()
    }

    #[test]
    fn c2_8_7_continuous_analytic_defocus_samples() {
        // Defocus Z̃_2^0 = √3 (2ρ² − 1), so the rim (ρ=1) is +√3 and the origin is −√3.
        // Compared to 1 ulp; this does not use a sampled pupil grid. (C2.8.7)
        let expected = 3.0_f64.sqrt();
        let at_rim = analytic_zernike(2, 0, 1.0, 0.0).unwrap();
        let at_origin = analytic_zernike(2, 0, 0.0, 0.0).unwrap();
        assert!(
            within_one_ulp(at_rim, expected),
            "Z̃_2^0(ρ=1, θ=0) = {at_rim}, expected √3 = {expected}"
        );
        assert!(
            within_one_ulp(at_origin, -expected),
            "Z̃_2^0(ρ=0) = {at_origin}, expected −√3 = {}",
            -expected
        );
    }

    #[test]
    fn c2_8_7_sampled_analytic_defocus_pre_normalization() {
        // On N_p = 256 the rim sample is not exactly ρ=1 (max on-axis ρ ≈ 0.996).
        // Do not stretch ξ, η to force a ρ=1 pixel. (C2.8.7)
        let pupil = contract_pupil();
        let screen = analytic_basis_screen(2, 0, &pupil).unwrap();
        let n = CONTRACT_N_PUPIL as usize;
        let expected = 3.0_f64.sqrt();

        let mut nearest_rim = (0, 0);
        let mut nearest_rim_d2 = f64::INFINITY;
        let mut min_rho = f64::INFINITY;
        let mut origin_pixels = Vec::new();
        for p in 0..n {
            for q in 0..n {
                let (rho, _) = rho_theta(p, q, n);
                let xi = xi_at_column(q, n);
                let eta = eta_at_row(p, n);
                let d2 = (xi - 1.0).powi(2) + eta.powi(2);
                if d2 < nearest_rim_d2 {
                    nearest_rim_d2 = d2;
                    nearest_rim = (p, q);
                }
                if rho < min_rho {
                    min_rho = rho;
                    origin_pixels.clear();
                    origin_pixels.push((p, q));
                } else if (rho - min_rho).abs() <= 1e-15 {
                    origin_pixels.push((p, q));
                }
            }
        }

        let z_rim = screen[[nearest_rim.0, nearest_rim.1]];
        assert!(
            (z_rim - expected).abs() < 0.03,
            "sampled Z̃ at pixel nearest (ρ=1, θ=0) = {z_rim}, |Z̃−√3| must be < 0.03"
        );
        for (p, q) in origin_pixels {
            let z0 = screen[[p, q]];
            assert!(
                (z0 + expected).abs() < 0.02,
                "sampled Z̃ at pixel nearest ρ=0 [{p},{q}] = {z0}, |Z̃+√3| must be < 0.02"
            );
        }
    }

    #[test]
    fn invalid_n_m_are_rejected() {
        assert!(validate_n_m(1, 0).is_err());
        assert!(validate_n_m(2, 3).is_err());
        assert!(validate_n_m(16, 0).is_err());
        assert!(analytic_zernike(1, 0, 0.5, 0.0).is_err());
        assert!(basis_screen(2, 3, &unit_pupil()).is_err());
        assert!(osa_index(16, 0).is_err());
    }

    fn assert_gram_diagonal_one(gram: &Array2<f64>, modes: &[(u32, i32)]) {
        for i in 0..modes.len() {
            assert!(
                (gram[[i, i]] - 1.0).abs() < 1e-12,
                "G[{i},{i}] = {} is not 1",
                gram[[i, i]]
            );
        }
    }

    #[test]
    fn gram_matrix_smoke_on_tiny_grid() {
        let modes = valid_modes_through(6);
        let gram = gram_matrix(&modes, &unit_pupil()).unwrap();
        assert_gram_diagonal_one(&gram, &modes);
        for i in 0..modes.len() {
            for j in 0..modes.len() {
                if i != j {
                    assert!(gram[[i, j]].is_finite());
                }
            }
        }
    }

    #[test]
    fn gram_matrix_circular_mask_n_le_6() {
        // Off-diagonal inner products on a circular mask should stay small; larger
        // values on a non-circular mask are expected and are reported, not treated as
        // a bug. v1 does not Gram–Schmidt. (C2.5)
        let modes = valid_modes_through(6);
        let gram = gram_matrix(&modes, &contract_pupil()).unwrap();
        assert_gram_diagonal_one(&gram, &modes);
        for i in 0..modes.len() {
            for j in 0..modes.len() {
                if i == j {
                    continue;
                }
                assert!(
                    gram[[i, j]].abs() < 0.05,
                    "off-diagonal G[{i},{j}] = {} exceeds 0.05",
                    gram[[i, j]]
                );
            }
        }
    }

    #[test]
    fn discrete_rms_after_normalization_is_one() {
        let pupil = unit_pupil();
        for &(n, m) in &valid_modes_through(4) {
            let z = basis_screen(n, m, &pupil).unwrap();
            let rms = discrete_rms(&z, &pupil.mask).unwrap();
            assert!(
                (rms - 1.0).abs() < 1e-12,
                "discrete RMS of Z_{n}^{m} = {rms}, expected 1"
            );
        }
    }

    #[test]
    fn pre_normalization_rms_near_one_on_contract_grid() {
        // On the v1 circular mask at N_p ≥ 256, sampled RMS of Z̃ should stay within
        // 2% of the continuous-disk value 1 for modes with n ≤ 8. Runtime still divides
        // by the sampled RMS. (C2.4)
        let pupil = contract_pupil();
        for &(n, m) in &valid_modes_through(8) {
            let z_tilde = analytic_basis_screen(n, m, &pupil).unwrap();
            let rms = discrete_rms(&z_tilde, &pupil.mask).unwrap();
            assert!(
                (rms - 1.0).abs() < 0.02,
                "pre-normalization rms of Z̃_{n}^{m} = {rms}, |rms−1| must be < 0.02"
            );
        }
    }

    #[test]
    fn osa_index_matches_formula() {
        assert_eq!(osa_index(0, 0).unwrap(), 0);
        assert_eq!(osa_index(1, -1).unwrap(), 1);
        assert_eq!(osa_index(1, 1).unwrap(), 2);
        assert_eq!(osa_index(2, -2).unwrap(), 3);
        assert_eq!(osa_index(2, 0).unwrap(), 4);
        assert_eq!(osa_index(2, 2).unwrap(), 5);
        assert_eq!(osa_index(3, -3).unwrap(), 6);
        assert_eq!(osa_index(3, 1).unwrap(), 8);
    }

    #[test]
    fn tilt_modes_nonzero_on_disk() {
        let pupil = unit_pupil();
        for m in [1, -1] {
            let z = basis_screen(1, m, &pupil).unwrap();
            let peak = z.iter().fold(0.0_f64, |a, &v| a.max(v.abs()));
            assert!(peak > 0.0, "Z_1^{m} is identically zero on the disk");
        }
    }

    #[test]
    fn radial_polynomial_matches_c2_2_identity() {
        for n in 0..=15 {
            for m_abs in 0..=n {
                if (n - m_abs) % 2 != 0 {
                    continue;
                }
                for &rho in &[0.0, 0.3, 0.7, 1.0] {
                    let got = radial_polynomial(n, m_abs, rho).unwrap();
                    let expected = radial_identity(n, m_abs, rho);
                    assert!(
                        (got - expected).abs() < 1e-10,
                        "R_{n}^{m_abs}({rho}) Horner={got} identity={expected}"
                    );
                }
            }
        }
    }

    #[test]
    fn phase_derivative_is_two_pi_times_basis() {
        let pupil = unit_pupil();
        let z = basis_screen(2, 0, &pupil).unwrap();
        let dphi = phase_derivative(2, 0, &pupil).unwrap();
        let expected = &z * (2.0 * std::f64::consts::PI);
        let max_diff = dphi
            .iter()
            .zip(expected.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(max_diff < 1e-15);
        let phi = phase_screen(
            &[PhaseCoefficient {
                n: 2,
                m: 0,
                waves_rms: 1.0,
            }],
            &pupil,
        )
        .unwrap();
        let max_phi_diff = phi
            .iter()
            .zip(dphi.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(max_phi_diff < 1e-15);
    }

    #[test]
    fn empty_mask_is_input_error() {
        let mut pupil = unit_pupil();
        for row in &mut pupil.mask {
            for v in row {
                *v = 0.0;
            }
        }
        let err = basis_screen(0, 0, &pupil).unwrap_err();
        assert_eq!(err.code, crate::error::ErrorCode::Input);
        assert_eq!(err.module, ErrorModule::Zernike);
    }

    /// Piston is a global phase; intensity |FT{A e^{iΦ}}|² must not change when a_{0,0}
    /// goes from 0 to 1 wave. (C2.8.1)
    #[test]
    #[ignore = "requires C9 PSF pipeline"]
    fn c2_8_1_piston_independence() {
        let intensity_piston = temporary_stamp_stub();
        let intensity_zero = temporary_stamp_stub();
        let max_i0 = stamp_max(&intensity_zero);
        let max_abs_diff = intensity_piston
            .iter()
            .zip(&intensity_zero)
            .flat_map(|(a, b)| a.iter().zip(b).map(|(x, y)| (x - y).abs()))
            .fold(0.0_f64, f64::max);
        assert!(max_abs_diff / max_i0 < 1e-10);
    }

    /// Zero aberrations must yield a centered Airy pattern with no azimuthal m=1 leak. (C2.8.2)
    #[test]
    #[ignore = "requires C9 PSF pipeline"]
    fn c2_8_2_zero_coefficients_airy_centered_azimuthal() {
        let intensity = temporary_stamp_stub();
        let c_star = (STAMP_SIZE as f64 - 1.0) / 2.0;
        let peak = unique_peak(&intensity).expect("peak pixel must be unique");
        assert_eq!(peak, (c_star as usize, c_star as usize));
        assert!(max_azimuthal_relative_rms(&intensity, c_star) < 1e-4);
    }

    /// Defocus is even in a_{2,0}: I(+α) must match I(−α). (C2.8.3)
    #[test]
    #[ignore = "requires C9 PSF pipeline"]
    fn c2_8_3_defocus_evenness() {
        let intensity_pos = temporary_stamp_stub();
        let intensity_neg = temporary_stamp_stub();
        let intensity_zero = temporary_stamp_stub();
        let c_star = (STAMP_SIZE as f64 - 1.0) / 2.0;
        let max_i0 = stamp_max(&intensity_zero);
        let max_abs_diff = intensity_pos
            .iter()
            .zip(&intensity_neg)
            .flat_map(|(a, b)| a.iter().zip(b).map(|(x, y)| (x - y).abs()))
            .fold(0.0_f64, f64::max);
        assert!(max_abs_diff / max_i0 < 1e-8);
        assert!(max_azimuthal_relative_rms(&intensity_pos, c_star) < 1e-4);
        assert!(max_azimuthal_relative_rms(&intensity_neg, c_star) < 1e-4);
    }

    /// Even function of defocus ⇒ Jacobian column for a_{2,0} vanishes at a=0. (C2.8.4)
    #[test]
    #[ignore = "requires C9 PSF pipeline"]
    fn c2_8_4_defocus_at_zero_vanishing_jacobian_column() {
        let jacobian_at_zero = temporary_jacobian_column_stub();
        let jacobian_at_defocus = temporary_jacobian_column_stub();
        let n0 = l2_norm(&jacobian_at_zero);
        let n_defocus = l2_norm(&jacobian_at_defocus);
        assert!(n0 < 1e-8 * n_defocus);
    }

    /// Orthonormal Zernike coma has zero Z-tilt but nonzero G-tilt, so the first
    /// moment follows ⟨∂Φ/∂ξ⟩ along the mode axis. The 5×10^{-3} px bound is the
    /// cross-axis residual, where G-tilt is identically zero. (C2.8.5)
    #[test]
    #[ignore = "requires C9 PSF pipeline"]
    fn c2_8_5_zernike_coma_g_tilt() {
        let c_star = (STAMP_SIZE as f64 - 1.0) / 2.0;
        for m in [1_i32, -1] {
            let intensity = temporary_stamp_stub();
            let (x_bar, y_bar) = first_moment(&intensity);
            if m > 0 {
                assert!((y_bar - c_star).abs() < 5e-3);
            } else {
                assert!((x_bar - c_star).abs() < 5e-3);
            }
        }
    }

    /// Flipping the sign of a_{3,1} must reverse the coma flare across the stamp center. (C2.8.6)
    #[test]
    #[ignore = "requires C9 PSF pipeline"]
    fn c2_8_6_coma_sign_flips_the_flare() {
        let c_star = (STAMP_SIZE as f64 - 1.0) / 2.0;
        let intensity_pos = temporary_stamp_stub();
        let intensity_neg = temporary_stamp_stub();
        let x_pos = intensity_squared_centroid_x(&intensity_pos, 0.05);
        let x_neg = intensity_squared_centroid_x(&intensity_neg, 0.05);
        assert!((x_pos - c_star).abs() >= 0.02);
        assert!((x_neg - c_star).abs() >= 0.02);
        assert!((x_pos - c_star).signum() != (x_neg - c_star).signum());
    }

    fn intensity_squared_centroid_x(intensity: &[Vec<f64>], fraction_of_max: f64) -> f64 {
        let max_i = stamp_max(intensity);
        let threshold = fraction_of_max * max_i;
        let mut wx = 0.0;
        let mut w = 0.0;
        for row in intensity {
            for (i, &v) in row.iter().enumerate() {
                if v > threshold {
                    let i2 = v * v;
                    wx += i as f64 * i2;
                    w += i2;
                }
            }
        }
        wx / w
    }
}
