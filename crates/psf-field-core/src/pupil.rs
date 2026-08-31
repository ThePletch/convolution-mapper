//! Pupil ξ, η grid (C9.2) and circular mask (C9.4).

use crate::types::{PupilSpec, SCHEMA_VERSION};

/// ξ of column `q` on an `n_pupil` × `n_pupil` grid (C9.2).
#[must_use]
pub fn xi_at_column(column: usize, n_pupil: usize) -> f64 {
    let n = n_pupil as f64;
    (column as f64 - (n - 1.0) / 2.0) / (n / 2.0)
}

/// η of row `p` on an `n_pupil` × `n_pupil` grid (C9.2).
#[must_use]
pub fn eta_at_row(row: usize, n_pupil: usize) -> f64 {
    let n = n_pupil as f64;
    (row as f64 - (n - 1.0) / 2.0) / (n / 2.0)
}

/// Polar `(ρ, θ)` at pixel `[row, column]` (NOR.10).
#[must_use]
pub fn rho_theta(row: usize, column: usize, n_pupil: usize) -> (f64, f64) {
    let xi = xi_at_column(column, n_pupil);
    let eta = eta_at_row(row, n_pupil);
    (xi.hypot(eta), eta.atan2(xi))
}

/// Circular unobstructed mask: `1` iff `ρ ≤ 1`, boundary included (C9.4).
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

/// `PupilSpec` whose mask is [`circular_mask`]. Does not run ingest (so unit
/// tests may use `n_pupil` outside `{128, 256, 512}`).
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
}
