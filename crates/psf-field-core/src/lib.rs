//! PSF field modeler core.
//!
//! Numeric modules are placeholders until later PRs. Types, ingest, and
//! errors land in PR-B.

pub mod catalog;
pub mod diag;
pub mod error;
pub mod eval;
pub mod fftutil;
pub mod forward;
pub mod jacobian;
pub mod kernels;
pub mod lm;
pub mod pupil;
pub mod resample;
pub mod stage2;
pub mod theta;
pub mod types;
pub mod zernike;

#[cfg(feature = "python")]
mod python_module;

#[cfg(test)]
mod tests {
    #[test]
    fn crate_name_matches_extension_module() {
        assert_eq!(env!("CARGO_PKG_NAME"), "psf-field-core");
    }
}
