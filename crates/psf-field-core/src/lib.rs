//! PSF field modeler core.
//!
//! Types, ingest, errors, and the generic `(n, m)` Zernike engine live here.
//! Other numeric modules are unimplemented until their contracts are built.

pub mod catalog;
pub mod diag;
pub mod error;
pub mod eval;
pub mod fftutil;
pub mod forward;
pub mod ingest;
pub mod jacobian;
pub mod kernels;
pub mod lm;
pub mod pupil;
pub mod resample;
pub mod stage2;
pub mod theta;
pub mod types;
pub mod zernike;

pub use error::{ErrorCode, ErrorModule, PsfFieldError};

#[cfg(feature = "python")]
mod python_module;

#[cfg(test)]
mod tests {
    #[test]
    fn crate_name_matches_extension_module() {
        assert_eq!(env!("CARGO_PKG_NAME"), "psf-field-core");
    }
}
