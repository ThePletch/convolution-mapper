//! PSF field modeler core.

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

pub use catalog::{kernel_parameters, KernelParameter};
pub use error::{ErrorCode, ErrorModule, PsfFieldError};
pub use theta::{
    assemble_layout, assemble_theta, defocus_moment_init, dtheta_du, evaluate_gaussian_prior,
    initialize_theta, moffat_alpha_from_fwhm, prior_residual, schedule_steps, theta_from_unbounded,
    DefocusMomentInputs, EvaluatedGaussianPrior, ScheduleStep, Stage1Options, ThetaAssembly,
    ThetaLayout, MOFFAT_BETA,
};

#[cfg(feature = "python")]
mod python_module;

#[cfg(test)]
mod tests {
    #[test]
    fn crate_name_matches_extension_module() {
        assert_eq!(env!("CARGO_PKG_NAME"), "psf-field-core");
    }
}
