//! Flat θ assembly, freeze masks, init, priors, and schedule (C4).

use std::collections::HashSet;

use crate::catalog::kernel_parameters;
use crate::error::{ErrorModule, PsfFieldError};
use crate::types::{
    check_finite, check_positive, Catalog, ErrorTerm, InitMethod, ParamMeta, PriorMean, PriorSpec,
    Scope, StarRecord,
};

pub const MOFFAT_BETA: f64 = 2.5;
const PRIOR_SIGMA_FLOOR: f64 = 1e-3;
const DEFOCUS_MOMENT_FACTOR: f64 = 0.35;
const DEFOCUS_D_MAX: f64 = 2.0;
const DEFOCUS_A_ABS_MAX: f64 = 2.0;
const SIGMA0_SQ_FLOOR: f64 = 1e-12;

fn input(message: impl Into<String>) -> PsfFieldError {
    PsfFieldError::input(ErrorModule::Boundary, message)
}

/// Stage-1 options that C3/C4 require before Levenberg–Marquardt.
#[derive(Debug, Clone, PartialEq)]
pub struct Stage1Options {
    pub freeze_mask: Option<Vec<bool>>,
    pub use_schedule: bool,
    pub expected_fwhm_px: f64,
}

impl Stage1Options {
    #[must_use]
    pub fn new(expected_fwhm_px: f64) -> Self {
        Self {
            freeze_mask: None,
            use_schedule: true,
            expected_fwhm_px,
        }
    }
}

/// Precomputed second-moment inputs for C3.5.1 (no forward PSF).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DefocusMomentInputs {
    pub sigma_meas_sq: f64,
    pub sigma0_sq: f64,
    pub known_defocus_waves: f64,
}

/// Evaluated Gaussian prior (μ, σ) after C3.5 init.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EvaluatedGaussianPrior {
    pub mu: f64,
    pub sigma: f64,
}

/// Enabled-term θ layout: sidecar, bounds, and free indices into the full vector.
#[derive(Debug, Clone, PartialEq)]
pub struct ThetaLayout {
    pub param_meta: Vec<ParamMeta>,
    pub bounds: Vec<Option<[f64; 2]>>,
    pub free_index: Vec<usize>,
}

impl ThetaLayout {
    #[must_use]
    pub fn n_all(&self) -> usize {
        self.param_meta.len()
    }

    #[must_use]
    pub fn n_free(&self) -> usize {
        self.free_index.len()
    }
}

/// Initialized θ together with layout, bounds, and per-slot Gaussian priors.
#[derive(Debug, Clone, PartialEq)]
pub struct ThetaAssembly {
    pub theta: Vec<f64>,
    pub param_meta: Vec<ParamMeta>,
    pub bounds: Vec<Option<[f64; 2]>>,
    pub free_index: Vec<usize>,
    pub priors: Vec<Option<EvaluatedGaussianPrior>>,
}

impl ThetaAssembly {
    #[must_use]
    pub fn n_all(&self) -> usize {
        self.param_meta.len()
    }

    #[must_use]
    pub fn n_free(&self) -> usize {
        self.free_index.len()
    }
}

/// One freeze/unfreeze step with free indices into the full θ vector.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScheduleStep {
    pub name: String,
    pub unfrozen_term_ids: Vec<String>,
    pub free_index: Vec<usize>,
}

struct Slot {
    term_id: String,
    role: String,
    scope: Scope,
    frozen: bool,
    unit: String,
    bounds: Option<[f64; 2]>,
}

fn slots_for_catalog(catalog: &Catalog) -> Vec<Slot> {
    let mut slots = Vec::new();
    for term in &catalog.terms {
        if !term.enabled() {
            continue;
        }
        match term {
            ErrorTerm::Phase { .. } => slots.push(Slot {
                term_id: term.term_id().to_string(),
                role: "local_value".to_string(),
                scope: term.scope(),
                frozen: term.frozen(),
                unit: term.units().to_string(),
                bounds: term.bounds(),
            }),
            ErrorTerm::Photometric { .. } => slots.push(Slot {
                term_id: term.term_id().to_string(),
                role: term.term_id().to_string(),
                scope: term.scope(),
                frozen: term.frozen(),
                unit: term.units().to_string(),
                bounds: term.bounds(),
            }),
            ErrorTerm::Kernel { kernel, .. } => {
                for (k, param) in kernel_parameters(kernel.id).iter().enumerate() {
                    slots.push(Slot {
                        term_id: term.term_id().to_string(),
                        role: param.role.to_string(),
                        scope: term.scope(),
                        frozen: term.frozen() || param.always_frozen,
                        unit: param.unit.to_string(),
                        bounds: if k == 0 { term.bounds() } else { None },
                    });
                }
            }
        }
    }
    slots
}

/// Assemble `param_meta` / `free_index` from enabled catalog terms (C4.1, C4.4).
pub fn assemble_layout(
    catalog: &Catalog,
    options: &Stage1Options,
) -> Result<ThetaLayout, PsfFieldError> {
    check_positive("expected_fwhm_px", options.expected_fwhm_px)?;
    let slots = slots_for_catalog(catalog);
    let mut param_meta: Vec<ParamMeta> = slots
        .iter()
        .map(|slot| ParamMeta {
            term_id: slot.term_id.clone(),
            role: slot.role.clone(),
            scope: slot.scope,
            frozen: slot.frozen,
            unit: slot.unit.clone(),
        })
        .collect();
    let bounds: Vec<Option<[f64; 2]>> = slots.iter().map(|slot| slot.bounds).collect();

    let catalog_free: Vec<usize> = param_meta
        .iter()
        .enumerate()
        .filter(|(_, meta)| !meta.frozen)
        .map(|(i, _)| i)
        .collect();
    let n_free = catalog_free.len();
    let free_index = match &options.freeze_mask {
        None => catalog_free,
        Some(mask) => {
            if mask.len() != n_free {
                return Err(input(format!(
                    "freeze_mask length {} does not equal n_free {n_free}",
                    mask.len()
                )));
            }
            let mut kept = Vec::new();
            for (&idx, &freeze) in catalog_free.iter().zip(mask.iter()) {
                if freeze {
                    param_meta[idx].frozen = true;
                } else {
                    kept.push(idx);
                }
            }
            kept
        }
    };
    Ok(ThetaLayout {
        param_meta,
        bounds,
        free_index,
    })
}

/// Fill θ and Gaussian priors for an assembled layout (C3.4, C3.5).
pub fn initialize_theta(
    catalog: &Catalog,
    layout: &ThetaLayout,
    options: &Stage1Options,
    star: Option<&StarRecord>,
    defocus: Option<DefocusMomentInputs>,
) -> Result<(Vec<f64>, Vec<Option<EvaluatedGaussianPrior>>), PsfFieldError> {
    check_positive("expected_fwhm_px", options.expected_fwhm_px)?;
    let mut theta = Vec::with_capacity(layout.n_all());
    let mut priors = Vec::with_capacity(layout.n_all());
    let mut slot_i = 0usize;
    for term in &catalog.terms {
        if !term.enabled() {
            continue;
        }
        let n_slots = match term {
            ErrorTerm::Kernel { kernel, .. } => kernel_parameters(kernel.id).len(),
            _ => 1,
        };
        for k in 0..n_slots {
            let meta = &layout.param_meta[slot_i];
            let value = init_slot_value(term, &meta.role, k, options, star, defocus)?;
            let prior = if k == 0 {
                evaluate_gaussian_prior(term.prior(), value)?
            } else {
                None
            };
            theta.push(value);
            priors.push(prior);
            slot_i += 1;
        }
    }
    if slot_i != layout.n_all() {
        return Err(PsfFieldError::internal(
            ErrorModule::Boundary,
            "θ slot count does not match catalog layout",
        ));
    }
    Ok((theta, priors))
}

fn init_slot_value(
    term: &ErrorTerm,
    role: &str,
    slot_in_term: usize,
    options: &Stage1Options,
    star: Option<&StarRecord>,
    defocus: Option<DefocusMomentInputs>,
) -> Result<f64, PsfFieldError> {
    if role == "beta" {
        return Ok(MOFFAT_BETA);
    }
    if slot_in_term > 0 {
        return Ok(0.0);
    }
    match term.init().method {
        InitMethod::Zero => zero_init(term.prior()),
        InitMethod::FluxSum => {
            let Some(star) = star else {
                return Err(input("flux_sum init requires a StarRecord"));
            };
            Ok(star.flux_sum_adu.max(0.0))
        }
        InitMethod::MoffatFwhm => moffat_alpha_from_fwhm(options.expected_fwhm_px),
        InitMethod::DefocusMoment => {
            let Some(defocus) = defocus else {
                return Err(input(
                    "defocus_moment init requires sigma_meas_sq, sigma0_sq, and known_defocus_waves",
                ));
            };
            defocus_moment_init(
                defocus.sigma_meas_sq,
                defocus.sigma0_sq,
                defocus.known_defocus_waves,
            )
        }
    }
}

fn zero_init(prior: &PriorSpec) -> Result<f64, PsfFieldError> {
    match prior {
        PriorSpec::None { .. } => Ok(0.0),
        PriorSpec::Gaussian {
            mean: PriorMean::Number(mu),
            ..
        } => {
            check_finite("prior.mean", *mu)?;
            Ok(*mu)
        }
        PriorSpec::Gaussian {
            mean: PriorMean::Init(_),
            ..
        } => Err(input(
            "zero init cannot use gaussian mean \"init\" (rejected at ingest)",
        )),
    }
}

/// Closed-form C3.5.1 defocus init from second moments, in waves.
pub fn defocus_moment_init(
    sigma_meas_sq: f64,
    sigma0_sq: f64,
    known_defocus_waves: f64,
) -> Result<f64, PsfFieldError> {
    check_finite("sigma_meas_sq", sigma_meas_sq)?;
    check_finite("sigma0_sq", sigma0_sq)?;
    check_finite("known_defocus_waves", known_defocus_waves)?;
    let extra = (sigma_meas_sq - sigma0_sq).max(0.0);
    let denom = sigma0_sq.max(SIGMA0_SQ_FLOOR);
    let d = (DEFOCUS_MOMENT_FACTOR * extra / denom).clamp(0.0, DEFOCUS_D_MAX);
    Ok((d - known_defocus_waves).clamp(-DEFOCUS_A_ABS_MAX, DEFOCUS_A_ABS_MAX))
}

/// Moffat α from FWHM with frozen β = 2.5 (C3.5.2).
pub fn moffat_alpha_from_fwhm(expected_fwhm_px: f64) -> Result<f64, PsfFieldError> {
    check_positive("expected_fwhm_px", expected_fwhm_px)?;
    let half_width_factor = (2.0_f64.powf(1.0 / MOFFAT_BETA) - 1.0).sqrt();
    Ok(expected_fwhm_px / (2.0 * half_width_factor))
}

/// Evaluate a term prior after init. `family: none` yields `None`.
pub fn evaluate_gaussian_prior(
    prior: &PriorSpec,
    a0: f64,
) -> Result<Option<EvaluatedGaussianPrior>, PsfFieldError> {
    match prior {
        PriorSpec::None { .. } => Ok(None),
        PriorSpec::Gaussian {
            mean: PriorMean::Number(mu),
            sigma,
            ..
        } => {
            let sigma = sigma.ok_or_else(|| input("gaussian prior is missing sigma"))?;
            check_finite("prior.mean", *mu)?;
            check_positive("prior.sigma", sigma)?;
            Ok(Some(EvaluatedGaussianPrior { mu: *mu, sigma }))
        }
        PriorSpec::Gaussian {
            mean: PriorMean::Init(_),
            sigma_rel,
            ..
        } => {
            let sigma_rel =
                sigma_rel.ok_or_else(|| input("gaussian prior is missing sigma_rel"))?;
            check_positive("prior.sigma_rel", sigma_rel)?;
            check_finite("a0", a0)?;
            Ok(Some(EvaluatedGaussianPrior {
                mu: a0,
                sigma: (sigma_rel * a0.abs()).max(PRIOR_SIGMA_FLOOR),
            }))
        }
    }
}

#[must_use]
pub fn prior_residual(a: f64, prior: &EvaluatedGaussianPrior) -> f64 {
    (a - prior.mu) / prior.sigma
}

fn logistic(u: f64) -> f64 {
    if u >= 0.0 {
        let e = (-u).exp();
        1.0 / (1.0 + e)
    } else {
        let e = u.exp();
        e / (1.0 + e)
    }
}

/// C4.6 map \(θ(u) = lo + (hi-lo)\,σ(u)\).
#[must_use]
pub fn theta_from_unbounded(u: f64, lo: f64, hi: f64) -> f64 {
    lo + (hi - lo) * logistic(u)
}

/// C4.6 derivative \(dθ/du = (hi-lo)\,σ(1-σ)\).
#[must_use]
pub fn dtheta_du(u: f64, lo: f64, hi: f64) -> f64 {
    let s = logistic(u);
    (hi - lo) * s * (1.0 - s)
}

pub fn assemble_theta(
    catalog: &Catalog,
    options: &Stage1Options,
    star: Option<&StarRecord>,
    defocus: Option<DefocusMomentInputs>,
) -> Result<ThetaAssembly, PsfFieldError> {
    let layout = assemble_layout(catalog, options)?;
    let (theta, priors) = initialize_theta(catalog, &layout, options, star, defocus)?;
    Ok(ThetaAssembly {
        theta,
        param_meta: layout.param_meta,
        bounds: layout.bounds,
        free_index: layout.free_index,
        priors,
    })
}

/// Freeze/unfreeze steps. After the last step the covariance free set is `assembly.free_index`.
pub fn schedule_steps(
    catalog: &Catalog,
    assembly: &ThetaAssembly,
    options: &Stage1Options,
) -> Vec<ScheduleStep> {
    if !options.use_schedule {
        let mut ids = Vec::new();
        let mut seen = HashSet::new();
        for &idx in &assembly.free_index {
            let term_id = &assembly.param_meta[idx].term_id;
            if seen.insert(term_id.clone()) {
                ids.push(term_id.clone());
            }
        }
        return vec![ScheduleStep {
            name: "all".to_string(),
            unfrozen_term_ids: ids,
            free_index: assembly.free_index.clone(),
        }];
    }

    let n_all = assembly.n_all();
    let mut is_catalog_free = vec![false; n_all];
    for &idx in &assembly.free_index {
        is_catalog_free[idx] = true;
    }

    catalog
        .fit_schedule
        .iter()
        .map(|step| {
            let unfrozen: HashSet<&str> =
                step.unfrozen_term_ids.iter().map(String::as_str).collect();
            let free_index = assembly
                .param_meta
                .iter()
                .enumerate()
                .filter(|(i, meta)| is_catalog_free[*i] && unfrozen.contains(meta.term_id.as_str()))
                .map(|(i, _)| i)
                .collect();
            ScheduleStep {
                name: step.name.clone(),
                unfrozen_term_ids: step.unfrozen_term_ids.clone(),
                free_index,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog::ingested_default_catalog;
    use crate::error::ErrorCode;
    use crate::types::Flag;

    fn options() -> Stage1Options {
        Stage1Options::new(2.8)
    }

    fn dummy_star() -> StarRecord {
        let size = 15;
        let center = (size as f64 - 1.0) / 2.0;
        StarRecord {
            schema_version: "1.0.0".into(),
            star_id: "s1".into(),
            exposure_id: "exp1".into(),
            session_id: "sess1".into(),
            field_xy_mm: [0.0, 0.0],
            source_xy_px: [511.5, 511.5],
            stamp: vec![vec![1.0; size]; size],
            variance: vec![vec![1.0; size]; size],
            centroid_xy_px: [center, center],
            pixel_mask: vec![vec![0; size]; size],
            flags: vec![Flag::Selected],
            flux_sum_adu: 42.0,
        }
    }

    fn dummy_defocus() -> DefocusMomentInputs {
        DefocusMomentInputs {
            sigma_meas_sq: 1.0,
            sigma0_sq: 1.0,
            known_defocus_waves: 0.0,
        }
    }

    #[test]
    fn default_catalog_free_count_matches_enabled_unfrozen_slots() {
        let catalog = ingested_default_catalog();
        let layout = assemble_layout(&catalog, &options()).unwrap();
        let free: Vec<(&str, &str)> = layout
            .free_index
            .iter()
            .map(|&i| {
                (
                    layout.param_meta[i].term_id.as_str(),
                    layout.param_meta[i].role.as_str(),
                )
            })
            .collect();
        assert_eq!(
            free,
            vec![
                ("zernike_2_0", "local_value"),
                ("zernike_2_2", "local_value"),
                ("zernike_2_m2", "local_value"),
                ("zernike_3_1", "local_value"),
                ("zernike_3_m1", "local_value"),
                ("zernike_3_3", "local_value"),
                ("zernike_3_m3", "local_value"),
                ("zernike_4_0", "local_value"),
                ("zernike_4_2", "local_value"),
                ("zernike_4_m2", "local_value"),
                ("flux", "flux"),
                ("moffat_seeing", "alpha_px"),
                ("gaussian_jitter", "sigma_px"),
                ("charge_diffusion", "sigma_px"),
            ]
        );
        assert_eq!(layout.n_free(), free.len());
        assert_eq!(layout.param_meta.len(), layout.n_all());
    }

    #[test]
    fn disabled_terms_absent_frozen_enabled_present() {
        let catalog = ingested_default_catalog();
        let layout = assemble_layout(&catalog, &options()).unwrap();
        let ids: HashSet<&str> = layout
            .param_meta
            .iter()
            .map(|m| m.term_id.as_str())
            .collect();
        for disabled in [
            "zernike_1_1",
            "zernike_1_m1",
            "linear_drift",
            "field_rotation",
            "periodic_error",
        ] {
            assert!(!ids.contains(disabled), "{disabled} should be omitted");
        }

        let piston = layout
            .param_meta
            .iter()
            .find(|m| m.term_id == "zernike_0_0")
            .unwrap();
        assert!(piston.frozen);
        assert_eq!(piston.role, "local_value");

        let sky = layout
            .param_meta
            .iter()
            .find(|m| m.term_id == "sky")
            .unwrap();
        assert!(sky.frozen);
        assert_eq!(sky.role, "sky");

        let beta = layout
            .param_meta
            .iter()
            .find(|m| m.term_id == "moffat_seeing" && m.role == "beta")
            .unwrap();
        assert!(beta.frozen);

        let alpha = layout
            .param_meta
            .iter()
            .find(|m| m.term_id == "moffat_seeing" && m.role == "alpha_px")
            .unwrap();
        assert!(!alpha.frozen);
    }

    #[test]
    fn freeze_mask_length_must_equal_n_free() {
        let catalog = ingested_default_catalog();
        let n_free = assemble_layout(&catalog, &options()).unwrap().n_free();
        assert_eq!(n_free, 14);

        let mut ok = options();
        ok.freeze_mask = Some(vec![false; n_free]);
        assemble_layout(&catalog, &ok).unwrap();

        let mut too_long = options();
        too_long.freeze_mask = Some(vec![false; 15]);
        let err = assemble_layout(&catalog, &too_long).unwrap_err();
        assert_eq!(err.code, ErrorCode::Input);
        assert!(err.message.contains("freeze_mask length"));
        assert!(err.message.contains("does not equal n_free"));

        let mut too_short = options();
        too_short.freeze_mask = Some(vec![false; n_free - 1]);
        let err = assemble_layout(&catalog, &too_short).unwrap_err();
        assert_eq!(err.code, ErrorCode::Input);
        assert!(err.message.contains("freeze_mask length"));
    }

    #[test]
    fn freeze_mask_does_not_zip_truncate() {
        let catalog = ingested_default_catalog();
        let mut options = options();
        options.freeze_mask = Some(vec![true; 20]);
        let err = assemble_layout(&catalog, &options).unwrap_err();
        assert_eq!(err.code, ErrorCode::Input);
        assert!(err.message.contains("does not equal n_free"));
    }

    #[test]
    fn moffat_mean_init_prior_uses_alpha0() {
        let catalog = ingested_default_catalog();
        let assembly = assemble_theta(
            &catalog,
            &options(),
            Some(&dummy_star()),
            Some(dummy_defocus()),
        )
        .unwrap();
        let alpha_i = assembly
            .param_meta
            .iter()
            .position(|m| m.term_id == "moffat_seeing" && m.role == "alpha_px")
            .unwrap();
        let alpha0 = assembly.theta[alpha_i];
        let expected = moffat_alpha_from_fwhm(2.8).unwrap();
        assert!((alpha0 - expected).abs() < 1e-12);

        let prior = assembly.priors[alpha_i].unwrap();
        assert!((prior.mu - alpha0).abs() < 1e-15);
        let expected_sigma = (0.5 * alpha0.abs()).max(1e-3);
        assert!((prior.sigma - expected_sigma).abs() < 1e-15);

        let beta_i = assembly
            .param_meta
            .iter()
            .position(|m| m.term_id == "moffat_seeing" && m.role == "beta")
            .unwrap();
        assert!((assembly.theta[beta_i] - MOFFAT_BETA).abs() < 1e-15);
        assert!(assembly.priors[beta_i].is_none());
    }

    #[test]
    fn default_schedule_has_three_steps_and_moffat_beta_stays_frozen() {
        let catalog = ingested_default_catalog();
        let assembly = assemble_theta(
            &catalog,
            &options(),
            Some(&dummy_star()),
            Some(dummy_defocus()),
        )
        .unwrap();
        let steps = schedule_steps(&catalog, &assembly, &options());
        assert_eq!(steps.len(), 3);
        assert_eq!(steps[0].name, "coarse");
        assert_eq!(steps[1].name, "mid");
        assert_eq!(steps[2].name, "full");
        assert_eq!(steps[0].free_index.len(), 5);
        assert_eq!(steps[1].free_index.len(), 10);
        assert_eq!(steps[2].free_index.len(), 14);
        assert_eq!(steps[2].free_index, assembly.free_index);

        let beta_i = assembly
            .param_meta
            .iter()
            .position(|m| m.term_id == "moffat_seeing" && m.role == "beta")
            .unwrap();
        for step in &steps {
            assert!(
                !step.free_index.contains(&beta_i),
                "moffat beta must stay frozen in {}",
                step.name
            );
        }
    }

    #[test]
    fn use_schedule_false_frees_all_unfrozen_in_one_step() {
        let catalog = ingested_default_catalog();
        let mut options = options();
        options.use_schedule = false;
        let assembly = assemble_theta(
            &catalog,
            &options,
            Some(&dummy_star()),
            Some(dummy_defocus()),
        )
        .unwrap();
        let steps = schedule_steps(&catalog, &assembly, &options);
        assert_eq!(steps.len(), 1);
        assert_eq!(steps[0].name, "all");
        assert_eq!(steps[0].free_index, assembly.free_index);
        assert_eq!(steps[0].free_index.len(), 14);
        assert!(!steps[0].unfrozen_term_ids.iter().any(|id| id == "sky"));
        assert!(!steps[0]
            .unfrozen_term_ids
            .iter()
            .any(|id| id == "zernike_0_0"));
    }

    #[test]
    fn defocus_moment_closed_form() {
        let a = defocus_moment_init(2.0, 1.0, 0.0).unwrap();
        assert!((a - 0.35).abs() < 1e-15);

        let a = defocus_moment_init(1.0, 1.0, 0.3).unwrap();
        assert!((a + 0.3).abs() < 1e-15);

        let a = defocus_moment_init(100.0, 1.0, 0.0).unwrap();
        assert!((a - 2.0).abs() < 1e-15);
    }

    #[test]
    fn logistic_bounds_map_endpoints() {
        let lo = 0.0;
        let hi = 20.0;
        assert!((theta_from_unbounded(0.0, lo, hi) - 10.0).abs() < 1e-12);
        assert!(dtheta_du(0.0, lo, hi) > 0.0);
        assert!(theta_from_unbounded(40.0, lo, hi) > 19.999);
        assert!(theta_from_unbounded(-40.0, lo, hi) < 0.001);
    }

    #[test]
    fn flux_sum_init_clips_at_zero() {
        let catalog = ingested_default_catalog();
        let mut star = dummy_star();
        star.flux_sum_adu = -5.0;
        let assembly =
            assemble_theta(&catalog, &options(), Some(&star), Some(dummy_defocus())).unwrap();
        let flux_i = assembly
            .param_meta
            .iter()
            .position(|m| m.role == "flux")
            .unwrap();
        assert_eq!(assembly.theta[flux_i], 0.0);
    }

    #[test]
    fn jitter_zero_init_uses_numeric_prior_mean() {
        let catalog = ingested_default_catalog();
        let assembly = assemble_theta(
            &catalog,
            &options(),
            Some(&dummy_star()),
            Some(dummy_defocus()),
        )
        .unwrap();
        let i = assembly
            .param_meta
            .iter()
            .position(|m| m.term_id == "gaussian_jitter")
            .unwrap();
        assert!((assembly.theta[i] - 0.1).abs() < 1e-15);
        let prior = assembly.priors[i].unwrap();
        assert!((prior.mu - 0.1).abs() < 1e-15);
        assert!((prior.sigma - 0.3).abs() < 1e-15);
        assert!((prior_residual(0.4, &prior) - 1.0).abs() < 1e-12);
    }
}
