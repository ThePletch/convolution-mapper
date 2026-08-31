//! Output / report types 1:1 with the remaining JSON schemas.

use serde::{Deserialize, Serialize};

use crate::types::common::{DegeneratePair, ParamMeta};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Termination {
    ConvergedFtol,
    ConvergedXtol,
    ConvergedGtol,
    ConvergedZeroResidual,
    MaxEval,
    Numerical,
    LostPatience,
    User,
    Unknown,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Stage1ErrorPayload {
    pub code: String,
    pub message: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Stage1Result {
    pub schema_version: String,
    pub star_id: String,
    pub theta: Vec<f64>,
    pub theta_init: Vec<f64>,
    pub param_meta: Vec<ParamMeta>,
    pub free_index: Vec<i64>,
    pub covariance: Vec<Vec<f64>>,
    pub covariance_chi2_scaled: Vec<Vec<f64>>,
    pub covariance_ok: bool,
    pub correlation: Vec<Vec<f64>>,
    pub degenerate_pairs: Vec<DegeneratePair>,
    pub residual_image: Vec<Vec<f64>>,
    pub weighted_residual: Vec<f64>,
    pub chi2: f64,
    pub chi2_reduced: f64,
    pub n_iter: i64,
    pub n_fev: i64,
    pub termination: Termination,
    pub converged: bool,
    pub defocus_sign_ambiguous: bool,
    pub flag_at_bound: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub error: Option<Stage1ErrorPayload>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FieldMap {
    pub coefficients: Vec<f64>,
    pub covariance: Vec<Vec<f64>>,
    pub cond2: f64,
    pub ill_conditioned: bool,
    pub independence_assumption: bool,
    pub residuals: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Stage2Result {
    pub schema_version: String,
    pub maps: std::collections::BTreeMap<String, FieldMap>,
    pub kernel_globals: serde_json::Value,
    pub dropped_star_ids: Vec<String>,
    pub n_stars_used: i64,
    pub param_meta: Vec<serde_json::Value>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub star_ids_used: Option<Vec<String>>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PsfEval {
    pub schema_version: String,
    pub psf: Vec<Vec<f64>>,
    pub zernike_vector: std::collections::BTreeMap<String, f64>,
    pub kernel_vector: std::collections::BTreeMap<String, f64>,
    pub field_xy_mm: [f64; 2],
    pub u_v: [f64; 2],
    pub image_meta_digest: String,
    pub catalog_id: String,
    pub stage2_schema_version: String,
    pub extrapolated: bool,
    pub outside_unit_square: bool,
    pub outside_hull: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FdReport {
    pub schema_version: String,
    pub column_errors: Vec<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub column_errors_unconstrained: Option<Vec<f64>>,
    pub passed: bool,
    pub passed_unconstrained: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScoreEntry {
    pub term_id: String,
    pub score: f64,
    pub suggest_add: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub undefined: Option<bool>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StarScore {
    pub star_id: String,
    pub chi2_reduced: f64,
    pub structured_residual: bool,
    pub centroid_leak_suspect: bool,
    pub scores: Vec<ScoreEntry>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScoreReport {
    pub schema_version: String,
    pub per_star: Vec<StarScore>,
    pub session_degeneracies: Vec<serde_json::Value>,
    pub stage2_maps: serde_json::Value,
    pub weak_phase_all_fraction: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CoverageReport {
    pub schema_version: String,
    pub n_detected: i64,
    pub n_candidate: i64,
    pub n_selected: i64,
    pub frac_selected_of_detected: f64,
    pub grid_3x3: [i64; 9],
    pub empty_cells: i64,
    pub convex_hull_area_mm2: f64,
    pub detector_area_mm2: f64,
    pub hull_fill: f64,
    pub design_cond_plane: f64,
    pub design_cond_quad: f64,
}
