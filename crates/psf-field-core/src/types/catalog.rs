//! Catalog JSON types (C3). Semantic ingest (θ, freeze) is a later PR.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Scope {
    PerStar,
    PerExposure,
    PerSession,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum InitMethod {
    Zero,
    FluxSum,
    DefocusMoment,
    MoffatFwhm,
}

/// How to initialize a term's local coefficient before Levenberg–Marquardt.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct InitSpec {
    pub method: InitMethod,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum KernelId {
    GaussianIso,
    MoffatIso,
    LinearDrift,
    FieldRotation,
    PeriodicError,
}

/// Which closed-form image-space kernel a catalog term implements.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct KernelSpec {
    pub id: KernelId,
}

/// Which monomials in normalized field (u, v) a term's stage-2 map may use.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FieldBasis {
    pub family: FieldFamily,
    pub degree: u32,
    pub terms: Vec<[u32; 2]>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FieldFamily {
    Monomial,
}

/// Optional Gaussian prior on stage-2 field-map coefficients, aligned with [`FieldBasis::terms`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Stage2Prior {
    pub mean: Vec<f64>,
    pub sigma: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum PriorMean {
    Number(f64),
    Init(InitMean),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum InitMean {
    Init,
}

/// Stage-1 (and optional stage-2) prior on a catalog term's local coefficient.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "family", rename_all = "snake_case")]
pub enum PriorSpec {
    None {
        #[serde(default, skip_serializing_if = "Option::is_none")]
        stage2: Option<Stage2Prior>,
    },
    Gaussian {
        mean: PriorMean,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        sigma: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        sigma_rel: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        stage2: Option<Stage2Prior>,
    },
}

/// One catalog term: Zernike phase, image-space kernel, or photometric scalar.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum ErrorTerm {
    Phase {
        term_id: String,
        name: String,
        scope: Scope,
        frozen: bool,
        enabled: bool,
        #[serde(default)]
        bounds: Option<[f64; 2]>,
        init: InitSpec,
        prior: PriorSpec,
        units: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        report: Option<serde_json::Value>,
        n: u32,
        m: i32,
        field_basis: FieldBasis,
    },
    Kernel {
        term_id: String,
        name: String,
        scope: Scope,
        frozen: bool,
        enabled: bool,
        #[serde(default)]
        bounds: Option<[f64; 2]>,
        init: InitSpec,
        prior: PriorSpec,
        units: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        report: Option<serde_json::Value>,
        kernel: KernelSpec,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        field_basis: Option<FieldBasis>,
    },
    Photometric {
        term_id: String,
        name: String,
        scope: Scope,
        frozen: bool,
        enabled: bool,
        #[serde(default)]
        bounds: Option<[f64; 2]>,
        init: InitSpec,
        prior: PriorSpec,
        units: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        report: Option<serde_json::Value>,
    },
}

impl ErrorTerm {
    #[must_use]
    pub fn term_id(&self) -> &str {
        match self {
            Self::Phase { term_id, .. }
            | Self::Kernel { term_id, .. }
            | Self::Photometric { term_id, .. } => term_id,
        }
    }
}

/// Named linear grouping of catalog terms; v1 ships with a null sensitivity matrix.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Bundle {
    pub bundle_id: String,
    pub name: String,
    pub term_ids: Vec<String>,
    pub matrix: Option<serde_json::Value>,
}

/// One freeze/unfreeze step in the staged stage-1 fit.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FitScheduleStep {
    pub name: String,
    pub unfrozen_term_ids: Vec<String>,
}

/// Named aberration and kernel terms, optional bundles, and the default fit schedule.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Catalog {
    pub schema_version: String,
    pub catalog_id: String,
    pub terms: Vec<ErrorTerm>,
    pub bundles: Vec<Bundle>,
    pub fit_schedule: Vec<FitScheduleStep>,
}

impl Catalog {
    /// Structural checks only (schema version, ids). Semantic θ ingest is C3/C4.
    pub fn ingest(self) -> Result<Self, crate::error::PsfFieldError> {
        use crate::types::common::{check_schema_version, check_term_id};
        check_schema_version(&self.schema_version)?;
        if self.catalog_id.is_empty() {
            return Err(crate::error::PsfFieldError::input(
                crate::error::ErrorModule::Boundary,
                "catalog_id is empty",
            ));
        }
        for term in &self.terms {
            check_term_id("term_id", term.term_id())?;
        }
        for bundle in &self.bundles {
            check_term_id("bundle_id", &bundle.bundle_id)?;
            if bundle.matrix.is_some() {
                return Err(crate::error::PsfFieldError::input(
                    crate::error::ErrorModule::Boundary,
                    "bundle.matrix must be null in v1",
                ));
            }
        }
        Ok(self)
    }
}
