//! Catalog JSON types: named aberration, kernel, and photometric terms that the
//! engine interprets without switching on a human label such as `"coma"`. (C3)

use serde::{Deserialize, Serialize};

/// How widely a coefficient is shared. Stage-1 still fits a local copy per star;
/// the scope decides how stage-2 groups those copies into a field map. (C4.2)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Scope {
    PerStar,
    PerExposure,
    PerSession,
}

/// How to choose the starting value of a term's local coefficient before fitting.
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

/// Closed set of image-space convolution kernels in v1. Unknown identifiers are rejected.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum KernelId {
    GaussianIso,
    MoffatIso,
    LinearDrift,
    FieldRotation,
}

/// Which closed-form image-space kernel a catalog term implements.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct KernelSpec {
    pub id: KernelId,
}

/// Which monomials in normalized field (u, v) a term's stage-2 map may use.
/// The local stage-1 fit still recovers a single scalar a; the monomials appear in stage-2.
/// (C3.3)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FieldBasis {
    /// Polynomial family; v1 allows only monomials in (u, v). Other families are rejected.
    pub family: FieldFamily,
    /// Highest total degree max(i+j) among `terms`. Must match that maximum exactly.
    pub degree: u32,
    /// Included monomials as `[i, j]` pairs meaning the term `u^i v^j`.
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
    /// Prior mean for each monomial coefficient c_{ij}, in the same order as `FieldBasis.terms`.
    pub mean: Vec<f64>,
    /// Prior standard deviation for each monomial coefficient, same length as `mean`.
    pub sigma: Vec<f64>,
}

/// Mean of a Gaussian prior: a catalog number, or the string `"init"` meaning "use a₀".
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
        /// Absolute prior σ when `mean` is numeric. Omitted when `mean` is `"init"`.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        sigma: Option<f64>,
        /// Relative prior σ when `mean` is `"init"`: σ = max(σ_rel |a₀|, 10⁻³).
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
        /// When true, Levenberg–Marquardt holds the coefficient at initialization.
        frozen: bool,
        /// When false, the term is omitted from θ entirely (not frozen at zero).
        enabled: bool,
        #[serde(default)]
        bounds: Option<[f64; 2]>,
        init: InitSpec,
        prior: PriorSpec,
        units: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        report: Option<serde_json::Value>,
        /// ANSI/OSA radial order n.
        n: u32,
        /// ANSI/OSA azimuthal index m: positive cosine, negative sine, zero if rotationally symmetric.
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
        /// `None` means the kernel is spatially uniform across the field. (C3.2.2)
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

    #[must_use]
    pub fn enabled(&self) -> bool {
        match self {
            Self::Phase { enabled, .. }
            | Self::Kernel { enabled, .. }
            | Self::Photometric { enabled, .. } => *enabled,
        }
    }

    #[must_use]
    pub fn frozen(&self) -> bool {
        match self {
            Self::Phase { frozen, .. }
            | Self::Kernel { frozen, .. }
            | Self::Photometric { frozen, .. } => *frozen,
        }
    }

    #[must_use]
    pub fn scope(&self) -> Scope {
        match self {
            Self::Phase { scope, .. }
            | Self::Kernel { scope, .. }
            | Self::Photometric { scope, .. } => *scope,
        }
    }

    #[must_use]
    pub fn units(&self) -> &str {
        match self {
            Self::Phase { units, .. }
            | Self::Kernel { units, .. }
            | Self::Photometric { units, .. } => units,
        }
    }

    #[must_use]
    pub fn bounds(&self) -> Option<[f64; 2]> {
        match self {
            Self::Phase { bounds, .. }
            | Self::Kernel { bounds, .. }
            | Self::Photometric { bounds, .. } => *bounds,
        }
    }

    #[must_use]
    pub fn init(&self) -> &InitSpec {
        match self {
            Self::Phase { init, .. }
            | Self::Kernel { init, .. }
            | Self::Photometric { init, .. } => init,
        }
    }

    #[must_use]
    pub fn prior(&self) -> &PriorSpec {
        match self {
            Self::Phase { prior, .. }
            | Self::Kernel { prior, .. }
            | Self::Photometric { prior, .. } => prior,
        }
    }

    #[must_use]
    pub fn field_basis(&self) -> Option<&FieldBasis> {
        match self {
            Self::Phase { field_basis, .. } => Some(field_basis),
            Self::Kernel { field_basis, .. } => field_basis.as_ref(),
            Self::Photometric { .. } => None,
        }
    }
}

impl InitMethod {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Zero => "zero",
            Self::FluxSum => "flux_sum",
            Self::DefocusMoment => "defocus_moment",
            Self::MoffatFwhm => "moffat_fwhm",
        }
    }
}

/// Named linear grouping of catalog terms. v1 ships with a null mechanical-sensitivity
/// matrix; a non-null matrix is reserved and rejected. (C3.7.1)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Bundle {
    pub bundle_id: String,
    pub name: String,
    pub term_ids: Vec<String>,
    /// Mechanical sensitivity matrix, shape (n_mech, n_field_coeffs). Must be null in v1.
    pub matrix: Option<serde_json::Value>,
}

/// One freeze/unfreeze step in the staged stage-1 fit. Each step starts from the
/// previous step's θ; coefficients not listed stay held. (C4.4)
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
