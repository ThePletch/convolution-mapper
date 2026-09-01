//! Catalog JSON types: named aberration, kernel, and photometric terms that the
//! engine interprets without switching on a human label such as `"coma"`. (C3)
//!
//! Illegal field combinations are rejected while converting the JSON shape into
//! these types, so later stages never see a half-valid term.

use serde::{de::Error as DeError, ser::SerializeStruct, Deserialize, Serialize};

/// How widely a coefficient is shared. Stage-1 still fits a local copy per star;
/// the scope decides how stage-2 groups those copies into a field map. (C4.2)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Scope {
    PerStar,
    PerExposure,
    PerSession,
}

/// How to choose the starting value of a Zernike coefficient before fitting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseInit {
    Zero,
    DefocusMoment,
}

/// How to choose the starting value of per-star flux before fitting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FluxInit {
    Zero,
    FluxSum,
}

/// How to choose the starting Moffat scale α before fitting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MoffatInit {
    Zero,
    Fwhm,
}

/// Closed set of image-space convolution kernels in v1. Unknown identifiers are rejected.
/// Stage-1 field rotation is a local trail `(length_px, angle_rad)`, not the three globals. (C3.6)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KernelKind {
    GaussianIso,
    MoffatIso { init: MoffatInit },
    LinearDrift,
    FieldRotation,
}

impl KernelKind {
    #[must_use]
    pub const fn id_str(self) -> &'static str {
        match self {
            Self::GaussianIso => "gaussian_iso",
            Self::MoffatIso { .. } => "moffat_iso",
            Self::LinearDrift => "linear_drift",
            Self::FieldRotation => "field_rotation",
        }
    }
}

/// Photometric scalars are only flux and residual sky; any other identifier is rejected. (C3.2)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhotometricId {
    Flux,
    Sky,
}

impl PhotometricId {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Flux => "flux",
            Self::Sky => "sky",
        }
    }
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

/// Stage-1 (and optional stage-2) prior on a catalog term's local coefficient.
/// The Gaussian arms are split so an absolute σ and a relative σ_rel cannot coexist. (C3.4)
#[derive(Debug, Clone, PartialEq)]
pub enum PriorSpec {
    None {
        stage2: Option<Stage2Prior>,
    },
    Gaussian {
        mean: f64,
        sigma: f64,
        stage2: Option<Stage2Prior>,
    },
    GaussianFromInit {
        sigma_rel: f64,
        stage2: Option<Stage2Prior>,
    },
}

impl PriorSpec {
    #[must_use]
    pub fn stage2(&self) -> Option<&Stage2Prior> {
        match self {
            Self::None { stage2 }
            | Self::Gaussian { stage2, .. }
            | Self::GaussianFromInit { stage2, .. } => stage2.as_ref(),
        }
    }

    #[must_use]
    pub fn is_from_init(&self) -> bool {
        matches!(self, Self::GaussianFromInit { .. })
    }
}

impl Serialize for PriorSpec {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        match self {
            Self::None { stage2 } => {
                let mut state =
                    serializer.serialize_struct("PriorSpec", 1 + usize::from(stage2.is_some()))?;
                state.serialize_field("family", "none")?;
                if let Some(stage2) = stage2 {
                    state.serialize_field("stage2", stage2)?;
                }
                state.end()
            }
            Self::Gaussian {
                mean,
                sigma,
                stage2,
            } => {
                let mut state =
                    serializer.serialize_struct("PriorSpec", 3 + usize::from(stage2.is_some()))?;
                state.serialize_field("family", "gaussian")?;
                state.serialize_field("mean", mean)?;
                state.serialize_field("sigma", sigma)?;
                if let Some(stage2) = stage2 {
                    state.serialize_field("stage2", stage2)?;
                }
                state.end()
            }
            Self::GaussianFromInit { sigma_rel, stage2 } => {
                let mut state =
                    serializer.serialize_struct("PriorSpec", 3 + usize::from(stage2.is_some()))?;
                state.serialize_field("family", "gaussian")?;
                state.serialize_field("mean", "init")?;
                state.serialize_field("sigma_rel", sigma_rel)?;
                if let Some(stage2) = stage2 {
                    state.serialize_field("stage2", stage2)?;
                }
                state.end()
            }
        }
    }
}

impl<'de> Deserialize<'de> for PriorSpec {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(rename_all = "snake_case")]
        enum Family {
            None,
            Gaussian,
        }

        #[derive(Deserialize)]
        #[serde(untagged)]
        enum Mean {
            Number(f64),
            Init(InitMean),
        }

        #[derive(Deserialize)]
        #[serde(rename_all = "lowercase")]
        enum InitMean {
            Init,
        }

        #[derive(Deserialize)]
        #[serde(deny_unknown_fields)]
        struct Raw {
            family: Family,
            #[serde(default)]
            mean: Option<Mean>,
            #[serde(default)]
            sigma: Option<f64>,
            #[serde(default)]
            sigma_rel: Option<f64>,
            #[serde(default)]
            stage2: Option<Stage2Prior>,
        }

        let raw = Raw::deserialize(deserializer)?;
        match raw.family {
            Family::None => {
                if raw.mean.is_some() || raw.sigma.is_some() || raw.sigma_rel.is_some() {
                    return Err(DeError::custom(
                        "none prior must omit mean, sigma, and sigma_rel",
                    ));
                }
                Ok(Self::None { stage2: raw.stage2 })
            }
            Family::Gaussian => match raw.mean {
                Some(Mean::Number(mean)) => {
                    if raw.sigma_rel.is_some() {
                        return Err(DeError::custom(
                            "gaussian prior with numeric mean requires sigma and must omit sigma_rel",
                        ));
                    }
                    let Some(sigma) = raw.sigma else {
                        return Err(DeError::custom(
                            "gaussian prior with numeric mean requires sigma and must omit sigma_rel",
                        ));
                    };
                    Ok(Self::Gaussian {
                        mean,
                        sigma,
                        stage2: raw.stage2,
                    })
                }
                Some(Mean::Init(_)) => {
                    if raw.sigma.is_some() {
                        return Err(DeError::custom(
                            "gaussian prior with mean \"init\" requires sigma_rel and must omit sigma",
                        ));
                    }
                    let Some(sigma_rel) = raw.sigma_rel else {
                        return Err(DeError::custom(
                            "gaussian prior with mean \"init\" requires sigma_rel and must omit sigma",
                        ));
                    };
                    Ok(Self::GaussianFromInit {
                        sigma_rel,
                        stage2: raw.stage2,
                    })
                }
                None => Err(DeError::custom("gaussian prior requires mean")),
            },
        }
    }
}

/// One catalog term: Zernike phase, image-space kernel, or photometric scalar.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "wire::ErrorTerm", into = "wire::ErrorTerm")]
pub enum ErrorTerm {
    Phase {
        term_id: String,
        name: String,
        scope: Scope,
        /// When true, Levenberg–Marquardt holds the coefficient at initialization.
        frozen: bool,
        /// When false, the term is omitted from θ entirely (not frozen at zero).
        enabled: bool,
        bounds: Option<[f64; 2]>,
        init: PhaseInit,
        prior: PriorSpec,
        units: String,
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
        bounds: Option<[f64; 2]>,
        prior: PriorSpec,
        units: String,
        report: Option<serde_json::Value>,
        kernel: KernelKind,
        /// `None` means the kernel is spatially uniform across the field. (C3.2.2)
        field_basis: Option<FieldBasis>,
    },
    Flux {
        name: String,
        scope: Scope,
        frozen: bool,
        enabled: bool,
        bounds: Option<[f64; 2]>,
        init: FluxInit,
        prior: PriorSpec,
        units: String,
        report: Option<serde_json::Value>,
    },
    Sky {
        name: String,
        scope: Scope,
        frozen: bool,
        enabled: bool,
        bounds: Option<[f64; 2]>,
        prior: PriorSpec,
        units: String,
        report: Option<serde_json::Value>,
    },
}

impl ErrorTerm {
    #[must_use]
    pub fn term_id(&self) -> &str {
        match self {
            Self::Phase { term_id, .. } | Self::Kernel { term_id, .. } => term_id,
            Self::Flux { .. } => PhotometricId::Flux.as_str(),
            Self::Sky { .. } => PhotometricId::Sky.as_str(),
        }
    }

    #[must_use]
    pub fn enabled(&self) -> bool {
        match self {
            Self::Phase { enabled, .. }
            | Self::Kernel { enabled, .. }
            | Self::Flux { enabled, .. }
            | Self::Sky { enabled, .. } => *enabled,
        }
    }

    #[must_use]
    pub fn frozen(&self) -> bool {
        match self {
            Self::Phase { frozen, .. }
            | Self::Kernel { frozen, .. }
            | Self::Flux { frozen, .. }
            | Self::Sky { frozen, .. } => *frozen,
        }
    }

    #[must_use]
    pub fn scope(&self) -> Scope {
        match self {
            Self::Phase { scope, .. }
            | Self::Kernel { scope, .. }
            | Self::Flux { scope, .. }
            | Self::Sky { scope, .. } => *scope,
        }
    }

    #[must_use]
    pub fn units(&self) -> &str {
        match self {
            Self::Phase { units, .. }
            | Self::Kernel { units, .. }
            | Self::Flux { units, .. }
            | Self::Sky { units, .. } => units,
        }
    }

    #[must_use]
    pub fn bounds(&self) -> Option<[f64; 2]> {
        match self {
            Self::Phase { bounds, .. }
            | Self::Kernel { bounds, .. }
            | Self::Flux { bounds, .. }
            | Self::Sky { bounds, .. } => *bounds,
        }
    }

    #[must_use]
    pub fn prior(&self) -> &PriorSpec {
        match self {
            Self::Phase { prior, .. }
            | Self::Kernel { prior, .. }
            | Self::Flux { prior, .. }
            | Self::Sky { prior, .. } => prior,
        }
    }

    #[must_use]
    pub fn field_basis(&self) -> Option<&FieldBasis> {
        match self {
            Self::Phase { field_basis, .. } => Some(field_basis),
            Self::Kernel { field_basis, .. } => field_basis.as_ref(),
            Self::Flux { .. } | Self::Sky { .. } => None,
        }
    }

    #[must_use]
    pub fn init_is_zero(&self) -> bool {
        match self {
            Self::Phase {
                init: PhaseInit::Zero,
                ..
            }
            | Self::Flux {
                init: FluxInit::Zero,
                ..
            }
            | Self::Sky { .. }
            | Self::Kernel {
                kernel:
                    KernelKind::GaussianIso
                    | KernelKind::LinearDrift
                    | KernelKind::FieldRotation
                    | KernelKind::MoffatIso {
                        init: MoffatInit::Zero,
                    },
                ..
            } => true,
            Self::Phase {
                init: PhaseInit::DefocusMoment,
                ..
            }
            | Self::Flux {
                init: FluxInit::FluxSum,
                ..
            }
            | Self::Kernel {
                kernel:
                    KernelKind::MoffatIso {
                        init: MoffatInit::Fwhm,
                    },
                ..
            } => false,
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
    pub matrix: NullMatrix,
}

/// JSON `null` for [`Bundle::matrix`]. Any other value fails to deserialize.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct NullMatrix;

impl Serialize for NullMatrix {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_none()
    }
}

impl<'de> Deserialize<'de> for NullMatrix {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        match Option::<serde_json::Value>::deserialize(deserializer)? {
            None => Ok(Self),
            Some(_) => Err(DeError::custom("bundle.matrix must be null in v1")),
        }
    }
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

mod wire {
    use super::{
        FieldBasis, FluxInit, KernelKind, MoffatInit, PhaseInit, PhotometricId, PriorSpec, Scope,
    };
    use serde::{Deserialize, Serialize};

    #[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
    #[serde(rename_all = "snake_case")]
    pub enum InitMethod {
        Zero,
        FluxSum,
        DefocusMoment,
        MoffatFwhm,
    }

    impl InitMethod {
        pub const fn as_str(self) -> &'static str {
            match self {
                Self::Zero => "zero",
                Self::FluxSum => "flux_sum",
                Self::DefocusMoment => "defocus_moment",
                Self::MoffatFwhm => "moffat_fwhm",
            }
        }
    }

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
    }

    #[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
    pub struct KernelSpec {
        pub id: KernelId,
    }

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

    fn reject_init_mean_with_zero(
        prior: &PriorSpec,
        init_is_zero: bool,
        term_id: &str,
    ) -> Result<(), String> {
        if prior.is_from_init() && init_is_zero {
            return Err(format!(
                "term {term_id} pairs mean \"init\" with init method zero"
            ));
        }
        Ok(())
    }

    impl TryFrom<ErrorTerm> for super::ErrorTerm {
        type Error = String;

        fn try_from(term: ErrorTerm) -> Result<Self, Self::Error> {
            match term {
                ErrorTerm::Phase {
                    term_id,
                    name,
                    scope,
                    frozen,
                    enabled,
                    bounds,
                    init,
                    prior,
                    units,
                    report,
                    n,
                    m,
                    field_basis,
                } => {
                    let init = match init.method {
                        InitMethod::Zero => PhaseInit::Zero,
                        InitMethod::DefocusMoment if n == 2 && m == 0 => PhaseInit::DefocusMoment,
                        other => {
                            return Err(format!(
                                "init method {} is not valid for term {term_id}",
                                other.as_str()
                            ));
                        }
                    };
                    reject_init_mean_with_zero(&prior, matches!(init, PhaseInit::Zero), &term_id)?;
                    Ok(Self::Phase {
                        term_id,
                        name,
                        scope,
                        frozen,
                        enabled,
                        bounds,
                        init,
                        prior,
                        units,
                        report,
                        n,
                        m,
                        field_basis,
                    })
                }
                ErrorTerm::Kernel {
                    term_id,
                    name,
                    scope,
                    frozen,
                    enabled,
                    bounds,
                    init,
                    prior,
                    units,
                    report,
                    kernel,
                    field_basis,
                } => {
                    let kernel = match (kernel.id, init.method) {
                        (KernelId::GaussianIso, InitMethod::Zero) => KernelKind::GaussianIso,
                        (KernelId::MoffatIso, InitMethod::Zero) => KernelKind::MoffatIso {
                            init: MoffatInit::Zero,
                        },
                        (KernelId::MoffatIso, InitMethod::MoffatFwhm) => KernelKind::MoffatIso {
                            init: MoffatInit::Fwhm,
                        },
                        (KernelId::LinearDrift, InitMethod::Zero) => KernelKind::LinearDrift,
                        (KernelId::FieldRotation, InitMethod::Zero) => KernelKind::FieldRotation,
                        (_, method) => {
                            return Err(format!(
                                "init method {} is not valid for term {term_id}",
                                method.as_str()
                            ));
                        }
                    };
                    reject_init_mean_with_zero(
                        &prior,
                        matches!(
                            kernel,
                            KernelKind::MoffatIso {
                                init: MoffatInit::Zero
                            } | KernelKind::GaussianIso
                                | KernelKind::LinearDrift
                                | KernelKind::FieldRotation
                        ),
                        &term_id,
                    )?;
                    Ok(Self::Kernel {
                        term_id,
                        name,
                        scope,
                        frozen,
                        enabled,
                        bounds,
                        prior,
                        units,
                        report,
                        kernel,
                        field_basis,
                    })
                }
                ErrorTerm::Photometric {
                    term_id,
                    name,
                    scope,
                    frozen,
                    enabled,
                    bounds,
                    init,
                    prior,
                    units,
                    report,
                } => match term_id.as_str() {
                    "flux" => {
                        let init = match init.method {
                            InitMethod::Zero => FluxInit::Zero,
                            InitMethod::FluxSum => FluxInit::FluxSum,
                            other => {
                                return Err(format!(
                                    "init method {} is not valid for term {term_id}",
                                    other.as_str()
                                ));
                            }
                        };
                        reject_init_mean_with_zero(
                            &prior,
                            matches!(init, FluxInit::Zero),
                            &term_id,
                        )?;
                        Ok(Self::Flux {
                            name,
                            scope,
                            frozen,
                            enabled,
                            bounds,
                            init,
                            prior,
                            units,
                            report,
                        })
                    }
                    "sky" => {
                        if init.method != InitMethod::Zero {
                            return Err(format!(
                                "init method {} is not valid for term {term_id}",
                                init.method.as_str()
                            ));
                        }
                        reject_init_mean_with_zero(&prior, true, &term_id)?;
                        Ok(Self::Sky {
                            name,
                            scope,
                            frozen,
                            enabled,
                            bounds,
                            prior,
                            units,
                            report,
                        })
                    }
                    other => Err(format!("photometric term_id {other} is not flux or sky")),
                },
            }
        }
    }

    impl From<super::ErrorTerm> for ErrorTerm {
        fn from(term: super::ErrorTerm) -> Self {
            match term {
                super::ErrorTerm::Phase {
                    term_id,
                    name,
                    scope,
                    frozen,
                    enabled,
                    bounds,
                    init,
                    prior,
                    units,
                    report,
                    n,
                    m,
                    field_basis,
                } => Self::Phase {
                    term_id,
                    name,
                    scope,
                    frozen,
                    enabled,
                    bounds,
                    init: InitSpec {
                        method: match init {
                            PhaseInit::Zero => InitMethod::Zero,
                            PhaseInit::DefocusMoment => InitMethod::DefocusMoment,
                        },
                    },
                    prior,
                    units,
                    report,
                    n,
                    m,
                    field_basis,
                },
                super::ErrorTerm::Kernel {
                    term_id,
                    name,
                    scope,
                    frozen,
                    enabled,
                    bounds,
                    prior,
                    units,
                    report,
                    kernel,
                    field_basis,
                } => {
                    let (id, method) = match kernel {
                        KernelKind::GaussianIso => (KernelId::GaussianIso, InitMethod::Zero),
                        KernelKind::MoffatIso {
                            init: MoffatInit::Zero,
                        } => (KernelId::MoffatIso, InitMethod::Zero),
                        KernelKind::MoffatIso {
                            init: MoffatInit::Fwhm,
                        } => (KernelId::MoffatIso, InitMethod::MoffatFwhm),
                        KernelKind::LinearDrift => (KernelId::LinearDrift, InitMethod::Zero),
                        KernelKind::FieldRotation => (KernelId::FieldRotation, InitMethod::Zero),
                    };
                    Self::Kernel {
                        term_id,
                        name,
                        scope,
                        frozen,
                        enabled,
                        bounds,
                        init: InitSpec { method },
                        prior,
                        units,
                        report,
                        kernel: KernelSpec { id },
                        field_basis,
                    }
                }
                super::ErrorTerm::Flux {
                    name,
                    scope,
                    frozen,
                    enabled,
                    bounds,
                    init,
                    prior,
                    units,
                    report,
                } => Self::Photometric {
                    term_id: PhotometricId::Flux.as_str().to_string(),
                    name,
                    scope,
                    frozen,
                    enabled,
                    bounds,
                    init: InitSpec {
                        method: match init {
                            FluxInit::Zero => InitMethod::Zero,
                            FluxInit::FluxSum => InitMethod::FluxSum,
                        },
                    },
                    prior,
                    units,
                    report,
                },
                super::ErrorTerm::Sky {
                    name,
                    scope,
                    frozen,
                    enabled,
                    bounds,
                    prior,
                    units,
                    report,
                } => Self::Photometric {
                    term_id: PhotometricId::Sky.as_str().to_string(),
                    name,
                    scope,
                    frozen,
                    enabled,
                    bounds,
                    init: InitSpec {
                        method: InitMethod::Zero,
                    },
                    prior,
                    units,
                    report,
                },
            }
        }
    }
}
