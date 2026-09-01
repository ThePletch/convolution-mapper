//! Semantic catalog ingest: validating aberration terms, field bases, priors, and
//! initialization methods before they are assembled into the stage-1 parameter vector.
//! (C3)

use std::collections::HashSet;

use crate::error::{ErrorModule, PsfFieldError};
use crate::types::{
    check_finite, check_positive, check_schema_version, check_term_id, Catalog, ErrorTerm,
    FieldBasis, InitMethod, KernelId, PriorMean, PriorSpec, Stage2Prior,
};

/// Highest Zernike radial order `n` accepted in v1. Higher orders are rejected here,
/// not later during evaluation. (C2.2)
const MAX_ZERNIKE_N: u32 = 15;

/// Highest polynomial degree of a monomial field map in v1: constant, linear, or quadratic.
/// (C3.3)
const MAX_FIELD_DEGREE: u32 = 2;

fn input(message: impl Into<String>) -> PsfFieldError {
    PsfFieldError::input(ErrorModule::Boundary, message)
}

/// One local kernel coefficient, in the order the catalog lists kernel parameters.
/// Stage-1 concatenates these after photometric slots. (C3.6, C4.1)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KernelParameter {
    /// Name stored on `ParamMeta.role` (for example `sigma_px` or `alpha_px`).
    pub role: &'static str,
    /// Physical unit of this coefficient (pixels, radians, seconds, or dimensionless).
    pub unit: &'static str,
    /// When true, this coefficient stays frozen even if the parent catalog term is free.
    /// Moffat `beta` is frozen at 2.5 in v1 so it is never a Levenberg–Marquardt unknown.
    pub always_frozen: bool,
}

/// Kernel coefficients in catalog listing order, used when flattening a term into θ.
#[must_use]
pub fn kernel_parameters(id: KernelId) -> &'static [KernelParameter] {
    match id {
        KernelId::GaussianIso => &[KernelParameter {
            // Isotropic Gaussian width on the detector, in pixels. (C3.6.1)
            role: "sigma_px",
            unit: "px",
            always_frozen: false,
        }],
        KernelId::MoffatIso => &[
            KernelParameter {
                // Moffat scale α on the detector, in pixels. Related to FWHM by
                // HWHM = α √(2^{1/β} − 1). (C3.6.2, C3.5.2)
                role: "alpha_px",
                unit: "px",
                always_frozen: false,
            },
            KernelParameter {
                // Moffat exponent β; dimensionless and frozen in v1. (C3.6.2)
                role: "beta",
                unit: "1",
                always_frozen: true,
            },
        ],
        KernelId::LinearDrift => &[
            KernelParameter {
                // Trail length of a unit-sum Gaussian line segment, in detector pixels. (C3.6.3)
                role: "length_px",
                unit: "px",
                always_frozen: false,
            },
            KernelParameter {
                // Trail direction from detector +x toward +y, in radians. (C3.6.3)
                role: "angle_rad",
                unit: "rad",
                always_frozen: false,
            },
        ],
        KernelId::FieldRotation => &[
            KernelParameter {
                // Stage-1 fits a local trail equivalent to field rotation, not the three
                // global (center, ω) parameters. Length in detector pixels. (C3.6.4, C4.1)
                role: "length_px",
                unit: "px",
                always_frozen: false,
            },
            KernelParameter {
                role: "angle_rad",
                unit: "rad",
                always_frozen: false,
            },
        ],
        KernelId::PeriodicError => &[
            KernelParameter {
                // Peak-to-center amplitude of the sinusoidal RA trail, in detector pixels. (C3.6.5)
                role: "amp_px",
                unit: "px",
                always_frozen: false,
            },
            KernelParameter {
                // Period of the mount periodic error, in seconds. Often frozen from metadata. (C3.6.5)
                role: "period_s",
                unit: "s",
                always_frozen: false,
            },
            KernelParameter {
                // Phase of the sinusoid at exposure midpoint, in radians. (C3.6.5, NOR.11)
                role: "phase_rad",
                unit: "rad",
                always_frozen: false,
            },
        ],
    }
}

impl Catalog {
    /// Structural and semantic checks: unique term identifiers, valid Zernike indices,
    /// field-basis monomials, prior/initialization pairings, and a null bundle matrix in v1.
    pub fn ingest(self) -> Result<Self, PsfFieldError> {
        check_schema_version(&self.schema_version)?;
        if self.catalog_id.is_empty() {
            return Err(input("catalog_id is empty"));
        }

        let mut seen_ids = HashSet::new();
        for term in &self.terms {
            check_term_id("term_id", term.term_id())?;
            if !seen_ids.insert(term.term_id().to_string()) {
                return Err(input(format!("duplicate term_id {}", term.term_id())));
            }
            validate_term(term)?;
        }

        for bundle in &self.bundles {
            check_term_id("bundle_id", &bundle.bundle_id)?;
            if bundle.matrix.is_some() {
                return Err(input("bundle.matrix must be null in v1"));
            }
        }
        Ok(self)
    }
}

fn validate_term(term: &ErrorTerm) -> Result<(), PsfFieldError> {
    if let Some(basis) = term.field_basis() {
        validate_field_basis(basis)?;
    }
    match term {
        ErrorTerm::Phase { n, m, .. } => validate_zernike_nm(*n, *m)?,
        ErrorTerm::Photometric { term_id, .. } => {
            if term_id != "flux" && term_id != "sky" {
                return Err(input(format!(
                    "photometric term_id {term_id} is not flux or sky"
                )));
            }
        }
        ErrorTerm::Kernel { .. } => {}
    }
    validate_bounds(term.bounds())?;
    validate_prior(term.prior())?;
    validate_init_pairing(term)?;
    if prior_mean_is_init(term.prior()) && term.init().method == InitMethod::Zero {
        return Err(input(format!(
            "term {} pairs mean \"init\" with init method zero",
            term.term_id()
        )));
    }
    Ok(())
}

/// Reject ANSI/OSA indices that are not a Zernike mode: radial order `n` and
/// azimuthal frequency `m` must satisfy `|m| ≤ n` with `n − |m|` even, and `n`
/// must not exceed the v1 cap. (C2.1, C2.2)
pub fn validate_zernike_nm(n: u32, m: i32) -> Result<(), PsfFieldError> {
    if n > MAX_ZERNIKE_N {
        return Err(input(format!(
            "Zernike n={n} exceeds frozen max {MAX_ZERNIKE_N}"
        )));
    }
    let m_abs = m.unsigned_abs();
    if m_abs > n {
        return Err(input(format!("invalid Zernike (n,m)=({n},{m}): |m| > n")));
    }
    if (n - m_abs) % 2 == 1 {
        return Err(input(format!(
            "invalid Zernike (n,m)=({n},{m}): n - |m| is odd"
        )));
    }
    Ok(())
}

/// A field basis is a list of monomials `u^i v^j` in normalized field coordinates.
/// `degree` must equal `max(i+j)`, pairs must be unique, and they must be sorted
/// by total degree then by the `v` exponent so the default catalog JSON ingests.
/// (C3.3)
pub fn validate_field_basis(basis: &FieldBasis) -> Result<(), PsfFieldError> {
    if basis.degree > MAX_FIELD_DEGREE {
        return Err(input(format!(
            "field basis degree {} exceeds v1 max {MAX_FIELD_DEGREE}",
            basis.degree
        )));
    }
    if basis.terms.is_empty() {
        return Err(input("field basis terms must not be empty"));
    }
    let mut seen = HashSet::new();
    let mut max_degree = 0_u32;
    let mut previous: Option<(u32, u32)> = None;
    for pair in &basis.terms {
        let [i, j] = *pair;
        let degree = i + j;
        if degree > basis.degree {
            return Err(input(format!(
                "field basis term [{i}, {j}] has i+j greater than degree {}",
                basis.degree
            )));
        }
        if !seen.insert((i, j)) {
            return Err(input(format!("duplicate field basis term [{i}, {j}]")));
        }
        max_degree = max_degree.max(degree);
        // The catalog JSON lists monomials by (i+j, j). The prose sort key is (i+j, i);
        // ingest follows the JSON so the shipped default catalog is accepted. (C3.3)
        let key = (degree, j);
        if let Some(prev) = previous {
            if key <= prev {
                return Err(input(format!(
                    "field basis terms are not sorted by (i+j, j); saw [{i}, {j}] out of order"
                )));
            }
        }
        previous = Some(key);
    }
    if max_degree != basis.degree {
        return Err(input(format!(
            "field basis degree {} does not equal max(i+j)={max_degree}",
            basis.degree
        )));
    }
    Ok(())
}

fn validate_bounds(bounds: Option<[f64; 2]>) -> Result<(), PsfFieldError> {
    let Some([lo, hi]) = bounds else {
        return Ok(());
    };
    check_finite("bounds[0]", lo)?;
    check_finite("bounds[1]", hi)?;
    if lo >= hi {
        return Err(input("bounds lo must be < hi"));
    }
    Ok(())
}

fn prior_mean_is_init(prior: &PriorSpec) -> bool {
    matches!(
        prior,
        PriorSpec::Gaussian {
            mean: PriorMean::Init(_),
            ..
        }
    )
}

fn validate_prior(prior: &PriorSpec) -> Result<(), PsfFieldError> {
    match prior {
        PriorSpec::None { stage2 } => validate_stage2_prior(stage2.as_ref()),
        PriorSpec::Gaussian {
            mean,
            sigma,
            sigma_rel,
            stage2,
        } => {
            match mean {
                PriorMean::Number(mu) => {
                    check_finite("prior.mean", *mu)?;
                    let Some(sigma) = *sigma else {
                        return Err(input(
                            "gaussian prior with numeric mean requires sigma and must omit sigma_rel",
                        ));
                    };
                    if sigma_rel.is_some() {
                        return Err(input(
                            "gaussian prior with numeric mean requires sigma and must omit sigma_rel",
                        ));
                    }
                    check_positive("prior.sigma", sigma)?;
                }
                PriorMean::Init(_) => {
                    let Some(sigma_rel) = *sigma_rel else {
                        return Err(input(
                            "gaussian prior with mean \"init\" requires sigma_rel and must omit sigma",
                        ));
                    };
                    if sigma.is_some() {
                        return Err(input(
                            "gaussian prior with mean \"init\" requires sigma_rel and must omit sigma",
                        ));
                    }
                    check_positive("prior.sigma_rel", sigma_rel)?;
                }
            }
            validate_stage2_prior(stage2.as_ref())
        }
    }
}

fn validate_stage2_prior(stage2: Option<&Stage2Prior>) -> Result<(), PsfFieldError> {
    let Some(stage2) = stage2 else {
        return Ok(());
    };
    if stage2.mean.len() != stage2.sigma.len() {
        return Err(input(
            "stage2 prior mean and sigma must have the same length",
        ));
    }
    for (i, mu) in stage2.mean.iter().enumerate() {
        check_finite(&format!("stage2.mean[{i}]"), *mu)?;
    }
    for (i, sigma) in stage2.sigma.iter().enumerate() {
        check_positive(&format!("stage2.sigma[{i}]"), *sigma)?;
    }
    Ok(())
}

/// Each initialization method applies only to specific terms: `flux_sum` to flux,
/// `defocus_moment` to Zernike defocus `(n, m) = (2, 0)`, `moffat_fwhm` to a
/// Moffat kernel. Other pairings are rejected. (C3.5)
fn validate_init_pairing(term: &ErrorTerm) -> Result<(), PsfFieldError> {
    let method = term.init().method;
    let ok = match (method, term) {
        (InitMethod::Zero, _) => true,
        (InitMethod::FluxSum, ErrorTerm::Photometric { term_id, .. }) => term_id == "flux",
        (InitMethod::DefocusMoment, ErrorTerm::Phase { n: 2, m: 0, .. }) => true,
        (InitMethod::MoffatFwhm, ErrorTerm::Kernel { kernel, .. }) => {
            kernel.id == KernelId::MoffatIso
        }
        _ => false,
    };
    if ok {
        Ok(())
    } else {
        Err(input(format!(
            "init method {} is not valid for term {}",
            method.as_str(),
            term.term_id()
        )))
    }
}

#[cfg(test)]
pub(crate) const DEFAULT_CATALOG_JSON: &str =
    include_str!("../../../docs/contracts/schemas/psf_field_v1_default.catalog.json");

#[cfg(test)]
pub(crate) fn ingested_default_catalog() -> Catalog {
    let catalog: Catalog =
        serde_json::from_str(DEFAULT_CATALOG_JSON).expect("default catalog JSON parses");
    catalog.ingest().expect("default catalog ingest")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCode;
    use crate::ingest::ingest_catalog;
    use serde_json::{json, Value};

    fn default_value() -> Value {
        serde_json::from_str(DEFAULT_CATALOG_JSON).unwrap()
    }

    fn term_index(value: &Value, term_id: &str) -> usize {
        value["terms"]
            .as_array()
            .unwrap()
            .iter()
            .position(|term| term["term_id"] == term_id)
            .unwrap_or_else(|| panic!("missing {term_id}"))
    }

    fn ingest_err(value: &Value) -> PsfFieldError {
        ingest_catalog(value).expect_err("expected ingest rejection")
    }

    #[test]
    fn default_catalog_ingests() {
        let catalog = ingested_default_catalog();
        assert_eq!(catalog.catalog_id, "psf_field_v1_default");
        assert_eq!(catalog.fit_schedule.len(), 3);
    }

    #[test]
    fn ingest_rejects_semantic_violations() {
        struct Case {
            name: &'static str,
            patch: fn(&mut Value),
            needle: &'static str,
        }
        let cases = [
            Case {
                name: "duplicate term_id",
                patch: |v| {
                    let terms = v["terms"].as_array_mut().unwrap();
                    let clone = terms[0].clone();
                    terms.push(clone);
                },
                needle: "duplicate term_id",
            },
            Case {
                name: "n greater than 15",
                patch: |v| {
                    let i = term_index(v, "zernike_0_0");
                    v["terms"][i]["n"] = json!(16);
                },
                needle: "exceeds frozen max 15",
            },
            Case {
                name: "|m| greater than n",
                patch: |v| {
                    let i = term_index(v, "zernike_2_2");
                    v["terms"][i]["n"] = json!(1);
                    v["terms"][i]["m"] = json!(2);
                },
                needle: "|m| > n",
            },
            Case {
                name: "n minus |m| odd",
                patch: |v| {
                    let i = term_index(v, "zernike_2_2");
                    v["terms"][i]["n"] = json!(2);
                    v["terms"][i]["m"] = json!(1);
                },
                needle: "n - |m| is odd",
            },
            Case {
                name: "unsorted field basis",
                patch: |v| {
                    let i = term_index(v, "zernike_2_2");
                    v["terms"][i]["field_basis"]["terms"] = json!([[1, 0], [0, 0], [0, 1]]);
                },
                needle: "not sorted",
            },
            Case {
                name: "duplicate field basis term",
                patch: |v| {
                    let i = term_index(v, "zernike_4_0");
                    v["terms"][i]["field_basis"]["terms"] = json!([[0, 0], [0, 0]]);
                },
                needle: "duplicate field basis term",
            },
            Case {
                name: "inconsistent field degree",
                patch: |v| {
                    let i = term_index(v, "zernike_4_0");
                    v["terms"][i]["field_basis"]["degree"] = json!(1);
                },
                needle: "does not equal max(i+j)",
            },
            Case {
                name: "defocus_moment on non-defocus phase",
                patch: |v| {
                    let i = term_index(v, "zernike_2_2");
                    v["terms"][i]["init"]["method"] = json!("defocus_moment");
                },
                needle: "init method defocus_moment is not valid",
            },
            Case {
                name: "flux_sum on sky",
                patch: |v| {
                    let i = term_index(v, "sky");
                    v["terms"][i]["init"]["method"] = json!("flux_sum");
                },
                needle: "init method flux_sum is not valid",
            },
            Case {
                name: "moffat_fwhm on gaussian_iso",
                patch: |v| {
                    let i = term_index(v, "gaussian_jitter");
                    v["terms"][i]["init"]["method"] = json!("moffat_fwhm");
                },
                needle: "init method moffat_fwhm is not valid",
            },
            Case {
                name: "mean init with zero method",
                patch: |v| {
                    let i = term_index(v, "moffat_seeing");
                    v["terms"][i]["init"]["method"] = json!("zero");
                },
                needle: "mean \"init\" with init method zero",
            },
            Case {
                name: "photometric term_id not flux or sky",
                patch: |v| {
                    let i = term_index(v, "flux");
                    v["terms"][i]["term_id"] = json!("background");
                },
                needle: "not flux or sky",
            },
            Case {
                name: "gaussian numeric mean with sigma_rel",
                patch: |v| {
                    let i = term_index(v, "gaussian_jitter");
                    v["terms"][i]["prior"] = json!({
                        "family": "gaussian",
                        "mean": 0.1,
                        "sigma": 0.3,
                        "sigma_rel": 0.5
                    });
                },
                needle: "must omit sigma_rel",
            },
        ];

        for case in cases {
            let mut value = default_value();
            (case.patch)(&mut value);
            let err = ingest_err(&value);
            assert_eq!(err.code, ErrorCode::Input, "{}", case.name);
            assert!(
                err.message.contains(case.needle),
                "{}: expected {:?} in {}",
                case.name,
                case.needle,
                err.message
            );
        }
    }

    #[test]
    fn reserved_field_family_rejected() {
        let mut value = default_value();
        let i = term_index(&value, "zernike_0_0");
        value["terms"][i]["field_basis"]["family"] = json!("zernike_field");
        let err = ingest_err(&value);
        assert_eq!(err.code, ErrorCode::Input);
    }

    #[test]
    fn unknown_kernel_id_rejected() {
        let mut value = default_value();
        let i = term_index(&value, "gaussian_jitter");
        value["terms"][i]["kernel"]["id"] = json!("seidel_blob");
        let err = ingest_err(&value);
        assert_eq!(err.code, ErrorCode::Input);
    }
}
