//! Shared identifiers, schema versioning (NOR.2), and small value types.

use serde::{Deserialize, Serialize};

use crate::error::{ErrorModule, PsfFieldError};

/// Frozen v1 `schema_version` (NOR.2).
pub const SCHEMA_VERSION: &str = "1.0.0";

/// Frozen v1 MAJOR. Unknown MAJOR is rejected at ingest (NOR.2).
pub const SCHEMA_MAJOR: u64 = 1;

const ID_MAX: usize = 128;

/// Reject unknown MAJOR and any version other than the frozen v1 string.
pub fn check_schema_version(version: &str) -> Result<(), PsfFieldError> {
    let mut parts = version.split('.');
    let (Some(major_s), Some(minor_s), Some(patch_s), None) =
        (parts.next(), parts.next(), parts.next(), parts.next())
    else {
        return Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("invalid schema_version {version:?}"),
        ));
    };
    if !is_uint(major_s) || !is_uint(minor_s) || !is_uint(patch_s) {
        return Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("invalid schema_version {version:?}"),
        ));
    }
    let major: u64 = major_s.parse().expect("digits");
    if major != SCHEMA_MAJOR {
        return Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("unknown schema MAJOR {major}"),
        ));
    }
    if version != SCHEMA_VERSION {
        return Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("unsupported schema_version {version:?}"),
        ));
    }
    Ok(())
}

fn is_uint(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

/// NOR.5 `star_id` / `exposure_id` / `session_id`.
pub fn check_id(kind: &str, value: &str) -> Result<(), PsfFieldError> {
    if value.is_empty() || value.len() > ID_MAX {
        return Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("{kind} length must be in 1..=128"),
        ));
    }
    let ok = value
        .bytes()
        .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'.' | b'_' | b':' | b'-'));
    if !ok {
        return Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("{kind} {value:?} does not match NOR.5"),
        ));
    }
    Ok(())
}

/// NOR.5 `term_id` / `bundle_id`.
pub fn check_term_id(kind: &str, value: &str) -> Result<(), PsfFieldError> {
    let mut chars = value.chars();
    let Some(first) = chars.next() else {
        return Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("{kind} is empty"),
        ));
    };
    if !first.is_ascii_lowercase()
        || !chars.all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_')
    {
        return Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("{kind} {value:?} does not match NOR.5"),
        ));
    }
    Ok(())
}

/// C1.2.1 stamp size.
pub fn check_stamp_size(s: usize) -> Result<(), PsfFieldError> {
    if s % 2 == 0 {
        return Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("even stamp size {s} is rejected"),
        ));
    }
    if !(15..=63).contains(&s) {
        return Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("stamp size {s} is outside the allowed odd set [15, 63]"),
        ));
    }
    Ok(())
}

pub fn check_finite(name: &str, x: f64) -> Result<(), PsfFieldError> {
    if x.is_finite() {
        Ok(())
    } else {
        Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("{name} must be finite"),
        ))
    }
}

pub fn check_positive(name: &str, x: f64) -> Result<(), PsfFieldError> {
    check_finite(name, x)?;
    if x > 0.0 {
        Ok(())
    } else {
        Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("{name} must be > 0"),
        ))
    }
}

pub fn check_non_negative(name: &str, x: f64) -> Result<(), PsfFieldError> {
    check_finite(name, x)?;
    if x >= 0.0 {
        Ok(())
    } else {
        Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("{name} must be >= 0"),
        ))
    }
}

pub fn check_vec2(name: &str, v: [f64; 2]) -> Result<(), PsfFieldError> {
    check_finite(&format!("{name}[0]"), v[0])?;
    check_finite(&format!("{name}[1]"), v[1])?;
    Ok(())
}

/// C1.7 closed flag vocabulary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum Flag {
    Saturated,
    Blended,
    Edge,
    Shape,
    Underdetermined,
    UserExclude,
    Selected,
}

/// C4.7 / C5.7 parameter annotation sidecar.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ParamMeta {
    pub term_id: String,
    pub role: String,
    pub scope: super::catalog::Scope,
    pub frozen: bool,
    pub unit: String,
}

/// C5.6 degenerate pair serialized as `[term_a, term_b, rho]`.
#[derive(Debug, Clone, PartialEq)]
pub struct DegeneratePair {
    pub term_a: String,
    pub term_b: String,
    pub rho: f64,
}

impl Serialize for DegeneratePair {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeSeq;
        let mut seq = serializer.serialize_seq(Some(3))?;
        seq.serialize_element(&self.term_a)?;
        seq.serialize_element(&self.term_b)?;
        seq.serialize_element(&self.rho)?;
        seq.end()
    }
}

impl<'de> Deserialize<'de> for DegeneratePair {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let raw = <(String, String, f64)>::deserialize(deserializer)?;
        Ok(Self {
            term_a: raw.0,
            term_b: raw.1,
            rho: raw.2,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_unknown_major() {
        let err = check_schema_version("2.0.0").unwrap_err();
        assert!(err.message.contains("unknown schema MAJOR 2"));
        assert_eq!(err.code, crate::error::ErrorCode::Input);
    }

    #[test]
    fn accepts_frozen_v1() {
        check_schema_version(SCHEMA_VERSION).unwrap();
    }
}
