//! Typed errors that map to the Python hierarchy in C1B.6.

use std::fmt;

/// Short `code` string carried on the Python exception.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorCode {
    Input,
    Convergence,
    Numerics,
    Internal,
}

impl ErrorCode {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Input => "input",
            Self::Convergence => "convergence",
            Self::Numerics => "numerics",
            Self::Internal => "internal",
        }
    }
}

impl fmt::Display for ErrorCode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// `module` string carried on the Python exception (C1B.6).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorModule {
    Zernike,
    Pipeline,
    Lm,
    Stage2,
    Eval,
    Boundary,
}

impl ErrorModule {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Zernike => "zernike",
            Self::Pipeline => "pipeline",
            Self::Lm => "lm",
            Self::Stage2 => "stage2",
            Self::Eval => "eval",
            Self::Boundary => "boundary",
        }
    }
}

impl fmt::Display for ErrorModule {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Core error type. Panics never cross the PyO3 boundary (C1B.6).
#[derive(Debug, Clone, thiserror::Error)]
#[error("{message}")]
pub struct PsfFieldError {
    pub code: ErrorCode,
    pub module: ErrorModule,
    pub message: String,
    pub star_id: Option<String>,
}

impl PsfFieldError {
    #[must_use]
    pub fn new(
        code: ErrorCode,
        module: ErrorModule,
        message: impl Into<String>,
        star_id: Option<String>,
    ) -> Self {
        Self {
            code,
            module,
            message: message.into(),
            star_id,
        }
    }

    #[must_use]
    pub fn input(module: ErrorModule, message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Input, module, message, None)
    }

    #[must_use]
    pub fn input_star(
        module: ErrorModule,
        message: impl Into<String>,
        star_id: impl Into<String>,
    ) -> Self {
        Self::new(ErrorCode::Input, module, message, Some(star_id.into()))
    }

    #[must_use]
    pub fn numerics(module: ErrorModule, message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Numerics, module, message, None)
    }

    #[must_use]
    pub fn convergence(module: ErrorModule, message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Convergence, module, message, None)
    }

    #[must_use]
    pub fn internal(module: ErrorModule, message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Internal, module, message, None)
    }

    #[must_use]
    pub fn with_star_id(mut self, star_id: impl Into<String>) -> Self {
        self.star_id = Some(star_id.into());
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_match_c1b6() {
        assert_eq!(ErrorCode::Input.as_str(), "input");
        assert_eq!(ErrorCode::Convergence.as_str(), "convergence");
        assert_eq!(ErrorCode::Numerics.as_str(), "numerics");
        assert_eq!(ErrorCode::Internal.as_str(), "internal");
    }

    #[test]
    fn modules_match_c1b6() {
        assert_eq!(ErrorModule::Zernike.as_str(), "zernike");
        assert_eq!(ErrorModule::Pipeline.as_str(), "pipeline");
        assert_eq!(ErrorModule::Lm.as_str(), "lm");
        assert_eq!(ErrorModule::Stage2.as_str(), "stage2");
        assert_eq!(ErrorModule::Eval.as_str(), "eval");
        assert_eq!(ErrorModule::Boundary.as_str(), "boundary");
    }
}
