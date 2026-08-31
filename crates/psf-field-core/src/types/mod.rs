//! serde types 1:1 with `docs/contracts/schemas/` (C1.9.2).

mod catalog;
mod common;
mod inputs;
mod outputs;

pub use catalog::{
    Bundle, Catalog, ErrorTerm, FieldBasis, FitScheduleStep, InitMethod, InitSpec, KernelId,
    KernelSpec, PriorMean, PriorSpec, Scope, Stage2Prior,
};
pub use common::{
    check_id, check_schema_version, check_stamp_size, check_term_id, DegeneratePair, Flag,
    ParamMeta, SCHEMA_MAJOR, SCHEMA_VERSION,
};
pub use inputs::{
    ExtractionConfig, ImageMeta, PixelMaskBits, PupilSpec, Stage1InputConfig, StarRecord,
};
pub use outputs::{
    CoverageReport, FdReport, FieldMap, PsfEval, ScoreEntry, ScoreReport, Stage1Result,
    Stage2Result, StarScore, Termination,
};
