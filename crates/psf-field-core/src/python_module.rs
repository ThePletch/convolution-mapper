use pyo3::prelude::*;

/// Compiled submodule `psf_field_core._core` (C1B.2). Numeric entry points land in later PRs.
#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
