//! Synthetic corpus generator and scorer (C10). Implemented in a later PR.

#[cfg(test)]
mod tests {
    #[test]
    fn crate_links() {
        assert_eq!(env!("CARGO_PKG_NAME"), "psf-field-corpus");
    }
}
