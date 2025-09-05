use thiserror::Error;

/// Crate-wide result type alias.
pub type Result<T> = std::result::Result<T, Error>;

/// Comprehensive error type for the simulation core.
///
/// This error enum is used throughout the application to avoid `.unwrap()`/`.expect()`
/// in accordance with AGENT.md. Each variant carries enough context to be actionable.
#[derive(Debug, Error)]
pub enum Error {
    /// Invalid user or API parameter.
    #[error("invalid parameter: {0}")]
    InvalidParam(String),

    /// Numerical or geometric issue (e.g., degenerate quadratic, near-zero normal, catastrophic cancellation).
    #[error("numerical error: {0}")]
    MathError(String),

    /// No collision detected for a given pair within tolerances.
    #[error("no collision detected within tolerance")]
    NoCollision,

    /// State drift or positions out of bounds after advancement.
    #[error("out of bounds: {0}")]
    OutOfBounds(String),

    /// Functionality not yet implemented (used to stage phases; should be eliminated over time).
    #[error("not implemented: {0}")]
    NotImplemented(&'static str),

    /// Propagated I/O errors (e.g., for future data exports).
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_display_is_informative() {
        let e = Error::InvalidParam("radius must be > 0".to_string());
        let msg = format!("{e}");
        assert!(msg.contains("invalid parameter"));
        assert!(msg.contains("radius"));
    }

    #[test]
    fn result_type_alias_compiles() -> Result<()> {
        // Simple smoke test for the alias
        Ok(())
    }
}
