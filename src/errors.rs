use std::path::PathBuf;
use thiserror::Error;

/// Generic result type for functions in this crate with
/// a global error class.
pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Error, Debug, PartialEq)]
pub enum Error {
    #[error("file not found: {path}")]
    FileNotFound { path: PathBuf },
    #[error("invalid (non-unicode) characters in path")]
    NonUnicodePath,
    #[error("The given position is too large to be converted to i64")]
    PositionTooLarge,
    #[error("sequence {sequence} not found in index")]
    UnknownSequence { sequence: String },
    #[error("previous iterator generation failed")]
    NoIter,
    #[error("truncated tabix record")]
    TabixTruncatedRecord,
    #[error("invTabixalid tabix index")]
    TabixInvalidIndex,
    #[error("failed to fetch region in index")]
    Fetch,
    #[error("error setting threads for for file reading")]
    SetThreads,
}
