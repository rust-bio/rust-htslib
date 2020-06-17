use std::path::PathBuf;
use thiserror::Error;

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Error, Debug)]
pub enum Error {
    #[error("previous iterator generation failed")]
    NoIter,
    #[error("truncated tabix record")]
    TruncatedRecord,
    #[error("invalid tabix index")]
    InvalidIndex,
    #[error("file not found: {path}")]
    FileNotFound { path: PathBuf },
    #[error("invalid (non-unique) characters in path")]
    NonUnicodePath,
    #[error("sequence {sequence} not found in tabix index")]
    UnknownSequence { sequence: String },
    #[error("failed to fetch region in tabix index")]
    Fetch,
    #[error("error setting threads for for tabix file reading")]
    SetThreads,
}
