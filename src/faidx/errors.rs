use std::path::PathBuf;
use thiserror::Error;

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Error, Debug, PartialEq)]
pub enum Error {
    #[error("file not found: {path}")]
    FileNotFound { path: PathBuf },
    #[error("invalid (non-unicode) characters in path")]
    NonUnicodePath,
    #[error("The given position is to large to be converted to i64")]
    PositionTooLarge,
}
