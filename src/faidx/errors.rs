use thiserror::Error;
use std::path::PathBuf;

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Error, Debug, PartialEq)]
pub enum Error {
    #[error("file not found: {path}")]
    FileNotFound { path: PathBuf },
    #[error("invalid (non-unicode) characters in path")]
    NonUnicodePath,
}
