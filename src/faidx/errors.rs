use snafu::Snafu;
use std::path::PathBuf;

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Snafu, Debug, PartialEq)]
#[snafu(visibility = "pub")]
pub enum Error {
    #[snafu(display("file not found: {}", path.display()))]
    FileNotFound { path: PathBuf },
    #[snafu(display("invalid (non-unique) characters in path"))]
    NonUnicodePath,
}
