use snafu::Snafu;
use std::path::PathBuf;

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Snafu, Debug, PartialEq)]
#[snafu(visibility = "pub")]
pub enum Error {
    #[snafu(display("previous iterator generation failed"))]
    NoIter,
    #[snafu(display("truncated tabix record"))]
    TruncatedRecord,
    #[snafu(display("invalid tabix index"))]
    InvalidIndex,
    #[snafu(display("file not found: {}", path.display()))]
    FileNotFound { path: PathBuf },
    #[snafu(display("invalid (non-unique) characters in path"))]
    NonUnicodePath,
    #[snafu(display("sequence {} not found in tabix index", sequence))]
    UnknownSequence { sequence: String },
    #[snafu(display("failed to fetch region in tabix index"))]
    Fetch,
    #[snafu(display("error setting threads for for tabix file reading"))]
    SetThreads,
}
