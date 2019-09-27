use snafu::Snafu;
use std::path::PathBuf;

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Snafu, Debug, PartialEq)]
#[snafu(visibility = "pub")]
pub enum Error {
    #[snafu(display("error parsing CIGAR string: {}", msg))]
    ParseCigar { msg: String },
    #[snafu(display("unexpected CIGAR operation: {}", msg))]
    UnexpectedCigarOperation { msg: String },
    #[snafu(display("error parsing SAM record: {}", rec))]
    ParseSAM { rec: String },
    #[snafu(display("error setting threads for writing SAM/BAM/CRAM file(s)"))]
    SetThreads,
    #[snafu(display("invalid reference path {}", path.display()))]
    InvalidReferencePath { path: PathBuf },
    #[snafu(display("invalid compression level {}", level))]
    InvalidCompressionLevel { level: u32 },
    #[snafu(display("file not found: {}", path.display()))]
    FileNotFound { path: PathBuf },
    #[snafu(display("invalid (non-unique) characters in path"))]
    NonUnicodePath,
    #[snafu(display("unable to open SAM/BAM/CRAM file at {}", target))]
    Open { target: String },
    #[snafu(display("unable to open SAM/BAM/CRAM index for {}", target))]
    InvalidIndex { target: String },
    #[snafu(display("failed to write SAM/BAM/CRAM record (out of disk space?)"))]
    Write,
    #[snafu(display("invalid record in SAM/BAM/CRAM file"))]
    InvalidRecord,
    #[snafu(display("truncated record in SAM/BAM/CRAM file"))]
    TruncatedRecord,
    #[snafu(display("error fetching region in SAM/BAM/CRAM file"))]
    Fetch,
    #[snafu(display("error seeking to offset in SAM/BAM/CRAM file"))]
    Seek,
    #[snafu(display("sequence {} not found in SAM/BAM/CRAM file header", sequence))]
    UnknownSequence { sequence: String },
    #[snafu(display("format of SAM files are not indexable"))]
    NotIndexable,
    #[snafu(display("failed to write BAM/CRAM index (out of disk space?)"))]
    WriteIndex,
    #[snafu(display("failed to build BAM/CRAM index"))]
    BuildIndex,
    #[snafu(display("failed to create SAM/BAM/CRAM pileup"))]
    Pileup,
}
