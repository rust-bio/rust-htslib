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

    #[error("failed to fetch region")]
    Fetch,
    #[error("error seeking to offset")]
    Seek,
    #[error("sequence {sequence} not found in index")]
    UnknownSequence { sequence: String },

    #[error("previous iterator generation failed")]
    NoIter,
    #[error("truncated tabix record")]
    TabixTruncatedRecord,
    #[error("invTabixalid tabix index")]
    TabixInvalidIndex,
    #[error("error setting threads for for file reading")]
    SetThreads,

    #[error("error parsing CIGAR string: {msg}")]
    ParseCigar { msg: String },
    #[error("unexpected CIGAR operation: {msg}")]
    UnexpectedCigarOperation { msg: String },
    #[error("error parsing SAM record: {rec}")]
    ParseSAM { rec: String },
    #[error("invalid reference path {path}")]
    InvalidReferencePath { path: PathBuf },
    #[error("invalid compression level {level}")]
    InvalidCompressionLevel { level: u32 },
    #[error("unable to open SAM/BAM/CRAM file at {target}")]
    BamOpen { target: String },
    #[error("unable to open SAM/BAM/CRAM index for {target}")]
    BamInvalidIndex { target: String },
    #[error("failed to write record (out of disk space?)")]
    Write,
    #[error("invalid record in SAM/BAM/CRAM file")]
    BamInvalidRecord,
    #[error("truncated record in SAM/BAM/CRAM file")]
    BamTruncatedRecord,

    #[error("format of SAM files are not indexable")]
    NotIndexable,
    #[error("failed to write BAM/CRAM index (out of disk space?)")]
    WriteIndex,
    #[error("failed to build BAM/CRAM index")]
    BuildIndex,
    #[error("failed to create SAM/BAM/CRAM pileup")]
    Pileup,
}
