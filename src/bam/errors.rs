use std::path::PathBuf;

use thiserror::Error;

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Error, Debug)]
pub enum Error {
    #[error("error parsing CIGAR string: {msg:?}")]
    ParseCigar { msg: String },
    #[error("unexpected CIGAR operation: {msg:?}")]
    UnexpectedCigarOperation { msg: String },
    #[error("error parsing SAM record: {rec:?}")]
    ParseSAM { rec: String },
    #[error("error setting threads for writing SAM/BAM/CRAM file(s)")]
    SetThreads,
    #[error("invalid reference path {path:?}")]
    InvalidReferencePath { path: PathBuf },
    #[error("invalid compression level {level:?}")]
    InvalidCompressionLevel { level: u32 },
    #[error("file not found: {path:?}")]
    FileNotFound { path: PathBuf },
    #[error("invalid (non-unique) characters in path")]
    NonUnicodePath,
    #[error("unable to open SAM/BAM/CRAM file at {target:?}")]
    Open { target: String },
    #[error("unable to open SAM/BAM/CRAM index for {target:?}")]
    InvalidIndex { target: String },
    #[error("failed to write SAM/BAM/CRAM record (out of disk space?)")]
    Write,
    #[error("invalid record in SAM/BAM/CRAM file")]
    InvalidRecord,
    #[error("truncated record in SAM/BAM/CRAM file")]
    TruncatedRecord,
    #[error("error fetching region in SAM/BAM/CRAM file")]
    Fetch,
    #[error("error seeking to offset in SAM/BAM/CRAM file")]
    Seek,
    #[error("sequence {sequence:?} not found in SAM/BAM/CRAM file header")]
    UnknownSequence { sequence: String },
    #[error("format of SAM files are not indexable")]
    NotIndexable,
    #[error("failed to write BAM/CRAM index (out of disk space?)")]
    WriteIndex,
    #[error("failed to build BAM/CRAM index")]
    BuildIndex,
    #[error("failed to create SAM/BAM/CRAM pileup")]
    Pileup,
}
