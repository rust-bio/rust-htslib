use std::path::PathBuf;
use thiserror::Error;

/// Generic result type for functions in this crate with
/// a global error class.
pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Error, Debug, PartialEq)]
pub enum Error {
    // General errors
    #[error("file not found: {path}")]
    FileNotFound { path: PathBuf },
    #[error("invalid (non-unicode) characters in path")]
    NonUnicodePath,
    #[error("failed to fetch region")]
    Fetch,
    #[error("error seeking to file offset")]
    FileSeek,
    #[error("error seeking to {contig:?}:{start} in indexed file")]
    GenomicSeek { contig: String, start: u64 },
    #[error("sequence {sequence} not found in index")]
    UnknownSequence { sequence: String },
    #[error("error setting threads for file reading")]
    SetThreads,
    #[error("failed to create htslib thread pool")]
    ThreadPool,

    #[error("failed to write BAM/BCF record (out of disk space?)")]
    WriteRecord,

    // Errors for faidx
    #[error("The given position is too large to be converted to i64")]
    PositionTooLarge,

    // Errors for Tbx
    #[error("previous iterator generation failed")]
    NoIter,
    #[error("truncated tabix record")]
    TabixTruncatedRecord,
    #[error("invalid tabix index")]
    TabixInvalidIndex,

    // Errors for BAM
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
    #[error("invalid record in SAM/BAM/CRAM file")]
    BamInvalidRecord,
    #[error("truncated record in SAM/BAM/CRAM file")]
    BamTruncatedRecord,
    #[error("format not indexable by htslib (format is detected as something else than SAM/BAM/CRAM)")]
    NotIndexable,
    #[error("failed to write BAM/CRAM index (out of disk space?)")]
    WriteIndex,
    #[error("failed to build BAM/CRAM index")]
    BuildIndex,
    #[error("failed to create SAM/BAM/CRAM pileup")]
    Pileup,

    // Errors for BCF
    #[error("error allocating internal data structure for BCF/VCF reader (out of memory?)")]
    AllocationError,
    #[error("failed to open BCF/VCF from {target:?}")]
    Open { target: String },
    #[error("invalid record in BCF/VCF file")]
    InvalidRecord,
    #[error("tag {tag} undefined in BCF/VCF header")]
    UndefinedTag { tag: String },
    #[error("unexpected type for tag {tag} in BCF/VCF file")]
    UnexpectedType { tag: String },
    #[error("tag {tag} missing from record {record} in BCF/VCF file")]
    MissingTag { tag: String, record: String },
    #[error("error setting tag {tag} in BCF/VCF record (out of memory?)")]
    SetTag { tag: String },
    #[error("ID {rid} not found in BCF/VCF header")]
    UnknownRID { rid: u32 },
    #[error("contig {contig} not found in BCF/VCF header")]
    UnknownContig { contig: String },
    #[error("ID {id} not found in BCF/VCF header")]
    UnknownID { id: String },
    #[error("sample {name} not found in BCF/VCF header")]
    UnknownSample { name: String },
    #[error("duplicate sample names given for subsetting BCF/VCF")]
    DuplicateSampleNames,
    #[error("failed to set values in BCF/VCF record (out of memory?)")]
    SetValues,
    #[error("failed to remove alleles in BCF/VCF record")]
    RemoveAlleles,
}
