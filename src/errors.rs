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
    FaidxPositionTooLarge,
    #[error("bad conversion of sequence name")]
    FaidxBadSeqName,

    // Errors for Tbx
    #[error("previous iterator generation failed")]
    TabixNoIter,
    #[error("truncated tabix record")]
    TabixTruncatedRecord,
    #[error("invalid tabix index")]
    TabixInvalidIndex,

    // Errors for BAM
    #[error("error parsing CIGAR string: {msg}")]
    BamParseCigar { msg: String },
    #[error("unexpected CIGAR operation: {msg}")]
    BamUnexpectedCigarOperation { msg: String },
    #[error("error parsing SAM record: {rec}")]
    BamParseSAM { rec: String },
    #[error("invalid path to CRAM-reference {path}")]
    BamInvalidReferencePath { path: PathBuf },
    #[error("invalid compression level {level}")]
    BamInvalidCompressionLevel { level: u32 },
    #[error("unable to open SAM/BAM/CRAM file at {target}")]
    BamOpen { target: String },
    #[error("unable to open SAM/BAM/CRAM index for {target}; please create an index")]
    BamInvalidIndex { target: String },
    #[error("invalid record in SAM/BAM/CRAM file")]
    BamInvalidRecord,
    #[error("truncated record in SAM/BAM/CRAM file")]
    BamTruncatedRecord,
    #[error(
        "format not indexable by htslib (format is detected as something else than SAM/BAM/CRAM)"
    )]
    BamNotIndexable,
    #[error("failed to write BAM/CRAM index (out of disk space?)")]
    BamWriteIndex,
    #[error("failed to build BAM/CRAM index")]
    BamBuildIndex,
    #[error("failed to create SAM/BAM/CRAM pileup")]
    BamPileup,
    #[error("file is not sorted by position")]
    BamUnsorted,

    // Errors for BAM auxiliary fields
    #[error("failed to add aux field (out of memory?)")]
    BamAux,
    #[error("provided string contains internal 0 byte(s)")]
    BamAuxStringError,
    #[error("failed to parse aux data")]
    BamAuxParsingError,
    #[error("the specified tag does could not be found")]
    BamAuxTagNotFound,
    #[error("data type of aux field is not known")]
    BamAuxUnknownType,
    #[error("failed to add aux field, tag is already present")]
    BamAuxTagAlreadyPresent,

    // Errors for base modification fields
    #[error("no base modification tag found for record")]
    BamBaseModificationTagNotFound,
    #[error("no base modification with the specified code found in record")]
    BamBaseModificationTypeNotFound,
    #[error("base modification iteration failed")]
    BamBaseModificationIterationFailed,
    #[error("base modification found too many modifications")]
    BamBaseModificationTooManyMods,

    // Errors for BCF
    #[error("error allocating internal data structure for BCF/VCF reader (out of memory?)")]
    BcfAllocationError,
    #[error("failed to open BCF/VCF from {target:?}")]
    BcfOpen { target: String },
    #[error("invalid record in BCF/VCF file")]
    BcfInvalidRecord,
    #[error("tag {tag} undefined in BCF/VCF header")]
    BcfUndefinedTag { tag: String },
    #[error("unexpected type for tag {tag} in BCF/VCF file")]
    BcfUnexpectedType { tag: String },
    #[error("tag {tag} missing from record {record} in BCF/VCF file")]
    BcfMissingTag { tag: String, record: String },
    #[error("error setting tag {tag} in BCF/VCF record (out of memory?)")]
    BcfSetTag { tag: String },
    #[error("ID {rid} not found in BCF/VCF header")]
    BcfUnknownRID { rid: u32 },
    #[error("contig {contig} not found in BCF/VCF header")]
    BcfUnknownContig { contig: String },
    #[error("ID {id} not found in BCF/VCF header")]
    BcfUnknownID { id: String },
    #[error("sample {name} not found in BCF/VCF header")]
    BcfUnknownSample { name: String },
    #[error("duplicate sample names given for subsetting BCF/VCF")]
    BcfDuplicateSampleNames,
    #[error("failed to set values in BCF/VCF record (out of memory?)")]
    BcfSetValues,
    #[error("failed to remove alleles in BCF/VCF record")]
    BcfRemoveAlleles,

    #[error("invalid compression level {level}")]
    BgzfInvalidCompressionLevel { level: i8 },
    #[error("failed setting hts reading options")]
    HtsSetOpt,
    #[error("failed calculating slow index statistics")]
    SlowIdxStats,
    #[error("invalid tid {tid}")]
    InvalidTid { tid: i32 },
    #[error("No sequences in the reference")]
    NoSequencesInReference,
}
