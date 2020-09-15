use std::path::PathBuf;

use thiserror::Error;

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Error, Debug)]
pub enum Error {
    #[error("error allocating internal data structure for BCF/VCF reader (out of memory?)")]
    AllocationError,
    #[error("failed to open BCF/VCF from {target:?}")]
    Open { target: String },
    #[error("invalid record in BCF/VCF file")]
    InvalidRecord,
    #[error("error setting threads for writing BCF/VCF file(s)")]
    SetThreads,
    #[error("error seeking to {contig:?}:{start} in indexed BCF/VCF file")]
    Seek { contig: String, start: u64 },
    #[error("error writing record to BCF/VCF file")]
    Write,
    #[error("tag {tag} undefined in BCF/VCF header")]
    UndefinedTag { tag: String },
    #[error("unexpected type for tag {tag} BCF/VCF file")]
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
    #[error("invalid (non-unique) characters in path")]
    NonUnicodePath,
    #[error("file not found: {path}")]
    FileNotFound { path: String },
    #[error("failed to set values in BCF/VCF record (out of memory?)")]
    SetValues,
    #[error("failed to remove alleles in BCF/VCF record")]
    RemoveAlleles,
}
