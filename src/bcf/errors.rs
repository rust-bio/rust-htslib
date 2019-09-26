use snafu::Snafu;

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Snafu, Debug)]
#[snafu(visibility = "pub")]
pub enum Error {
    #[snafu(display(
        "Error allocating internal data structure for BCF/VCF reader (out of memory?)."
    ))]
    AllocationError,
    #[snafu(display("Error opening BCF/VCF from {}.", target))]
    Open { target: String },
    #[snafu(display("Invalid record in BCF/VCF file."))]
    InvalidRecord,
    #[snafu(display("Error setting threads for writing BCF/VCF file(s)."))]
    SetThreads,
    #[snafu(display("Error seeking to {}:{} in indexed BCF/VCF file.", contig, start))]
    Seek { contig: String, start: u32 },
    #[snafu(display("Error writing record to BCF/VCF file."))]
    Write,
    #[snafu(display("Tag {} undefined in BCF/VCF header.", tag))]
    UndefinedTag { tag: String },
    #[snafu(display("Unexpected type for tag {} BCF/VCF file.", tag))]
    UnexpectedType { tag: String },
    #[snafu(display("Tag {} missing from record {} in BCF/VCF file.", tag, record))]
    MissingTag { tag: String, record: String },
    #[snafu(display("Error setting tag {} in BCF/VCF record (out of memory?).", tag))]
    SetTag { tag: String },
    #[snafu(display("ID {} not found in BCF/VCF header.", rid))]
    UnknownRID { rid: u32 },
    #[snafu(display("Contig {} not found in BCF/VCF header.", contig))]
    UnknownContig { contig: String },
    #[snafu(display("ID {} not found in BCF/VCF header.", id))]
    UnknownID { id: String },
    #[snafu(display("Sample {} not found in BCF/VCF header.", name))]
    UnknownSample { name: String },
    #[snafu(display("Duplicate sample names given for subsetting BCF/VCF."))]
    DuplicateSampleNames,
    #[snafu(display("invalid (non-unique) characters in path"))]
    NonUnicodePath,
    #[snafu(display("failed to set values in BCF/VCF record (out of memory?)"))]
    SetValues,
    #[snafu(display("failed to remove alleles in BCF/VCF record"))]
    RemoveAlleles,
}
