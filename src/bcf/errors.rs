use snafu::Snafu;

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Snafu, Debug, PartialEq)]
#[snafu(visibility = "pub")]
pub enum Error {
    #[snafu(display(
        "error allocating internal data structure for BCF/VCF reader (out of memory?)"
    ))]
    AllocationError,
    #[snafu(display("failed to open BCF/VCF from {}", target))]
    Open { target: String },
    #[snafu(display("invalid record in BCF/VCF file"))]
    InvalidRecord,
    #[snafu(display("error setting threads for writing BCF/VCF file(s)"))]
    SetThreads,
    #[snafu(display("error seeking to {}:{} in indexed BCF/VCF file", contig, start))]
    Seek { contig: String, start: u64 },
    #[snafu(display("error writing record to BCF/VCF file"))]
    Write,
    #[snafu(display("tag {} undefined in BCF/VCF header", tag))]
    UndefinedTag { tag: String },
    #[snafu(display("unexpected type for tag {} BCF/VCF file", tag))]
    UnexpectedType { tag: String },
    #[snafu(display("tag {} missing from record {} in BCF/VCF file", tag, record))]
    MissingTag { tag: String, record: String },
    #[snafu(display("error setting tag {} in BCF/VCF record (out of memory?)", tag))]
    SetTag { tag: String },
    #[snafu(display("ID {} not found in BCF/VCF header", rid))]
    UnknownRID { rid: u32 },
    #[snafu(display("contig {} not found in BCF/VCF header", contig))]
    UnknownContig { contig: String },
    #[snafu(display("ID {} not found in BCF/VCF header", id))]
    UnknownID { id: String },
    #[snafu(display("sample {} not found in BCF/VCF header", name))]
    UnknownSample { name: String },
    #[snafu(display("duplicate sample names given for subsetting BCF/VCF"))]
    DuplicateSampleNames,
    #[snafu(display("invalid (non-unique) characters in path"))]
    NonUnicodePath,
    #[snafu(display("failed to set values in BCF/VCF record (out of memory?)"))]
    SetValues,
    #[snafu(display("failed to remove alleles in BCF/VCF record"))]
    RemoveAlleles,
}
