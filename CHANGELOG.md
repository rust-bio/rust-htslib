# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.35.1] - 2020-11-23
### Changes
- Fixed wrongly define missing value constants in bcf::record (@johanneskoester).
- Bump hts-sys depedency to the latest version, containing build fixes for macOS (@johanneskoester).


## [0.35.0] - 2020-11-19
### Changes
- BREAKING: info and format field access in BCF records now allocates a separate buffer each time. In addition, it is also possible to pass a buffer that has been created manually before (@johanneskoester)
- Fixes for building on macOS (@brainstorm)

### Added
- ability to push genotypes into BCF record (@MaltheSR, @tedil).

## [0.34.0] - 2020-11-13
### Added
- Ability to set minimum refetch distance in `bam::RecordBuffer`.

## [0.33.0] - 2020-11-04
### Changes
- BREAKING: Rename feature 'serde' as 'serde_feature' (for technical reasons)
- BREAKING: Consolidate module-wide errors into a crate-wide error module
- Making `bcf::IndexedReader` always unpack records to reflect the behaviour of `bcf::Reader`.
- Adding `bcf::errors::Error::FileNotFound` and using it.
- Fixes for musl compilation (@brainstorm).
- Improved BCF constants handling (@juliangehring)
- Fixes for tabix reader (@felix-clark, @brainstorm).
- Fixes for BCF handling (@holtgrewe, @tedil).
- Documentation improvements (@vsoch, @brainstorm, @edmundlth).
- BREAKING: Improved, more ergonomic BAM fetch API (@TiberiusPrime, @brainstorm, @tedil).
- BREAKING: Let BamRecordExtensions return iterators instead of vectors (@TiberiusPrime).
- Handle all errors via a unified single thiserror based enum (@landesfeind).
- BREAKING: BAM read API now returns Option<Result> (@slazicoicr).
### Added
- Support for reading indexed FASTA files (@landesfeind, @pmarks, @brainstorm).
- Support for shared threadpools when reading and writing BAM (@pmarks, @nlhepler).
- Serde support for Cigar strings (@FelixMoelder, @pmarks, @johanneskoester).
- Expose bgzf functionality (@landesfeind).
- Iterator over BAM records using Rc-pointers (@TiberiusPrime, @tedil).
- Ability to obtain pairs of read and genome intervals from BAM (aligned_block_pairs) (@TiberiusPrime, @brainstorm).

## [0.32.0] - 2020-07-09
### Changes
- Method `seq_len()` of `bam::Record` is now public.
- Speedup when parsing BAM headers (thanks to @juliangehring).
- Compatibility fixes for older rust versions (thanks to @pmarks and @brainstorm).

## [0.31.0] - 2020-06-22
### Changes
- Bam record buffer now returns reference counted (Rc) objects. This makes the API more ergonomic to use.
- Switched to thiserror instead of snafu for error handling.
- Various cleanups and little fixes.

## [0.30.0] - 2020-04-03
### Changes
- Removed `fn header_mut()` from `bam::Read` trait.
- Fixed a major performance regression when reading bam files (issue #195).

## [0.29.0] - 2020-03-26
### Changes
- Migrate buffer intervals to u64. 

## [0.28.0] - 2020-03-26
### Changes
- Return u64 wherever htslib has migrated to using 64 bit.
- Implement more bio-types (Interval, Locus, Strand).

## [0.27.0] - 2020-03-17
### Changes
- Updated to Htslib 1.10.2.
- bam::Record.set() will panic if seq.len() != qual.len(). Previously, mismatched length would cause
  uninitialized memory to be written into the BAM file.
- use `serde_bytes` to serialize .data section of bam::Record when using serde - large speed improvement.
- change build.rs to avoid re-running when htslib or wrapper.h haven't changed.
- update some dependencies.
- refactor native dependency into htslib-sys crate, for greater versioning flexibility
- Record::from_sam require `&mut HeaderView`. Provide the appropriate accessor.
- set() no longer invalidates tag data.
- Various minor improvements.

## [0.26.1] - 2019-12-03
### Changes
- Various bug fixes in CIGAR string handling, INFO tag reading and FORMAT tag reading.

## [0.26.0] - 2019-09-27
### Changes
- Allow caching of CIGAR in bam::RecordBuffer.

## [0.25.0] - 2019-09-27
### Changes
- Migrated error handling to the snafu crate: https://docs.rs/snafu.
- Cleaned up and simplified API (including breaking changes).
- Allow writing SAM files from the bam::Writer.

## [0.24.0] - 2019-05-31
### Added
- Allow setting unmapped BAM record (without Cigar string).
- Various bug fixes and error handling improvements.
- Various Pysam-derived methods for interpreting Cigar strings.

## [0.23.0] - 2019-05-02
### Added
- Support for BAM indices that are not placed beside a file.bam as file.bam.bai
- Implement SequenceRead trait for BAM records.
- Add function to build an index for a BAM file.
- CRAM support for BAM reader and writer.
### Changes
- Allow to specify particular index filename when instantiating a BAM reader.
- Various bug and API fixes.

## [0.22.0] - 2018-11-02
### Changes
- Support compilation against musl.
- Support for removing alleles.
- Improvements to SyncedReader API.

## [0.21.0] - 2018-08-01
### Changes
- Adding `bcf::synced::SyncedReader::fetch()`, changing error type for `bcf::synced::SyncedReader::read_next()`.
- Adding `bcf::Record::unpack()` for explicitely unpacking BCF records.
- Fixed `bcf::synced::SyncedReader::record()`.
- `bam::Record::cigar()` now returns a reference (in constant time) and needs `bam::Record::unpack_cigar()` to be called first.
- Allow to create Cigar string from `bio_types::Alignment`.
- Provide a cached variant of obtaining cigar string.
- Lots of small usability improvements.

## [0.20.0] - 2018-06-18
### Added
- Initial implementation of synced BCF reader interface.
- Several small helper methods for BAM readers.
### Changes
- Not skipping `fileformat=` header any more.
- BCF records are always unpacked when reading.

## [0.19.1] - 2018-06-07
### Changed
- Moved unpacking of BCF records into constructor to prevent race conditions.
- Fixed bug in retrieving BCF record id.
- Fixed bug in the filter iterator of BCF.

## [0.19.0] - 2018-06-01
### Added
- more push functions for BCF.

## [0.18.0] - 2018-05-04
### Added
- bcf::IndexedReader
- support for writing bcf FILTER field
- setting thread count in all readers and writers
- setting ID and alleles in bcf records
- support for tabix indexes
- convert CIGAR to and from strings

## [0.17.0] - 2018-02-22
### Added
- Serde support for bam records.
### Changed
- Various convenience improvements in the API.

## [0.16.0] - 2018-01-05
### Changed
- Raw Htslib bindings are now generated on the fly.
- Switched to Htslib 1.6.
- Fixed a potential dangling pointer to the header in bcf records.
- Various small API improvements.

## [0.15.0] - 2017-12-05
### Changed
- HeaderView of bam and bcf can now be cloned.


## [0.14.0] - 2017-12-03
### Added
- An efficient ringbuffer for accessing BCF regions
- An efficient ringbuffer for accessing BAM regions
### Changed
- Improved mutability annotation for readers.

## [0.13.0] - 2017-09-22
### Added
- Ability to clone bam records.
- Ability to set only qname.
### Changed
- Further improved CIGAR string API.
- Improved documentation.


## [0.12.1] - 2017-06-12
### Changed
- Adapt to changes in Rust 1.18 that caused compilation issues.


## [0.12.0] - 2017-06-01
### Added
- Support seek and tell to handle virtual offsets.
### Changed
- Renamed previous seek method into fetch (in line with pysam).
- Improved CIGAR API.
- Updated dependencies.

## [0.11.0] - 2017-05-01
### Added
- A SAM writer.
### Changed
- Improved CIGAR string API using a newtype wrapper.
- Improved pileup API.
- Support threaded writing for BAM files.


## [0.10.0] - 2016-11-10
### Added
- Prelude module to easily import all relevant traits.
### Changed
- fine-grained constructors for STDIN/STDOUT, paths and URLs
- better template handling with bam files


## [0.9.0] - 2016-11-02
### Changed
- improved genotype handling
- improved error handling
- improved missing value handling

## [0.8.1] - 2016-08-17
### Changed
- Finally converted the last unit error types to real error types (really now!).

## [0.8.0] - 2016-08-17
### Changed
- More error types.

## [0.7.0] - 2016-08-16
### Changed
- Error types now properly implement the Display and the Error trait.

## [0.6.2] - 2016-07-22
### Changed
- Mark all records as Send and Sync.

## [0.6.1] - 2016-07-20
### Changed
- Improved error messages.
- Check existence of BAM when instantiating Readers.

## [0.6.0] - 2016-06-01
### Changed
- Improved handling of memory allocation for BCF and BAM records.
- Fixed a memory leak occuring when creating a new BAM record (thanks to @andrelmartins).
