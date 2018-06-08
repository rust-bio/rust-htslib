# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.19.1] - 2018-06-07
### Changed
- Fixed bug in missing unpack when retrieving BCF record id.
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
