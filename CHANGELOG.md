# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

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
