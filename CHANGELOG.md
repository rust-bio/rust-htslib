# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.46.0](https://github.com/rust-bio/rust-htslib/compare/v0.45.0...v0.46.0) (2024-02-22)


### Features

* making several RecordBuffer methods public ([6757f52](https://github.com/rust-bio/rust-htslib/commit/6757f5219955fd4edba4f61e62978ce1e001068e))


### Bug Fixes

* fix building libz-sys ([#420](https://github.com/rust-bio/rust-htslib/issues/420)) ([01c8849](https://github.com/rust-bio/rust-htslib/commit/01c884945686e7a6756406b579fde28657f70b36))

## [0.45.0](https://github.com/rust-bio/rust-htslib/compare/v0.44.1...v0.45.0) (2024-02-07)


### Features

* adding function to get sequence length to faidx mod ([#410](https://github.com/rust-bio/rust-htslib/issues/410)) ([ae79eba](https://github.com/rust-bio/rust-htslib/commit/ae79eba82ef6929105bdbe08246a8e973660899e))


### Bug Fixes

* Loosen acceptable types to support current linux build on aarch64 ([#415](https://github.com/rust-bio/rust-htslib/issues/415)) ([1d78d12](https://github.com/rust-bio/rust-htslib/commit/1d78d1251a052461605d28cd8cf832ccad93ef73))

## [0.44.1](https://github.com/rust-bio/rust-htslib/compare/v0.44.0...v0.44.1) (2023-06-21)


### Bug Fixes

* use correct return value in bcf_get_format and bcf_get_info_values ([#398](https://github.com/rust-bio/rust-htslib/issues/398)) ([f9a1981](https://github.com/rust-bio/rust-htslib/commit/f9a1981fa84eef39e35f868ddfc773ea265b94b3))

## [0.44.0](https://github.com/rust-bio/rust-htslib/compare/v0.43.1...v0.44.0) (2023-06-20)


### Features

* implement Clone for bcf::Record ([#394](https://github.com/rust-bio/rust-htslib/issues/394)) ([e89538d](https://github.com/rust-bio/rust-htslib/commit/e89538d5a9971c6508ac38d92ac468f3d70241aa))
* implement htslib basemod api ([#385](https://github.com/rust-bio/rust-htslib/issues/385)) ([8beee14](https://github.com/rust-bio/rust-htslib/commit/8beee145a116f7ae936f1b6e36d876116dca18f1))


### Bug Fixes

* include doctests in test coverage calculations ([#397](https://github.com/rust-bio/rust-htslib/issues/397)) ([8ed0837](https://github.com/rust-bio/rust-htslib/commit/8ed083783fa1dce09535564a090d37f687fc832f))

## [0.43.1](https://github.com/rust-bio/rust-htslib/compare/v0.43.0...v0.43.1) (2023-05-16)


### Bug Fixes

* implement Drop for faidx::Reader, destroying the fai handle ([#391](https://github.com/rust-bio/rust-htslib/issues/391)) ([0e6d6ac](https://github.com/rust-bio/rust-htslib/commit/0e6d6acec9a6d24ed6baf810e56f02394737a046))

## [0.43.0](https://github.com/rust-bio/rust-htslib/compare/v0.42.0...v0.43.0) (2023-05-12)


### Features

* HeaderRecord::push_tag: Value may be owned ([#388](https://github.com/rust-bio/rust-htslib/issues/388)) ([b64537d](https://github.com/rust-bio/rust-htslib/commit/b64537db011b76e5ace73a5e74c598a62a0a018b))
* Index for `bam::IndexedReader` ([#387](https://github.com/rust-bio/rust-htslib/issues/387)) ([fb74387](https://github.com/rust-bio/rust-htslib/commit/fb743875182c9232894e07007367f09f05d6e275))

## [0.42.0](https://github.com/rust-bio/rust-htslib/compare/v0.41.1...v0.42.0) (2023-03-30)


### Features

* Add ability to fetch number of sequences and I-th sequence from FAI index ([#377](https://github.com/rust-bio/rust-htslib/issues/377)) ([6ecc4bd](https://github.com/rust-bio/rust-htslib/commit/6ecc4bd1f88110da278c7f934453024e4e64ac74))

## [0.41.1](https://github.com/rust-bio/rust-htslib/compare/v0.41.0...v0.41.1) (2023-03-03)

### Features

* Revised calculation of leading- and trailing-softclips ([#375](https://github.com/rust-bio/rust-htslib/issues/375)) ([b61dd2c](https://github.com/rust-bio/rust-htslib/commit/b61dd2cfb2b74c0180f2d76bbd4ed4eb14fa09b3))


### Performance Improvements

* update htslib and corresponding bindings to 1.16 ([#366](https://github.com/rust-bio/rust-htslib/issues/366)) ([f597ce0](https://github.com/rust-bio/rust-htslib/commit/f597ce0451e3f3c393166a7291486bbc2bde4c39))


## [0.40.2](https://github.com/rust-bio/rust-htslib/compare/rust-htslib-v0.40.1...rust-htslib-v0.40.2) (2022-10-13)


### Performance Improvements

* update htslib and corresponding bindings to 1.16 ([#366](https://github.com/rust-bio/rust-htslib/issues/366)) ([f597ce0](https://github.com/rust-bio/rust-htslib/commit/f597ce0451e3f3c393166a7291486bbc2bde4c39))

## [0.40.1](https://github.com/rust-bio/rust-htslib/compare/rust-htslib-v0.40.0...rust-htslib-v0.40.1) (2022-08-24)


### Bug Fixes

* Header::to_hashmap skips `@CO` tags, add `comments()`  method ([#363](https://github.com/rust-bio/rust-htslib/issues/363)) ([c24a7f6](https://github.com/rust-bio/rust-htslib/commit/c24a7f69fbe5d2db4a6f1fbd6eda3922fe7f1c18))

## [0.40.0](https://github.com/rust-bio/rust-htslib/compare/rust-htslib-v0.39.5...rust-htslib-v0.40.0) (2022-07-05)


### Features

* Add wrapper of BGZF writer ([#349](https://github.com/rust-bio/rust-htslib/issues/349)) ([965ed88](https://github.com/rust-bio/rust-htslib/commit/965ed886a0c24ee3070e48cc192c0772ac5cbaf4))


### Bug Fixes

* update to latest release-please ([b130634](https://github.com/rust-bio/rust-htslib/commit/b130634b3d096e620dcfe59acae2200df3e4d847))

### [0.39.5](https://www.github.com/rust-bio/rust-htslib/compare/rust-htslib-v0.39.4...rust-htslib-v0.39.5) (2022-05-09)


### Bug Fixes

* set path in release-please config ([d8f7c6e](https://www.github.com/rust-bio/rust-htslib/commit/d8f7c6e8f31accb7576e150fa1439e177f7816cb))

### [0.39.4](https://www.github.com/rust-bio/rust-htslib/compare/rust-htslib-v0.39.3...rust-htslib-v0.39.4) (2022-05-09)


### Bug Fixes

* perform checkout before running release please ([cbc6a0a](https://www.github.com/rust-bio/rust-htslib/commit/cbc6a0ad37c5623d14f2ed0bcbb4c5289d012fcb))

### [0.39.3](https://www.github.com/rust-bio/rust-htslib/compare/rust-htslib-v0.39.2...rust-htslib-v0.39.3) (2021-11-20)


### Bug Fixes

* change the type to c_char so it can be compiled for aarch64 ([#337](https://www.github.com/rust-bio/rust-htslib/issues/337)) ([a21aff2](https://www.github.com/rust-bio/rust-htslib/commit/a21aff289bc03c7549afc7a958084ed57e8c93f2))

### [0.39.2](https://www.github.com/rust-bio/rust-htslib/compare/rust-htslib-v0.39.1...rust-htslib-v0.39.2) (2021-08-23)


### Bug Fixes

* Configuration when cross-compiling. Even when cross-compiling, build.rs runs on the build host. Hence within build.rs `#[cfg(target_os)]` always reflects the host, not the target. Use $CARGO_CFG_TARGET_OS instead to query target properties. ([#329](https://www.github.com/rust-bio/rust-htslib/issues/329)) ([d5198e6](https://www.github.com/rust-bio/rust-htslib/commit/d5198e6c777fdbbfdd9c73a820f1be983a458ce2))

### [0.39.1](https://www.github.com/rust-bio/rust-htslib/compare/rust-htslib-v0.39.0...rust-htslib-v0.39.1) (2021-07-06)


### Bug Fixes

* bump hts-sys version to 2.0.1 ([336c8b8](https://www.github.com/rust-bio/rust-htslib/commit/336c8b8a1779422afea1065e37bcc44f54abac42))

## [0.39.0](https://www.github.com/rust-bio/rust-htslib/compare/rust-htslib-vrust-htslib-0.38.3...rust-htslib-v0.39.0) (2021-07-06)


### ⚠ BREAKING CHANGES

* dummy major version bump to move away from previous versions that were following htslib versions.
* bump to new major version (for technical reasons).
* dummy breaking change to increase hts-sys major version.

### Bug Fixes

* bump to new major version (for technical reasons). ([9c6db30](https://www.github.com/rust-bio/rust-htslib/commit/9c6db3060818692070db1411d63e113dc7effd64))
* dummy breaking change to increase hts-sys major version. ([93415cb](https://www.github.com/rust-bio/rust-htslib/commit/93415cbb82e4f11d257a2b2cedba2664f86a034d))
* dummy changes ([3af5ede](https://www.github.com/rust-bio/rust-htslib/commit/3af5ede13a6b44ce5d1e7f0eb90836a692e711ec))
* dummy major version bump to move away from previous versions that were following htslib versions. ([aaa70a8](https://www.github.com/rust-bio/rust-htslib/commit/aaa70a85ef9a908d3b101f23879189e84a15d23f))
* dummy release ([74d1565](https://www.github.com/rust-bio/rust-htslib/commit/74d1565329fc862f1172c0925c7b66ceb8bcf988))
* dummy release ([af2f84e](https://www.github.com/rust-bio/rust-htslib/commit/af2f84eb0411507f8866b3cc05e9a6ba9d81d172))
* dummy release ([b97915f](https://www.github.com/rust-bio/rust-htslib/commit/b97915f2c70da4c914f2e69861bf78eec5979baf))
* handle subcrate with release-please ([0a4605f](https://www.github.com/rust-bio/rust-htslib/commit/0a4605f165cb2edf4428d8fb39f7e4787585f4e1))
* trigger dummy release ([7c5a7de](https://www.github.com/rust-bio/rust-htslib/commit/7c5a7de33e2a92052126e5f44389d421974d1e02))
* update changelog ([deef08f](https://www.github.com/rust-bio/rust-htslib/commit/deef08feb0b5ba2d8abf98f2cc6d327236da8aef))

### [0.38.3](https://www.github.com/rust-bio/rust-htslib/compare/v0.38.2...v0.38.3) (2021-07-06)


### Bug Fixes

* dummy fix for triggering release ([e92e6b1](https://www.github.com/rust-bio/rust-htslib/commit/e92e6b10b0a7e5db50b12e2fbe2c42b467eb369e))

### [0.38.2](https://www.github.com/rust-bio/rust-htslib/compare/v0.38.1...v0.38.2) (2021-07-06)


### Bug Fixes

* add ID to automatic release handling ([1244393](https://www.github.com/rust-bio/rust-htslib/commit/124439300e1e3e01e1d847f7549747d560c01989))

### [0.38.1](https://www.github.com/rust-bio/rust-htslib/compare/v0.38.0...v0.38.1) (2021-07-06)


### Bug Fixes

* improved documentation ([cb0b66c](https://www.github.com/rust-bio/rust-htslib/commit/cb0b66c4a92d4f03debe38dfb2a014b154c7dd96))

## [0.38.0](https://www.github.com/rust-bio/rust-htslib/compare/v0.37.0...v0.38.0) (2021-07-06)


### ⚠ BREAKING CHANGES

* Improve bcf Record filter interface and improve docs (#306)

### Features

* Improve bcf Record filter interface and improve docs ([#306](https://www.github.com/rust-bio/rust-htslib/issues/306)) ([f45e91d](https://www.github.com/rust-bio/rust-htslib/commit/f45e91dfdc64ecb662d676f2996ed4f14c079995))


### Bug Fixes

* enum name usage in doc example ([#311](https://www.github.com/rust-bio/rust-htslib/issues/311)) ([6e9ba49](https://www.github.com/rust-bio/rust-htslib/commit/6e9ba4928b60c3105490a8179d074c705ea06fd7))

## [Unreleased]
### Changes
- `bcf::Record` methods `has_filter`, `remove_filter`, `push_filter`, and `set_filter`
  all now take a byte slice (`&[u8]`) or an `Id`. 

[Unreleased]: https://github.com/rust-bio/rust-htslib/compare/v0.37.0...HEAD

## [0.37.0] - 2021-07-05
### Added
- `bcf::Record` methods `end`, `clear`, and `rlen` (@mbhall88).

### Changes
- `bcf::IndexReader::fetch` parameter `end` is now an `Option<u64>`. This is inline with
  htslib regions, which do not require an end position (@mbhall88)
- Removed unused dependencies (@sreenathkrishnan).
- Improved documentation (@mbhall88).
- Improved error message when failing to load index files (@johanneskoester).
- Improved API for accessing AUX fields in BAM records (@jch-13).
- Fixed compiler warnings (@fxwiegand).
- BAM header representation is now always kept in sync between textual and binary (@jch-13).

## [0.36.0] - 2020-11-23
### Changes
- Improved genotype API in VCF/BCF records (@MaltheSR).
- Read pair orientation inference for BAM records (@johanneskoester).

## [0.35.2] - 2020-11-23
### Changes
- let hts-sys dependency comply to semver.

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
