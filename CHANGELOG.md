# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

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
