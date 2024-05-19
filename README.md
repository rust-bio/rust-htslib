[![Crates.io](https://img.shields.io/crates/d/rust-htslib.svg)](https://crates.io/crates/rust-htslib)
[![Crates.io](https://img.shields.io/crates/v/rust-htslib.svg)](https://crates.io/crates/rust-htslib)
[![Crates.io](https://img.shields.io/crates/l/rust-htslib.svg)](https://crates.io/crates/rust-htslib)
[![docs.rs](https://docs.rs/rust-htslib/badge.svg)](https://docs.rs/rust-htslib)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/rust-bio/rust-htslib/rust.yml?branch=master&label=tests)
[![Coverage Status](https://coveralls.io/repos/github/rust-bio/rust-htslib/badge.svg?branch=master)](https://coveralls.io/github/rust-bio/rust-htslib?branch=master)

# HTSlib bindings for Rust

This library provides HTSlib bindings and a high level Rust API for reading and writing BAM files.

To clone this repository, issue

```shell
$ git clone --recursive https://github.com/rust-bio/rust-htslib.git
```

ensuring that the HTSlib submodule is fetched, too.
If you only want to use the library, there is no need to clone the repository. Go on to the **Usage** section in this case.

## Requirements

rust-htslib comes with pre-built bindings to htslib for Mac and Linux. You will need a C toolchain compatible with the `cc` crate. The build script for this crate will automatically build a link htslib.

## Usage

Add this to your `Cargo.toml`:
```toml
[dependencies]
rust-htslib = "*"
```

By default `rust-htslib` links to `bzip2-sys` and `lzma-sys` for full CRAM support. If you do not need CRAM support, or you do need to support CRAM files
with these compression methods, you can deactivate these features to reduce you dependency count:

```toml
[dependencies]
rust-htslib = { version = "*", default-features = false }
```

`rust-htslib` has optional support for `serde`, to allow (de)serialization of `bam::Record` via any serde-supported format.

Http access to files is available with the `curl` feature.

Beta-level S3 and Google Cloud Storge support is available with the `s3` and `gcs` features.

`rust-htslib` can optionally use `bindgen` to generate bindings to htslib. This can slow down the build substantially. Enabling the `bindgen` feature will 
cause `hts-sys` to use a create a binding file for your architecture. Pre-built bindings are supplied for Mac and Linux. The `bindgen` feature on Windows is untested - please file a bug if you need help.



```toml
[dependencies]
rust-htslib = { version = "*", features = ["serde_feature"] }
```

For more information, please see the [docs](https://docs.rs/rust-htslib).

# Alternatives

There's [noodles](https://github.com/zaeleus/noodles) by [Michael Macias](https://github.com/zaeleus) which implements a large part of htslib's C functionality in pure Rust (still experimental though).

# Authors

* [Johannes Köster](https://github.com/johanneskoester)
* [Christopher Schröder](https://github.com/christopher-schroeder)
* [Patrick Marks](https://github.com/pmarks)
* [David Lähnemann](https://github.com/dlaehnemann)
* [Manuel Holtgrewe](https://github.com/holtgrewe)
* [Julian Gehring](https://github.com/juliangehring)

For other contributors, see [here](https://github.com/rust-bio/rust-htslib/graphs/contributors).

## License

Licensed under the MIT license https://opensource.org/licenses/MIT. This project may not be copied, modified, or distributed except according to those terms.
Some test files are taken from https://github.com/samtools/htslib.
