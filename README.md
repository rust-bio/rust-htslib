[![Crates.io](https://img.shields.io/crates/d/rust-htslib.svg)](https://crates.io/crates/rust-htslib)
[![Crates.io](https://img.shields.io/crates/v/rust-htslib.svg)](https://crates.io/crates/rust-htslib)
[![Crates.io](https://img.shields.io/crates/l/rust-htslib.svg)](https://crates.io/crates/rust-htslib)
[![docs.rs](https://docs.rs/rust-htslib/badge.svg)](https://docs.rs/rust-htslib)
![GitHub Workflow Status](https://img.shields.io/github/workflow/status/rust-bio/rust-htslib/CI/master?label=tests)

# HTSlib bindings for Rust

This library provides HTSlib bindings and a high level Rust API for reading and writing BAM files.

To clone this repository, issue

```shell
$ git clone --recursive https://github.com/rust-bio/rust-htslib.git
```

ensuring that the HTSlib submodule is fetched, too.
If you only want to use the library, there is no need to clone the repository. Go on to the **Usage** section in this case.

## Requirements

To compile this crate you need docker and cross:

```shell
$ cargo install cross
$ cross build 				              # will build with GNU toolchain
```

If you want to run rust-htslib code on AWS lambda, you'll need to statically compile it with MUSL as follows:

```shell
$ export CFLAGS="-I/usr/local/musl/include"
$ cross build --target x86_64-unknown-linux-musl      # will build with MUSL toolchain
```

Alternatively, you can also install it locally by installing the development headers of zlib, bzip2 and xz. For instance, in Debian systems one needs the following dependencies:

```shell
$ sudo apt-get install zlib1g-dev libbz2-dev liblzma-dev clang pkg-config
```

On OSX, this will take a significant amount of time due to musl cross compiling toolchain:

```shell
$ brew install FiloSottile/musl-cross/musl-cross
$ brew install bzip2 zlib xz curl-openssl
```

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

By default `rust-htslib` uses `bindgen` to generate bindings to htslib. This can slow down the build substantially. Disabling the `bindgen` feature will 
cause `hts-sys` to use a prebuilt binding.rs file for your architecture.  Note
Windows is not currently supported & this approach is experimental -- please 
exercise caution.

`rust-htslib` also has optional support for `serde`, to allow (de)serialization of `bam::Record` via any serde-supported format:

```toml
[dependencies]
rust-htslib = { version = "*", features = ["serde"] }
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

Licensed under the MIT license http://opensource.org/licenses/MIT. This project may not be copied, modified, or distributed except according to those terms.
