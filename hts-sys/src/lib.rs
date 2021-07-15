#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(clippy::all)]
#![allow(improper_ctypes)]
//! This module exposes the raw [HTSlib](https://github.com/samtools/htslib) bindings for reading and writing
//! genomics file formats like SAM, BAM, CRAM, VCF, BCF and tabix.
//! Instead of using this crate directly, it is recommended to use [rust-htslib](https://docs.rs/rust-htslib)
//! instead, which has a more idiomatic and high-level API and builds on top of this crate.

#[cfg(feature = "bzip2")]
extern crate bzip2_sys;
#[cfg(feature = "curl")]
extern crate curl_sys;
#[cfg(feature = "libdeflate")]
extern crate libdeflate_sys;
extern crate libz_sys;
#[cfg(feature = "lzma")]
extern crate lzma_sys;

// include on-the-fly generated bindings
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
