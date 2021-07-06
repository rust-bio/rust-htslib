#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(clippy::all)]
#![allow(improper_ctypes)]
//! This module exposes the raw HTSlib bindings.

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
