#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

//! This module exposes the raw Htslib bindings.

// include on-the-fly generated bindings
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
