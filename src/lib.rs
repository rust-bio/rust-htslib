// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Rust-HTSlib provides a high level BAM API.
//! Reading and writing BAM files is as easy as
//!
//! ```
//! use rust_htslib::bam;
//! use rust_htslib::bam::Read;
//!
//! let bam = bam::Reader::new(&"test/test.bam").ok().expect("Error opening bam.");
//! let mut out = bam::Writer::with_template(&"test/test.bam", &"test/out.bam").ok().expect("Error opening bam.");
//!
//! // copy reverse reads to new BAM file
//! for r in bam.records() {
//!     let record = r.ok().expect("Error reading BAM file.");
//!     if record.is_reverse() {
//!         out.write(&record).ok().expect("Error writing BAM file.");
//!     }
//! }
//! ```
//!
//! Pileups can be performed with
//!
//! ```
//! use rust_htslib::bam;
//! use rust_htslib::bam::Read;
//!
//! let bam = bam::Reader::new(&"test/test.bam").ok().expect("Error opening bam.");
//!
//! // pileup over all covered sites
//! for p in bam.pileup() {
//!     let pileup = p.ok().expect("Error reading BAM file.");
//!     println!("{}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());
//!
//!     for alignment in pileup.alignments() {
//!         match alignment.indel() {
//!             bam::pileup::Indel::Ins(len) => println!("Insertion of length {}", len),
//!             bam::pileup::Indel::Del(len) => println!("Deletion of length {}", len),
//!             _ => println!("Base {}", alignment.record().seq()[alignment.qpos()])
//!         }
//!     }
//! }
//! ```
//!
//! In both cases, indexed BAM files can be seeked for specific regions, constraining either the record iterator or the pileups:
//!
//! ```
//! use rust_htslib::bam;
//!
//! let mut bam = bam::IndexedReader::new(&"test/test.bam").ok().expect("Error opening indexed BAM.");
//!
//! // seek to chr1:50000-100000
//! let tid = bam.header.tid(b"CHROMOSOME_I").unwrap();
//! bam.seek(tid, 0, 20).ok().expect("Error seeking BAM file.");
//! // afterwards, read or pileup in this region
//! ```


extern crate libc;
extern crate itertools;
#[macro_use]
extern crate quick_error;
#[macro_use]
extern crate newtype_derive;

pub mod htslib;
pub mod bam;
pub mod bcf;
pub mod utils;
