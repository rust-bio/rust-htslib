// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Rust-HTSlib provides a high level BAM API.
//! Reading and writing BAM files is as easy as
//!
//! ```
//! use rust_htslib::bam;
//! use rust_htslib::prelude::*;
//!
//! let mut bam = bam::Reader::from_path(&"test/test.bam").unwrap();
//! let header = bam::Header::from_template(bam.header());
//! let mut out = bam::Writer::from_path(&"test/out.bam", &header).unwrap();
//!
//! // copy reverse reads to new BAM file
//! for r in bam.records() {
//!     let record = r.unwrap();
//!     if record.is_reverse() {
//!         out.write(&record).unwrap();
//!     }
//! }
//! ```
//!
//! Pileups can be performed with
//!
//! ```
//! use rust_htslib::bam;
//! use rust_htslib::prelude::*;
//!
//! let mut bam = bam::Reader::from_path(&"test/test.bam").unwrap();
//!
//! // pileup over all covered sites
//! for p in bam.pileup() {
//!     let pileup = p.unwrap();
//!     println!("{}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());
//!
//!     for alignment in pileup.alignments() {
//!         if !alignment.is_del() && !alignment.is_refskip() {
//!             println!("Base {}", alignment.record().seq()[alignment.qpos().unwrap()]);
//!         }
//!         // mark indel start
//!         match alignment.indel() {
//!             bam::pileup::Indel::Ins(len) => println!("Insertion of length {} between this and next position.", len),
//!             bam::pileup::Indel::Del(len) => println!("Deletion of length {} between this and next position.", len),
//!             bam::pileup::Indel::None => ()
//!         }
//!     }
//! }
//! ```
//!
//! In both cases, indexed BAM files can be seeked for specific regions, constraining either the record iterator or the pileups:
//!
//! ```
//! use rust_htslib::bam;
//! use rust_htslib::prelude::*;
//!
//! let mut bam = bam::IndexedReader::from_path(&"test/test.bam").unwrap();
//!
//! // seek to chr1:50000-100000
//! let tid = bam.header().tid(b"CHROMOSOME_I").unwrap();
//! bam.fetch(tid, 0, 20).unwrap();
//! // afterwards, read or pileup in this region
//! ```


extern crate libc;
extern crate itertools;
#[macro_use]
extern crate quick_error;
#[macro_use]
extern crate newtype_derive;
#[macro_use]
extern crate custom_derive;
extern crate url;
extern crate ieee754;
#[macro_use]
extern crate lazy_static;
extern crate bitflags;

pub mod htslib;
pub mod bam;
pub mod bcf;
pub mod utils;
pub mod prelude;
pub mod sam;
