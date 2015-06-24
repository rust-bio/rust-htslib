// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


#![feature(libc)]
#![feature(step_by)]
#![feature(convert)]
#![feature(vec_push_all)]
#![feature(copy_lifetime)]
#![feature(slice_bytes)]

extern crate libc;
pub mod htslib;
pub mod bam;
pub mod bcf;
