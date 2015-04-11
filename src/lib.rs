// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


#![feature(libc)]
#![feature(std_misc)]
#![feature(core)]
#![feature(collections)]
#![feature(step_by)]
#![feature(convert)]


extern crate libc;
pub mod htslib;
pub mod bam;
pub mod bcf;
