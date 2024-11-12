// Copyright 2019 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Module for working with BAM or CRAM indices.

use std::path::Path;
use std::ptr;

use crate::errors::{Error, Result};
use crate::htslib;
use crate::utils;

/// Index type to build.
pub enum Type {
    /// BAI index
    Bai,
    /// CSI index, with given minimum shift
    Csi(u32),
}

/// Build a BAM index.
pub fn build<P: AsRef<Path>>(
    bam_path: P,
    idx_path: Option<P>,
    idx_type: Type,
    n_threads: u32,
) -> Result<()> {
    let min_shift = match idx_type {
        Type::Bai => 0,
        Type::Csi(min_shift) => min_shift as i32,
    };
    let idx_path_cstr;
    let idx_path_ptr = if let Some(p) = idx_path {
        idx_path_cstr =
            utils::path_to_cstring(&p).expect("path_to_cstring unexpectedly returned with Err");
        idx_path_cstr.as_ptr()
    } else {
        ptr::null()
    };
    let ret = unsafe {
        htslib::sam_index_build3(
            utils::path_to_cstring(&bam_path).unwrap().as_ptr(),
            idx_path_ptr,
            min_shift,
            n_threads as i32,
        )
    };
    match ret {
        0 => Ok(()),
        -1 => Err(Error::BamBuildIndex),
        -2 => Err(Error::BamOpen {
            target: bam_path.as_ref().to_str().unwrap().to_owned(),
        }),
        -3 => Err(Error::BamNotIndexable),
        -4 => Err(Error::BamWriteIndex),
        e => panic!("unexpected error code from sam_index_build3: {}", e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index_build() {
        let test_bam = "test/test_index_build.bam";

        // test BAI index creation with 1 thread
        let idx1 = "test/results/test1.bam.bai";
        build(test_bam, Some(idx1), Type::Bai, 1).unwrap();
        assert!(Path::new(idx1).exists());

        // test BAI index creation with 2 threads
        let idx2 = "test/results/test2.bam.bai";
        build(test_bam, Some(idx2), Type::Bai, 2).unwrap();
        assert!(Path::new(idx2).exists());

        // test CSI index creation with 2 threads
        let idx3 = "test/results/test3.bam.csi";
        build(test_bam, Some(idx3), Type::Csi(2), 2).unwrap();
        assert!(Path::new(idx3).exists());

        // test CSI index creation with 2 threads and default file name
        build(test_bam, None, Type::Csi(5), 2).unwrap();
        assert!(Path::new("test/test_index_build.bam.csi").exists());
    }
}
