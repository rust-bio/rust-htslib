// Copyright 2019 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Module for working with BAM or CRAM indices.

use std::path::Path;
use std::ptr;

use htslib;
use utils;

/// Index type to build.
pub enum Type {
    /// BAI index
    BAI,
    /// CSI index, with given minimum shift
    CSI(u32),
}

/// Build a BAM index.
pub fn build<P: AsRef<Path>>(
    bam_path: P,
    idx_path: Option<P>,
    idx_type: Type,
    n_threads: u32,
) -> Result<(), IndexBuildError> {
    let min_shift = match idx_type {
        Type::BAI => 0,
        Type::CSI(min_shift) => min_shift as i32,
    };
    let idx_path = idx_path.map(|p| utils::path_to_cstring(&p).unwrap());
    let ret = unsafe {
        htslib::sam_index_build3(
            utils::path_to_cstring(&bam_path).unwrap().as_ptr(),
            idx_path.map(|p| p.as_ptr()).unwrap_or(ptr::null()),
            min_shift,
            n_threads as i32,
        )
    };
    match ret {
        0 => Ok(()),
        -1 => Err(IndexBuildError::Some),
        -2 => Err(IndexBuildError::Opening),
        -3 => Err(IndexBuildError::InvalidFormat),
        -4 => Err(IndexBuildError::Saving),
        e => panic!("unexpected error code from sam_index_build3: {}", e),
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum IndexBuildError {
        Some {
            description("error building index")
        }
        Opening {
            description("error opening input file for building index")
        }
        InvalidFormat {
            description("format is not indexable")
        }
        Saving {
            description("could not save index")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index_build() {
        let idx = "test/results/test.bam.bai";
        build("test/test.bam", Some(idx), Type::BAI, 1).unwrap();
        assert!(Path::new(idx).exists());
        build("test/test.bam", Some(idx), Type::BAI, 2).unwrap();
        build(
            "test/test.bam",
            Some("test/results/test.bam.csi"),
            Type::CSI(2),
            1,
        )
        .unwrap();
    }
}
