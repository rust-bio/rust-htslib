// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Module with utility code.

use crate::errors::{Error, Result};
use std::ffi;
use std::path::Path;
use std::ptr;

/// Copies data from `src` to `dst`
/// TODO remove once stable in standard library.
///
/// Panics if the length of `dst` is less than the length of `src`.
#[inline]
pub fn copy_memory(src: &[u8], dst: &mut [u8]) {
    let len_src = src.len();
    assert!(
        dst.len() >= len_src,
        "dst len {} < src len {}",
        dst.len(),
        src.len()
    );
    // `dst` is unaliasable, so we know statically it doesn't overlap
    // with `src`.
    unsafe {
        ptr::copy_nonoverlapping(src.as_ptr(), dst.as_mut_ptr(), len_src);
    }
}

pub fn path_to_cstring<P: AsRef<Path>>(path: &P) -> Option<ffi::CString> {
    path.as_ref()
        .to_str()
        .and_then(|p| ffi::CString::new(p).ok())
}

/// Convert a path into a byte-vector
pub fn path_as_bytes<'a, P: 'a + AsRef<Path>>(path: P, must_exist: bool) -> Result<Vec<u8>> {
    if path.as_ref().exists() || !must_exist {
        Ok(path
            .as_ref()
            .to_str()
            .ok_or(Error::NonUnicodePath)?
            .as_bytes()
            .to_owned())
    } else {
        Err(Error::FileNotFound {
            path: path.as_ref().to_owned(),
        })
    }
}
