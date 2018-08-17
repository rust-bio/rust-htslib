// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Module with utility code.

use bv::BitVec;
use bv::{BitSlice, Bits};
use htslib;
use std::ffi;
use std::mem;
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
        format!("dst len {} < src len {}", dst.len(), src.len())
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

//
// pub fn kbitset(values: &[bool]) -> *mut htslib::kbitset_t {
//     let block_size = mem::size_of::<u64>();
//     let n = (values.len() + block_size - 1) / block_size;
//     let kbs = unsafe {
// 		::libc::malloc(mem::size_of::<usize>() + block_size * (n + 1)) as *mut htslib::kbitset_t
//     };
//     unsafe {
//         (*kbs).n = n;
//     }
//     for (i, &b) in values.iter().enumerate() {
//         let j = i / block_size;
//         unsafe {
//             ptr::write(&mut (*kbs).b[j], (*kbs).b[j] | 1u64 << (i % block_size));
//         }
//     }
//     unsafe {
//         ptr::write(&mut (*kbs).b[n], !0u64);
//     }
//
//     kbs
// }
