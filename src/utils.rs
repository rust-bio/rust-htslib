
use std::ptr;
use std::ffi;
use std::path::Path;

/// Copies data from `src` to `dst`
/// TODO remove once stable in standard library.
///
/// Panics if the length of `dst` is less than the length of `src`.
#[inline]
pub fn copy_memory(src: &[u8], dst: &mut [u8]) {
    let len_src = src.len();
    assert!(dst.len() >= len_src, format!("dst len {} < src len {}", dst.len(), src.len()));
    // `dst` is unaliasable, so we know statically it doesn't overlap
    // with `src`.
    unsafe {
        ptr::copy_nonoverlapping(src.as_ptr(),
                                 dst.as_mut_ptr(),
                                 len_src);
    }
}


pub fn path_to_cstring<P: AsRef<Path>>(path: &P) -> Option<ffi::CString> {
    path.as_ref().to_str().and_then(|p| ffi::CString::new(p).ok())
}
