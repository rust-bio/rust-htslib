// Copyright 2019 Manuel Landesfeind, Evotec International GmbH
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//!
//! Module for working with faidx-indexed FASTA files.
//!

use std::ffi;
use std::path::Path;
use url::Url;

use crate::htslib;

pub mod errors;
pub use errors::{Error, Result};

fn path_as_bytes<'a, P: 'a + AsRef<Path>>(path: P, must_exist: bool) -> Result<Vec<u8>> {
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
/// A Fasta reader.
#[derive(Debug)]
pub struct Reader {
    inner: *mut htslib::faidx_t,
}

impl Reader {
    /// Create a new Reader from a path.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, Error> {
        Self::new(&path_as_bytes(path, true)?)
    }

    /// Create a new Reader from an URL.
    ///
    /// # Arguments
    ///
    /// * `url` - the url to open
    pub fn from_url(url: &Url) -> Result<Self, Error> {
        Self::new(url.as_str().as_bytes())
    }

    /// Internal function to create a Reader from some sort of path (could be file path but also URL).
    /// The path or URL will be handled by the c-implementation transparently.
    ///
    /// # Arguments
    ///
    /// * `path` - the path or URL to open
    fn new(path: &[u8]) -> Result<Self, Error> {
        let cpath = ffi::CString::new(path).unwrap();
        let inner = unsafe { htslib::fai_load(cpath.as_ptr()) };
        Ok(Self { inner })
    }

    /// Fetches the sequence and returns it.
    ///
    /// # Arguments
    ///
    /// * `name` - the name of the template sequence (e.g., "chr1")
    /// * `begin` - the offset within the template sequence (starting with 0)
    /// * `end` - the end position to return (if smaller than `begin`, the behavior is undefined).
    pub fn fetch_seq<N: AsRef<str>>(&self, name: N, begin: usize, end: usize) -> String {
        let cname = ffi::CString::new(name.as_ref().as_bytes()).unwrap();
        let len_out: i32 = 0;
        let cseq = unsafe {
            let ptr = htslib::faidx_fetch_seq(
                self.inner,                              //*const faidx_t,
                cname.as_ptr(),                          // c_name
                begin as ::std::os::raw::c_int,          // p_beg_i
                end as ::std::os::raw::c_int,            // p_end_i
                &mut (len_out as ::std::os::raw::c_int), //len
            );
            ffi::CStr::from_ptr(ptr)
        };

        cseq.to_str().unwrap().to_owned()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn open_reader() -> Reader {
        Reader::from_path(format!("{}/test/test_cram.fa", env!("CARGO_MANIFEST_DIR")))
            .ok()
            .unwrap()
    }

    #[test]
    fn faidx_read_chr_first_base() {
        let r = open_reader();
        let seq = r.fetch_seq("chr1", 0, 0);
        assert_eq!(seq.len(), 1);
        assert_eq!(seq, "G");
    }

    #[test]
    fn faidx_read_chr_start() {
        let r = open_reader();
        let seq = r.fetch_seq("chr1", 0, 9);
        assert_eq!(seq.len(), 10);
        assert_eq!(seq, "GGGCACAGCC");
    }

    #[test]
    fn faidx_read_chr_between() {
        let r = open_reader();
        let seq = r.fetch_seq("chr1", 4, 14);
        assert_eq!(seq.len(), 11);
        assert_eq!(seq, "ACAGCCTCACC");
    }

    #[test]
    fn faidx_read_chr_end() {
        let r = open_reader();
        let seq = r.fetch_seq("chr1", 110, 120);
        assert_eq!(seq.len(), 10);
        assert_eq!(seq, "CCCCTCCGTG");
    }

    #[test]
    fn faidx_read_twice() {
        let r = open_reader();
        let seq = r.fetch_seq("chr1", 110, 120);
        assert_eq!(seq.len(), 10);
        assert_eq!(seq, "CCCCTCCGTG");

        let seq = r.fetch_seq("chr1", 5, 9);
        assert_eq!(seq.len(), 5);
        assert_eq!(seq, "CAGCC");
    }
}
