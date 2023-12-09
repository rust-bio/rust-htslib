// Copyright 2020 Manuel Landesfeind, Evotec International GmbH
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

use crate::errors::{Error, Result};
use crate::utils::path_as_bytes;

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

    /// Fetch the sequence as a byte array.
    ///
    /// # Arguments
    ///
    /// * `name` - the name of the template sequence (e.g., "chr1")
    /// * `begin` - the offset within the template sequence (starting with 0)
    /// * `end` - the end position to return (if smaller than `begin`, the behavior is undefined).
    pub fn fetch_seq<N: AsRef<str>>(&self, name: N, begin: usize, end: usize) -> Result<&[u8]> {
        if begin > std::i64::MAX as usize {
            return Err(Error::FaidxPositionTooLarge);
        }
        if end > std::i64::MAX as usize {
            return Err(Error::FaidxPositionTooLarge);
        }
        let cname = ffi::CString::new(name.as_ref().as_bytes()).unwrap();
        let len_out: i64 = 0;
        let cseq = unsafe {
            let ptr = htslib::faidx_fetch_seq64(
                self.inner,                          //*const faidx_t,
                cname.as_ptr(),                      // c_name
                begin as htslib::hts_pos_t,          // p_beg_i
                end as htslib::hts_pos_t,            // p_end_i
                &mut (len_out as htslib::hts_pos_t), //len
            );
            ffi::CStr::from_ptr(ptr)
        };

        Ok(cseq.to_bytes())
    }

    /// Fetches the sequence and returns it as string.
    ///
    /// # Arguments
    ///
    /// * `name` - the name of the template sequence (e.g., "chr1")
    /// * `begin` - the offset within the template sequence (starting with 0)
    /// * `end` - the end position to return (if smaller than `begin`, the behavior is undefined).
    pub fn fetch_seq_string<N: AsRef<str>>(
        &self,
        name: N,
        begin: usize,
        end: usize,
    ) -> Result<String> {
        let bytes = self.fetch_seq(name, begin, end)?;
        Ok(std::str::from_utf8(bytes).unwrap().to_owned())
    }

    /// Fetches the number of sequences in the fai index
    pub fn n_seqs(&self) -> u64 {
        let n = unsafe { htslib::faidx_nseq(self.inner) };
        n as u64
    }

    /// Fetches the i-th sequence name
    ///
    /// # Arguments
    ///
    /// * `i` - index to query
    pub fn seq_name(&self, i: i32) -> Result<String> {
        let cname = unsafe {
            let ptr = htslib::faidx_iseq(self.inner, i);
            ffi::CStr::from_ptr(ptr)
        };

        let out = match cname.to_str() {
            Ok(s) => s.to_string(),
            Err(_) => {
                return Err(Error::FaidxBadSeqName);
            }
        };

        Ok(out)
    }

    /// Fetches the length of the given sequence name.
    ///
    /// # Arguments
    ///
    /// * `name` - the name of the template sequence (e.g., "chr1")
    pub fn fetch_seq_len<N: AsRef<str>>(&self, name: N) -> u64 {
        let cname = ffi::CString::new(name.as_ref().as_bytes()).unwrap();
        let seq_len = unsafe { htslib::faidx_seq_len(self.inner, cname.as_ptr()) };
        seq_len as u64
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        unsafe {
            htslib::fai_destroy(self.inner);
        }
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
    fn faidx_open() {
        open_reader();
    }

    #[test]
    fn faidx_read_chr_first_base() {
        let r = open_reader();

        let bseq = r.fetch_seq("chr1", 0, 0).unwrap();
        assert_eq!(bseq.len(), 1);
        assert_eq!(bseq, b"G");

        let seq = r.fetch_seq_string("chr1", 0, 0).unwrap();
        assert_eq!(seq.len(), 1);
        assert_eq!(seq, "G");
    }

    #[test]
    fn faidx_read_chr_start() {
        let r = open_reader();

        let bseq = r.fetch_seq("chr1", 0, 9).unwrap();
        assert_eq!(bseq.len(), 10);
        assert_eq!(bseq, b"GGGCACAGCC");

        let seq = r.fetch_seq_string("chr1", 0, 9).unwrap();
        assert_eq!(seq.len(), 10);
        assert_eq!(seq, "GGGCACAGCC");
    }

    #[test]
    fn faidx_read_chr_between() {
        let r = open_reader();

        let bseq = r.fetch_seq("chr1", 4, 14).unwrap();
        assert_eq!(bseq.len(), 11);
        assert_eq!(bseq, b"ACAGCCTCACC");

        let seq = r.fetch_seq_string("chr1", 4, 14).unwrap();
        assert_eq!(seq.len(), 11);
        assert_eq!(seq, "ACAGCCTCACC");
    }

    #[test]
    fn faidx_read_chr_end() {
        let r = open_reader();

        let bseq = r.fetch_seq("chr1", 110, 120).unwrap();
        assert_eq!(bseq.len(), 10);
        assert_eq!(bseq, b"CCCCTCCGTG");

        let seq = r.fetch_seq_string("chr1", 110, 120).unwrap();
        assert_eq!(seq.len(), 10);
        assert_eq!(seq, "CCCCTCCGTG");
    }

    #[test]
    fn faidx_read_twice_string() {
        let r = open_reader();
        let seq = r.fetch_seq_string("chr1", 110, 120).unwrap();
        assert_eq!(seq.len(), 10);
        assert_eq!(seq, "CCCCTCCGTG");

        let seq = r.fetch_seq_string("chr1", 5, 9).unwrap();
        assert_eq!(seq.len(), 5);
        assert_eq!(seq, "CAGCC");
    }

    #[test]
    fn faidx_read_twice_bytes() {
        let r = open_reader();
        let seq = r.fetch_seq("chr1", 110, 120).unwrap();
        assert_eq!(seq.len(), 10);
        assert_eq!(seq, b"CCCCTCCGTG");

        let seq = r.fetch_seq("chr1", 5, 9).unwrap();
        assert_eq!(seq.len(), 5);
        assert_eq!(seq, b"CAGCC");
    }

    #[test]
    fn faidx_position_too_large() {
        let r = open_reader();
        let position_too_large = i64::MAX as usize;
        let res = r.fetch_seq("chr1", position_too_large, position_too_large + 1);
        assert_eq!(res, Err(Error::FaidxPositionTooLarge));
    }

    #[test]
    fn faidx_n_seqs() {
        let r = open_reader();
        assert_eq!(r.n_seqs(), 3);
    }

    #[test]
    fn faidx_seq_name() {
        let r = open_reader();
        let n = r.seq_name(1).unwrap();
        assert_eq!(n, "chr2");
    }

    #[test]
    fn faidx_get_seq_len() {
        let r = open_reader();
        let chr1_len = r.fetch_seq_len("chr1");
        let chr2_len = r.fetch_seq_len("chr2");
        assert_eq!(chr1_len, 120u64);
        assert_eq!(chr2_len, 120u64);
    }

    #[test]
    fn open_many_readers() {
        for _ in 0..500_000 {
            let reader = open_reader();
            drop(reader);
        }
    }
}
