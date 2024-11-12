// Copyright 2018 Manuel Holtgrewe, Berlin Institute of Health.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Module for working with tabix-indexed text files.
//!
//! This module allows to read tabix-indexed text files (such as BED) in a convenient but in a
//! line-based (and thus format-agnostic way). For accessing tabix-inxed VCF files, using the
//! `bcf` module is probably a better choice as this module gives you lines from the text files
//! which you then have to take care of parsing.
//!
//! In general, for reading tabix-indexed files, first to open the file by creating a `tbx::Reader`
//! objects, possibly translate the chromosome name to its numeric ID in the file, fetch the region
//! of interest using `fetch()`, and finally iterate over the records using `records()`.
//!
//! # Examples
//!
//! ```rust,no_run
//! use rust_htslib::tbx::{self, Read};
//!
//! // Create a tabix reader for reading a tabix-indexed BED file.
//! let path_bed = "file.bed.gz";
//! let mut tbx_reader = tbx::Reader::from_path(&path_bed)
//!     .expect(&format!("Could not open {}", path_bed));
//!
//! // Resolve chromosome name to numeric ID.
//! let tid = match tbx_reader.tid("chr1") {
//!     Ok(tid) => tid,
//!     Err(_) => panic!("Could not resolve 'chr1' to contig ID"),
//! };
//!
//! // Set region to fetch.
//! tbx_reader
//!     .fetch(tid, 0, 100_000)
//!     .expect("Could not seek to chr1:1-100,000");
//!
//! // Read through all records in region.
//! for record in tbx_reader.records() {
//!     // ... actually do some work
//! }
//! ```

use std::ffi;
use std::path::Path;
use std::ptr;
use url::Url;

use crate::errors::{Error, Result};
use crate::htslib;
use crate::utils::path_as_bytes;

/// A trait for a Tabix reader with a read method.
pub trait Read: Sized {
    /// Read next line into the given `Vec<u8>` (i.e., ASCII string).
    ///
    /// Use this method in combination with a single allocated record to avoid the reallocations
    /// occurring with the iterator.
    ///
    /// # Arguments
    ///
    /// * `record` - the `Vec<u8>` to be filled
    ///
    /// # Returns
    /// Ok(true) if record was read, Ok(false) if no more record in file
    fn read(&mut self, record: &mut Vec<u8>) -> Result<bool>;

    /// Iterator over the lines/records of the seeked region.
    ///
    /// Note that, while being convenient, this is less efficient than pre-allocating a
    /// `Vec<u8>` and reading into it with the `read()` method, since every iteration involves
    /// the allocation of a new `Vec<u8>`.
    fn records(&mut self) -> Records<'_, Self>;

    /// Return the text headers, split by line.
    fn header(&self) -> &Vec<String>;
}

/// A Tabix file reader.
///
/// This struct and its associated functions are meant for reading plain-text tabix indexed
/// by `tabix`.
///
/// Note that the `tabix` command from `htslib` can actually several more things, including
/// building indices and converting BCF to VCF text output.  Both is out of scope here.
#[derive(Debug)]
pub struct Reader {
    /// The header lines (if any).
    header: Vec<String>,

    /// The file to read from.
    hts_file: *mut htslib::htsFile,
    /// The file format information.
    hts_format: htslib::htsExactFormat,
    /// The tbx_t structure to read from.
    tbx: *mut htslib::tbx_t,
    /// The current buffer.
    buf: htslib::kstring_t,
    /// Iterator over the buffer.
    itr: Option<*mut htslib::hts_itr_t>,

    /// The currently fetch region's tid.
    tid: i64,
    /// The currently fetch region's 0-based begin pos.
    start: i64,
    /// The currently fetch region's 0-based end pos.
    end: i64,
}

unsafe impl Send for Reader {}

/// Redefinition of `KS_SEP_LINE` from `htslib/kseq.h`.
const KS_SEP_LINE: i32 = 2;

impl Reader {
    /// Create a new Reader from path.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::new(&path_as_bytes(path, true)?)
    }

    pub fn from_url(url: &Url) -> Result<Self> {
        Self::new(url.as_str().as_bytes())
    }

    /// Create a new Reader.
    ///
    /// # Arguments
    ///
    /// * `path` - the path.
    fn new(path: &[u8]) -> Result<Self> {
        let path = ffi::CString::new(path).unwrap();
        let c_str = ffi::CString::new("r").unwrap();
        let hts_file = unsafe { htslib::hts_open(path.as_ptr(), c_str.as_ptr()) };
        let hts_format: u32 = unsafe {
            let file_format: *const hts_sys::htsFormat = htslib::hts_get_format(hts_file);
            (*file_format).format
        };

        let tbx = unsafe { htslib::tbx_index_load(path.as_ptr()) };
        if tbx.is_null() {
            return Err(Error::TabixInvalidIndex);
        }
        let mut header = Vec::new();
        let mut buf = htslib::kstring_t {
            l: 0,
            m: 0,
            s: ptr::null_mut(),
        };
        unsafe {
            while htslib::hts_getline(hts_file, KS_SEP_LINE, &mut buf) >= 0 {
                if buf.l > 0 && i32::from(*buf.s) == (*tbx).conf.meta_char {
                    header.push(String::from(ffi::CStr::from_ptr(buf.s).to_str().unwrap()));
                } else {
                    break;
                }
            }
        }

        Ok(Reader {
            header,
            hts_file,
            hts_format,
            tbx,
            buf,
            itr: None,
            tid: -1,
            start: -1,
            end: -1,
        })
    }

    /// Get sequence/target ID from sequence name.
    pub fn tid(&self, name: &str) -> Result<u64> {
        let name_cstr = ffi::CString::new(name.as_bytes()).unwrap();
        let res = unsafe { htslib::tbx_name2id(self.tbx, name_cstr.as_ptr()) };
        if res < 0 {
            Err(Error::UnknownSequence {
                sequence: name.to_owned(),
            })
        } else {
            Ok(res as u64)
        }
    }

    /// Fetch region given by numeric sequence number and 0-based begin and end position.
    pub fn fetch(&mut self, tid: u64, start: u64, end: u64) -> Result<()> {
        self.tid = tid as i64;
        self.start = start as i64;
        self.end = end as i64;

        if let Some(itr) = self.itr {
            unsafe {
                htslib::hts_itr_destroy(itr);
            }
        }
        let itr = unsafe {
            htslib::hts_itr_query(
                (*self.tbx).idx,
                tid as i32,
                start as i64,
                end as i64,
                Some(htslib::tbx_readrec),
            )
        };
        if itr.is_null() {
            self.itr = None;
            Err(Error::Fetch)
        } else {
            self.itr = Some(itr);
            Ok(())
        }
    }

    /// Return the sequence contig names.
    pub fn seqnames(&self) -> Vec<String> {
        let mut result = Vec::new();

        let mut nseq: i32 = 0;
        let seqs = unsafe { htslib::tbx_seqnames(self.tbx, &mut nseq) };
        for i in 0..nseq {
            unsafe {
                result.push(String::from(
                    ffi::CStr::from_ptr(*seqs.offset(i as isize))
                        .to_str()
                        .unwrap(),
                ));
            }
        }
        unsafe {
            libc::free(seqs as *mut libc::c_void);
        };

        result
    }

    /// Activate multi-threaded BGZF read support in htslib. This should permit faster
    /// reading of large BGZF files.
    ///
    /// # Arguments
    ///
    /// * `n_threads` - number of extra background reader threads to use
    pub fn set_threads(&mut self, n_threads: usize) -> Result<()> {
        assert!(n_threads > 0, "n_threads must be > 0");

        let r = unsafe { htslib::hts_set_threads(self.hts_file, n_threads as i32) };
        if r != 0 {
            Err(Error::SetThreads)
        } else {
            Ok(())
        }
    }
}

/// Return whether the two given genomic intervals overlap.
fn overlap(tid1: i64, begin1: i64, end1: i64, tid2: i64, begin2: i64, end2: i64) -> bool {
    (tid1 == tid2) && (begin1 < end2) && (begin2 < end1)
}

impl Read for Reader {
    fn read(&mut self, record: &mut Vec<u8>) -> Result<bool> {
        match self.itr {
            Some(itr) => {
                loop {
                    // Try to read next line.
                    let ret = unsafe {
                        htslib::hts_itr_next(
                            htslib::hts_get_bgzfp(self.hts_file),
                            itr,
                            //mem::transmute(&mut self.buf),
                            &mut self.buf as *mut htslib::kstring_t as *mut libc::c_void,
                            //mem::transmute(self.tbx),
                            self.tbx as *mut libc::c_void,
                        )
                    };
                    // Handle errors first.
                    if ret == -1 {
                        return Ok(false);
                    } else if ret == -2 {
                        return Err(Error::TabixTruncatedRecord);
                    } else if ret < 0 {
                        panic!("Return value should not be <0 but was: {}", ret);
                    }
                    // Return first overlapping record (loop will stop when `hts_itr_next(...)`
                    // returns `< 0`).
                    let (tid, start, end) =
                        unsafe { ((*itr).curr_tid, (*itr).curr_beg, (*itr).curr_end) };
                    // XXX: Careful with this tid conversion!!!
                    if overlap(self.tid, self.start, self.end, tid as i64, start, end) {
                        *record =
                            unsafe { Vec::from(ffi::CStr::from_ptr(self.buf.s).to_str().unwrap()) };
                        return Ok(true);
                    }
                }
            }
            _ => Err(Error::TabixNoIter),
        }
    }

    fn records(&mut self) -> Records<'_, Self> {
        Records { reader: self }
    }

    fn header(&self) -> &Vec<String> {
        &self.header
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        unsafe {
            if self.itr.is_some() {
                htslib::hts_itr_destroy(self.itr.unwrap());
            }
            htslib::tbx_destroy(self.tbx);
            htslib::hts_close(self.hts_file);
        }
    }
}

/// Iterator over the lines of a tabix file.
#[derive(Debug)]
pub struct Records<'a, R: Read> {
    reader: &'a mut R,
}

impl<'a, R: Read> Iterator for Records<'a, R> {
    type Item = Result<Vec<u8>>;

    #[allow(clippy::read_zero_byte_vec)]
    fn next(&mut self) -> Option<Result<Vec<u8>>> {
        let mut record = Vec::new();
        match self.reader.read(&mut record) {
            Ok(false) => None,
            Ok(true) => Some(Ok(record)),
            Err(err) => Some(Err(err)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bed_basic() {
        let reader =
            Reader::from_path("test/tabix_reader/test_bed3.bed.gz").expect("Error opening file.");

        // Check sequence name vector.
        assert_eq!(
            reader.seqnames(),
            vec![String::from("chr1"), String::from("chr2")]
        );

        // Check mapping between name and idx.
        assert_eq!(reader.tid("chr1").unwrap(), 0);
        assert_eq!(reader.tid("chr2").unwrap(), 1);
        assert!(reader.tid("chr3").is_err());
    }

    #[test]
    fn bed_fetch_from_chr1_read_api() {
        let mut reader =
            Reader::from_path("test/tabix_reader/test_bed3.bed.gz").expect("Error opening file.");

        let chr1_id = reader.tid("chr1").unwrap();
        assert!(reader.fetch(chr1_id, 1000, 1003).is_ok());

        let mut record = Vec::new();
        assert!(reader.read(&mut record).is_ok());
        assert_eq!(record, Vec::from("chr1\t1001\t1002"));
        assert_eq!(reader.read(&mut record), Ok(false)); // EOF
    }

    #[test]
    fn bed_fetch_from_chr1_iterator_api() {
        let mut reader =
            Reader::from_path("test/tabix_reader/test_bed3.bed.gz").expect("Error opening file.");

        let chr1_id = reader.tid("chr1").unwrap();
        assert!(reader.fetch(chr1_id, 1000, 1003).is_ok());

        let records: Vec<Vec<u8>> = reader.records().map(|r| r.unwrap()).collect();
        assert_eq!(records, vec![Vec::from("chr1\t1001\t1002")]);
    }

    #[test]
    fn test_fails_on_bam() {
        let reader = Reader::from_path("test/test.bam");
        assert!(reader.is_err());
    }

    #[test]
    fn test_fails_on_non_existiant() {
        let reader = Reader::from_path("test/no_such_file");
        assert!(reader.is_err());
    }

    #[test]
    fn test_fails_on_vcf() {
        let reader = Reader::from_path("test/test_left.vcf");
        assert!(reader.is_err());
    }

    #[test]
    fn test_text_header_regions() {
        // This file has chromosome, start, and end positions with a header line.
        Reader::from_path("test/tabix_reader/genomic_regions_header.txt.gz")
            .expect("Error opening file.");
    }

    #[test]
    fn test_text_header_positions() {
        // This file has chromosome and position with a header line, indexed with
        // `tabix -b2 -e2 <file>`.
        Reader::from_path("test/tabix_reader/genomic_positions_header.txt.gz")
            .expect("Error opening file.");
    }

    #[test]
    fn test_text_bad_header() {
        // This is a duplicate of the above file but the index file is nonsense text.
        Reader::from_path("test/tabix_reader/bad_header.txt.gz")
            .expect_err("Invalid index file should fail.");
    }
}
