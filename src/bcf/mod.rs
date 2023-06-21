// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Module for working with VCF and BCF files.
//!
//! # Performance Remarks
//!
//! Note that BCF corresponds to the in-memory representation of BCF/VCF records in Htslib
//! itself. Thus, it comes without a runtime penalty for parsing, in contrast to reading VCF
//! files.
//!
//! # Example (reading)
//!
//!   - Obtaining 0-based locus index of the VCF record.
//!   - Obtaining alleles of the VCF record.
//!   - calculate alt-allele dosage in a mutli-sample VCF / BCF
//!
//! ```
//! use crate::rust_htslib::bcf::{Reader, Read};
//! use std::convert::TryFrom;
//!
//! let path = &"test/test_string.vcf";
//! let mut bcf = Reader::from_path(path).expect("Error opening file.");
//! // iterate through each row of the vcf body.
//! for (i, record_result) in bcf.records().enumerate() {
//!     let mut record = record_result.expect("Fail to read record");
//!     let mut s = String::new();
//!      for allele in record.alleles() {
//!          for c in allele {
//!              s.push(char::from(*c))
//!          }
//!          s.push(' ')
//!      }
//!     // 0-based position and the list of alleles
//!     println!("Locus: {}, Alleles: {}", record.pos(), s);
//!     // number of sample in the vcf
//!     let sample_count = usize::try_from(record.sample_count()).unwrap();
//!
//!     // Counting ref, alt and missing alleles for each sample
//!     let mut n_ref = vec![0; sample_count];
//!     let mut n_alt = vec![0; sample_count];
//!     let mut n_missing = vec![0; sample_count];
//!     let gts = record.genotypes().expect("Error reading genotypes");
//!     for sample_index in 0..sample_count {
//!         // for each sample
//!         for gta in gts.get(sample_index).iter() {
//!             // for each allele
//!             match gta.index() {
//!                 Some(0) => n_ref[sample_index] += 1,  // reference allele
//!                 Some(_) => n_alt[sample_index] += 1,  // alt allele
//!                 None => n_missing[sample_index] += 1, // missing allele
//!             }
//!         }
//!     }
//! }
//! ```
//!
//! # Example (writing)
//!
//!   - Setting up a VCF writer from scratch (including a simple header)
//!   - Creating a VCF record and writing it to the VCF file
//!
//! ```
//! use rust_htslib::bcf::{Format, Writer};
//! use rust_htslib::bcf::header::Header;
//! use rust_htslib::bcf::record::GenotypeAllele;
//!
//! // Create minimal VCF header with a single contig and a single sample
//! let mut header = Header::new();
//! let header_contig_line = r#"##contig=<ID=1,length=10>"#;
//! header.push_record(header_contig_line.as_bytes());
//! let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
//! header.push_record(header_gt_line.as_bytes());
//! header.push_sample("test_sample".as_bytes());
//!
//! // Write uncompressed VCF to stdout with above header and get an empty record
//! let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
//! let mut record = vcf.empty_record();
//!
//! // Set chrom and pos to 1 and 7, respectively - note the 0-based positions
//! let rid = vcf.header().name2rid(b"1").unwrap();
//! record.set_rid(Some(rid));
//! record.set_pos(6);
//!
//! // Set record genotype to 0|1 - note first allele is always unphased
//! let alleles = &[GenotypeAllele::Unphased(0), GenotypeAllele::Phased(1)];
//! record.push_genotypes(alleles).unwrap();
//!
//! // Write record
//! vcf.write(&record).unwrap()
//! ```
//!
//! This will print the following VCF to stdout:
//!
//! ```lang-none
//! ##fileformat=VCFv4.2
//! ##FILTER=<ID=PASS,Description="All filters passed">
//! ##contig=<ID=1,length=10>
//! ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
//! #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  test_sample
//! 1       7       .       .       .       0       .       .       GT      0|1
//! ```

use std::ffi;
use std::path::Path;
use std::rc::Rc;
use std::str;

use url::Url;

pub mod buffer;
pub mod header;
pub mod record;

use crate::bcf::header::{HeaderView, SampleSubset};
use crate::errors::{Error, Result};
use crate::htslib;

pub use crate::bcf::header::{Header, HeaderRecord};
pub use crate::bcf::record::Record;

/// A trait for a BCF reader with a read method.
pub trait Read: Sized {
    /// Read the next record.
    ///
    /// # Arguments
    /// * record - an empty record, that can be created with `bcf::Reader::empty_record`.
    ///
    /// # Returns
    /// None if end of file was reached, otherwise Some will contain
    /// a result with an error in case of failure.
    fn read(&mut self, record: &mut record::Record) -> Option<Result<()>>;

    /// Return an iterator over all records of the VCF/BCF file.
    fn records(&mut self) -> Records<'_, Self>;

    /// Return the header.
    fn header(&self) -> &HeaderView;

    /// Return empty record.  Can be reused multiple times.
    fn empty_record(&self) -> Record;

    /// Activate multi-threaded BCF/VCF read support in htslib. This should permit faster
    /// reading of large VCF files.
    ///
    /// Setting `nthreads` to `0` does not change the current state.  Note that it is not
    /// possible to set the number of background threads below `1` once it has been set.
    ///
    /// # Arguments
    ///
    /// * `n_threads` - number of extra background writer threads to use, must be `> 0`.
    fn set_threads(&mut self, n_threads: usize) -> Result<()>;
}

/// A VCF/BCF reader.
#[derive(Debug)]
pub struct Reader {
    inner: *mut htslib::htsFile,
    header: Rc<HeaderView>,
}

unsafe impl Send for Reader {}
/// # Safety
///
/// Implementation for `Reader::set_threads()` and `Writer::set_threads`.
pub unsafe fn set_threads(hts_file: *mut htslib::htsFile, n_threads: usize) -> Result<()> {
    assert!(n_threads > 0, "n_threads must be > 0");

    let r = htslib::hts_set_threads(hts_file, n_threads as i32);
    if r != 0 {
        Err(Error::SetThreads)
    } else {
        Ok(())
    }
}

impl Reader {
    /// Create a new reader from a given path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        match path.as_ref().to_str() {
            Some(p) if !path.as_ref().exists() => Err(Error::FileNotFound { path: p.into() }),
            Some(p) => Self::new(p.as_bytes()),
            _ => Err(Error::NonUnicodePath),
        }
    }

    /// Create a new reader from a given URL.
    pub fn from_url(url: &Url) -> Result<Self> {
        Self::new(url.as_str().as_bytes())
    }

    /// Create a new reader from standard input.
    pub fn from_stdin() -> Result<Self> {
        Self::new(b"-")
    }

    fn new(path: &[u8]) -> Result<Self> {
        let htsfile = bcf_open(path, b"r")?;
        let header = unsafe { htslib::bcf_hdr_read(htsfile) };
        Ok(Reader {
            inner: htsfile,
            header: Rc::new(HeaderView::new(header)),
        })
    }
}

impl Read for Reader {
    fn read(&mut self, record: &mut record::Record) -> Option<Result<()>> {
        match unsafe { htslib::bcf_read(self.inner, self.header.inner, record.inner) } {
            0 => {
                unsafe {
                    // Always unpack record.
                    htslib::bcf_unpack(record.inner_mut(), htslib::BCF_UN_ALL as i32);
                }
                record.set_header(Rc::clone(&self.header));
                Some(Ok(()))
            }
            -1 => None,
            _ => Some(Err(Error::BcfInvalidRecord)),
        }
    }

    fn records(&mut self) -> Records<'_, Self> {
        Records { reader: self }
    }

    fn set_threads(&mut self, n_threads: usize) -> Result<()> {
        unsafe { set_threads(self.inner, n_threads) }
    }

    fn header(&self) -> &HeaderView {
        &self.header
    }

    /// Return empty record.  Can be reused multiple times.
    fn empty_record(&self) -> Record {
        Record::new(Rc::clone(&self.header))
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        unsafe {
            htslib::hts_close(self.inner);
        }
    }
}

/// An indexed VCF/BCF reader.
#[derive(Debug)]
pub struct IndexedReader {
    /// The synced VCF/BCF reader to use internally.
    inner: *mut htslib::bcf_srs_t,
    /// The header.
    header: Rc<HeaderView>,

    /// The position of the previous fetch, if any.
    current_region: Option<(u32, u64, Option<u64>)>,
}

unsafe impl Send for IndexedReader {}

impl IndexedReader {
    /// Create a new `IndexedReader` from path.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        match path.to_str() {
            Some(p) if path.exists() => {
                Self::new(&ffi::CString::new(p).map_err(|_| Error::NonUnicodePath)?)
            }
            Some(p) => Err(Error::FileNotFound { path: p.into() }),
            None => Err(Error::NonUnicodePath),
        }
    }

    /// Create a new `IndexedReader` from an URL.
    pub fn from_url(url: &Url) -> Result<Self> {
        Self::new(&ffi::CString::new(url.as_str()).unwrap())
    }

    /// Create a new `IndexedReader`.
    ///
    /// # Arguments
    ///
    /// * `path` - the path. Use "-" for stdin.
    fn new(path: &ffi::CStr) -> Result<Self> {
        // Create reader and require existence of index file.
        let ser_reader = unsafe { htslib::bcf_sr_init() };
        unsafe {
            htslib::bcf_sr_set_opt(ser_reader, 0);
        } // 0: BCF_SR_REQUIRE_IDX
          // Attach a file with the path from the arguments.
        if unsafe { htslib::bcf_sr_add_reader(ser_reader, path.as_ptr()) } >= 0 {
            let header = Rc::new(HeaderView::new(unsafe {
                htslib::bcf_hdr_dup((*(*ser_reader).readers.offset(0)).header)
            }));
            Ok(IndexedReader {
                inner: ser_reader,
                header,
                current_region: None,
            })
        } else {
            Err(Error::BcfOpen {
                target: path.to_str().unwrap().to_owned(),
            })
        }
    }

    /// Jump to the given region.
    ///
    /// # Arguments
    ///
    /// * `rid` - numeric ID of the reference to jump to; use `HeaderView::name2rid` for resolving
    ///           contig name to ID.
    /// * `start` - `0`-based **inclusive** start coordinate of region on reference.
    /// * `end` - Optional `0`-based **inclusive** end coordinate of region on reference. If `None`
    /// is given, records are fetched from `start` until the end of the contig.
    ///
    /// # Note
    /// The entire contig can be fetched by setting `start` to `0` and `end` to `None`.
    pub fn fetch(&mut self, rid: u32, start: u64, end: Option<u64>) -> Result<()> {
        let contig = self.header.rid2name(rid)?;
        let contig = ffi::CString::new(contig).unwrap();
        if unsafe { htslib::bcf_sr_seek(self.inner, contig.as_ptr(), start as i64) } != 0 {
            Err(Error::GenomicSeek {
                contig: contig.to_str().unwrap().to_owned(),
                start,
            })
        } else {
            self.current_region = Some((rid, start, end));
            Ok(())
        }
    }
}

impl Read for IndexedReader {
    fn read(&mut self, record: &mut record::Record) -> Option<Result<()>> {
        match unsafe { htslib::bcf_sr_next_line(self.inner) } {
            0 => {
                if unsafe { (*self.inner).errnum } != 0 {
                    Some(Err(Error::BcfInvalidRecord))
                } else {
                    None
                }
            }
            i => {
                assert!(i > 0, "Must not be negative");
                // Note that the sync BCF reader has a different interface than the others
                // as it keeps its own buffer already for each record.  An alternative here
                // would be to replace the `inner` value by an enum that can be a pointer
                // into a synced reader or an owning popinter to an allocated record.
                unsafe {
                    htslib::bcf_copy(
                        record.inner,
                        *(*(*self.inner).readers.offset(0)).buffer.offset(0),
                    );
                }

                unsafe {
                    // Always unpack record.
                    htslib::bcf_unpack(record.inner_mut(), htslib::BCF_UN_ALL as i32);
                }

                record.set_header(Rc::clone(&self.header));

                match self.current_region {
                    Some((rid, _start, end)) => {
                        let endpos = match end {
                            Some(e) => e,
                            None => u64::MAX,
                        };
                        if Some(rid) == record.rid() && record.pos() as u64 <= endpos {
                            Some(Ok(()))
                        } else {
                            None
                        }
                    }
                    None => Some(Ok(())),
                }
            }
        }
    }

    fn records(&mut self) -> Records<'_, Self> {
        Records { reader: self }
    }

    fn set_threads(&mut self, n_threads: usize) -> Result<()> {
        assert!(n_threads > 0, "n_threads must be > 0");

        let r = unsafe { htslib::bcf_sr_set_threads(self.inner, n_threads as i32) };
        if r != 0 {
            Err(Error::SetThreads)
        } else {
            Ok(())
        }
    }

    fn header(&self) -> &HeaderView {
        &self.header
    }

    fn empty_record(&self) -> Record {
        Record::new(Rc::clone(&self.header))
    }
}

impl Drop for IndexedReader {
    fn drop(&mut self) {
        unsafe { htslib::bcf_sr_destroy(self.inner) };
    }
}

/// This module contains the `SyncedReader` class and related code.
pub mod synced {

    use super::*;

    /// This module contains bitmask constants for `SyncedReader`.
    pub mod pairing {
        /// Allow different alleles, as long as they all are SNPs.
        pub const SNPS: u32 = crate::htslib::BCF_SR_PAIR_SNPS;
        /// The same as above, but with indels.
        pub const INDELS: u32 = crate::htslib::BCF_SR_PAIR_INDELS;
        /// Any combination of alleles can be returned by `bcf_sr_next_line()`.
        pub const ANY: u32 = crate::htslib::BCF_SR_PAIR_ANY;
        /// At least some of multiallelic ALTs must match.  Implied by all the others with the exception of `EXACT`.
        pub const SOME: u32 = crate::htslib::BCF_SR_PAIR_SOME;
        /// Allow REF-only records with SNPs.
        pub const SNP_REF: u32 = crate::htslib::BCF_SR_PAIR_SNP_REF;
        /// Allow REF-only records with indels.
        pub const INDEL_REF: u32 = crate::htslib::BCF_SR_PAIR_INDEL_REF;
        /// Require the exact same set of alleles in all files.
        pub const EXACT: u32 = crate::htslib::BCF_SR_PAIR_EXACT;
        /// `SNPS | INDELS`.
        pub const BOTH: u32 = crate::htslib::BCF_SR_PAIR_BOTH;
        /// `SNPS | INDELS | SNP_REF | INDEL_REF`.
        pub const BOTH_REF: u32 = crate::htslib::BCF_SR_PAIR_BOTH_REF;
    }

    /// A wrapper for `bcf_srs_t`; allows joint traversal of multiple VCF and/or BCF files.
    #[derive(Debug)]
    pub struct SyncedReader {
        /// Internal handle for the synced reader.
        inner: *mut crate::htslib::bcf_srs_t,

        /// RC's of `HeaderView`s of the readers.
        headers: Vec<Rc<HeaderView>>,

        /// The position of the previous fetch, if any.
        current_region: Option<(u32, u64, u64)>,
    }

    // TODO: add interface for setting threads, ensure that the pool is freed properly
    impl SyncedReader {
        pub fn new() -> Result<Self> {
            let inner = unsafe { crate::htslib::bcf_sr_init() };
            if inner.is_null() {
                return Err(Error::BcfAllocationError);
            }

            Ok(SyncedReader {
                inner,
                headers: Vec::new(),
                current_region: None,
            })
        }

        /// Enable or disable requiring of index
        pub fn set_require_index(&mut self, do_require: bool) {
            unsafe {
                (*self.inner).require_index = if do_require { 1 } else { 0 };
            }
        }

        /// Set the given bitmask of values from `sr_pairing` module.
        pub fn set_pairing(&mut self, bitmask: u32) {
            unsafe {
                // TODO: 1 actually is BCF_SR_PAIR_LOGIC but is not available here?
                crate::htslib::bcf_sr_set_opt(self.inner, 1, bitmask);
            }
        }

        /// Add new reader with the path to the file.
        pub fn add_reader<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
            match path.as_ref().to_str() {
                Some(p) if path.as_ref().exists() => {
                    let p_cstring = ffi::CString::new(p).unwrap();
                    let res =
                        unsafe { crate::htslib::bcf_sr_add_reader(self.inner, p_cstring.as_ptr()) };

                    if res == 0 {
                        return Err(Error::BcfOpen {
                            target: p.to_owned(),
                        });
                    }

                    let i = (self.reader_count() - 1) as isize;
                    let header = Rc::new(HeaderView::new(unsafe {
                        crate::htslib::bcf_hdr_dup((*(*self.inner).readers.offset(i)).header)
                    }));
                    self.headers.push(header);
                    Ok(())
                }
                _ => Err(Error::NonUnicodePath),
            }
        }

        /// Remove reader with the given index.
        pub fn remove_reader(&mut self, idx: u32) {
            if idx >= self.reader_count() {
                panic!("Invalid reader!");
            } else {
                unsafe {
                    crate::htslib::bcf_sr_remove_reader(self.inner, idx as i32);
                }
                self.headers.remove(idx as usize);
            }
        }

        /// Return number of open files/readers.
        pub fn reader_count(&self) -> u32 {
            unsafe { (*self.inner).nreaders as u32 }
        }

        /// Read next line and return number of readers that have the given line (0 if end of all files is reached).
        pub fn read_next(&mut self) -> Result<u32> {
            let num = unsafe { crate::htslib::bcf_sr_next_line(self.inner) as u32 };

            if num == 0 {
                if unsafe { (*self.inner).errnum } != 0 {
                    return Err(Error::BcfInvalidRecord);
                }
                Ok(0)
            } else {
                assert!(num > 0, "num returned by htslib must not be negative");
                match self.current_region {
                    Some((rid, _start, end)) => {
                        for idx in 0..self.reader_count() {
                            if !self.has_line(idx) {
                                continue;
                            }
                            unsafe {
                                let record = *(*(*self.inner).readers.offset(idx as isize))
                                    .buffer
                                    .offset(0);
                                if (*record).rid != (rid as i32) || (*record).pos >= (end as i64) {
                                    return Ok(0);
                                }
                            }
                        }
                        Ok(num)
                    }
                    None => Ok(num),
                }
            }
        }

        /// Return whether the given reader has the line.
        pub fn has_line(&self, idx: u32) -> bool {
            if idx >= self.reader_count() {
                panic!("Invalid reader!");
            } else {
                unsafe { (*(*self.inner).has_line.offset(idx as isize)) != 0 }
            }
        }

        /// Return record from the given reader, if any.
        pub fn record(&self, idx: u32) -> Option<Record> {
            if self.has_line(idx) {
                let record = Record::new(self.headers[idx as usize].clone());
                unsafe {
                    crate::htslib::bcf_copy(
                        record.inner,
                        *(*(*self.inner).readers.offset(idx as isize))
                            .buffer
                            .offset(0),
                    );
                }
                Some(record)
            } else {
                None
            }
        }

        /// Return header from the given reader.
        pub fn header(&self, idx: u32) -> &HeaderView {
            // TODO: is the mutability here correct?
            if idx >= self.reader_count() {
                panic!("Invalid reader!");
            } else {
                &self.headers[idx as usize]
            }
        }

        /// Jump to the given region.
        ///
        /// # Arguments
        ///
        /// * `rid` - numeric ID of the reference to jump to; use `HeaderView::name2rid` for resolving
        ///           contig name to ID.
        /// * `start` - `0`-based start coordinate of region on reference.
        /// * `end` - `0`-based end coordinate of region on reference.
        pub fn fetch(&mut self, rid: u32, start: u64, end: u64) -> Result<()> {
            let contig = {
                let contig = self.header(0).rid2name(rid).unwrap(); //.clone();
                ffi::CString::new(contig).unwrap()
            };
            if unsafe { htslib::bcf_sr_seek(self.inner, contig.as_ptr(), start as i64) } != 0 {
                Err(Error::GenomicSeek {
                    contig: contig.to_str().unwrap().to_owned(),
                    start,
                })
            } else {
                self.current_region = Some((rid, start, end));
                Ok(())
            }
        }
    }

    impl Drop for SyncedReader {
        fn drop(&mut self) {
            unsafe { crate::htslib::bcf_sr_destroy(self.inner) };
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub enum Format {
    Vcf,
    Bcf,
}

/// A VCF/BCF writer.
#[derive(Debug)]
pub struct Writer {
    inner: *mut htslib::htsFile,
    header: Rc<HeaderView>,
    subset: Option<SampleSubset>,
}

unsafe impl Send for Writer {}

impl Writer {
    /// Create a new writer that writes to the given path.
    ///
    /// # Arguments
    ///
    /// * `path` - the path
    /// * `header` - header definition to use
    /// * `uncompressed` - disable compression
    /// * `vcf` - write VCF instead of BCF
    pub fn from_path<P: AsRef<Path>>(
        path: P,
        header: &Header,
        uncompressed: bool,
        format: Format,
    ) -> Result<Self> {
        if let Some(p) = path.as_ref().to_str() {
            Ok(Self::new(p.as_bytes(), header, uncompressed, format)?)
        } else {
            Err(Error::NonUnicodePath)
        }
    }

    /// Create a new writer from a URL.
    ///
    /// # Arguments
    ///
    /// * `url` - the URL
    /// * `header` - header definition to use
    /// * `uncompressed` - disable compression
    /// * `vcf` - write VCF instead of BCF
    pub fn from_url(
        url: &Url,
        header: &Header,
        uncompressed: bool,
        format: Format,
    ) -> Result<Self> {
        Self::new(url.as_str().as_bytes(), header, uncompressed, format)
    }

    /// Create a new writer to stdout.
    ///
    /// # Arguments
    ///
    /// * `header` - header definition to use
    /// * `uncompressed` - disable compression
    /// * `vcf` - write VCF instead of BCF
    pub fn from_stdout(header: &Header, uncompressed: bool, format: Format) -> Result<Self> {
        Self::new(b"-", header, uncompressed, format)
    }

    fn new(path: &[u8], header: &Header, uncompressed: bool, format: Format) -> Result<Self> {
        let mode: &[u8] = match (uncompressed, format) {
            (true, Format::Vcf) => b"w",
            (false, Format::Vcf) => b"wz",
            (true, Format::Bcf) => b"wbu",
            (false, Format::Bcf) => b"wb",
        };

        let htsfile = bcf_open(path, mode)?;
        unsafe { htslib::bcf_hdr_write(htsfile, header.inner) };
        Ok(Writer {
            inner: htsfile,
            header: Rc::new(HeaderView::new(unsafe {
                htslib::bcf_hdr_dup(header.inner)
            })),
            subset: header.subset.clone(),
        })
    }

    /// Obtain reference to the lightweight `HeaderView` of the BCF header.
    pub fn header(&self) -> &HeaderView {
        &self.header
    }

    /// Create empty record for writing to this writer.
    ///
    /// This record can then be reused multiple times.
    pub fn empty_record(&self) -> Record {
        record::Record::new(Rc::clone(&self.header))
    }

    /// Translate record to header of this writer.
    ///
    /// # Arguments
    ///
    /// - `record` - The `Record` to translate.
    pub fn translate(&mut self, record: &mut record::Record) {
        unsafe {
            htslib::bcf_translate(self.header.inner, record.header().inner, record.inner);
        }
        record.set_header(Rc::clone(&self.header));
    }

    /// Subset samples of record to match header of this writer.
    ///
    /// # Arguments
    ///
    /// - `record` - The `Record` to modify.
    pub fn subset(&mut self, record: &mut record::Record) {
        if let Some(ref mut subset) = self.subset {
            unsafe {
                htslib::bcf_subset(
                    self.header.inner,
                    record.inner,
                    subset.len() as i32,
                    subset.as_mut_ptr(),
                );
            }
        }
    }

    /// Write `record` to the Writer.
    ///
    /// # Arguments
    ///
    /// - `record` - The `Record` to write.
    pub fn write(&mut self, record: &record::Record) -> Result<()> {
        if unsafe { htslib::bcf_write(self.inner, self.header.inner, record.inner) } == -1 {
            Err(Error::WriteRecord)
        } else {
            Ok(())
        }
    }

    /// Activate multi-threaded BCF write support in htslib. This should permit faster
    /// writing of large BCF files.
    ///
    /// # Arguments
    ///
    /// * `n_threads` - number of extra background writer threads to use, must be `> 0`.
    pub fn set_threads(&mut self, n_threads: usize) -> Result<()> {
        unsafe { set_threads(self.inner, n_threads) }
    }
}

impl Drop for Writer {
    fn drop(&mut self) {
        unsafe {
            htslib::hts_close(self.inner);
        }
    }
}

#[derive(Debug)]
pub struct Records<'a, R: Read> {
    reader: &'a mut R,
}

impl<'a, R: Read> Iterator for Records<'a, R> {
    type Item = Result<record::Record>;

    fn next(&mut self) -> Option<Result<record::Record>> {
        let mut record = self.reader.empty_record();
        match self.reader.read(&mut record) {
            Some(Err(e)) => Some(Err(e)),
            Some(Ok(_)) => Some(Ok(record)),
            None => None,
        }
    }
}

/// Wrapper for opening a BCF file.
fn bcf_open(target: &[u8], mode: &[u8]) -> Result<*mut htslib::htsFile> {
    let p = ffi::CString::new(target).unwrap();
    let c_str = ffi::CString::new(mode).unwrap();
    let ret = unsafe { htslib::hts_open(p.as_ptr(), c_str.as_ptr()) };

    if ret.is_null() {
        return Err(Error::BcfOpen {
            target: str::from_utf8(target).unwrap().to_owned(),
        });
    }

    unsafe {
        if !(mode.contains(&b'w')
            || (*ret).format.category == htslib::htsFormatCategory_variant_data)
        {
            return Err(Error::BcfOpen {
                target: str::from_utf8(target).unwrap().to_owned(),
            });
        }
    }
    Ok(ret)
}

#[cfg(test)]
mod tests {
    use super::record::Buffer;
    use super::*;
    use crate::bcf::header::Id;
    use crate::bcf::record::GenotypeAllele;
    use crate::bcf::record::Numeric;
    use crate::bcf::Reader;
    use std::convert::TryFrom;
    use std::fs::File;
    use std::io::prelude::Read as IoRead;
    use std::path::Path;
    use std::str;

    fn _test_read<P: AsRef<Path>>(path: &P) {
        let mut bcf = Reader::from_path(path).expect("Error opening file.");
        assert_eq!(bcf.header.samples(), [b"NA12878.subsample-0.25-0"]);

        for (i, rec) in bcf.records().enumerate() {
            let record = rec.expect("Error reading record.");
            assert_eq!(record.sample_count(), 1);

            assert_eq!(record.rid().expect("Error reading rid."), 0);
            assert_eq!(record.pos(), 10021 + i as i64);
            assert!((record.qual() - 0f32).abs() < std::f32::EPSILON);
            let mut buffer = Buffer::new();
            assert!(
                (record
                    .info_shared_buffer(b"MQ0F", &mut buffer)
                    .float()
                    .expect("Error reading info.")
                    .expect("Missing tag")[0]
                    - 1.0)
                    .abs()
                    < std::f32::EPSILON
            );
            if i == 59 {
                assert!(
                    (record
                        .info_shared_buffer(b"SGB", &mut buffer)
                        .float()
                        .expect("Error reading info.")
                        .expect("Missing tag")[0]
                        - -0.379885)
                        .abs()
                        < std::f32::EPSILON
                );
            }
            // the artificial "not observed" allele is present in each record.
            assert_eq!(record.alleles().iter().last().unwrap(), b"<X>");

            let fmt = record.format(b"PL");
            let pl = fmt.integer().expect("Error reading format.");
            assert_eq!(pl.len(), 1);
            if i == 59 {
                assert_eq!(pl[0].len(), 6);
            } else {
                assert_eq!(pl[0].len(), 3);
            }
        }
    }

    #[test]
    fn test_read() {
        _test_read(&"test/test.bcf");
    }

    #[test]
    fn test_reader_set_threads() {
        let path = &"test/test.bcf";
        let mut bcf = Reader::from_path(path).expect("Error opening file.");
        bcf.set_threads(2).unwrap();
    }

    #[test]
    fn test_writer_set_threads() {
        let path = &"test/test.bcf";
        let tmp = tempfile::Builder::new()
            .prefix("rust-htslib")
            .tempdir()
            .expect("Cannot create temp dir");
        let bcfpath = tmp.path().join("test.bcf");
        let bcf = Reader::from_path(path).expect("Error opening file.");
        let header = Header::from_template_subset(&bcf.header, &[b"NA12878.subsample-0.25-0"])
            .expect("Error subsetting samples.");
        let mut writer =
            Writer::from_path(&bcfpath, &header, false, Format::Bcf).expect("Error opening file.");
        writer.set_threads(2).unwrap();
    }

    #[test]
    fn test_fetch() {
        let mut bcf = IndexedReader::from_path(&"test/test.bcf").expect("Error opening file.");
        bcf.set_threads(2).unwrap();
        let rid = bcf
            .header()
            .name2rid(b"1")
            .expect("Translating from contig '1' to ID failed.");
        bcf.fetch(rid, 10_033, Some(10_060))
            .expect("Fetching failed");
        assert_eq!(bcf.records().count(), 28);
    }

    #[test]
    fn test_fetch_all() {
        let mut bcf = IndexedReader::from_path(&"test/test.bcf").expect("Error opening file.");
        bcf.set_threads(2).unwrap();
        let rid = bcf
            .header()
            .name2rid(b"1")
            .expect("Translating from contig '1' to ID failed.");
        bcf.fetch(rid, 0, None).expect("Fetching failed");
        assert_eq!(bcf.records().count(), 62);
    }

    #[test]
    fn test_fetch_open_ended() {
        let mut bcf = IndexedReader::from_path(&"test/test.bcf").expect("Error opening file.");
        bcf.set_threads(2).unwrap();
        let rid = bcf
            .header()
            .name2rid(b"1")
            .expect("Translating from contig '1' to ID failed.");
        bcf.fetch(rid, 10077, None).expect("Fetching failed");
        assert_eq!(bcf.records().count(), 6);
    }

    #[test]
    fn test_fetch_start_greater_than_last_vcf_pos() {
        let mut bcf = IndexedReader::from_path(&"test/test.bcf").expect("Error opening file.");
        bcf.set_threads(2).unwrap();
        let rid = bcf
            .header()
            .name2rid(b"1")
            .expect("Translating from contig '1' to ID failed.");
        bcf.fetch(rid, 20077, None).expect("Fetching failed");
        assert_eq!(bcf.records().count(), 0);
    }

    #[test]
    fn test_write() {
        let mut bcf = Reader::from_path(&"test/test_multi.bcf").expect("Error opening file.");
        let tmp = tempfile::Builder::new()
            .prefix("rust-htslib")
            .tempdir()
            .expect("Cannot create temp dir");
        let bcfpath = tmp.path().join("test.bcf");
        println!("{:?}", bcfpath);
        {
            let header = Header::from_template_subset(&bcf.header, &[b"NA12878.subsample-0.25-0"])
                .expect("Error subsetting samples.");
            let mut writer = Writer::from_path(&bcfpath, &header, false, Format::Bcf)
                .expect("Error opening file.");
            for rec in bcf.records() {
                let mut record = rec.expect("Error reading record.");
                writer.translate(&mut record);
                writer.subset(&mut record);
                record.trim_alleles().expect("Error trimming alleles.");
                writer.write(&record).expect("Error writing record");
            }
        }
        {
            _test_read(&bcfpath);
        }
        tmp.close().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_strings() {
        let mut vcf = Reader::from_path(&"test/test_string.vcf").expect("Error opening file.");
        let fs1 = [
            &b"LongString1"[..],
            &b"LongString2"[..],
            &b"."[..],
            &b"LongString4"[..],
            &b"evenlength"[..],
            &b"ss6"[..],
        ];
        let mut buffer = Buffer::new();
        for (i, rec) in vcf.records().enumerate() {
            println!("record {}", i);
            let record = rec.expect("Error reading record.");
            assert_eq!(
                record
                    .info_shared_buffer(b"S1", &mut buffer)
                    .string()
                    .expect("Error reading string.")
                    .expect("Missing tag")[0],
                format!("string{}", i + 1).as_bytes()
            );
            let fs1_str_vec = record
                .format_shared_buffer(b"FS1", &mut buffer)
                .string()
                .expect("Error reading string.");
            assert_eq!(fs1_str_vec.len(), 2);
            println!("{}", String::from_utf8_lossy(fs1_str_vec[0]));
            assert_eq!(
                record
                    .format(b"FS1")
                    .string()
                    .expect("Error reading string.")[0],
                fs1[i]
            );
        }
    }

    #[test]
    fn test_missing() {
        let mut vcf = Reader::from_path(&"test/test_missing.vcf").expect("Error opening file.");
        let fn4 = [
            &[
                i32::missing(),
                i32::missing(),
                i32::missing(),
                i32::missing(),
            ][..],
            &[i32::missing()][..],
        ];
        let f1 = [false, true];
        let mut buffer = Buffer::new();
        for (i, rec) in vcf.records().enumerate() {
            let record = rec.expect("Error reading record.");
            assert_eq!(
                record
                    .info_shared_buffer(b"F1", &mut buffer)
                    .float()
                    .expect("Error reading float.")
                    .expect("Missing tag")[0]
                    .is_nan(),
                f1[i]
            );
            assert_eq!(
                record
                    .format(b"FN4")
                    .integer()
                    .expect("Error reading integer.")[1],
                fn4[i]
            );
            assert!(
                record.format(b"FF4").float().expect("Error reading float.")[1]
                    .iter()
                    .all(|&v| v.is_missing())
            );
        }
    }

    #[test]
    fn test_genotypes() {
        let mut vcf = Reader::from_path(&"test/test_string.vcf").expect("Error opening file.");
        let expected = ["./1", "1|1", "0/1", "0|1", "1|.", "1/1"];
        for (rec, exp_gt) in vcf.records().zip(expected.iter()) {
            let rec = rec.expect("Error reading record.");
            let genotypes = rec.genotypes().expect("Error reading genotypes");
            assert_eq!(&format!("{}", genotypes.get(0)), exp_gt);
        }
    }

    #[test]
    fn test_header_ids() {
        let vcf = Reader::from_path(&"test/test_string.vcf").expect("Error opening file.");
        let header = &vcf.header();
        use crate::bcf::header::Id;

        assert_eq!(header.id_to_name(Id(4)), b"GT");
        assert_eq!(header.name_to_id(b"GT").unwrap(), Id(4));
        assert!(header.name_to_id(b"XX").is_err());
    }

    #[test]
    fn test_header_samples() {
        let vcf = Reader::from_path(&"test/test_string.vcf").expect("Error opening file.");
        let header = &vcf.header();

        assert_eq!(header.id_to_sample(Id(0)), b"one");
        assert_eq!(header.id_to_sample(Id(1)), b"two");
        assert_eq!(header.sample_to_id(b"one").unwrap(), Id(0));
        assert_eq!(header.sample_to_id(b"two").unwrap(), Id(1));
        assert!(header.sample_to_id(b"three").is_err());
    }

    #[test]
    fn test_header_contigs() {
        let vcf = Reader::from_path(&"test/test_multi.bcf").expect("Error opening file.");
        let header = &vcf.header();

        assert_eq!(header.contig_count(), 86);

        // test existing contig names and IDs
        assert_eq!(header.rid2name(0).unwrap(), b"1");
        assert_eq!(header.name2rid(b"1").unwrap(), 0);

        assert_eq!(header.rid2name(85).unwrap(), b"hs37d5");
        assert_eq!(header.name2rid(b"hs37d5").unwrap(), 85);

        // test nonexistent contig names and IDs
        assert!(header.name2rid(b"nonexistent_contig").is_err());
        assert!(header.rid2name(100).is_err());
    }

    #[test]
    fn test_header_records() {
        let vcf = Reader::from_path(&"test/test_string.vcf").expect("Error opening file.");
        let records = vcf.header().header_records();
        assert_eq!(records.len(), 10);

        match records[1] {
            HeaderRecord::Filter {
                ref key,
                ref values,
            } => {
                assert_eq!(key, "FILTER");
                assert_eq!(values["ID"], "PASS");
            }
            _ => {
                panic!("Invalid HeaderRecord");
            }
        }
    }

    #[test]
    fn test_header_info_types() {
        let vcf = Reader::from_path(&"test/test.bcf").unwrap();
        let header = vcf.header();
        let truth = vec![
            (
                // INFO=<ID=INDEL,Number=0,Type=Flag>
                "INDEL",
                header::TagType::Flag,
                header::TagLength::Fixed(0),
            ),
            (
                // INFO=<ID=DP,Number=1,Type=Integer>
                "DP",
                header::TagType::Integer,
                header::TagLength::Fixed(1),
            ),
            (
                // INFO=<ID=QS,Number=R,Type=Float>
                "QS",
                header::TagType::Float,
                header::TagLength::Alleles,
            ),
            (
                // INFO=<ID=I16,Number=16,Type=Float>
                "I16",
                header::TagType::Float,
                header::TagLength::Fixed(16),
            ),
        ];
        for (ref_name, ref_type, ref_length) in truth {
            let (tag_type, tag_length) = header.info_type(ref_name.as_bytes()).unwrap();
            assert_eq!(tag_type, ref_type);
            assert_eq!(tag_length, ref_length);
        }

        let vcf = Reader::from_path(&"test/test_svlen.vcf").unwrap();
        let header = vcf.header();
        let truth = vec![
            (
                // INFO=<ID=IMPRECISE,Number=0,Type=Flag>
                "IMPRECISE",
                header::TagType::Flag,
                header::TagLength::Fixed(0),
            ),
            (
                // INFO=<ID=SVTYPE,Number=1,Type=String>
                "SVTYPE",
                header::TagType::String,
                header::TagLength::Fixed(1),
            ),
            (
                // INFO=<ID=SVLEN,Number=.,Type=Integer>
                "SVLEN",
                header::TagType::Integer,
                header::TagLength::Variable,
            ),
            (
                // INFO<ID=CIGAR,Number=A,Type=String>
                "CIGAR",
                header::TagType::String,
                header::TagLength::AltAlleles,
            ),
        ];
        for (ref_name, ref_type, ref_length) in truth {
            let (tag_type, tag_length) = header.info_type(ref_name.as_bytes()).unwrap();
            assert_eq!(tag_type, ref_type);
            assert_eq!(tag_length, ref_length);
        }

        assert!(header.info_type(b"NOT_THERE").is_err());
    }

    #[test]
    fn test_remove_alleles() {
        let mut bcf = Reader::from_path(&"test/test_multi.bcf").unwrap();
        for res in bcf.records() {
            let mut record = res.unwrap();
            if record.pos() == 10080 {
                record.remove_alleles(&[false, false, true]).unwrap();
                assert_eq!(record.alleles(), [b"A", b"C"]);
            }
        }
    }

    // Helper function reading full file into string.
    fn read_all<P: AsRef<Path>>(path: P) -> String {
        let mut file = File::open(path.as_ref())
            .unwrap_or_else(|_| panic!("Unable to open the file: {:?}", path.as_ref()));
        let mut contents = String::new();
        file.read_to_string(&mut contents)
            .unwrap_or_else(|_| panic!("Unable to read the file: {:?}", path.as_ref()));
        contents
    }

    // Open `test_various.vcf`, add a record from scratch to it and write it out again.
    //
    // This exercises the full functionality of updating information in a `record::Record`.
    #[test]
    fn test_write_various() {
        // Open reader, then create writer.
        let tmp = tempfile::Builder::new()
            .prefix("rust-htslib")
            .tempdir()
            .expect("Cannot create temp dir");
        let out_path = tmp.path().join("test_various.out.vcf");

        let vcf = Reader::from_path(&"test/test_various.vcf").expect("Error opening file.");
        // The writer goes into its own block so we can ensure that the file is closed and
        // all data is written below.
        {
            let mut writer = Writer::from_path(
                &out_path,
                &Header::from_template(&vcf.header()),
                true,
                Format::Vcf,
            )
            .expect("Error opening file.");

            // Setup empty record, filled below.
            let mut record = writer.empty_record();

            record.set_rid(Some(0));
            assert_eq!(record.rid().unwrap(), 0);

            record.set_pos(12);
            assert_eq!(record.pos(), 12);

            assert_eq!(str::from_utf8(record.id().as_ref()).unwrap(), ".");
            record.set_id(b"to_be_cleared").unwrap();
            assert_eq!(
                str::from_utf8(record.id().as_ref()).unwrap(),
                "to_be_cleared"
            );
            record.clear_id().unwrap();
            assert_eq!(str::from_utf8(record.id().as_ref()).unwrap(), ".");
            record.set_id(b"first_id").unwrap();
            record.push_id(b"second_id").unwrap();
            record.push_id(b"first_id").unwrap();

            assert!(record.filters().next().is_none());
            record.set_filters(&["q10".as_bytes()]).unwrap();
            record.push_filter("s50".as_bytes()).unwrap();
            record.remove_filter("q10".as_bytes(), true).unwrap();
            record.push_filter("q10".as_bytes()).unwrap();

            record.set_alleles(&[b"C", b"T", b"G"]).unwrap();

            record.set_qual(10.0);

            record.push_info_integer(b"N1", &[32]).unwrap();
            record.push_info_float(b"F1", &[33.0]).unwrap();
            record.push_info_string(b"S1", &[b"fourtytwo"]).unwrap();
            record.push_info_flag(b"X1").unwrap();

            record
                .push_genotypes(&[
                    GenotypeAllele::Unphased(0),
                    GenotypeAllele::Unphased(1),
                    GenotypeAllele::Unphased(1),
                    GenotypeAllele::Phased(1),
                ])
                .unwrap();

            record
                .push_format_string(b"FS1", &[&b"yes"[..], &b"no"[..]])
                .unwrap();
            record.push_format_integer(b"FF1", &[43, 11]).unwrap();
            record.push_format_float(b"FN1", &[42.0, 10.0]).unwrap();
            record
                .push_format_char(b"CH1", &[b"A"[0], b"B"[0]])
                .unwrap();

            // Finally, write out the record.
            writer.write(&record).unwrap();
        }

        // Now, compare expected and real output.
        let expected = read_all("test/test_various.out.vcf");
        let actual = read_all(&out_path);
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_remove_headers() {
        let vcf = Reader::from_path(&"test/test_headers.vcf").expect("Error opening file.");
        let tmp = tempfile::Builder::new()
            .prefix("rust-htslib")
            .tempdir()
            .expect("Cannot create temp dir");
        let vcfpath = tmp.path().join("test.vcf");
        let mut header = Header::from_template(&vcf.header);
        header
            .remove_contig(b"contig2")
            .remove_info(b"INFO2")
            .remove_format(b"FORMAT2")
            .remove_filter(b"FILTER2")
            .remove_structured(b"Foo2")
            .remove_generic(b"Bar2");
        {
            let mut _writer = Writer::from_path(&vcfpath, &header, true, Format::Vcf)
                .expect("Error opening output file.");
            // Note that we don't need to write anything, we are just looking at the header.
        }

        let expected = read_all("test/test_headers.out.vcf");
        let actual = read_all(&vcfpath);
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_synced_reader() {
        let mut reader = synced::SyncedReader::new().unwrap();
        reader.set_require_index(true);
        reader.set_pairing(synced::pairing::SNPS);

        assert_eq!(reader.reader_count(), 0);
        reader.add_reader(&"test/test_left.vcf.gz").unwrap();
        reader.add_reader(&"test/test_right.vcf.gz").unwrap();
        assert_eq!(reader.reader_count(), 2);

        let res1 = reader.read_next();
        assert_eq!(res1.unwrap(), 2);
        assert!(reader.has_line(0));
        assert!(reader.has_line(1));

        let res2 = reader.read_next();
        assert_eq!(res2.unwrap(), 1);
        assert!(reader.has_line(0));
        assert!(!reader.has_line(1));

        let res3 = reader.read_next();
        assert_eq!(res3.unwrap(), 1);
        assert!(!reader.has_line(0));
        assert!(reader.has_line(1));

        let res4 = reader.read_next();
        assert_eq!(res4.unwrap(), 0);
    }

    #[test]
    fn test_synced_reader_fetch() {
        let mut reader = synced::SyncedReader::new().unwrap();
        reader.set_require_index(true);
        reader.set_pairing(synced::pairing::SNPS);

        assert_eq!(reader.reader_count(), 0);
        reader.add_reader(&"test/test_left.vcf.gz").unwrap();
        reader.add_reader(&"test/test_right.vcf.gz").unwrap();
        assert_eq!(reader.reader_count(), 2);

        reader.fetch(0, 0, 1000).unwrap();
        let res1 = reader.read_next();
        assert_eq!(res1.unwrap(), 2);
        assert!(reader.has_line(0));
        assert!(reader.has_line(1));

        let res2 = reader.read_next();
        assert_eq!(res2.unwrap(), 1);
        assert!(reader.has_line(0));
        assert!(!reader.has_line(1));

        let res3 = reader.read_next();
        assert_eq!(res3.unwrap(), 1);
        assert!(!reader.has_line(0));
        assert!(reader.has_line(1));

        let res4 = reader.read_next();
        assert_eq!(res4.unwrap(), 0);
    }

    #[test]
    fn test_svlen() {
        let mut reader = Reader::from_path("test/test_svlen.vcf").unwrap();

        let mut record = reader.empty_record();
        reader.read(&mut record).unwrap().unwrap();

        assert_eq!(
            *record.info(b"SVLEN").integer().unwrap().unwrap(),
            &[-127][..]
        );
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
    fn test_multi_string_info_tag() {
        let mut reader = Reader::from_path("test/test-info-multi-string.vcf").unwrap();
        let mut rec = reader.empty_record();
        let _ = reader.read(&mut rec);

        assert_eq!(
            rec.info_shared_buffer(b"ANN", Buffer::new())
                .string()
                .unwrap()
                .unwrap()
                .len(),
            14
        );
    }

    #[test]
    fn test_multi_string_info_tag_number_a() {
        let mut reader = Reader::from_path("test/test-info-multi-string-number=A.vcf").unwrap();
        let mut rec = reader.empty_record();
        let _ = reader.read(&mut rec);

        assert_eq!(
            rec.info_shared_buffer(b"X", Buffer::new())
                .string()
                .unwrap()
                .unwrap()
                .len(),
            2
        );
    }

    #[test]
    fn test_genotype_allele_conversion() {
        let allele = GenotypeAllele::Unphased(1);
        let converted: i32 = allele.into();
        let expected = 4;
        assert_eq!(converted, expected);
        let reverse_conversion = GenotypeAllele::from(expected);
        assert_eq!(allele, reverse_conversion);
    }

    #[test]
    fn test_genotype_missing_allele_conversion() {
        let allele = GenotypeAllele::PhasedMissing;
        let converted: i32 = allele.into();
        let expected = 1;
        assert_eq!(converted, expected);
        let reverse_conversion = GenotypeAllele::from(expected);
        assert_eq!(allele, reverse_conversion);
    }

    #[test]
    fn test_alt_allele_dosage() {
        let path = &"test/test_string.vcf";
        let mut bcf = Reader::from_path(path).expect("Error opening file.");
        let _header = bcf.header();
        // FORMAT fields of first record of the vcf should look like:
        // GT:FS1:FN1	./1:LongString1:1	1/1:ss1:2
        let first_record = bcf.records().next().unwrap().expect("Fail to read record");
        let sample_count = usize::try_from(first_record.sample_count()).unwrap();
        assert_eq!(sample_count, 2);
        let mut n_ref = vec![0; sample_count];
        let mut n_alt = vec![0; sample_count];
        let mut n_missing = vec![0; sample_count];
        let gts = first_record.genotypes().expect("Error reading genotypes");
        for sample_index in 0..sample_count {
            // for each sample
            for gta in gts.get(sample_index).iter() {
                // for each allele
                match gta.index() {
                    Some(0) => n_ref[sample_index] += 1,  // reference allele
                    Some(_) => n_alt[sample_index] += 1,  // alt allele
                    None => n_missing[sample_index] += 1, // missing allele
                }
            }
        }
        assert_eq!(n_ref, [0, 0]);
        assert_eq!(n_alt, [1, 2]);
        assert_eq!(n_missing, [1, 0]);
    }

    #[test]
    fn test_obs_cornercase() {
        let mut reader = Reader::from_path("test/obs-cornercase.vcf").unwrap();
        let first_record = reader
            .records()
            .next()
            .unwrap()
            .expect("Fail to read record");

        assert_eq!(
            *first_record.info(b"EVENT").string().unwrap().unwrap(),
            [b"gridss33fb_1085"]
        );
        assert_eq!(
            *first_record.info(b"MATEID").string().unwrap().unwrap(),
            [b"gridss33fb_1085h"]
        );
    }

    // #[test]
    // fn test_buffer_lifetime() {
    //     let mut reader = Reader::from_path("test/obs-cornercase.vcf").unwrap();
    //     let first_record = reader
    //         .records()
    //         .next()
    //         .unwrap()
    //         .expect("Fail to read record");

    //     fn get_value<'a, 'b>(record: &'a Record) -> &'b [u8] {
    //         // FIXME: this should not be possible, because the slice outlives the buffer.
    //         let buffer: BufferBacked<'b, _, _> = record.info(b"EVENT").string().unwrap().unwrap();
    //         let value: &'b [u8] = buffer[0];
    //         value
    //     }

    //     let buffered = first_record.info(b"EVENT").string().unwrap().unwrap();
    //     assert_eq!(get_value(&first_record), buffered[0]);
    // }
}
