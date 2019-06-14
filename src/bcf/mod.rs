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

use std::ffi;
use std::path::Path;
use std::rc::Rc;

use url::Url;

pub mod buffer;
pub mod header;
pub mod record;

use crate::bcf::header::{HeaderView, SampleSubset};
use crate::htslib;

pub use crate::bcf::header::{Header, HeaderRecord};
pub use crate::bcf::record::Record;

/// Redefinition of corresponding `#define` in `vcf.h.`.
pub const GT_MISSING: i32 = 0;

/// A trait for a BCF reader with a read method.
pub trait Read: Sized {
    /// Read the next record.
    ///
    /// # Arguments
    /// * record - an empty record, that can be created with `bcf::Reader::empty_record`.
    fn read(&mut self, record: &mut record::Record) -> Result<(), ReadError>;

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
    fn set_threads(&mut self, n_threads: usize) -> Result<(), ThreadingError>;
}

/// A VCF/BCF reader.
#[derive(Debug)]
pub struct Reader {
    inner: *mut htslib::htsFile,
    header: Rc<HeaderView>,
}

unsafe impl Send for Reader {}

/// Implementation for `Reader::set_threads()` and `Writer::set_threads`.
pub fn set_threads(hts_file: *mut htslib::htsFile, n_threads: usize) -> Result<(), ThreadingError> {
    assert!(n_threads > 0, "n_threads must be > 0");

    let r = unsafe { htslib::hts_set_threads(hts_file, n_threads as i32) };
    if r != 0 {
        Err(ThreadingError::Some)
    } else {
        Ok(())
    }
}

impl Reader {
    /// Create a new reader from a given path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, BCFPathError> {
        match path.as_ref().to_str() {
            Some(p) if path.as_ref().exists() => Ok(r#try!(Self::new(p.as_bytes()))),
            _ => Err(BCFPathError::InvalidPath),
        }
    }

    /// Create a new reader from a given URL.
    pub fn from_url(url: &Url) -> Result<Self, BCFError> {
        Self::new(url.as_str().as_bytes())
    }

    /// Create a new reader from standard input.
    pub fn from_stdin() -> Result<Self, BCFError> {
        Self::new(b"-")
    }

    fn new(path: &[u8]) -> Result<Self, BCFError> {
        let htsfile = r#try!(bcf_open(path, b"r"));
        let header = unsafe { htslib::bcf_hdr_read(htsfile) };
        Ok(Reader {
            inner: htsfile,
            header: Rc::new(HeaderView::new(header)),
        })
    }
}

impl Read for Reader {
    fn read(&mut self, record: &mut record::Record) -> Result<(), ReadError> {
        match unsafe { htslib::bcf_read(self.inner, self.header.inner, record.inner) } {
            0 => {
                unsafe {
                    // Always unpack record.
                    htslib::bcf_unpack(record.inner_mut(), htslib::BCF_UN_ALL as i32);
                }
                record.set_header(self.header.clone());
                Ok(())
            }
            -1 => Err(ReadError::NoMoreRecord),
            _ => Err(ReadError::Invalid),
        }
    }

    fn records(&mut self) -> Records<'_, Self> {
        Records { reader: self }
    }

    fn set_threads(&mut self, n_threads: usize) -> Result<(), ThreadingError> {
        set_threads(self.inner, n_threads)
    }

    fn header(&self) -> &HeaderView {
        return &self.header;
    }

    /// Return empty record.  Can be reused multiple times.
    fn empty_record(&self) -> Record {
        return Record::new(self.header.clone());
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
    current_region: Option<(u32, u32, u32)>,
}

unsafe impl Send for IndexedReader {}

impl IndexedReader {
    /// Create a new `IndexedReader` from path.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, IndexedReaderPathError> {
        match path.as_ref().to_str() {
            Some(p) if path.as_ref().exists() => {
                Ok(r#try!(Self::new(&ffi::CString::new(p).unwrap())))
            }
            _ => Err(IndexedReaderPathError::InvalidPath),
        }
    }

    /// Create a new `IndexedReader` from an URL.
    pub fn from_url(url: &Url) -> Result<Self, IndexedReaderError> {
        Self::new(&ffi::CString::new(url.as_str()).unwrap())
    }

    /// Create a new `IndexedReader`.
    ///
    /// # Arguments
    ///
    /// * `path` - the path. Use "-" for stdin.
    fn new(path: &ffi::CStr) -> Result<Self, IndexedReaderError> {
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
                header: header,
                current_region: None,
            })
        } else {
            Err(IndexedReaderError::InvalidPath)
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
    pub fn fetch(&mut self, rid: u32, start: u32, end: u32) -> Result<(), FetchError> {
        let contig = self.header.rid2name(rid).unwrap();
        let contig = ffi::CString::new(contig).unwrap();
        let contig = contig.as_ptr();
        if unsafe { htslib::bcf_sr_seek(self.inner, contig, start as i32) } != 0 {
            Err(FetchError::Some)
        } else {
            self.current_region = Some((rid, start, end));
            Ok(())
        }
    }
}

impl Read for IndexedReader {
    fn read(&mut self, record: &mut record::Record) -> Result<(), ReadError> {
        match unsafe { htslib::bcf_sr_next_line(self.inner) } {
            0 => {
                if unsafe { (*self.inner).errnum } != 0 {
                    Err(ReadError::SyncedBcfReaderError)
                } else {
                    Err(ReadError::NoMoreRecord)
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

                record.set_header(self.header.clone());

                match self.current_region {
                    Some((rid, _start, end)) => {
                        if record.rid().is_some()
                            && rid == record.rid().unwrap()
                            && record.pos() <= end
                        {
                            Ok(())
                        } else {
                            Err(ReadError::NoMoreRecord)
                        }
                    }
                    None => Ok(()),
                }
            }
        }
    }

    fn records(&mut self) -> Records<'_, Self> {
        Records { reader: self }
    }

    fn set_threads(&mut self, n_threads: usize) -> Result<(), ThreadingError> {
        assert!(n_threads > 0, "n_threads must be > 0");

        let r = unsafe { htslib::bcf_sr_set_threads(self.inner, n_threads as i32) };
        if r != 0 {
            Err(ThreadingError::Some)
        } else {
            Ok(())
        }
    }

    fn header(&self) -> &HeaderView {
        return &self.header;
    }

    fn empty_record(&self) -> Record {
        return Record::new(self.header.clone());
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
        current_region: Option<(u32, u32, u32)>,
    }

    // TODO: add interface for setting threads, ensure that the pool is freed properly

    impl SyncedReader {
        pub fn new() -> Result<Self, AllocationError> {
            let inner = unsafe { crate::htslib::bcf_sr_init() };
            if inner.is_null() {
                Err(AllocationError::Some)
            } else {
                Ok(SyncedReader {
                    inner: inner,
                    headers: Vec::new(),
                    current_region: None,
                })
            }
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
        pub fn add_reader<P: AsRef<Path>>(&mut self, path: P) -> Result<(), BCFError> {
            match path.as_ref().to_str() {
                Some(p) if path.as_ref().exists() => {
                    let p_cstring = ffi::CString::new(p).unwrap();
                    let res =
                        unsafe { crate::htslib::bcf_sr_add_reader(self.inner, p_cstring.as_ptr()) };
                    if res != 0 {
                        let i = (self.reader_count() - 1) as isize;
                        let header = Rc::new(HeaderView::new(unsafe {
                            crate::htslib::bcf_hdr_dup((*(*self.inner).readers.offset(i)).header)
                        }));
                        self.headers.push(header);
                        Ok(())
                    } else {
                        Err(BCFError::Some)
                    }
                }
                _ => Err(BCFError::Some),
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

        /// Read next line and return number of readers that have the given line.
        pub fn read_next(&mut self) -> Result<u32, ReadError> {
            let num = unsafe { crate::htslib::bcf_sr_next_line(self.inner) as u32 };

            if num == 0 {
                if unsafe { (*self.inner).errnum } != 0 {
                    Err(ReadError::Invalid)
                } else {
                    Ok(0)
                }
            } else {
                assert!(num > 0, "Must not be negative");
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
                                if (*record).rid != (rid as i32) || (*record).pos >= (end as i32) {
                                    return Err(ReadError::NoMoreRecord);
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
        pub fn fetch(&mut self, rid: u32, start: u32, end: u32) -> Result<(), FetchError> {
            let contig = {
                let contig = self.header(0).rid2name(rid).unwrap().clone();
                ffi::CString::new(contig).unwrap()
            };
            let contig = contig.as_ptr();
            if unsafe { htslib::bcf_sr_seek(self.inner, contig, start as i32) } != 0 {
                Err(FetchError::Some)
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
        vcf: bool,
    ) -> Result<Self, BCFPathError> {
        if let Some(p) = path.as_ref().to_str() {
            Ok(r#try!(Self::new(p.as_bytes(), header, uncompressed, vcf)))
        } else {
            Err(BCFPathError::InvalidPath)
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
        vcf: bool,
    ) -> Result<Self, BCFError> {
        Self::new(url.as_str().as_bytes(), header, uncompressed, vcf)
    }

    /// Create a new writer to stdout.
    ///
    /// # Arguments
    ///
    /// * `header` - header definition to use
    /// * `uncompressed` - disable compression
    /// * `vcf` - write VCF instead of BCF
    pub fn from_stdout(header: &Header, uncompressed: bool, vcf: bool) -> Result<Self, BCFError> {
        Self::new(b"-", header, uncompressed, vcf)
    }

    fn new(path: &[u8], header: &Header, uncompressed: bool, vcf: bool) -> Result<Self, BCFError> {
        let mode: &[u8] = match (uncompressed, vcf) {
            (true, true) => b"w",
            (false, true) => b"wz",
            (true, false) => b"wu",
            (false, false) => b"wb",
        };

        let htsfile = r#try!(bcf_open(path, mode));
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
        record::Record::new(self.header.clone())
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
        record.set_header(self.header.clone());
    }

    /// Subset samples of record to match header of this writer.
    ///
    /// # Arguments
    ///
    /// - `record` - The `Record` to modify.
    pub fn subset(&mut self, record: &mut record::Record) {
        match self.subset {
            Some(ref mut subset) => unsafe {
                htslib::bcf_subset(
                    self.header.inner,
                    record.inner,
                    subset.len() as i32,
                    subset.as_mut_ptr(),
                );
            },
            None => (),
        }
    }

    /// Write `record` to the Writer.
    ///
    /// # Arguments
    ///
    /// - `record` - The `Record` to write.
    pub fn write(&mut self, record: &record::Record) -> Result<(), WriteError> {
        if unsafe { htslib::bcf_write(self.inner, self.header.inner, record.inner) } == -1 {
            Err(WriteError::Some)
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
    pub fn set_threads(&mut self, n_threads: usize) -> Result<(), ThreadingError> {
        set_threads(self.inner, n_threads)
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
    type Item = Result<record::Record, ReadError>;

    fn next(&mut self) -> Option<Result<record::Record, ReadError>> {
        let mut record = self.reader.empty_record();
        match self.reader.read(&mut record) {
            Err(ReadError::NoMoreRecord) => None,
            Err(e) => Some(Err(e)),
            Ok(()) => Some(Ok(record)),
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum AllocationError {
        Some {
            description("Error allocating internal data structure")
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum BCFError {
        Some {
            description("error reading BCF/VCF file")
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum BCFPathError {
        InvalidPath {
            description("invalid path")
        }
        BCFError(err: BCFError) {
            from()
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum ThreadingError {
        Some {
            description("error setting threads for multi-threaded I/O")
        }
    }
}

/// Wrapper for opening a BCF file.
fn bcf_open(path: &[u8], mode: &[u8]) -> Result<*mut htslib::htsFile, BCFError> {
    let p = ffi::CString::new(path).unwrap();
    let ret = unsafe { htslib::hts_open(p.as_ptr(), ffi::CString::new(mode).unwrap().as_ptr()) };
    if ret.is_null() {
        Err(BCFError::Some)
    } else {
        unsafe {
            if !mode.contains(&b'w')
                && (*ret).format.category != htslib::htsFormatCategory_variant_data
            {
                Err(BCFError::Some)
            } else {
                Ok(ret)
            }
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum ReadError {
        Invalid {
            description("invalid record")
        }
        NoMoreRecord {
            description("no more record")
        }
        SyncedBcfReaderError {
            description("problem reading from synced bcf reader")
        }
    }
}

impl ReadError {
    /// Returns true if no record has been read because the end of the file was reached.
    pub fn is_eof(&self) -> bool {
        match self {
            &ReadError::NoMoreRecord => true,
            _ => false,
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum IndexedReaderPathError {
        InvalidPath {
            description("invalid path")
        }
        IndexedReaderError(err: IndexedReaderError) {
            from()
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum IndexedReaderError {
        InvalidPath {
            description("invalid path")
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum WriteError {
        Some {
            description("failed to write record")
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum FetchError {
        Some {
            description("error fetching a locus")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bcf::header::Id;
    use crate::bcf::record::Numeric;
    use std::fs::File;
    use std::io::prelude::Read as IoRead;
    use std::path::Path;
    use std::str;
    use tempdir;

    fn _test_read<P: AsRef<Path>>(path: &P) {
        let mut bcf = Reader::from_path(path).ok().expect("Error opening file.");
        assert_eq!(bcf.header.samples(), [b"NA12878.subsample-0.25-0"]);

        for (i, rec) in bcf.records().enumerate() {
            let mut record = rec.ok().expect("Error reading record.");
            assert_eq!(record.sample_count(), 1);

            assert_eq!(record.rid().expect("Error reading rid."), 0);
            assert_eq!(record.pos(), 10021 + i as u32);
            assert_eq!(record.qual(), 0f32);
            assert_eq!(
                record
                    .info(b"MQ0F")
                    .float()
                    .ok()
                    .expect("Error reading info.")
                    .expect("Missing tag"),
                [1.0]
            );
            if i == 59 {
                assert_eq!(
                    record
                        .info(b"SGB")
                        .float()
                        .ok()
                        .expect("Error reading info.")
                        .expect("Missing tag"),
                    [-0.379885]
                );
            }
            // the artificial "not observed" allele is present in each record.
            assert_eq!(record.alleles().iter().last().unwrap(), b"<X>");

            let mut fmt = record.format(b"PL");
            let pl = fmt.integer().ok().expect("Error reading format.");
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
        let mut bcf = Reader::from_path(path).ok().expect("Error opening file.");
        bcf.set_threads(2).unwrap();
    }

    #[test]
    fn test_writer_set_threads() {
        let path = &"test/test.bcf";
        let tmp = tempdir::TempDir::new("rust-htslib")
            .ok()
            .expect("Cannot create temp dir");
        let bcfpath = tmp.path().join("test.bcf");
        let bcf = Reader::from_path(path).ok().expect("Error opening file.");
        let header = Header::from_template_subset(&bcf.header, &[b"NA12878.subsample-0.25-0"])
            .ok()
            .expect("Error subsetting samples.");
        let mut writer = Writer::from_path(&bcfpath, &header, false, false)
            .ok()
            .expect("Error opening file.");
        writer.set_threads(2).unwrap();
    }

    #[test]
    fn test_fetch() {
        let mut bcf = IndexedReader::from_path(&"test/test.bcf")
            .ok()
            .expect("Error opening file.");
        bcf.set_threads(2).unwrap();
        let rid = bcf
            .header()
            .name2rid(b"1")
            .expect("Translating from contig '1' to ID failed.");
        bcf.fetch(rid, 10_033, 10_060).expect("Fetching failed");
        assert_eq!(bcf.records().count(), 28);
    }

    #[test]
    fn test_write() {
        let mut bcf = Reader::from_path(&"test/test_multi.bcf")
            .ok()
            .expect("Error opening file.");
        let tmp = tempdir::TempDir::new("rust-htslib")
            .ok()
            .expect("Cannot create temp dir");
        let bcfpath = tmp.path().join("test.bcf");
        println!("{:?}", bcfpath);
        {
            let header = Header::from_template_subset(&bcf.header, &[b"NA12878.subsample-0.25-0"])
                .ok()
                .expect("Error subsetting samples.");
            let mut writer = Writer::from_path(&bcfpath, &header, false, false)
                .ok()
                .expect("Error opening file.");
            for rec in bcf.records() {
                let mut record = rec.ok().expect("Error reading record.");
                writer.translate(&mut record);
                writer.subset(&mut record);
                record.trim_alleles().ok().expect("Error trimming alleles.");
                writer.write(&record).ok().expect("Error writing record");
            }
        }
        {
            _test_read(&bcfpath);
        }
        tmp.close().ok().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_strings() {
        let mut vcf = Reader::from_path(&"test/test_string.vcf")
            .ok()
            .expect("Error opening file.");
        let fs1 = [
            &b"LongString1"[..],
            &b"LongString2"[..],
            &b"."[..],
            &b"LongString4"[..],
            &b"evenlength"[..],
            &b"ss6"[..],
        ];
        for (i, rec) in vcf.records().enumerate() {
            println!("record {}", i);
            let mut record = rec.ok().expect("Error reading record.");
            assert_eq!(
                record
                    .info(b"S1")
                    .string()
                    .ok()
                    .expect("Error reading string.")
                    .expect("Missing tag")[0],
                format!("string{}", i + 1).as_bytes()
            );
            println!(
                "{}",
                String::from_utf8_lossy(
                    record
                        .format(b"FS1")
                        .string()
                        .ok()
                        .expect("Error reading string.")[0]
                )
            );
            assert_eq!(
                record
                    .format(b"FS1")
                    .string()
                    .ok()
                    .expect("Error reading string.")[0],
                fs1[i]
            );
        }
    }

    #[test]
    fn test_missing() {
        let mut vcf = Reader::from_path(&"test/test_missing.vcf")
            .ok()
            .expect("Error opening file.");
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
        for (i, rec) in vcf.records().enumerate() {
            let mut record = rec.ok().expect("Error reading record.");
            assert_eq!(
                record
                    .info(b"F1")
                    .float()
                    .ok()
                    .expect("Error reading float.")
                    .expect("Missing tag")[0]
                    .is_nan(),
                f1[i]
            );
            assert_eq!(
                record
                    .format(b"FN4")
                    .integer()
                    .ok()
                    .expect("Error reading integer.")[1],
                fn4[i]
            );
            assert!(record
                .format(b"FF4")
                .float()
                .ok()
                .expect("Error reading float.")[1]
                .iter()
                .all(|&v| v.is_missing()));
        }
    }

    #[test]
    fn test_genotypes() {
        let mut vcf = Reader::from_path(&"test/test_string.vcf")
            .ok()
            .expect("Error opening file.");
        let expected = ["./1", "1|1", "0/1", "0|1", "1|.", "1/1"];
        for (rec, exp_gt) in vcf.records().zip(expected.into_iter()) {
            let mut rec = rec.ok().expect("Error reading record.");
            let genotypes = rec.genotypes().expect("Error reading genotypes");
            assert_eq!(&format!("{}", genotypes.get(0)), exp_gt);
        }
    }

    #[test]
    fn test_header_ids() {
        let vcf = Reader::from_path(&"test/test_string.vcf")
            .ok()
            .expect("Error opening file.");
        let header = &vcf.header();
        use crate::bcf::header::Id;

        assert_eq!(header.id_to_name(Id(4)), b"GT");
        assert_eq!(header.name_to_id(b"GT").unwrap(), Id(4));
        assert!(header.name_to_id(b"XX").is_err());
    }

    #[test]
    fn test_header_samples() {
        let vcf = Reader::from_path(&"test/test_string.vcf")
            .ok()
            .expect("Error opening file.");
        let header = &vcf.header();

        assert_eq!(header.id_to_sample(Id(0)), b"one");
        assert_eq!(header.id_to_sample(Id(1)), b"two");
        assert_eq!(header.sample_to_id(b"one").unwrap(), Id(0));
        assert_eq!(header.sample_to_id(b"two").unwrap(), Id(1));
        assert!(header.sample_to_id(b"three").is_err());
    }

    #[test]
    fn test_header_contigs() {
        let vcf = Reader::from_path(&"test/test_multi.bcf")
            .ok()
            .expect("Error opening file.");
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
        let vcf = Reader::from_path(&"test/test_string.vcf")
            .ok()
            .expect("Error opening file.");
        let records = vcf.header().header_records();
        assert_eq!(records.len(), 10);

        match &records[1] {
            &HeaderRecord::Filter {
                ref key,
                ref values,
            } => {
                assert_eq!(key, "FILTER");
                assert_eq!(values["ID"], "PASS");
            }
            _ => {
                assert!(false);
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
            .expect(&format!("Unable to open the file: {:?}", path.as_ref()));
        let mut contents = String::new();
        file.read_to_string(&mut contents)
            .expect(&format!("Unable to read the file: {:?}", path.as_ref()));
        contents
    }

    // Open `test_various.vcf`, add a record from scratch to it and write it out again.
    //
    // This exercises the full functionality of updating information in a `record::Record`.
    #[test]
    fn test_write_various() {
        // Open reader, then create writer.
        let tmp = tempdir::TempDir::new("rust-htslib")
            .ok()
            .expect("Cannot create temp dir");
        let out_path = tmp.path().join("test_various.out.vcf");

        let vcf = Reader::from_path(&"test/test_various.vcf")
            .ok()
            .expect("Error opening file.");
        // The writer goes into its own block so we can ensure that the file is closed and
        // all data is written below.
        {
            let mut writer =
                Writer::from_path(&out_path, &Header::from_template(&vcf.header()), true, true)
                    .ok()
                    .expect("Error opening file.");
            let header = writer.header().clone();

            // Setup empty record, filled below.
            let mut record = writer.empty_record();

            record.set_rid(&Some(0));
            assert_eq!(record.rid().unwrap(), 0);

            record.set_pos(12);
            assert_eq!(record.pos(), 12);

            assert_eq!(str::from_utf8(record.id().as_ref()).ok().unwrap(), ".");
            record.set_id("to_be_cleared".as_bytes()).unwrap();
            assert_eq!(
                str::from_utf8(record.id().as_ref()).ok().unwrap(),
                "to_be_cleared"
            );
            record.clear_id().unwrap();
            assert_eq!(str::from_utf8(record.id().as_ref()).ok().unwrap(), ".");
            record.set_id("first_id".as_bytes()).unwrap();
            record.push_id("second_id".as_bytes()).unwrap();
            record.push_id("first_id".as_bytes()).unwrap();

            assert!(record.filters().next().is_none());
            record.set_filters(&[header.name_to_id(b"q10").unwrap()]);
            record.push_filter(header.name_to_id(b"s50").unwrap());
            record.remove_filter(header.name_to_id(b"q10").unwrap(), true);
            record.push_filter(header.name_to_id(b"q10").unwrap());

            record
                .set_alleles(&["C".as_bytes(), "T".as_bytes(), "G".as_bytes()])
                .unwrap();

            record.set_qual(10.0);

            record.push_info_integer(b"N1", &[32]).unwrap();
            record.push_info_float(b"F1", &[33.0]).unwrap();
            record
                .push_info_string(b"S1", &["fourtytwo".as_bytes()])
                .unwrap();
            record.push_info_flag(b"X1").unwrap();

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
        let vcf = Reader::from_path(&"test/test_headers.vcf")
            .ok()
            .expect("Error opening file.");
        let tmp = tempdir::TempDir::new("rust-htslib")
            .ok()
            .expect("Cannot create temp dir");
        let bcfpath = tmp.path().join("test.vcf");
        let mut header = Header::from_template(&vcf.header);
        header
            .remove_contig(b"contig2")
            .remove_info(b"INFO2")
            .remove_format(b"FORMAT2")
            .remove_filter(b"FILTER2")
            .remove_structured(b"Foo2")
            .remove_generic(b"Bar2");
        {
            let mut _writer = Writer::from_path(&bcfpath, &header, true, true)
                .ok()
                .expect("Error opening output file.");
            // Note that we don't need to write anything, we are just looking at the header.
        }

        let expected = read_all("test/test_headers.out.vcf");
        let actual = read_all(&bcfpath);
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
        reader.read(&mut record).unwrap();

        assert_eq!(record.info(b"SVLEN").integer().unwrap(), Some(&[-127][..]));
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

}
