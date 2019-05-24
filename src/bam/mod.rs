// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Module for working with BAM and CRAM files.

pub mod buffer;
pub mod ext;
pub mod header;
pub mod index;
pub mod pileup;
pub mod record;

#[cfg(feature = "serde")]
pub mod record_serde;

use libc;
use std::ffi;
use std::path::Path;
use std::slice;
use url::Url;

use htslib;

pub use bam::buffer::RecordBuffer;
pub use bam::header::Header;
pub use bam::record::Record;

/// Implementation for `Read::set_threads` and `Writer::set_threads`.
pub fn set_threads(htsfile: *mut htslib::htsFile, n_threads: usize) -> Result<(), ThreadingError> {
    assert!(n_threads != 0, "n_threads must be > 0");

    if unsafe { htslib::hts_set_threads(htsfile, n_threads as i32) } != 0 {
        Err(ThreadingError::Some)
    } else {
        Ok(())
    }
}

/// Set the reference FAI index path in a `htslib::htsFile` struct for reading CRAM format.
pub fn set_fai_filename(
    htsfile: *mut htslib::htsFile,
    path: &[u8],
) -> Result<(), ReferencePathError> {
    if unsafe { htslib::hts_set_fai_filename(htsfile, ffi::CString::new(path).unwrap().as_ptr()) }
        == 0
    {
        Ok(())
    } else {
        Err(ReferencePathError::InvalidPath)
    }
}

/// A trait for a BAM reader with a read method.
pub trait Read: Sized {
    /// Read next BAM record into given record.
    /// Use this method in combination with a single allocated record to avoid the reallocations
    /// occurring with the iterator.
    ///
    /// # Arguments
    ///
    /// * `record` - the record to be filled
    fn read(&mut self, record: &mut record::Record) -> Result<(), ReadError>;

    /// Iterator over the records of the seeked region.
    /// Note that, while being convenient, this is less efficient than pre-allocating a
    /// `Record` and reading into it with the `read` method, since every iteration involves
    /// the allocation of a new `Record`.
    fn records(&mut self) -> Records<Self>;

    /// Iterator over pileups.
    fn pileup(&mut self) -> pileup::Pileups<Self>;

    /// Return the htsFile struct
    fn htsfile(&self) -> *mut htslib::htsFile;

    /// Return the header.
    fn header(&self) -> &HeaderView;

    /// Seek to the given virtual offset in the file
    fn seek(&mut self, offset: i64) -> Result<(), SeekError> {
        let htsfile = unsafe { self.htsfile().as_ref() }.expect("bug: null pointer to htsFile");
        let ret = match htsfile.format.format {
            htslib::htsExactFormat_cram => unsafe {
                htslib::cram_seek(htsfile.fp.cram, offset as libc::off_t, libc::SEEK_SET) as i64
            },
            _ => unsafe { htslib::bgzf_seek(htsfile.fp.bgzf, offset, libc::SEEK_SET) },
        };

        if ret == 0 {
            Ok(())
        } else {
            Err(SeekError::Some)
        }
    }

    /// Report the current virtual offset
    fn tell(&self) -> i64 {
        // this reimplements the bgzf_tell macro
        let htsfile = unsafe { self.htsfile().as_ref() }.expect("bug: null pointer to htsFile");
        let bgzf = unsafe { *htsfile.fp.bgzf };
        (bgzf.block_address << 16) | (bgzf.block_offset as i64 & 0xFFFF)
    }

    /// Activate multi-threaded BAM read support in htslib. This should permit faster
    /// reading of large BAM files.
    ///
    /// Setting `nthreads` to `0` does not change the current state.  Note that it is not
    /// possible to set the number of background threads below `1` once it has been set.
    ///
    /// # Arguments
    ///
    /// * `n_threads` - number of extra background writer threads to use, must be `> 0`.
    fn set_threads(&mut self, n_threads: usize) -> Result<(), ThreadingError> {
        set_threads(self.htsfile(), n_threads)
    }
}

/// A BAM reader.
#[derive(Debug)]
pub struct Reader {
    htsfile: *mut htslib::htsFile,
    header: HeaderView,
}

unsafe impl Send for Reader {}

impl Reader {
    /// Create a new Reader from path.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, ReaderPathError> {
        match path.as_ref().to_str() {
            Some(p) if path.as_ref().exists() => Ok(try!(Self::new(p.as_bytes()))),
            _ => Err(ReaderPathError::InvalidPath),
        }
    }

    /// Create a new Reader from STDIN.
    pub fn from_stdin() -> Result<Self, BGZFError> {
        Self::new(b"-")
    }

    /// Create a new Reader from URL.
    pub fn from_url(url: &Url) -> Result<Self, BGZFError> {
        Self::new(url.as_str().as_bytes())
    }

    /// Create a new Reader.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open. Use "-" for stdin.
    fn new(path: &[u8]) -> Result<Self, BGZFError> {
        let htsfile = try!(hts_open(&ffi::CString::new(path).unwrap(), b"r"));

        let header = unsafe { htslib::sam_hdr_read(htsfile) };
        Ok(Reader {
            htsfile: htsfile,
            header: HeaderView::new(header),
        })
    }

    extern "C" fn pileup_read(
        data: *mut ::std::os::raw::c_void,
        record: *mut htslib::bam1_t,
    ) -> i32 {
        let _self = unsafe { &*(data as *mut Self) };
        unsafe { htslib::sam_read1(_self.htsfile(), &mut _self.header().inner(), record) }
    }

    /// Iterator over the records between the (optional) virtual offsets `start` and `end`
    ///
    /// # Arguments
    ///
    /// * `start` - Optional starting virtual offset to seek to. Throws an error if it is not
    /// a valid virtual offset.
    ///
    /// * `end` - Read until the virtual offset is less than `end`
    pub fn iter_chunk(&mut self, start: Option<i64>, end: Option<i64>) -> ChunkIterator<Self> {
        if let Some(pos) = start {
            self.seek(pos)
                .expect("Failed to seek to the starting position");
        };

        ChunkIterator {
            reader: self,
            end: end,
        }
    }

    /// Set the reference path for reading CRAM files.
    ///
    /// # Arguments
    ///
    /// * `path` - path to the FASTA reference
    pub fn set_reference<P: AsRef<Path>>(&mut self, path: P) -> Result<(), ReferencePathError> {
        let fai_path_str = format!("{}.fai", path.as_ref().to_str().unwrap());
        let fai_path = Path::new(&fai_path_str);

        match fai_path.to_str() {
            Some(p) if fai_path.exists() => set_fai_filename(self.htsfile, p.as_bytes()),
            _ => Err(ReferencePathError::InvalidPath),
        }
    }
}

impl Read for Reader {
    fn read(&mut self, record: &mut record::Record) -> Result<(), ReadError> {
        match unsafe { htslib::sam_read1(self.htsfile, &mut self.header.inner(), record.inner) } {
            -1 => Err(ReadError::NoMoreRecord),
            -2 => Err(ReadError::Truncated),
            -4 => Err(ReadError::Invalid),
            _ => Ok(()),
        }
    }

    /// Iterator over the records of the fetched region.
    /// Note that, while being convenient, this is less efficient than pre-allocating a
    /// `Record` and reading into it with the `read` method, since every iteration involves
    /// the allocation of a new `Record`.
    fn records(&mut self) -> Records<Self> {
        Records { reader: self }
    }

    fn pileup(&mut self) -> pileup::Pileups<Self> {
        let _self = self as *const Self;
        let itr = unsafe {
            htslib::bam_plp_init(
                Some(Reader::pileup_read),
                _self as *mut ::std::os::raw::c_void,
            )
        };
        pileup::Pileups::new(self, itr)
    }

    fn htsfile(&self) -> *mut htslib::htsFile {
        self.htsfile
    }

    fn header(&self) -> &HeaderView {
        &self.header
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        unsafe {
            htslib::hts_close(self.htsfile);
        }
    }
}

#[derive(Debug)]
pub struct IndexedReader {
    htsfile: *mut htslib::htsFile,
    header: HeaderView,
    idx: *mut htslib::hts_idx_t,
    itr: Option<*mut htslib::hts_itr_t>,
}

unsafe impl Send for IndexedReader {}

impl IndexedReader {
    /// Create a new Reader from path.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, IndexedReaderPathError> {
        match path.as_ref().to_str() {
            Some(p) if path.as_ref().exists() => {
                Ok(try!(Self::new(&ffi::CString::new(p).unwrap())))
            }
            _ => Err(IndexedReaderPathError::InvalidPath),
        }
    }
    pub fn from_path_and_index<P: AsRef<Path>>(
        path: P,
        index_path: P,
    ) -> Result<Self, IndexedReaderPathError> {
        if !path.as_ref().exists() || !index_path.as_ref().exists() {
            return Err(IndexedReaderPathError::InvalidPath);
        }
        let p = path
            .as_ref()
            .to_str()
            .ok_or(IndexedReaderPathError::InvalidPath)?;
        let idx_p = index_path
            .as_ref()
            .to_str()
            .ok_or(IndexedReaderPathError::InvalidPath)?;
        Ok(try!(Self::new_with_index_path(
            &ffi::CString::new(p).unwrap(),
            &ffi::CString::new(idx_p).unwrap(),
        )))
    }

    pub fn from_url(url: &Url) -> Result<Self, IndexedReaderError> {
        Self::new(&ffi::CString::new(url.as_str()).unwrap())
    }

    /// Create a new Reader.
    ///
    /// # Arguments
    ///
    /// * `path` - the path. Use "-" for stdin.
    fn new(path: &ffi::CStr) -> Result<Self, IndexedReaderError> {
        let htsfile = try!(hts_open(path, b"r"));
        let header = unsafe { htslib::sam_hdr_read(htsfile) };
        let idx = unsafe { htslib::sam_index_load(htsfile, path.as_ptr()) };
        if idx.is_null() {
            Err(IndexedReaderError::InvalidIndex)
        } else {
            Ok(IndexedReader {
                htsfile: htsfile,
                header: HeaderView::new(header),
                idx: idx,
                itr: None,
            })
        }
    }
    /// Create a new Reader.
    ///
    /// # Arguments
    ///
    /// * `path` - the path. Use "-" for stdin.
    /// * `index_path` - the index path to use
    fn new_with_index_path(
        path: &ffi::CStr,
        index_path: &ffi::CStr,
    ) -> Result<Self, IndexedReaderError> {
        let htsfile = try!(hts_open(path, b"r"));
        let header = unsafe { htslib::sam_hdr_read(htsfile) };
        let idx = unsafe { htslib::sam_index_load2(htsfile, path.as_ptr(), index_path.as_ptr()) };
        if idx.is_null() {
            Err(IndexedReaderError::InvalidIndex)
        } else {
            Ok(IndexedReader {
                htsfile: htsfile,
                header: HeaderView::new(header),
                idx: idx,
                itr: None,
            })
        }
    }

    pub fn fetch(&mut self, tid: u32, beg: u32, end: u32) -> Result<(), FetchError> {
        if let Some(itr) = self.itr {
            unsafe { htslib::hts_itr_destroy(itr) }
        }
        let itr = unsafe { htslib::sam_itr_queryi(self.idx, tid as i32, beg as i32, end as i32) };
        if itr.is_null() {
            self.itr = None;
            Err(FetchError::Some)
        } else {
            self.itr = Some(itr);
            Ok(())
        }
    }

    extern "C" fn pileup_read(
        data: *mut ::std::os::raw::c_void,
        record: *mut htslib::bam1_t,
    ) -> i32 {
        let _self = unsafe { &*(data as *mut Self) };
        match _self.itr {
            Some(itr) => itr_next(_self.htsfile, itr, record), // read fetched region
            None => unsafe { htslib::sam_read1(_self.htsfile, &mut _self.header.inner(), record) }, // ordinary reading
        }
    }

    /// Set the reference path for reading CRAM files.
    ///
    /// # Arguments
    ///
    /// * `path` - path to the FASTA reference
    pub fn set_reference<P: AsRef<Path>>(&mut self, path: P) -> Result<(), ReferencePathError> {
        let fai_path_str = format!("{}.fai", path.as_ref().to_str().unwrap());
        let fai_path = Path::new(&fai_path_str);

        match fai_path.to_str() {
            Some(p) if fai_path.exists() => set_fai_filename(self.htsfile, p.as_bytes()),
            _ => Err(ReferencePathError::InvalidPath),
        }
    }
}

impl Read for IndexedReader {
    fn read(&mut self, record: &mut record::Record) -> Result<(), ReadError> {
        match self.itr {
            Some(itr) => match itr_next(self.htsfile, itr, record.inner) {
                -1 => Err(ReadError::NoMoreRecord),
                -2 => Err(ReadError::Truncated),
                -4 => Err(ReadError::Invalid),
                _ => Ok(()),
            },
            None => Err(ReadError::NoMoreRecord),
        }
    }

    /// Iterator over the records of the fetched region.
    /// Note that, while being convenient, this is less efficient than pre-allocating a
    /// `Record` and reading into it with the `read` method, since every iteration involves
    /// the allocation of a new `Record`.
    fn records(&mut self) -> Records<Self> {
        Records { reader: self }
    }

    fn pileup(&mut self) -> pileup::Pileups<Self> {
        let _self = self as *const Self;
        let itr = unsafe {
            htslib::bam_plp_init(
                Some(IndexedReader::pileup_read),
                _self as *mut ::std::os::raw::c_void,
            )
        };
        pileup::Pileups::new(self, itr)
    }

    fn htsfile(&self) -> *mut htslib::htsFile {
        self.htsfile
    }

    fn header(&self) -> &HeaderView {
        &self.header
    }
}

impl Drop for IndexedReader {
    fn drop(&mut self) {
        unsafe {
            if self.itr.is_some() {
                htslib::hts_itr_destroy(self.itr.unwrap());
            }
            htslib::hts_idx_destroy(self.idx);
            htslib::hts_close(self.htsfile);
        }
    }
}

/// A BAM writer.
#[derive(Debug)]
pub struct Writer {
    f: *mut htslib::htsFile,
    header: HeaderView,
}

unsafe impl Send for Writer {}

impl Writer {
    /// Create a new BAM file.
    ///
    /// # Arguments
    ///
    /// * `path` - the path.
    /// * `header` - header definition to use
    pub fn from_path<P: AsRef<Path>>(
        path: P,
        header: &header::Header,
    ) -> Result<Self, WriterPathError> {
        if let Some(p) = path.as_ref().to_str() {
            Ok(try!(Self::new(p.as_bytes(), b"wb", header)))
        } else {
            Err(WriterPathError::InvalidPath)
        }
    }

    /// Create a new CRAM file.
    ///
    /// # Arguments
    ///
    /// * `path` - the path.
    /// * `header` - header definition to use
    pub fn from_cram_path<P: AsRef<Path>>(
        path: P,
        header: &header::Header,
    ) -> Result<Self, WriterPathError> {
        if let Some(p) = path.as_ref().to_str() {
            Ok(try!(Self::new(p.as_bytes(), b"wc", header)))
        } else {
            Err(WriterPathError::InvalidPath)
        }
    }

    /// Create a new BAM file at STDOUT.
    ///
    /// # Arguments
    ///
    /// * `header` - header definition to use
    pub fn from_stdout(header: &header::Header) -> Result<Self, BGZFError> {
        Self::new(b"-", b"wb", header)
    }

    /// Create a new BAM or CRAM file.
    ///
    /// # Arguments
    ///
    /// * `path` - the path. Use "-" for stdout.
    /// * `mode` - write mode, refer to htslib::hts_open()
    /// * `header` - header definition to use
    fn new(path: &[u8], mode: &[u8], header: &header::Header) -> Result<Self, BGZFError> {
        let f = try!(hts_open(&ffi::CString::new(path).unwrap(), mode));

        // sam_hdr_parse does not populate the text and l_text fields of the header_record.
        // This causes non-SQ headers to be dropped in the output BAM file.
        // To avoid this, we copy the full header to a new C-string that is allocated with malloc,
        // and set this into header_record manually.
        let header_record = unsafe {
            let mut header_string = header.to_bytes();
            if !header_string.is_empty() && header_string[header_string.len() - 1] != b'\n' {
                header_string.push(b'\n');
            }
            let l_text = header_string.len();
            let text = ::libc::malloc(l_text + 1);
            libc::memset(text, 0, l_text + 1);
            libc::memcpy(
                text,
                header_string.as_ptr() as *const ::libc::c_void,
                header_string.len(),
            );

            //println!("{}", str::from_utf8(&header_string).unwrap());
            let rec = htslib::sam_hdr_parse((l_text + 1) as i32, text as *const i8);

            (*rec).text = text as *mut i8;
            (*rec).l_text = l_text as u32;
            rec
        };

        unsafe {
            htslib::sam_hdr_write(f, header_record);
        }

        Ok(Writer {
            f: f,
            header: HeaderView::new(header_record),
        })
    }

    /// Activate multi-threaded BAM write support in htslib. This should permit faster
    /// writing of large BAM files.
    ///
    /// # Arguments
    ///
    /// * `n_threads` - number of extra background writer threads to use, must be `> 0`.
    pub fn set_threads(&mut self, n_threads: usize) -> Result<(), ThreadingError> {
        set_threads(self.f, n_threads)
    }

    /// Write record to BAM.
    ///
    /// # Arguments
    ///
    /// * `record` - the record to write
    pub fn write(&mut self, record: &record::Record) -> Result<(), WriteError> {
        if unsafe { htslib::sam_write1(self.f, &mut self.header.inner(), record.inner) } == -1 {
            Err(WriteError::Some)
        } else {
            Ok(())
        }
    }

    /// Return the header.
    pub fn header(&self) -> &HeaderView {
        &self.header
    }

    /// Set the reference path for reading CRAM files.
    ///
    /// # Arguments
    ///
    /// * `path` - path to the FASTA reference
    pub fn set_reference<P: AsRef<Path>>(&mut self, path: P) -> Result<(), ReferencePathError> {
        let fai_path_str = format!("{}.fai", path.as_ref().to_str().unwrap());
        let fai_path = Path::new(&fai_path_str);

        match fai_path.to_str() {
            Some(p) if fai_path.exists() => set_fai_filename(self.f, p.as_bytes()),
            _ => Err(ReferencePathError::InvalidPath),
        }
    }
}

impl Drop for Writer {
    fn drop(&mut self) {
        unsafe {
            htslib::hts_close(self.f);
        }
    }
}

/// Iterator over the records of a BAM.
#[derive(Debug)]
pub struct Records<'a, R: 'a + Read> {
    reader: &'a mut R,
}

impl<'a, R: Read> Iterator for Records<'a, R> {
    type Item = Result<record::Record, ReadError>;

    fn next(&mut self) -> Option<Result<record::Record, ReadError>> {
        let mut record = record::Record::new();
        match self.reader.read(&mut record) {
            Err(ReadError::NoMoreRecord) => None,
            Ok(()) => Some(Ok(record)),
            Err(err) => Some(Err(err)),
        }
    }
}

/// Iterator over the records of a BAM until the virtual offset is less than `end`
pub struct ChunkIterator<'a, R: 'a + Read> {
    reader: &'a mut R,
    end: Option<i64>,
}

impl<'a, R: Read> Iterator for ChunkIterator<'a, R> {
    type Item = Result<record::Record, ReadError>;
    fn next(&mut self) -> Option<Result<record::Record, ReadError>> {
        if let Some(pos) = self.end {
            if self.reader.tell() >= pos {
                return None;
            }
        }
        let mut record = record::Record::new();
        match self.reader.read(&mut record) {
            Err(ReadError::NoMoreRecord) => None,
            Ok(()) => Some(Ok(record)),
            Err(err) => Some(Err(err)),
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum ReadError {
        Truncated {
            description("truncated record")
        }
        Invalid {
            description("invalid record")
        }
        NoMoreRecord {
            description("no more record")
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
    pub enum IndexedReaderError {
        InvalidIndex {
            description("invalid index")
        }
        BGZFError(err: BGZFError) {
            from()
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum WriterPathError {
        InvalidPath {
            description("invalid path")
        }
        BGZFError(err: BGZFError) {
            from()
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
    pub enum BGZFError {
        Some {
            description("error reading BGZF file")
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum ReaderPathError {
        InvalidPath {
            description("invalid path")
        }
        BGZFError(err: BGZFError) {
            from()
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum ReferencePathError {
        InvalidPath {
            description("invalid path")
        }
        BGZFError(err: BGZFError) {
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

quick_error! {
    #[derive(Debug, Clone)]
    pub enum WriteError {
        Some {
            description("error writing record")
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

quick_error! {
    #[derive(Debug, Clone)]
    pub enum SeekError {
        Some {
            description("error seeking to voffset")
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum AuxWriteError {
        Some {
            description("error pushing aux data to record")
        }
    }
}

/// Wrapper for opening a BAM file.
fn hts_open(path: &ffi::CStr, mode: &[u8]) -> Result<*mut htslib::htsFile, BGZFError> {
    let ret = unsafe { htslib::hts_open(path.as_ptr(), ffi::CString::new(mode).unwrap().as_ptr()) };
    if ret.is_null() {
        Err(BGZFError::Some)
    } else {
        if !mode.contains(&b'w') {
            unsafe {
                if (*ret).format.category != htslib::htsFormatCategory_sequence_data {
                    return Err(BGZFError::Some);
                }
            }
        }
        Ok(ret)
    }
}

/// Wrapper for iterating an indexed BAM file.
fn itr_next(
    htsfile: *mut htslib::htsFile,
    itr: *mut htslib::hts_itr_t,
    record: *mut htslib::bam1_t,
) -> i32 {
    unsafe {
        htslib::hts_itr_next(
            (*htsfile).fp.bgzf,
            itr,
            record as *mut ::std::os::raw::c_void,
            htsfile as *mut ::std::os::raw::c_void,
        )
    }
}

#[derive(Debug)]
pub struct HeaderView {
    inner: *mut htslib::bam_hdr_t,
    owned: bool,
}

impl HeaderView {
    /// Create a new HeaderView from a pre-populated Header object
    pub fn from_header(header: &Header) -> Self {
        let mut header_string = header.to_bytes();
        if !header_string.is_empty() && header_string[header_string.len() - 1] != b'\n' {
            header_string.push(b'\n');
        }
        Self::from_bytes(&header_string)
    }

    /// Create a new HeaderView from bytes
    pub fn from_bytes(header_string: &[u8]) -> Self {
        let header_record = unsafe {
            let l_text = header_string.len();
            let text = ::libc::malloc(l_text + 1);
            ::libc::memset(text, 0, l_text + 1);
            ::libc::memcpy(
                text,
                header_string.as_ptr() as *const ::libc::c_void,
                header_string.len(),
            );

            let rec = htslib::sam_hdr_parse((l_text + 1) as i32, text as *const i8);
            (*rec).text = text as *mut i8;
            (*rec).l_text = l_text as u32;
            rec
        };

        HeaderView::new(header_record)
    }

    /// Create a new HeaderView from the underlying Htslib type, and own it.
    pub fn new(inner: *mut htslib::bam_hdr_t) -> Self {
        HeaderView {
            inner: inner,
            owned: true,
        }
    }

    #[inline]
    pub fn inner(&self) -> htslib::bam_hdr_t {
        unsafe { (*self.inner) }
    }

    #[inline]
    // Pointer to inner bam_hdr_t struct
    pub fn inner_ptr(&self) -> *const htslib::bam_hdr_t {
        self.inner
    }

    #[inline]
    // Mutable pointer to bam_hdr_t struct
    pub fn inner_ptr_mut(&self) -> *mut htslib::bam_hdr_t {
        self.inner
    }

    pub fn tid(&self, name: &[u8]) -> Option<u32> {
        let tid = unsafe {
            htslib::bam_name2id(
                self.inner,
                ffi::CString::new(name)
                    .ok()
                    .expect("Expected valid name.")
                    .as_ptr(),
            )
        };
        if tid < 0 {
            None
        } else {
            Some(tid as u32)
        }
    }

    pub fn target_count(&self) -> u32 {
        self.inner().n_targets as u32
    }

    pub fn target_names(&self) -> Vec<&[u8]> {
        let names = unsafe {
            slice::from_raw_parts(self.inner().target_name, self.target_count() as usize)
        };
        names
            .iter()
            .map(|name| unsafe { ffi::CStr::from_ptr(*name).to_bytes() })
            .collect()
    }

    pub fn target_len(&self, tid: u32) -> Option<u32> {
        let inner = unsafe { *self.inner };
        if (tid as i32) < inner.n_targets {
            let l: &[u32] =
                unsafe { slice::from_raw_parts(inner.target_len, inner.n_targets as usize) };
            Some(l[tid as usize])
        } else {
            None
        }
    }

    /// Retrieve the textual SAM header as bytes
    pub fn as_bytes<'a>(&'a self) -> &'a [u8] {
        unsafe { ffi::CStr::from_ptr((*self.inner).text).to_bytes() }
    }
}

impl Clone for HeaderView {
    fn clone(&self) -> Self {
        HeaderView {
            inner: unsafe { htslib::bam_hdr_dup(self.inner) },
            owned: true,
        }
    }
}

impl Drop for HeaderView {
    fn drop(&mut self) {
        if self.owned {
            unsafe {
                htslib::bam_hdr_destroy(self.inner);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    extern crate tempdir;
    use super::header::HeaderRecord;
    use super::record::{Aux, Cigar, CigarString};
    use super::*;
    use std::collections::HashMap;
    use std::path::Path;
    use std::str;

    fn gold() -> (
        [&'static [u8]; 6],
        [u16; 6],
        [&'static [u8]; 6],
        [&'static [u8]; 6],
        [CigarString; 6],
    ) {
        let names = [
            &b"I"[..],
            &b"II.14978392"[..],
            &b"III"[..],
            &b"IV"[..],
            &b"V"[..],
            &b"VI"[..],
        ];
        let flags = [16u16, 16u16, 16u16, 16u16, 16u16, 2048u16];
        let seqs = [
            &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCC\
TAAGCCTAAGCCTAAGCCTAA"[..],
            &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCC\
TAAGCCTAAGCCTAAGCCTAA"[..],
            &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCC\
TAAGCCTAAGCCTAAGCCTAA"[..],
            &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCC\
TAAGCCTAAGCCTAAGCCTAA"[..],
            &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCC\
TAAGCCTAAGCCTAAGCCTAA"[..],
            &b"ACTAAGCCTAAGCCTAAGCCTAAGCCAATTATCGATTTCTGAAAAAATTATCGAATTTTCTAGAAATTTTGCAAATTTT\
TTCATAAAATTATCGATTTTA"[..],
        ];
        let quals = [
            &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
CCCCCCCCCCCCCCCCCCC"[..],
            &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
CCCCCCCCCCCCCCCCCCC"[..],
            &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
CCCCCCCCCCCCCCCCCCC"[..],
            &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
CCCCCCCCCCCCCCCCCCC"[..],
            &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
CCCCCCCCCCCCCCCCCCC"[..],
            &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
CCCCCCCCCCCCCCCCCCC"[..],
        ];
        let cigars = [
            CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
            CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
            CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
            CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
            CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
            CigarString(vec![Cigar::Match(27), Cigar::Del(100000), Cigar::Match(73)]),
        ];
        (names, flags, seqs, quals, cigars)
    }

    #[test]
    fn test_read() {
        let (names, flags, seqs, quals, cigars) = gold();
        let mut bam = Reader::from_path(&Path::new("test/test.bam"))
            .ok()
            .expect("Error opening file.");
        let del_len = [1, 1, 1, 1, 1, 100000];

        for (i, record) in bam.records().enumerate() {
            let mut rec = record.ok().expect("Expected valid record");
            println!("{}", str::from_utf8(rec.qname()).ok().unwrap());
            assert_eq!(rec.qname(), names[i]);
            assert_eq!(rec.flags(), flags[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);

            let cigar = rec.cigar();
            assert_eq!(*cigar, cigars[i]);

            if let Ok(end_pos) = cigar.end_pos() {
                assert_eq!(end_pos, rec.pos() + 100 + del_len[i]);
                assert_eq!(
                    cigar
                        .read_pos(end_pos as u32 - 10, false, false)
                        .unwrap()
                        .unwrap(),
                    90
                );
            } else {
                panic!("bug: failed to fetch cigar.end_pos() in test_read()")
            }
            assert_eq!(
                cigar
                    .read_pos(rec.pos() as u32 + 20, false, false)
                    .unwrap()
                    .unwrap(),
                20
            );
            assert_eq!(cigar.read_pos(4000000, false, false).unwrap(), None);
            // fix qual offset
            let qual: Vec<u8> = quals[i].iter().map(|&q| q - 33).collect();
            assert_eq!(rec.qual(), &qual[..]);
        }
    }

    #[test]
    fn test_seek() {
        let mut bam = Reader::from_path(&Path::new("test/test.bam"))
            .ok()
            .expect("Error opening file.");

        let mut names_by_voffset = HashMap::new();

        let mut offset = bam.tell();
        let mut rec = Record::new();
        loop {
            if let Err(e) = bam.read(&mut rec) {
                if e.is_eof() {
                    break;
                } else {
                    panic!("error reading bam");
                }
            }
            let qname = str::from_utf8(rec.qname()).unwrap().to_string();
            println!("{} {}", offset, qname);
            names_by_voffset.insert(offset, qname);
            offset = bam.tell();
        }

        for (offset, qname) in names_by_voffset.iter() {
            println!("{} {}", offset, qname);
            bam.seek(*offset).unwrap();
            bam.read(&mut rec).unwrap();
            let rec_qname = str::from_utf8(rec.qname()).ok().unwrap().to_string();
            assert_eq!(qname, &rec_qname);
        }
    }

    #[test]
    fn test_read_sam_header() {
        let bam = Reader::from_path(&"test/test.bam")
            .ok()
            .expect("Error opening file.");

        let true_header =
            "@SQ\tSN:CHROMOSOME_I\tLN:15072423\n@SQ\tSN:CHROMOSOME_II\tLN:15279345\
             \n@SQ\tSN:CHROMOSOME_III\tLN:13783700\n@SQ\tSN:CHROMOSOME_IV\tLN:17493793\n@SQ\t\
             SN:CHROMOSOME_V\tLN:20924149\n"
                .to_string();
        let header_text = String::from_utf8(bam.header.as_bytes().to_owned()).unwrap();
        assert_eq!(header_text, true_header);
    }

    fn _test_read_indexed_common(mut bam: IndexedReader) {
        let (names, flags, seqs, quals, cigars) = gold();
        let tid = bam.header.tid(b"CHROMOSOME_I").expect("Expected tid.");
        assert!(bam.header.target_len(tid).expect("Expected target len.") == 15072423);

        // fetch to position containing reads
        bam.fetch(tid, 0, 2)
            .ok()
            .expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);

        // compare reads
        bam.fetch(tid, 0, 2)
            .ok()
            .expect("Expected successful fetch.");
        for (i, record) in bam.records().enumerate() {
            let mut rec = record.ok().expect("Expected valid record");

            println!("{}", str::from_utf8(rec.qname()).ok().unwrap());
            assert_eq!(rec.qname(), names[i]);
            assert_eq!(rec.flags(), flags[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);
            assert_eq!(*rec.cigar(), cigars[i]);
            // fix qual offset
            let qual: Vec<u8> = quals[i].iter().map(|&q| q - 33).collect();
            assert_eq!(rec.qual(), &qual[..]);
            assert_eq!(rec.aux(b"NotAvailableAux"), None);
        }

        // fetch to empty position
        bam.fetch(2, 1, 1).ok().expect("Expected successful fetch.");
        assert!(bam.records().count() == 0);
    }

    #[test]
    fn test_read_indexed() {
        let bam = IndexedReader::from_path(&"test/test.bam")
            .ok()
            .expect("Expected valid index.");
        _test_read_indexed_common(bam);
    }
    #[test]
    fn test_read_indexed_different_index_name() {
        let bam = IndexedReader::from_path_and_index(
            &"test/test_different_index_name.bam",
            &"test/test.bam.bai",
        )
        .ok()
        .expect("Expected valid index.");
        _test_read_indexed_common(bam);
    }

    #[test]
    fn test_set_record() {
        let (names, _, seqs, quals, cigars) = gold();

        let mut rec = record::Record::new();
        rec.set_reverse();
        rec.set(names[0], &cigars[0], seqs[0], quals[0]);
        // note: this segfaults if you push_aux() before set()
        //       because set() obliterates aux
        rec.push_aux(b"NM", &Aux::Integer(15)).unwrap();

        assert_eq!(rec.qname(), names[0]);
        assert_eq!(*rec.cigar(), cigars[0]);
        assert_eq!(rec.seq().as_bytes(), seqs[0]);
        assert_eq!(rec.qual(), quals[0]);
        assert!(rec.is_reverse());
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));
    }

    #[test]
    fn test_set_qname() {
        let (names, _, seqs, quals, cigars) = gold();

        assert!(names[0] != names[1]);

        for i in 0..names.len() {
            let mut rec = record::Record::new();
            rec.set(names[i], &cigars[i], seqs[i], quals[i]);
            rec.push_aux(b"NM", &Aux::Integer(15)).unwrap();

            assert_eq!(rec.qname(), names[i]);
            assert_eq!(*rec.cigar(), cigars[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);
            assert_eq!(rec.qual(), quals[i]);
            assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));

            // Equal length qname
            assert!(rec.qname()[0] != 'X' as u8);
            rec.set_qname(b"X");
            assert_eq!(rec.qname(), b"X");

            // Longer qname
            let mut longer_name = names[i].to_owned().clone();
            let extension = b"BuffaloBUffaloBUFFaloBUFFAloBUFFALoBUFFALO";
            longer_name.extend(extension.iter());
            rec.set_qname(&longer_name);

            assert_eq!(rec.qname(), longer_name.as_slice());
            assert_eq!(*rec.cigar(), cigars[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);
            assert_eq!(rec.qual(), quals[i]);
            assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));

            // Shorter qname
            let shorter_name = b"42";
            rec.set_qname(shorter_name);

            assert_eq!(rec.qname(), shorter_name);
            assert_eq!(*rec.cigar(), cigars[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);
            assert_eq!(rec.qual(), quals[i]);
            assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));

            // Zero-length qname
            rec.set_qname(b"");

            assert_eq!(rec.qname(), b"");
            assert_eq!(*rec.cigar(), cigars[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);
            assert_eq!(rec.qual(), quals[i]);
            assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));
        }
    }

    #[test]
    fn test_set_qname2() {
        let mut _header = Header::new();
        _header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"1")
                .push_tag(b"LN", &10000000),
        );
        let header = HeaderView::from_header(&_header);

        let line =
            b"blah1	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1\tli:i:0\ttf:Z:cC";

        let mut rec = Record::from_sam(&header, line).unwrap();
        assert_eq!(rec.qname(), b"blah1");
        rec.set_qname(b"r0");
        assert_eq!(rec.qname(), b"r0");
    }

    #[test]
    fn test_remove_aux() {
        let mut bam = Reader::from_path(&Path::new("test/test.bam"))
            .ok()
            .expect("Error opening file.");

        for record in bam.records() {
            let rec = record.ok().expect("Expected valid record");

            if rec.aux(b"XS").is_some() {
                rec.remove_aux(b"XS");
            }

            if rec.aux(b"YT").is_some() {
                rec.remove_aux(b"YT");
            }

            rec.remove_aux(b"ab");

            assert_eq!(rec.aux(b"XS"), None);
            assert_eq!(rec.aux(b"YT"), None);
        }
    }

    #[test]
    fn test_write() {
        let (names, _, seqs, quals, cigars) = gold();

        let tmp = tempdir::TempDir::new("rust-htslib")
            .ok()
            .expect("Cannot create temp dir");
        let bampath = tmp.path().join("test.bam");
        println!("{:?}", bampath);
        {
            let mut bam = Writer::from_path(
                &bampath,
                Header::new().push_record(
                    HeaderRecord::new(b"SQ")
                        .push_tag(b"SN", &"chr1")
                        .push_tag(b"LN", &15072423),
                ),
            )
            .ok()
            .expect("Error opening file.");

            for i in 0..names.len() {
                let mut rec = record::Record::new();
                rec.set(names[i], &cigars[i], seqs[i], quals[i]);
                rec.push_aux(b"NM", &Aux::Integer(15)).unwrap();

                bam.write(&mut rec).ok().expect("Failed to write record.");
            }
        }

        {
            let mut bam = Reader::from_path(&bampath)
                .ok()
                .expect("Error opening file.");

            for i in 0..names.len() {
                let mut rec = record::Record::new();
                bam.read(&mut rec).ok().expect("Failed to read record.");

                assert_eq!(rec.qname(), names[i]);
                assert_eq!(*rec.cigar(), cigars[i]);
                assert_eq!(rec.seq().as_bytes(), seqs[i]);
                assert_eq!(rec.qual(), quals[i]);
                assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));
            }
        }

        tmp.close().ok().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_write_threaded() {
        let (names, _, seqs, quals, cigars) = gold();

        let tmp = tempdir::TempDir::new("rust-htslib")
            .ok()
            .expect("Cannot create temp dir");
        let bampath = tmp.path().join("test.bam");
        println!("{:?}", bampath);
        {
            let mut bam = Writer::from_path(
                &bampath,
                Header::new().push_record(
                    HeaderRecord::new(b"SQ")
                        .push_tag(b"SN", &"chr1")
                        .push_tag(b"LN", &15072423),
                ),
            )
            .ok()
            .expect("Error opening file.");
            bam.set_threads(4).unwrap();

            for i in 0..10000 {
                let mut rec = record::Record::new();
                let idx = i % names.len();
                rec.set(names[idx], &cigars[idx], seqs[idx], quals[idx]);
                rec.push_aux(b"NM", &Aux::Integer(15)).unwrap();
                rec.set_pos(i as i32);

                bam.write(&mut rec).ok().expect("Failed to write record.");
            }
        }

        {
            let mut bam = Reader::from_path(&bampath)
                .ok()
                .expect("Error opening file.");

            for (i, _rec) in bam.records().enumerate() {
                let idx = i % names.len();

                let mut rec = _rec.expect("Failed to read record.");

                assert_eq!(rec.pos(), i as i32);
                assert_eq!(rec.qname(), names[idx]);
                assert_eq!(*rec.cigar(), cigars[idx]);
                assert_eq!(rec.seq().as_bytes(), seqs[idx]);
                assert_eq!(rec.qual(), quals[idx]);
                assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));
            }
        }

        tmp.close().ok().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_copy_template() {
        // Verify that BAM headers are transmitted correctly when using an existing BAM as a
        // template for headers.

        let tmp = tempdir::TempDir::new("rust-htslib")
            .ok()
            .expect("Cannot create temp dir");
        let bampath = tmp.path().join("test.bam");
        println!("{:?}", bampath);

        let mut input_bam = Reader::from_path(&"test/test.bam")
            .ok()
            .expect("Error opening file.");

        {
            let mut bam = Writer::from_path(&bampath, &Header::from_template(&input_bam.header()))
                .ok()
                .expect("Error opening file.");

            for rec in input_bam.records() {
                bam.write(&rec.unwrap())
                    .ok()
                    .expect("Failed to write record.");
            }
        }

        {
            let copy_bam = Reader::from_path(&bampath)
                .ok()
                .expect("Error opening file.");

            // Verify that the header came across correctly
            assert_eq!(input_bam.header().as_bytes(), copy_bam.header().as_bytes());
        }

        tmp.close().ok().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_pileup() {
        let (_, _, seqs, quals, _) = gold();

        let mut bam = Reader::from_path(&"test/test.bam")
            .ok()
            .expect("Error opening file.");
        let pileups = bam.pileup();
        for pileup in pileups.take(26) {
            let _pileup = pileup.ok().expect("Expected successful pileup.");
            let pos = _pileup.pos() as usize;
            assert_eq!(_pileup.depth(), 6);
            assert!(_pileup.tid() == 0);
            for (i, a) in _pileup.alignments().enumerate() {
                assert_eq!(a.indel(), pileup::Indel::None);
                let qpos = a.qpos().unwrap();
                assert_eq!(qpos, pos - 1);
                assert_eq!(a.record().seq()[qpos], seqs[i][qpos]);
                assert_eq!(a.record().qual()[qpos], quals[i][qpos] - 33);
            }
        }
    }

    #[test]
    fn test_idx_pileup() {
        let mut bam = IndexedReader::from_path(&"test/test.bam")
            .ok()
            .expect("Error opening file.");
        // read without fetch
        for pileup in bam.pileup() {
            pileup.unwrap();
        }
        // go back again
        let tid = bam.header().tid(b"CHROMOSOME_I").unwrap();
        bam.fetch(tid, 0, 5).unwrap();
        for p in bam.pileup() {
            println!("{}", p.unwrap().pos())
        }
    }

    #[test]
    fn parse_from_sam() {
        use std::fs::File;
        use std::io::Read;

        let bamfile = "./test/bam2sam_test.bam";
        let samfile = "./test/bam2sam_expected.sam";

        // Load BAM file:
        let mut rdr = Reader::from_path(bamfile).unwrap();
        let bam_recs: Vec<Record> = rdr.records().map(|v| v.unwrap()).collect();

        let mut sam = Vec::new();
        assert!(File::open(samfile).unwrap().read_to_end(&mut sam).is_ok());

        let sam_recs: Vec<Record> = sam
            .split(|x| *x == b'\n')
            .filter(|x| x.len() > 0 && x[0] != b'@')
            .map(|line| Record::from_sam(rdr.header(), line).unwrap())
            .collect();

        for (b1, s1) in bam_recs.iter().zip(sam_recs.iter()) {
            assert!(b1 == s1);
        }
    }

    #[test]
    fn test_cigar_modes() {
        // test the cached and uncached ways of getting the cigar string.

        let (_, _, _, _, cigars) = gold();
        let mut bam = Reader::from_path(&Path::new("test/test.bam"))
            .ok()
            .expect("Error opening file.");

        for (i, record) in bam.records().enumerate() {
            let mut rec = record.ok().expect("Expected valid record");

            let cigar = rec.cigar();
            assert_eq!(*cigar, cigars[i]);
        }

        for (i, record) in bam.records().enumerate() {
            let mut rec = record.ok().expect("Expected valid record");
            rec.cache_cigar();

            let cigar = rec.cigar_cached().unwrap();
            assert_eq!(**cigar, cigars[i]);

            let cigar = rec.cigar();
            assert_eq!(*cigar, cigars[i]);
        }
    }

    #[test]
    fn test_read_cram() {
        let cram_path = "./test/test_cram.cram";
        let bam_path = "./test/test_cram.bam";
        let ref_path = "./test/test_cram.fa";

        // Load CRAM file, records
        let mut cram_reader = Reader::from_path(cram_path).unwrap();
        cram_reader.set_reference(ref_path).unwrap();
        let cram_records: Vec<Record> = cram_reader.records().map(|v| v.unwrap()).collect();

        // Load BAM file, records
        let mut bam_reader = Reader::from_path(bam_path).unwrap();
        let bam_records: Vec<Record> = bam_reader.records().map(|v| v.unwrap()).collect();

        // Compare CRAM records to BAM records
        for (c1, b1) in cram_records.iter().zip(bam_records.iter()) {
            assert!(c1 == b1);
        }
    }

    #[test]
    fn test_write_cram() {
        // BAM file, records
        let bam_path = "./test/test_cram.bam";
        let ref_path = "./test/test_cram.fa";
        let mut bam_reader = Reader::from_path(bam_path).unwrap();
        let bam_records: Vec<Record> = bam_reader.records().map(|v| v.unwrap()).collect();

        // New CRAM file
        let tmp = tempdir::TempDir::new("rust-htslib")
            .ok()
            .expect("Cannot create temp dir");
        let cram_path = tmp.path().join("test.cram");

        // Write BAM records to new CRAM file
        {
            let mut header = Header::new();
            header.push_record(
                HeaderRecord::new(b"HD")
                    .push_tag(b"VN", &"1.5")
                    .push_tag(b"SO", &"coordinate"),
            );
            header.push_record(
                HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", &"chr1")
                    .push_tag(b"LN", &120)
                    .push_tag(b"M5", &"20a9a0fb770814e6c5e49946750f9724")
                    .push_tag(b"UR", &"test/test_cram.fa"),
            );
            header.push_record(
                HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", &"chr2")
                    .push_tag(b"LN", &120)
                    .push_tag(b"M5", &"7a2006ccca94ea92b6dae5997e1b0d70")
                    .push_tag(b"UR", &"test/test_cram.fa"),
            );
            header.push_record(
                HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", &"chr3")
                    .push_tag(b"LN", &120)
                    .push_tag(b"M5", &"a66b336bfe3ee8801c744c9545c87e24")
                    .push_tag(b"UR", &"test/test_cram.fa"),
            );

            let mut cram_writer = Writer::from_cram_path(&cram_path, &header)
                .ok()
                .expect("Error opening CRAM file.");
            cram_writer.set_reference(ref_path).unwrap();

            // Write BAM records to CRAM file
            for rec in bam_records.iter() {
                cram_writer
                    .write(&rec)
                    .ok()
                    .expect("Faied to write record to CRAM.");
            }
        }

        // Compare written CRAM records with BAM records
        {
            // Load written CRAM file
            let mut cram_reader = Reader::from_path(cram_path).unwrap();
            cram_reader.set_reference(ref_path).unwrap();
            let cram_records: Vec<Record> = cram_reader.records().map(|v| v.unwrap()).collect();

            // Compare CRAM records to BAM records
            for (c1, b1) in cram_records.iter().zip(bam_records.iter()) {
                assert!(c1 == b1);
            }
        }

        tmp.close().ok().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_bam_fails_on_vcf() {
        let bam_path = "./test/test_left.vcf";
        let bam_reader = Reader::from_path(bam_path);
        assert!(bam_reader.is_err());
    }

    #[test]
    fn test_indexde_bam_fails_on_vcf() {
        let bam_path = "./test/test_left.vcf";
        let bam_reader = IndexedReader::from_path(bam_path);
        assert!(bam_reader.is_err());
    }

    #[test]
    fn test_bam_fails_on_toml() {
        let bam_path = "./cargo.toml";
        let bam_reader = Reader::from_path(bam_path);
        assert!(bam_reader.is_err());
    }
}
