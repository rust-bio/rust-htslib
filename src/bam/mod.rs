// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Module for working with SAM, BAM, and CRAM files.

pub mod buffer;
pub mod ext;
pub mod header;
pub mod index;
pub mod pileup;
pub mod record;

#[cfg(feature = "serde_feature")]
pub mod record_serde;

use std::ffi;
use std::os::raw::c_char;
use std::path::Path;
use std::rc::Rc;
use std::slice;
use std::str;
use std::sync::Arc;

use url::Url;

use crate::errors::{Error, Result};
use crate::htslib;
use crate::tpool::ThreadPool;
use crate::utils::path_as_bytes;

pub use crate::bam::buffer::RecordBuffer;
pub use crate::bam::header::Header;
pub use crate::bam::record::Record;
use hts_sys::{hts_fmt_option, sam_fields};
use std::convert::{TryFrom, TryInto};
use std::mem::MaybeUninit;

/// # Safety
///
/// Implementation for `Read::set_threads` and `Writer::set_threads`.
unsafe fn set_threads(htsfile: *mut htslib::htsFile, n_threads: usize) -> Result<()> {
    assert!(n_threads != 0, "n_threads must be > 0");

    if htslib::hts_set_threads(htsfile, n_threads as i32) != 0 {
        Err(Error::SetThreads)
    } else {
        Ok(())
    }
}

unsafe fn set_thread_pool(htsfile: *mut htslib::htsFile, tpool: &ThreadPool) -> Result<()> {
    let mut b = tpool.handle.borrow_mut();

    if htslib::hts_set_thread_pool(htsfile, &mut b.inner as *mut _) != 0 {
        Err(Error::ThreadPool)
    } else {
        Ok(())
    }
}

/// # Safety
///
/// Set the reference FAI index path in a `htslib::htsFile` struct for reading CRAM format.
pub unsafe fn set_fai_filename<P: AsRef<Path>>(
    htsfile: *mut htslib::htsFile,
    fasta_path: P,
) -> Result<()> {
    let path = if let Some(ext) = fasta_path.as_ref().extension() {
        fasta_path
            .as_ref()
            .with_extension(format!("{}.fai", ext.to_str().unwrap()))
    } else {
        fasta_path.as_ref().with_extension(".fai")
    };
    let p: &Path = path.as_ref();
    let c_str = ffi::CString::new(p.to_str().unwrap().as_bytes()).unwrap();
    if htslib::hts_set_fai_filename(htsfile, c_str.as_ptr()) == 0 {
        Ok(())
    } else {
        Err(Error::BamInvalidReferencePath { path: p.to_owned() })
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
    ///
    /// # Returns
    ///
    /// Some(Ok(())) if the record was read and None if no more records to read
    ///
    /// Example:
    /// ```
    /// use rust_htslib::errors::Error;
    /// use rust_htslib::bam::{Read, IndexedReader, Record};
    ///
    /// let mut bam = IndexedReader::from_path("test/test.bam").unwrap();
    /// bam.fetch((0, 1000, 2000)); // reads on tid 0, from 1000bp to 2000bp
    /// let mut record = Record::new();
    /// while let Some(result) = bam.read(&mut record) {
    ///     match result {
    ///         Ok(_) => {
    ///             println!("Read sequence: {:?}", record.seq().as_bytes());
    ///         }
    ///         Err(_) => panic!("BAM parsing failed...")
    ///     }
    /// }
    /// ```
    ///
    /// Consider using [`rc_records`](#tymethod.rc_records) instead.
    fn read(&mut self, record: &mut record::Record) -> Option<Result<()>>;

    /// Iterator over the records of the seeked region.
    /// Note that, while being convenient, this is less efficient than pre-allocating a
    /// `Record` and reading into it with the `read` method, since every iteration involves
    /// the allocation of a new `Record`.
    ///
    /// This is about 10% slower than record in micro benchmarks.
    ///
    /// Consider using [`rc_records`](#tymethod.rc_records) instead.
    fn records(&mut self) -> Records<'_, Self>;

    /// Records iterator using an Rc to avoid allocating a Record each turn.
    /// This is about 1% slower than the [`read`](#tymethod.read) based API in micro benchmarks,
    /// but has nicer ergonomics (and might not actually be slower in your applications).
    ///
    /// Example:
    /// ```
    /// use rust_htslib::errors::Error;
    /// use rust_htslib::bam::{Read, Reader, Record};
    /// use rust_htslib::htslib; // for BAM_F*
    /// let mut bam = Reader::from_path("test/test.bam").unwrap();
    ///
    /// for read in
    ///     bam.rc_records()
    ///     .map(|x| x.expect("Failure parsing Bam file"))
    ///     .filter(|read|
    ///         read.flags()
    ///          & (htslib::BAM_FUNMAP
    ///              | htslib::BAM_FSECONDARY
    ///              | htslib::BAM_FQCFAIL
    ///              | htslib::BAM_FDUP) as u16
    ///          == 0
    ///     )
    ///     .filter(|read| !read.is_reverse()) {
    ///     println!("Found a forward read: {:?}", read.qname());
    /// }
    ///
    /// //or to add the read qnames into a Vec
    /// let collected: Vec<_> = bam.rc_records().map(|read| read.unwrap().qname().to_vec()).collect();
    ///
    ///
    /// ```
    fn rc_records(&mut self) -> RcRecords<'_, Self>;

    /// Iterator over pileups.
    fn pileup(&mut self) -> pileup::Pileups<'_, Self>;

    /// Return the htsFile struct
    fn htsfile(&self) -> *mut htslib::htsFile;

    /// Return the header.
    fn header(&self) -> &HeaderView;

    /// Seek to the given virtual offset in the file
    fn seek(&mut self, offset: i64) -> Result<()> {
        let htsfile = unsafe { self.htsfile().as_ref() }.expect("bug: null pointer to htsFile");
        let ret = match htsfile.format.format {
            htslib::htsExactFormat_cram => unsafe {
                i64::from(htslib::cram_seek(
                    htsfile.fp.cram,
                    offset as libc::off_t,
                    libc::SEEK_SET,
                ))
            },
            _ => unsafe { htslib::bgzf_seek(htsfile.fp.bgzf, offset, libc::SEEK_SET) },
        };

        if ret == 0 {
            Ok(())
        } else {
            Err(Error::FileSeek)
        }
    }

    /// Report the current virtual offset
    fn tell(&self) -> i64 {
        // this reimplements the bgzf_tell macro
        let htsfile = unsafe { self.htsfile().as_ref() }.expect("bug: null pointer to htsFile");
        let bgzf = unsafe { *htsfile.fp.bgzf };
        (bgzf.block_address << 16) | (i64::from(bgzf.block_offset) & 0xFFFF)
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
    fn set_threads(&mut self, n_threads: usize) -> Result<()> {
        unsafe { set_threads(self.htsfile(), n_threads) }
    }

    /// Use a shared thread-pool for writing. This permits controlling the total
    /// thread count when multiple readers and writers are working simultaneously.
    /// A thread pool can be created with `crate::tpool::ThreadPool::new(n_threads)`
    ///
    /// # Arguments
    ///
    /// * `tpool` - thread pool to use for compression work.
    fn set_thread_pool(&mut self, tpool: &ThreadPool) -> Result<()>;

    /// If the underlying file is in CRAM format, allows modifying CRAM options.
    /// Note that this method does *not* check that the underlying file actually is in CRAM format.
    ///
    /// # Examples
    ///
    /// Set the required fields to RNAME and FLAG,
    /// potentially allowing htslib to skip over the rest,
    /// resulting in faster iteration:
    /// ```
    /// use rust_htslib::bam::{Read, Reader};
    /// use hts_sys;
    /// let mut cram = Reader::from_path("test/test_cram.cram").unwrap();
    /// cram.set_cram_options(hts_sys::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
    ///             hts_sys::sam_fields_SAM_RNAME | hts_sys::sam_fields_SAM_FLAG).unwrap();
    /// ```
    fn set_cram_options(&mut self, fmt_opt: hts_fmt_option, fields: sam_fields) -> Result<()> {
        unsafe {
            if hts_sys::hts_set_opt(self.htsfile(), fmt_opt, fields) != 0 {
                Err(Error::HtsSetOpt)
            } else {
                Ok(())
            }
        }
    }
}

/// A BAM reader.
#[derive(Debug)]
pub struct Reader {
    htsfile: *mut htslib::htsFile,
    header: Arc<HeaderView>,
    tpool: Option<ThreadPool>,
}

unsafe impl Send for Reader {}

impl Reader {
    /// Create a new Reader from path.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::new(&path_as_bytes(path, true)?)
    }

    /// Create a new Reader from STDIN.
    pub fn from_stdin() -> Result<Self> {
        Self::new(b"-")
    }

    /// Create a new Reader from URL.
    pub fn from_url(url: &Url) -> Result<Self> {
        Self::new(url.as_str().as_bytes())
    }

    /// Create a new Reader.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open. Use "-" for stdin.
    fn new(path: &[u8]) -> Result<Self> {
        let htsfile = hts_open(path, b"r")?;

        let header = unsafe { htslib::sam_hdr_read(htsfile) };
        if header.is_null() {
            return Err(Error::BamOpen {
                target: String::from_utf8_lossy(path).to_string(),
            });
        }

        // Invalidate the `text` representation of the header
        unsafe {
            let _ = htslib::sam_hdr_line_name(header, b"SQ".as_ptr().cast::<c_char>(), 0);
        }

        Ok(Reader {
            htsfile,
            header: Arc::new(HeaderView::new(header)),
            tpool: None,
        })
    }

    extern "C" fn pileup_read(
        data: *mut ::std::os::raw::c_void,
        record: *mut htslib::bam1_t,
    ) -> i32 {
        let mut _self = unsafe { (data as *mut Self).as_mut().unwrap() };
        unsafe {
            htslib::sam_read1(
                _self.htsfile(),
                _self.header().inner_ptr() as *mut hts_sys::sam_hdr_t,
                record,
            )
        }
    }

    /// Iterator over the records between the (optional) virtual offsets `start` and `end`
    ///
    /// # Arguments
    ///
    /// * `start` - Optional starting virtual offset to seek to. Throws an error if it is not
    ///   a valid virtual offset.
    ///
    /// * `end` - Read until the virtual offset is less than `end`
    pub fn iter_chunk(&mut self, start: Option<i64>, end: Option<i64>) -> ChunkIterator<'_, Self> {
        if let Some(pos) = start {
            self.seek(pos)
                .expect("Failed to seek to the starting position");
        };

        ChunkIterator { reader: self, end }
    }

    /// Set the reference path for reading CRAM files.
    ///
    /// # Arguments
    ///
    /// * `path` - path to the FASTA reference
    pub fn set_reference<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        unsafe { set_fai_filename(self.htsfile, path) }
    }
}

impl Read for Reader {
    /// Read the next BAM record into the given `Record`.
    /// Returns `None` if there are no more records.
    ///
    /// This method is useful if you want to read records as fast as possible as the
    /// `Record` can be reused. A more ergonomic approach is to use the [records](Reader::records)
    /// iterator.
    ///
    /// # Errors
    /// If there are any issues with reading the next record an error will be returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use rust_htslib::errors::Error;
    /// use rust_htslib::bam::{Read, Reader, Record};
    ///
    /// let mut bam = Reader::from_path("test/test.bam")?;
    /// let mut record = Record::new();
    ///
    /// // Print the TID of each record
    /// while let Some(r) = bam.read(&mut record) {
    ///    r.expect("Failed to parse record");
    ///    println!("TID: {}", record.tid())
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    fn read(&mut self, record: &mut record::Record) -> Option<Result<()>> {
        match unsafe {
            htslib::sam_read1(
                self.htsfile,
                self.header().inner_ptr() as *mut hts_sys::sam_hdr_t,
                record.inner_ptr_mut(),
            )
        } {
            -1 => None,
            -2 => Some(Err(Error::BamTruncatedRecord)),
            -4 => Some(Err(Error::BamInvalidRecord)),
            _ => {
                record.set_header(Arc::clone(&self.header));

                Some(Ok(()))
            }
        }
    }

    /// Iterator over the records of the fetched region.
    /// Note that, while being convenient, this is less efficient than pre-allocating a
    /// `Record` and reading into it with the `read` method, since every iteration involves
    /// the allocation of a new `Record`.
    fn records(&mut self) -> Records<'_, Self> {
        Records { reader: self }
    }

    fn rc_records(&mut self) -> RcRecords<'_, Self> {
        RcRecords {
            reader: self,
            record: Rc::new(record::Record::new()),
        }
    }

    fn pileup(&mut self) -> pileup::Pileups<'_, Self> {
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

    fn set_thread_pool(&mut self, tpool: &ThreadPool) -> Result<()> {
        unsafe { set_thread_pool(self.htsfile(), tpool)? }
        self.tpool = Some(tpool.clone());
        Ok(())
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        unsafe {
            htslib::hts_close(self.htsfile);
        }
    }
}

/// Conversion type for start/stop coordinates
/// only public because it's leaked by the conversions
#[doc(hidden)]
pub struct FetchCoordinate(i64);

//the old sam spec
impl From<i32> for FetchCoordinate {
    fn from(coord: i32) -> FetchCoordinate {
        FetchCoordinate(coord as i64)
    }
}

// to support un-annotated literals (type interference fails on those)
impl From<u32> for FetchCoordinate {
    fn from(coord: u32) -> FetchCoordinate {
        FetchCoordinate(coord as i64)
    }
}

//the new sam spec
impl From<i64> for FetchCoordinate {
    fn from(coord: i64) -> FetchCoordinate {
        FetchCoordinate(coord)
    }
}

//what some of our header methods return
impl From<u64> for FetchCoordinate {
    fn from(coord: u64) -> FetchCoordinate {
        FetchCoordinate(coord.try_into().expect("Coordinate exceeded 2^^63-1"))
    }
}

/// Enum for [IndexdReader.fetch()](struct.IndexedReader.html#method.fetch) arguments.
///
/// tids may be converted From<>:
/// * i32 (correct as per spec)
/// * u32 (because of header.tid. Will panic if above 2^31-1).
///
///Coordinates may be (via FetchCoordinate)
/// * i32 (as of the sam v1 spec)
/// * i64 (as of the htslib 'large coordinate' extension (even though they are not supported in BAM)
/// * u32 (because that's what rust literals will default to)
/// * u64 (because of header.target_len(). Will panic if above 2^^63-1).
#[derive(Debug)]
pub enum FetchDefinition<'a> {
    /// tid, start, stop,
    Region(i32, i64, i64),
    /// 'named-reference', start, stop tuple.
    RegionString(&'a [u8], i64, i64),
    ///complete reference. May be i32 or u32 (which panics if above 2^31-')
    CompleteTid(i32),
    ///complete reference by name (&[u8] or &str)
    String(&'a [u8]),
    /// Every read
    All,
    /// Only reads with the BAM flag BAM_FUNMAP (which might not be all reads with reference = -1)
    Unmapped,
}

impl<'a, X: Into<FetchCoordinate>, Y: Into<FetchCoordinate>> From<(i32, X, Y)>
    for FetchDefinition<'a>
{
    fn from(tup: (i32, X, Y)) -> FetchDefinition<'a> {
        let start: FetchCoordinate = tup.1.into();
        let stop: FetchCoordinate = tup.2.into();
        FetchDefinition::Region(tup.0, start.0, stop.0)
    }
}

impl<'a, X: Into<FetchCoordinate>, Y: Into<FetchCoordinate>> From<(u32, X, Y)>
    for FetchDefinition<'a>
{
    fn from(tup: (u32, X, Y)) -> FetchDefinition<'a> {
        let start: FetchCoordinate = tup.1.into();
        let stop: FetchCoordinate = tup.2.into();
        FetchDefinition::Region(
            tup.0.try_into().expect("Tid exceeded 2^31-1"),
            start.0,
            stop.0,
        )
    }
}

//non tuple impls
impl<'a> From<i32> for FetchDefinition<'a> {
    fn from(tid: i32) -> FetchDefinition<'a> {
        FetchDefinition::CompleteTid(tid)
    }
}

impl<'a> From<u32> for FetchDefinition<'a> {
    fn from(tid: u32) -> FetchDefinition<'a> {
        let tid: i32 = tid.try_into().expect("tid exceeded 2^31-1");
        FetchDefinition::CompleteTid(tid)
    }
}

impl<'a> From<&'a str> for FetchDefinition<'a> {
    fn from(s: &'a str) -> FetchDefinition<'a> {
        FetchDefinition::String(s.as_bytes())
    }
}

//also accept &[u8;n] literals
impl<'a> From<&'a [u8]> for FetchDefinition<'a> {
    fn from(s: &'a [u8]) -> FetchDefinition<'a> {
        FetchDefinition::String(s)
    }
}

//also accept &[u8;n] literals
impl<'a, T: AsRef<[u8]>> From<&'a T> for FetchDefinition<'a> {
    fn from(s: &'a T) -> FetchDefinition<'a> {
        FetchDefinition::String(s.as_ref())
    }
}

impl<'a, X: Into<FetchCoordinate>, Y: Into<FetchCoordinate>> From<(&'a str, X, Y)>
    for FetchDefinition<'a>
{
    fn from(tup: (&'a str, X, Y)) -> FetchDefinition<'a> {
        let start: FetchCoordinate = tup.1.into();
        let stop: FetchCoordinate = tup.2.into();
        FetchDefinition::RegionString(tup.0.as_bytes(), start.0, stop.0)
    }
}

impl<'a, X: Into<FetchCoordinate>, Y: Into<FetchCoordinate>> From<(&'a [u8], X, Y)>
    for FetchDefinition<'a>
{
    fn from(tup: (&'a [u8], X, Y)) -> FetchDefinition<'a> {
        let start: FetchCoordinate = tup.1.into();
        let stop: FetchCoordinate = tup.2.into();
        FetchDefinition::RegionString(tup.0, start.0, stop.0)
    }
}

//also accept &[u8;n] literals
impl<'a, T: AsRef<[u8]>, X: Into<FetchCoordinate>, Y: Into<FetchCoordinate>> From<(&'a T, X, Y)>
    for FetchDefinition<'a>
{
    fn from(tup: (&'a T, X, Y)) -> FetchDefinition<'a> {
        let start: FetchCoordinate = tup.1.into();
        let stop: FetchCoordinate = tup.2.into();
        FetchDefinition::RegionString(tup.0.as_ref(), start.0, stop.0)
    }
}

#[derive(Debug)]
pub struct IndexedReader {
    htsfile: *mut htslib::htsFile,
    header: Arc<HeaderView>,
    idx: Option<IndexView>,
    itr: Option<*mut htslib::hts_itr_t>,
    tpool: Option<ThreadPool>,
}

unsafe impl Send for IndexedReader {}

impl IndexedReader {
    /// Create a new Reader from path.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::new(&path_as_bytes(path, true)?)
    }

    pub fn from_path_and_index<P: AsRef<Path>>(path: P, index_path: P) -> Result<Self> {
        Self::new_with_index_path(
            &path_as_bytes(path, true)?,
            &path_as_bytes(index_path, true)?,
        )
    }

    pub fn from_url(url: &Url) -> Result<Self> {
        Self::new(url.as_str().as_bytes())
    }

    /// Create a new Reader.
    ///
    /// # Arguments
    ///
    /// * `path` - the path. Use "-" for stdin.
    fn new(path: &[u8]) -> Result<Self> {
        let htsfile = hts_open(path, b"r")?;
        let header = unsafe { htslib::sam_hdr_read(htsfile) };
        let c_str = ffi::CString::new(path).unwrap();
        let idx = unsafe { htslib::sam_index_load(htsfile, c_str.as_ptr()) };
        if idx.is_null() {
            Err(Error::BamInvalidIndex {
                target: str::from_utf8(path).unwrap().to_owned(),
            })
        } else {
            Ok(IndexedReader {
                htsfile,
                header: Arc::new(HeaderView::new(header)),
                idx: Some(IndexView::new(idx)),
                itr: None,
                tpool: None,
            })
        }
    }
    /// Create a new Reader.
    ///
    /// # Arguments
    ///
    /// * `path` - the path. Use "-" for stdin.
    /// * `index_path` - the index path to use
    fn new_with_index_path(path: &[u8], index_path: &[u8]) -> Result<Self> {
        let htsfile = hts_open(path, b"r")?;
        let header = unsafe { htslib::sam_hdr_read(htsfile) };
        let c_str_path = ffi::CString::new(path).unwrap();
        let c_str_index_path = ffi::CString::new(index_path).unwrap();
        let idx = unsafe {
            htslib::sam_index_load2(htsfile, c_str_path.as_ptr(), c_str_index_path.as_ptr())
        };
        if idx.is_null() {
            Err(Error::BamInvalidIndex {
                target: str::from_utf8(path).unwrap().to_owned(),
            })
        } else {
            Ok(IndexedReader {
                htsfile,
                header: Arc::new(HeaderView::new(header)),
                idx: Some(IndexView::new(idx)),
                itr: None,
                tpool: None,
            })
        }
    }

    /// Define the region from which .read() or .records will retrieve reads.
    ///
    /// Both iterating (with [.records()](trait.Read.html#tymethod.records)) and looping without allocation (with [.read()](trait.Read.html#tymethod.read) are a two stage process:
    /// 1. 'fetch' the region of interest
    /// 2. iter/loop through the reads.
    ///
    /// Example:
    /// ```
    /// use rust_htslib::bam::{IndexedReader, Read};
    /// let mut bam = IndexedReader::from_path("test/test.bam").unwrap();
    /// bam.fetch(("chrX", 10000, 20000)); // coordinates 10000..20000 on reference named "chrX"
    /// for read in bam.records() {
    ///     println!("read name: {:?}", read.unwrap().qname());
    /// }
    /// ```
    ///
    /// The arguments may be anything that can be converted into a FetchDefinition
    /// such as
    ///
    /// * fetch(tid: u32) -> fetch everything on this reference
    /// * fetch(reference_name: &[u8] | &str) -> fetch everything on this reference
    /// * fetch((tid: i32, start: i64, stop: i64)): -> fetch in this region on this tid
    /// * fetch((reference_name: &[u8] | &str, start: i64, stop: i64) -> fetch in this region on this tid
    /// * fetch(FetchDefinition::All) or fetch(".") -> Fetch overything
    /// * fetch(FetchDefinition::Unmapped) or fetch("*") -> Fetch unmapped (as signified by the 'unmapped' flag in the BAM - might be unreliable with some aligners.
    ///
    /// The start / stop coordinates will take i64 (the correct type as of htslib's 'large
    /// coordinates' expansion), i32, u32, and u64 (with a possible panic! if the coordinate
    /// won't fit an i64).
    ///
    /// `start` and `stop` are zero-based. `start` is inclusive, `stop` is exclusive.
    ///
    /// This replaces the old fetch and fetch_str implementations.
    pub fn fetch<'a, T: Into<FetchDefinition<'a>>>(&mut self, fetch_definition: T) -> Result<()> {
        //this 'compile time redirect' safes us
        //from monomorphing the 'meat' of the fetch function
        self._inner_fetch(fetch_definition.into())
    }

    fn _inner_fetch(&mut self, fetch_definition: FetchDefinition) -> Result<()> {
        match fetch_definition {
            FetchDefinition::Region(tid, start, stop) => {
                self._fetch_by_coord_tuple(tid, start, stop)
            }
            FetchDefinition::RegionString(s, start, stop) => {
                let tid = self.header().tid(s);
                match tid {
                    Some(tid) => self._fetch_by_coord_tuple(tid as i32, start, stop),
                    None => Err(Error::Fetch),
                }
            }
            FetchDefinition::CompleteTid(tid) => {
                let len = self.header().target_len(tid as u32);
                match len {
                    Some(len) => self._fetch_by_coord_tuple(tid, 0, len as i64),
                    None => Err(Error::Fetch),
                }
            }
            FetchDefinition::String(s) => {
                // either a target-name or a samtools style definition
                let tid = self.header().tid(s);
                match tid {
                    Some(tid) => {
                        //'large position' spec says target len must will fit into an i64.
                        let len: i64 = self.header.target_len(tid).unwrap().try_into().unwrap();
                        self._fetch_by_coord_tuple(tid as i32, 0, len)
                    }
                    None => self._fetch_by_str(s),
                }
            }
            FetchDefinition::All => self._fetch_by_str(b"."),
            FetchDefinition::Unmapped => self._fetch_by_str(b"*"),
        }
    }

    fn _fetch_by_coord_tuple(&mut self, tid: i32, beg: i64, end: i64) -> Result<()> {
        if let Some(itr) = self.itr {
            unsafe { htslib::hts_itr_destroy(itr) }
        }
        let itr = unsafe { htslib::sam_itr_queryi(self.index().inner_ptr(), tid, beg, end) };
        if itr.is_null() {
            self.itr = None;
            Err(Error::Fetch)
        } else {
            self.itr = Some(itr);
            Ok(())
        }
    }

    fn _fetch_by_str(&mut self, region: &[u8]) -> Result<()> {
        if let Some(itr) = self.itr {
            unsafe { htslib::hts_itr_destroy(itr) }
        }
        let rstr = ffi::CString::new(region).unwrap();
        let rptr = rstr.as_ptr();
        let itr = unsafe {
            htslib::sam_itr_querys(
                self.index().inner_ptr(),
                self.header().inner_ptr() as *mut hts_sys::sam_hdr_t,
                rptr,
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

    extern "C" fn pileup_read(
        data: *mut ::std::os::raw::c_void,
        record: *mut htslib::bam1_t,
    ) -> i32 {
        let _self = unsafe { (data as *mut Self).as_mut().unwrap() };
        match _self.itr {
            Some(itr) => itr_next(_self.htsfile, itr, record), // read fetched region
            None => unsafe {
                htslib::sam_read1(
                    _self.htsfile,
                    _self.header().inner_ptr() as *mut hts_sys::sam_hdr_t,
                    record,
                )
            }, // ordinary reading
        }
    }

    /// Set the reference path for reading CRAM files.
    ///
    /// # Arguments
    ///
    /// * `path` - path to the FASTA reference
    pub fn set_reference<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        unsafe { set_fai_filename(self.htsfile, path) }
    }

    pub fn index(&self) -> &IndexView {
        self.idx.as_ref().unwrap()
    }

    // Analogous to slow_idxstats in samtools, see
    // https://github.com/samtools/samtools/blob/556c60fdff977c0e6cadc4c2581661f187098b4d/bam_index.c#L140-L199
    unsafe fn slow_idxstats(&mut self) -> Result<Vec<(i64, u64, u64, u64)>> {
        self.set_cram_options(
            hts_sys::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
            hts_sys::sam_fields_SAM_RNAME | hts_sys::sam_fields_SAM_FLAG,
        )?;
        let header = self.header();
        let h = header.inner;
        let mut ret;
        let mut last_tid = -2;
        let fp = self.htsfile();

        let nref =
            usize::try_from(hts_sys::sam_hdr_nref(h)).map_err(|_| Error::NoSequencesInReference)?;
        if nref == 0 {
            return Ok(vec![]);
        }
        let mut counts = vec![vec![0; 2]; nref + 1];
        let mut bb: hts_sys::bam1_t = MaybeUninit::zeroed().assume_init();
        let b = &mut bb as *mut hts_sys::bam1_t;
        loop {
            ret = hts_sys::sam_read1(fp, h, b);
            if ret < 0 {
                break;
            }
            let tid = (*b).core.tid;
            if tid >= nref as i32 || tid < -1 {
                return Err(Error::InvalidTid { tid });
            }

            // Map unmapped reads (tid == -1) to the last slot (nref) to avoid usize wrapping.
            let count_idx = if tid == -1 { nref } else { tid as usize };

            if tid != last_tid {
                if (last_tid >= -1) && (counts[count_idx][0] + counts[count_idx][1]) > 0 {
                    return Err(Error::BamUnsorted);
                }
                last_tid = tid;
            }

            let idx = if ((*b).core.flag as u32 & hts_sys::BAM_FUNMAP) > 0 {
                1
            } else {
                0
            };
            counts[count_idx][idx] += 1;
        }

        if ret == -1 {
            let res = (0..nref)
                .map(|i| {
                    (
                        i as i64,
                        header.target_len(i as u32).unwrap(),
                        counts[i][0],
                        counts[i][1],
                    )
                })
                .chain([(-1, 0, counts[nref][0], counts[nref][1])])
                .collect();
            Ok(res)
        } else {
            Err(Error::SlowIdxStats)
        }
    }

    /// Similar to samtools idxstats, this returns a vector of tuples
    /// containing the target id, length, number of mapped reads, and number of unmapped reads.
    /// The last entry in the vector corresponds to the unmapped reads for the entire file, with
    /// the tid set to -1.
    pub fn index_stats(&mut self) -> Result<Vec<(i64, u64, u64, u64)>> {
        let header = self.header();
        let index = self.index();
        if index.inner_ptr().is_null() {
            panic!("Index is null");
        }
        // the quick index stats method only works for BAM files, not SAM or CRAM
        unsafe {
            if (*self.htsfile()).format.format != htslib::htsExactFormat_bam {
                return self.slow_idxstats();
            }
        }
        Ok((0..header.target_count())
            .map(|tid| {
                let (mapped, unmapped) = index.number_mapped_unmapped(tid);
                let tlen = header.target_len(tid).unwrap();
                (tid as i64, tlen, mapped, unmapped)
            })
            .chain([(-1, 0, 0, index.number_unmapped())])
            .collect::<_>())
    }
}

#[derive(Debug)]
pub struct IndexView {
    inner: *mut hts_sys::hts_idx_t,
}

unsafe impl Send for IndexView {}
unsafe impl Sync for IndexView {}

impl IndexView {
    fn new(hts_idx: *mut hts_sys::hts_idx_t) -> Self {
        Self { inner: hts_idx }
    }

    #[inline]
    pub fn inner(&self) -> &hts_sys::hts_idx_t {
        unsafe { self.inner.as_ref().unwrap() }
    }

    #[inline]
    // Pointer to inner hts_idx_t struct
    pub fn inner_ptr(&self) -> *const hts_sys::hts_idx_t {
        self.inner
    }

    #[inline]
    pub fn inner_mut(&mut self) -> &mut hts_sys::hts_idx_t {
        unsafe { self.inner.as_mut().unwrap() }
    }

    #[inline]
    // Mutable pointer to hts_idx_t struct
    pub fn inner_ptr_mut(&mut self) -> *mut hts_sys::hts_idx_t {
        self.inner
    }

    /// Get the number of mapped and unmapped reads for a given target id
    /// FIXME only valid for BAM, not SAM/CRAM
    fn number_mapped_unmapped(&self, tid: u32) -> (u64, u64) {
        let (mut mapped, mut unmapped) = (0, 0);
        unsafe {
            hts_sys::hts_idx_get_stat(self.inner, tid as i32, &mut mapped, &mut unmapped);
        }
        (mapped, unmapped)
    }

    /// Get the total number of unmapped reads in the file
    /// FIXME only valid for BAM, not SAM/CRAM
    fn number_unmapped(&self) -> u64 {
        unsafe { hts_sys::hts_idx_get_n_no_coor(self.inner) }
    }
}

impl Drop for IndexView {
    fn drop(&mut self) {
        unsafe {
            htslib::hts_idx_destroy(self.inner);
        }
    }
}

impl Read for IndexedReader {
    fn read(&mut self, record: &mut record::Record) -> Option<Result<()>> {
        match self.itr {
            Some(itr) => {
                match itr_next(self.htsfile, itr, &mut record.inner as *mut htslib::bam1_t) {
                    -1 => None,
                    -2 => Some(Err(Error::BamTruncatedRecord)),
                    -4 => Some(Err(Error::BamInvalidRecord)),
                    _ => {
                        record.set_header(Arc::clone(&self.header));

                        Some(Ok(()))
                    }
                }
            }
            None => None,
        }
    }

    /// Iterator over the records of the fetched region.
    /// Note that, while being convenient, this is less efficient than pre-allocating a
    /// `Record` and reading into it with the `read` method, since every iteration involves
    /// the allocation of a new `Record`.
    fn records(&mut self) -> Records<'_, Self> {
        Records { reader: self }
    }

    fn rc_records(&mut self) -> RcRecords<'_, Self> {
        RcRecords {
            reader: self,
            record: Rc::new(record::Record::new()),
        }
    }

    fn pileup(&mut self) -> pileup::Pileups<'_, Self> {
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

    fn set_thread_pool(&mut self, tpool: &ThreadPool) -> Result<()> {
        unsafe { set_thread_pool(self.htsfile(), tpool)? }
        self.tpool = Some(tpool.clone());
        Ok(())
    }
}

impl Drop for IndexedReader {
    fn drop(&mut self) {
        unsafe {
            if let Some(itr) = self.itr.take() {
                htslib::hts_itr_destroy(itr);
            }

            // A CRAM index contains a pointer to the CRAM file handle.
            // Destroy the index before hts_close frees that handle.
            drop(self.idx.take());
            htslib::hts_close(self.htsfile);
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Sam,
    Bam,
    Cram,
}

impl Format {
    fn write_mode(self) -> &'static [u8] {
        match self {
            Format::Sam => b"w",
            Format::Bam => b"wb",
            Format::Cram => b"wc",
        }
    }
}

/// A BAM writer.
#[derive(Debug)]
pub struct Writer {
    f: *mut htslib::htsFile,
    header: Arc<HeaderView>,
    tpool: Option<ThreadPool>,
}

unsafe impl Send for Writer {}

impl Writer {
    /// Create a new SAM/BAM/CRAM file.
    ///
    /// # Arguments
    ///
    /// * `path` - the path.
    /// * `header` - header definition to use
    /// * `format` - the format to use (SAM/BAM/CRAM)
    pub fn from_path<P: AsRef<Path>>(
        path: P,
        header: &header::Header,
        format: Format,
    ) -> Result<Self> {
        Self::new(&path_as_bytes(path, false)?, format.write_mode(), header)
    }

    /// Create a new SAM/BAM/CRAM file at STDOUT.
    ///
    /// # Arguments
    ///
    /// * `header` - header definition to use
    /// * `format` - the format to use (SAM/BAM/CRAM)
    pub fn from_stdout(header: &header::Header, format: Format) -> Result<Self> {
        Self::new(b"-", format.write_mode(), header)
    }

    /// Create a new SAM/BAM/CRAM file.
    ///
    /// # Arguments
    ///
    /// * `path` - the path. Use "-" for stdout.
    /// * `mode` - write mode, refer to htslib::hts_open()
    /// * `header` - header definition to use
    fn new(path: &[u8], mode: &[u8], header: &header::Header) -> Result<Self> {
        let f = hts_open(path, mode)?;

        // sam_hdr_parse does not populate the text and l_text fields of the header_record.
        // This causes non-SQ headers to be dropped in the output BAM file.
        // To avoid this, we copy the All header to a new C-string that is allocated with malloc,
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
            let rec = htslib::sam_hdr_parse(l_text + 1, text as *const c_char);

            (*rec).text = text as *mut c_char;
            (*rec).l_text = l_text;
            rec
        };

        unsafe {
            htslib::sam_hdr_write(f, header_record);
        }

        Ok(Writer {
            f,
            header: Arc::new(HeaderView::new(header_record)),
            tpool: None,
        })
    }

    /// Activate multi-threaded BAM write support in htslib. This should permit faster
    /// writing of large BAM files.
    ///
    /// # Arguments
    ///
    /// * `n_threads` - number of extra background writer threads to use, must be `> 0`.
    pub fn set_threads(&mut self, n_threads: usize) -> Result<()> {
        unsafe { set_threads(self.f, n_threads) }
    }

    /// Use a shared thread-pool for writing. This permits controlling the total
    /// thread count when multiple readers and writers are working simultaneously.
    /// A thread pool can be created with `crate::tpool::ThreadPool::new(n_threads)`
    ///
    /// # Arguments
    ///
    /// * `tpool` - thread pool to use for compression work.
    pub fn set_thread_pool(&mut self, tpool: &ThreadPool) -> Result<()> {
        unsafe { set_thread_pool(self.f, tpool)? }
        self.tpool = Some(tpool.clone());
        Ok(())
    }

    /// Write record to BAM.
    ///
    /// # Arguments
    ///
    /// * `record` - the record to write
    pub fn write(&mut self, record: &record::Record) -> Result<()> {
        if unsafe { htslib::sam_write1(self.f, self.header.inner(), record.inner_ptr()) } == -1 {
            Err(Error::WriteRecord)
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
    pub fn set_reference<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        unsafe { set_fai_filename(self.f, path) }
    }

    /// Set the compression level for writing BAM/CRAM files.
    ///
    /// # Arguments
    ///
    /// * `compression_level` - `CompressionLevel` enum variant
    pub fn set_compression_level(&mut self, compression_level: CompressionLevel) -> Result<()> {
        let level = compression_level.convert()?;
        match unsafe {
            htslib::hts_set_opt(
                self.f,
                htslib::hts_fmt_option_HTS_OPT_COMPRESSION_LEVEL,
                level,
            )
        } {
            0 => Ok(()),
            _ => Err(Error::BamInvalidCompressionLevel { level }),
        }
    }
}

/// Compression levels in BAM/CRAM files
///
/// * Uncompressed: No compression, zlib level 0
/// * Fastest: Lowest compression level, zlib level 1
/// * Maximum: Highest compression level, zlib level 9
/// * Level(i): Custom compression level in the range [0, 9]
#[derive(Debug, Clone, Copy)]
pub enum CompressionLevel {
    Uncompressed,
    Fastest,
    Maximum,
    Level(u32),
}

impl CompressionLevel {
    // Convert and check the variants of the `CompressionLevel` enum to a numeric level
    fn convert(self) -> Result<u32> {
        match self {
            CompressionLevel::Uncompressed => Ok(0),
            CompressionLevel::Fastest => Ok(1),
            CompressionLevel::Maximum => Ok(9),
            CompressionLevel::Level(i @ 0..=9) => Ok(i),
            CompressionLevel::Level(i) => Err(Error::BamInvalidCompressionLevel { level: i }),
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
pub struct Records<'a, R: Read> {
    reader: &'a mut R,
}

impl<R: Read> Iterator for Records<'_, R> {
    type Item = Result<record::Record>;

    fn next(&mut self) -> Option<Result<record::Record>> {
        let mut record = record::Record::new();
        match self.reader.read(&mut record) {
            None => None,
            Some(Ok(_)) => Some(Ok(record)),
            Some(Err(err)) => Some(Err(err)),
        }
    }
}

/// Iterator over the records of a BAM, using an Rc.
///
/// See [rc_records](trait.Read.html#tymethod.rc_records).
#[derive(Debug)]
pub struct RcRecords<'a, R: Read> {
    reader: &'a mut R,
    record: Rc<record::Record>,
}

impl<R: Read> Iterator for RcRecords<'_, R> {
    type Item = Result<Rc<record::Record>>;

    fn next(&mut self) -> Option<Self::Item> {
        let record = match Rc::get_mut(&mut self.record) {
            //not make_mut, we don't need a clone
            Some(x) => x,
            None => {
                self.record = Rc::new(record::Record::new());
                Rc::get_mut(&mut self.record).unwrap()
            }
        };

        match self.reader.read(record) {
            None => None,
            Some(Ok(_)) => Some(Ok(Rc::clone(&self.record))),
            Some(Err(err)) => Some(Err(err)),
        }
    }
}

/// Iterator over the records of a BAM until the virtual offset is less than `end`
pub struct ChunkIterator<'a, R: Read> {
    reader: &'a mut R,
    end: Option<i64>,
}

impl<R: Read> Iterator for ChunkIterator<'_, R> {
    type Item = Result<record::Record>;
    fn next(&mut self) -> Option<Result<record::Record>> {
        if let Some(pos) = self.end {
            if self.reader.tell() >= pos {
                return None;
            }
        }
        let mut record = record::Record::new();
        match self.reader.read(&mut record) {
            None => None,
            Some(Ok(_)) => Some(Ok(record)),
            Some(Err(err)) => Some(Err(err)),
        }
    }
}

/// Wrapper for opening a BAM file.
fn hts_open(path: &[u8], mode: &[u8]) -> Result<*mut htslib::htsFile> {
    let cpath = ffi::CString::new(path).unwrap();
    let path = str::from_utf8(path).unwrap();
    let c_str = ffi::CString::new(mode).unwrap();
    let ret = unsafe { htslib::hts_open(cpath.as_ptr(), c_str.as_ptr()) };
    if ret.is_null() {
        Err(Error::BamOpen {
            target: path.to_owned(),
        })
    } else {
        if !mode.contains(&b'w') {
            unsafe {
                // Comparison against 'htsFormatCategory_sequence_data' doesn't handle text files correctly
                // hence the explicit checks against all supported exact formats
                if (*ret).format.format != htslib::htsExactFormat_sam
                    && (*ret).format.format != htslib::htsExactFormat_bam
                    && (*ret).format.format != htslib::htsExactFormat_cram
                {
                    return Err(Error::BamOpen {
                        target: path.to_owned(),
                    });
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
}

unsafe impl Send for HeaderView {}
unsafe impl Sync for HeaderView {}

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

            let rec = htslib::sam_hdr_parse(l_text + 1, text as *const c_char);
            (*rec).text = text as *mut c_char;
            (*rec).l_text = l_text;
            rec
        };

        HeaderView::new(header_record)
    }

    /// Create a new HeaderView from the underlying Htslib type, and own it.
    fn new(inner: *mut htslib::bam_hdr_t) -> Self {
        HeaderView { inner }
    }

    #[inline]
    pub fn inner(&self) -> &htslib::bam_hdr_t {
        unsafe { self.inner.as_ref().unwrap() }
    }

    #[inline]
    // Pointer to inner bam_hdr_t struct
    pub fn inner_ptr(&self) -> *const htslib::bam_hdr_t {
        self.inner
    }

    #[inline]
    pub fn inner_mut(&mut self) -> &mut htslib::bam_hdr_t {
        unsafe { self.inner.as_mut().unwrap() }
    }

    #[inline]
    // Mutable pointer to bam_hdr_t struct
    pub fn inner_ptr_mut(&mut self) -> *mut htslib::bam_hdr_t {
        self.inner
    }

    pub fn tid(&self, name: &[u8]) -> Option<u32> {
        let c_str = ffi::CString::new(name).expect("Expected valid name.");
        let tid = unsafe { htslib::sam_hdr_name2tid(self.inner, c_str.as_ptr()) };
        if tid < 0 {
            None
        } else {
            Some(tid as u32)
        }
    }

    pub fn tid2name(&self, tid: u32) -> &[u8] {
        let ptr = unsafe { htslib::sam_hdr_tid2name(self.inner, tid as i32) };
        if ptr.is_null() {
            b""
        } else {
            unsafe { ffi::CStr::from_ptr(ptr).to_bytes() }
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

    pub fn target_len(&self, tid: u32) -> Option<u64> {
        let inner = unsafe { *self.inner };
        if (tid as i32) < inner.n_targets {
            let l: &[u32] =
                unsafe { slice::from_raw_parts(inner.target_len, inner.n_targets as usize) };
            Some(l[tid as usize] as u64)
        } else {
            None
        }
    }

    /// Retrieve the textual SAM header as bytes
    pub fn as_bytes(&self) -> &[u8] {
        unsafe {
            let rebuilt_hdr = htslib::sam_hdr_str(self.inner);
            if rebuilt_hdr.is_null() {
                return b"";
            }
            ffi::CStr::from_ptr(rebuilt_hdr).to_bytes()
        }
    }
}

impl Clone for HeaderView {
    fn clone(&self) -> Self {
        HeaderView {
            inner: unsafe { htslib::sam_hdr_dup(self.inner) },
        }
    }
}

impl Drop for HeaderView {
    fn drop(&mut self) {
        unsafe {
            htslib::sam_hdr_destroy(self.inner);
        }
    }
}

#[cfg(test)]
mod mod_test;
