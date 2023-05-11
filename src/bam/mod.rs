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
    /// let mut bam = IndexedReader::from_path(&"test/test.bam").unwrap();
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
    /// let mut bam = Reader::from_path(&"test/test.bam").unwrap();
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
    /// let mut cram = Reader::from_path(&"test/test_cram.cram").unwrap();
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
    header: Rc<HeaderView>,
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
            header: Rc::new(HeaderView::new(header)),
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
    /// a valid virtual offset.
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
    /// let mut bam = Reader::from_path(&"test/test.bam")?;
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
                record.set_header(Rc::clone(&self.header));

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
    header: Rc<HeaderView>,
    idx: Rc<IndexView>,
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
                header: Rc::new(HeaderView::new(header)),
                idx: Rc::new(IndexView::new(idx)),
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
                header: Rc::new(HeaderView::new(header)),
                idx: Rc::new(IndexView::new(idx)),
                itr: None,
                tpool: None,
            })
        }
    }

    /// Define the region from which .read() or .records will retrieve reads.
    ///
    /// Both iterating (with [.records()](trait.Read.html#tymethod.records)) and looping without allocation (with [.read()](trait.Read.html#tymethod.read) are a two stage process:
    /// 1. 'fetch' the region of interest
    /// 2. iter/loop trough the reads.
    ///
    /// Example:
    /// ```
    /// use rust_htslib::bam::{IndexedReader, Read};
    /// let mut bam = IndexedReader::from_path(&"test/test.bam").unwrap();
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
        &self.idx
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

            if tid != last_tid {
                if (last_tid >= -1) && (counts[tid as usize][0] + counts[tid as usize][1]) > 0 {
                    return Err(Error::BamUnsorted);
                }
                last_tid = tid;
            }

            let idx = if ((*b).core.flag as u32 & hts_sys::BAM_FUNMAP) > 0 {
                1
            } else {
                0
            };
            counts[(*b).core.tid as usize][idx] += 1;
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
    owned: bool,
}

impl IndexView {
    fn new(hts_idx: *mut hts_sys::hts_idx_t) -> Self {
        Self {
            inner: hts_idx,
            owned: true,
        }
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
        if self.owned {
            unsafe {
                htslib::hts_idx_destroy(self.inner);
            }
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
                        record.set_header(Rc::clone(&self.header));

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
            if self.itr.is_some() {
                htslib::hts_itr_destroy(self.itr.unwrap());
            }
            htslib::hts_close(self.htsfile);
        }
    }
}

#[derive(Debug, Clone, Copy)]
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
    header: Rc<HeaderView>,
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
            let rec = htslib::sam_hdr_parse(
                ((l_text + 1) as usize).try_into().unwrap(),
                text as *const c_char,
            );

            (*rec).text = text as *mut c_char;
            (*rec).l_text = l_text as u64;
            rec
        };

        unsafe {
            htslib::sam_hdr_write(f, header_record);
        }

        Ok(Writer {
            f,
            header: Rc::new(HeaderView::new(header_record)),
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

impl<'a, R: Read> Iterator for Records<'a, R> {
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

impl<'a, R: Read> Iterator for RcRecords<'a, R> {
    type Item = Result<Rc<record::Record>>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = match Rc::get_mut(&mut self.record) {
            //not make_mut, we don't need a clone
            Some(x) => x,
            None => {
                self.record = Rc::new(record::Record::new());
                Rc::get_mut(&mut self.record).unwrap()
            }
        };

        match self.reader.read(&mut record) {
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

impl<'a, R: Read> Iterator for ChunkIterator<'a, R> {
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

            let rec = htslib::sam_hdr_parse((l_text + 1) as u64, text as *const c_char);
            (*rec).text = text as *mut c_char;
            (*rec).l_text = l_text as u64;
            rec
        };

        HeaderView::new(header_record)
    }

    /// Create a new HeaderView from the underlying Htslib type, and own it.
    pub fn new(inner: *mut htslib::bam_hdr_t) -> Self {
        HeaderView { inner, owned: true }
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
        unsafe { ffi::CStr::from_ptr(htslib::sam_hdr_tid2name(self.inner, tid as i32)).to_bytes() }
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
            owned: true,
        }
    }
}

impl Drop for HeaderView {
    fn drop(&mut self) {
        if self.owned {
            unsafe {
                htslib::sam_hdr_destroy(self.inner);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::header::HeaderRecord;
    use super::record::{Aux, Cigar, CigarString};
    use super::*;
    use std::collections::HashMap;
    use std::fs;
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

    fn compare_inner_bam_cram_records(cram_records: &[Record], bam_records: &[Record]) {
        // Selectively compares bam1_t struct fields from BAM and CRAM
        for (c1, b1) in cram_records.iter().zip(bam_records.iter()) {
            // CRAM vs BAM l_data is off by 3, see: https://github.com/rust-bio/rust-htslib/pull/184#issuecomment-590133544
            // The rest of the fields should be identical:
            assert_eq!(c1.cigar(), b1.cigar());
            assert_eq!(c1.inner().core.pos, b1.inner().core.pos);
            assert_eq!(c1.inner().core.mpos, b1.inner().core.mpos);
            assert_eq!(c1.inner().core.mtid, b1.inner().core.mtid);
            assert_eq!(c1.inner().core.tid, b1.inner().core.tid);
            assert_eq!(c1.inner().core.bin, b1.inner().core.bin);
            assert_eq!(c1.inner().core.qual, b1.inner().core.qual);
            assert_eq!(c1.inner().core.l_extranul, b1.inner().core.l_extranul);
            assert_eq!(c1.inner().core.flag, b1.inner().core.flag);
            assert_eq!(c1.inner().core.l_qname, b1.inner().core.l_qname);
            assert_eq!(c1.inner().core.n_cigar, b1.inner().core.n_cigar);
            assert_eq!(c1.inner().core.l_qseq, b1.inner().core.l_qseq);
            assert_eq!(c1.inner().core.isize, b1.inner().core.isize);
            //... except m_data
        }
    }

    #[test]
    fn test_read() {
        let (names, flags, seqs, quals, cigars) = gold();
        let mut bam = Reader::from_path(&Path::new("test/test.bam")).expect("Error opening file.");
        let del_len = [1, 1, 1, 1, 1, 100000];

        for (i, record) in bam.records().enumerate() {
            let rec = record.expect("Expected valid record");
            assert_eq!(rec.qname(), names[i]);
            assert_eq!(rec.flags(), flags[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);

            let cigar = rec.cigar();
            assert_eq!(*cigar, cigars[i]);

            let end_pos = cigar.end_pos();
            assert_eq!(end_pos, rec.pos() + 100 + del_len[i]);
            assert_eq!(
                cigar
                    .read_pos(end_pos as u32 - 10, false, false)
                    .unwrap()
                    .unwrap(),
                90
            );
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
        let mut bam = Reader::from_path(&Path::new("test/test.bam")).expect("Error opening file.");

        let mut names_by_voffset = HashMap::new();

        let mut offset = bam.tell();
        let mut rec = Record::new();
        loop {
            match bam.read(&mut rec) {
                Some(r) => r.expect("error reading bam"),
                None => break,
            };
            let qname = str::from_utf8(rec.qname()).unwrap().to_string();
            println!("{} {}", offset, qname);
            names_by_voffset.insert(offset, qname);
            offset = bam.tell();
        }

        for (offset, qname) in names_by_voffset.iter() {
            println!("{} {}", offset, qname);
            bam.seek(*offset).unwrap();
            match bam.read(&mut rec) {
                Some(r) => r.unwrap(),
                None => {}
            };
            let rec_qname = str::from_utf8(rec.qname()).unwrap().to_string();
            assert_eq!(qname, &rec_qname);
        }
    }

    #[test]
    fn test_read_sam_header() {
        let bam = Reader::from_path(&"test/test.bam").expect("Error opening file.");

        let true_header = "@SQ\tSN:CHROMOSOME_I\tLN:15072423\n@SQ\tSN:CHROMOSOME_II\tLN:15279345\
             \n@SQ\tSN:CHROMOSOME_III\tLN:13783700\n@SQ\tSN:CHROMOSOME_IV\tLN:17493793\n@SQ\t\
             SN:CHROMOSOME_V\tLN:20924149\n"
            .to_string();
        let header_text = String::from_utf8(bam.header.as_bytes().to_owned()).unwrap();
        assert_eq!(header_text, true_header);
    }

    #[test]
    fn test_read_against_sam() {
        let mut bam = Reader::from_path("./test/bam2sam_out.sam").unwrap();
        for read in bam.records() {
            let _read = read.unwrap();
        }
    }

    fn _test_read_indexed_common(mut bam: IndexedReader) {
        let (names, flags, seqs, quals, cigars) = gold();
        let sq_1 = b"CHROMOSOME_I";
        let sq_2 = b"CHROMOSOME_II";
        let tid_1 = bam.header.tid(sq_1).expect("Expected tid.");
        let tid_2 = bam.header.tid(sq_2).expect("Expected tid.");
        assert!(bam.header.target_len(tid_1).expect("Expected target len.") == 15072423);

        // fetch to position containing reads
        bam.fetch((tid_1, 0, 2))
            .expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);

        // compare reads
        bam.fetch((tid_1, 0, 2))
            .expect("Expected successful fetch.");
        for (i, record) in bam.records().enumerate() {
            let rec = record.expect("Expected valid record");

            println!("{}", str::from_utf8(rec.qname()).unwrap());
            assert_eq!(rec.qname(), names[i]);
            assert_eq!(rec.flags(), flags[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);
            assert_eq!(*rec.cigar(), cigars[i]);
            // fix qual offset
            let qual: Vec<u8> = quals[i].iter().map(|&q| q - 33).collect();
            assert_eq!(rec.qual(), &qual[..]);
            assert_eq!(rec.aux(b"NotAvailableAux"), Err(Error::BamAuxTagNotFound));
        }

        // fetch to empty position
        bam.fetch((tid_2, 1, 1))
            .expect("Expected successful fetch.");
        assert!(bam.records().count() == 0);

        // repeat with byte-string based fetch

        // fetch to position containing reads
        // using coordinate-string chr:start-stop
        bam.fetch(format!("{}:{}-{}", str::from_utf8(sq_1).unwrap(), 0, 2).as_bytes())
            .expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);
        // using &str and exercising some of the coordinate conversion funcs
        bam.fetch((str::from_utf8(sq_1).unwrap(), 0 as u32, 2 as u64))
            .expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);
        // using a slice
        bam.fetch((&sq_1[..], 0, 2))
            .expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);
        // using a literal
        bam.fetch((sq_1, 0, 2)).expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);

        // using a tid
        bam.fetch((0i32, 0u32, 2i64))
            .expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);
        // using a tid:u32
        bam.fetch((0u32, 0u32, 2i64))
            .expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);

        // compare reads
        bam.fetch(format!("{}:{}-{}", str::from_utf8(sq_1).unwrap(), 0, 2).as_bytes())
            .expect("Expected successful fetch.");
        for (i, record) in bam.records().enumerate() {
            let rec = record.expect("Expected valid record");

            println!("{}", str::from_utf8(rec.qname()).unwrap());
            assert_eq!(rec.qname(), names[i]);
            assert_eq!(rec.flags(), flags[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);
            assert_eq!(*rec.cigar(), cigars[i]);
            // fix qual offset
            let qual: Vec<u8> = quals[i].iter().map(|&q| q - 33).collect();
            assert_eq!(rec.qual(), &qual[..]);
            assert_eq!(rec.aux(b"NotAvailableAux"), Err(Error::BamAuxTagNotFound));
        }

        // fetch to empty position
        bam.fetch(format!("{}:{}-{}", str::from_utf8(sq_2).unwrap(), 1, 1).as_bytes())
            .expect("Expected successful fetch.");
        assert!(bam.records().count() == 0);

        //all on a tid
        bam.fetch(0).expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);
        //all on a tid:u32
        bam.fetch(0u32).expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);

        //all on a tid - by &[u8]
        bam.fetch(sq_1).expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);
        //all on a tid - by str
        bam.fetch(str::from_utf8(sq_1).unwrap())
            .expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);

        //all reads
        bam.fetch(FetchDefinition::All)
            .expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);

        //all reads
        bam.fetch(".").expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);

        //all unmapped
        bam.fetch(FetchDefinition::Unmapped)
            .expect("Expected successful fetch.");
        assert_eq!(bam.records().count(), 1); // expect one 'truncade record' Record.

        bam.fetch("*").expect("Expected successful fetch.");
        assert_eq!(bam.records().count(), 1); // expect one 'truncade record' Record.
    }

    #[test]
    fn test_read_indexed() {
        let bam = IndexedReader::from_path(&"test/test.bam").expect("Expected valid index.");
        _test_read_indexed_common(bam);
    }

    #[test]
    fn test_read_indexed_different_index_name() {
        let bam = IndexedReader::from_path_and_index(
            &"test/test_different_index_name.bam",
            &"test/test.bam.bai",
        )
        .expect("Expected valid index.");
        _test_read_indexed_common(bam);
    }

    #[test]
    fn test_set_record() {
        let (names, _, seqs, quals, cigars) = gold();

        let mut rec = record::Record::new();
        rec.set_reverse();
        rec.set(names[0], Some(&cigars[0]), seqs[0], quals[0]);
        // note: this segfaults if you push_aux() before set()
        //       because set() obliterates aux
        rec.push_aux(b"NM", Aux::I32(15)).unwrap();

        assert_eq!(rec.qname(), names[0]);
        assert_eq!(*rec.cigar(), cigars[0]);
        assert_eq!(rec.seq().as_bytes(), seqs[0]);
        assert_eq!(rec.qual(), quals[0]);
        assert!(rec.is_reverse());
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));
    }

    #[test]
    fn test_set_repeated() {
        let mut rec = Record::new();
        rec.set(
            b"123",
            Some(&CigarString(vec![Cigar::Match(3)])),
            b"AAA",
            b"III",
        );
        rec.push_aux(b"AS", Aux::I32(12345)).unwrap();
        assert_eq!(rec.qname(), b"123");
        assert_eq!(rec.seq().as_bytes(), b"AAA");
        assert_eq!(rec.qual(), b"III");
        assert_eq!(rec.aux(b"AS").unwrap(), Aux::I32(12345));

        rec.set(
            b"1234",
            Some(&CigarString(vec![Cigar::SoftClip(1), Cigar::Match(3)])),
            b"AAAA",
            b"IIII",
        );
        assert_eq!(rec.qname(), b"1234");
        assert_eq!(rec.seq().as_bytes(), b"AAAA");
        assert_eq!(rec.qual(), b"IIII");
        assert_eq!(rec.aux(b"AS").unwrap(), Aux::I32(12345));

        rec.set(
            b"12",
            Some(&CigarString(vec![Cigar::Match(2)])),
            b"AA",
            b"II",
        );
        assert_eq!(rec.qname(), b"12");
        assert_eq!(rec.seq().as_bytes(), b"AA");
        assert_eq!(rec.qual(), b"II");
        assert_eq!(rec.aux(b"AS").unwrap(), Aux::I32(12345));
    }

    #[test]
    fn test_set_qname() {
        let (names, _, seqs, quals, cigars) = gold();

        assert!(names[0] != names[1]);

        for i in 0..names.len() {
            let mut rec = record::Record::new();
            rec.set(names[i], Some(&cigars[i]), seqs[i], quals[i]);
            rec.push_aux(b"NM", Aux::I32(15)).unwrap();

            assert_eq!(rec.qname(), names[i]);
            assert_eq!(*rec.cigar(), cigars[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);
            assert_eq!(rec.qual(), quals[i]);
            assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));

            // Equal length qname
            assert!(rec.qname()[0] != b'X');
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
            assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));

            // Shorter qname
            let shorter_name = b"42";
            rec.set_qname(shorter_name);

            assert_eq!(rec.qname(), shorter_name);
            assert_eq!(*rec.cigar(), cigars[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);
            assert_eq!(rec.qual(), quals[i]);
            assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));

            // Zero-length qname
            rec.set_qname(b"");

            assert_eq!(rec.qname(), b"");
            assert_eq!(*rec.cigar(), cigars[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);
            assert_eq!(rec.qual(), quals[i]);
            assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));
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
        let mut bam = Reader::from_path(&Path::new("test/test.bam")).expect("Error opening file.");

        for record in bam.records() {
            let mut rec = record.expect("Expected valid record");

            if rec.aux(b"XS").is_ok() {
                rec.remove_aux(b"XS").unwrap();
            }

            if rec.aux(b"YT").is_ok() {
                rec.remove_aux(b"YT").unwrap();
            }

            assert!(rec.remove_aux(b"ab").is_err());

            assert_eq!(rec.aux(b"XS"), Err(Error::BamAuxTagNotFound));
            assert_eq!(rec.aux(b"YT"), Err(Error::BamAuxTagNotFound));
        }
    }

    #[test]
    fn test_write() {
        let (names, _, seqs, quals, cigars) = gold();

        let tmp = tempfile::Builder::new()
            .prefix("rust-htslib")
            .tempdir()
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
                Format::Bam,
            )
            .expect("Error opening file.");

            for i in 0..names.len() {
                let mut rec = record::Record::new();
                rec.set(names[i], Some(&cigars[i]), seqs[i], quals[i]);
                rec.push_aux(b"NM", Aux::I32(15)).unwrap();

                bam.write(&rec).expect("Failed to write record.");
            }
        }

        {
            let mut bam = Reader::from_path(&bampath).expect("Error opening file.");

            for i in 0..names.len() {
                let mut rec = record::Record::new();
                match bam.read(&mut rec) {
                    Some(r) => r.expect("Failed to read record."),
                    None => {}
                };

                assert_eq!(rec.qname(), names[i]);
                assert_eq!(*rec.cigar(), cigars[i]);
                assert_eq!(rec.seq().as_bytes(), seqs[i]);
                assert_eq!(rec.qual(), quals[i]);
                assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));
            }
        }

        tmp.close().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_write_threaded() {
        let (names, _, seqs, quals, cigars) = gold();

        let tmp = tempfile::Builder::new()
            .prefix("rust-htslib")
            .tempdir()
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
                Format::Bam,
            )
            .expect("Error opening file.");
            bam.set_threads(4).unwrap();

            for i in 0..10000 {
                let mut rec = record::Record::new();
                let idx = i % names.len();
                rec.set(names[idx], Some(&cigars[idx]), seqs[idx], quals[idx]);
                rec.push_aux(b"NM", Aux::I32(15)).unwrap();
                rec.set_pos(i as i64);

                bam.write(&rec).expect("Failed to write record.");
            }
        }

        {
            let mut bam = Reader::from_path(&bampath).expect("Error opening file.");

            for (i, _rec) in bam.records().enumerate() {
                let idx = i % names.len();

                let rec = _rec.expect("Failed to read record.");

                assert_eq!(rec.pos(), i as i64);
                assert_eq!(rec.qname(), names[idx]);
                assert_eq!(*rec.cigar(), cigars[idx]);
                assert_eq!(rec.seq().as_bytes(), seqs[idx]);
                assert_eq!(rec.qual(), quals[idx]);
                assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));
            }
        }

        tmp.close().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_write_shared_tpool() {
        let (names, _, seqs, quals, cigars) = gold();

        let tmp = tempfile::Builder::new()
            .prefix("rust-htslib")
            .tempdir()
            .expect("Cannot create temp dir");
        let bampath1 = tmp.path().join("test1.bam");
        let bampath2 = tmp.path().join("test2.bam");

        {
            let (mut bam1, mut bam2) = {
                let pool = crate::tpool::ThreadPool::new(4).unwrap();

                let mut bam1 = Writer::from_path(
                    &bampath1,
                    Header::new().push_record(
                        HeaderRecord::new(b"SQ")
                            .push_tag(b"SN", &"chr1")
                            .push_tag(b"LN", &15072423),
                    ),
                    Format::Bam,
                )
                .expect("Error opening file.");

                let mut bam2 = Writer::from_path(
                    &bampath2,
                    Header::new().push_record(
                        HeaderRecord::new(b"SQ")
                            .push_tag(b"SN", &"chr1")
                            .push_tag(b"LN", &15072423),
                    ),
                    Format::Bam,
                )
                .expect("Error opening file.");

                bam1.set_thread_pool(&pool).unwrap();
                bam2.set_thread_pool(&pool).unwrap();
                (bam1, bam2)
            };

            for i in 0..10000 {
                let mut rec = record::Record::new();
                let idx = i % names.len();
                rec.set(names[idx], Some(&cigars[idx]), seqs[idx], quals[idx]);
                rec.push_aux(b"NM", Aux::I32(15)).unwrap();
                rec.set_pos(i as i64);

                bam1.write(&rec).expect("Failed to write record.");
                bam2.write(&rec).expect("Failed to write record.");
            }
        }

        {
            let pool = crate::tpool::ThreadPool::new(2).unwrap();

            for p in vec![bampath1, bampath2] {
                let mut bam = Reader::from_path(&p).expect("Error opening file.");
                bam.set_thread_pool(&pool).unwrap();

                for (i, _rec) in bam.iter_chunk(None, None).enumerate() {
                    let idx = i % names.len();

                    let rec = _rec.expect("Failed to read record.");

                    assert_eq!(rec.pos(), i as i64);
                    assert_eq!(rec.qname(), names[idx]);
                    assert_eq!(*rec.cigar(), cigars[idx]);
                    assert_eq!(rec.seq().as_bytes(), seqs[idx]);
                    assert_eq!(rec.qual(), quals[idx]);
                    assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));
                }
            }
        }

        tmp.close().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_copy_template() {
        // Verify that BAM headers are transmitted correctly when using an existing BAM as a
        // template for headers.

        let tmp = tempfile::Builder::new()
            .prefix("rust-htslib")
            .tempdir()
            .expect("Cannot create temp dir");
        let bampath = tmp.path().join("test.bam");
        println!("{:?}", bampath);

        let mut input_bam = Reader::from_path(&"test/test.bam").expect("Error opening file.");

        {
            let mut bam = Writer::from_path(
                &bampath,
                &Header::from_template(&input_bam.header()),
                Format::Bam,
            )
            .expect("Error opening file.");

            for rec in input_bam.records() {
                bam.write(&rec.unwrap()).expect("Failed to write record.");
            }
        }

        {
            let copy_bam = Reader::from_path(&bampath).expect("Error opening file.");

            // Verify that the header came across correctly
            assert_eq!(input_bam.header().as_bytes(), copy_bam.header().as_bytes());
        }

        tmp.close().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_pileup() {
        let (_, _, seqs, quals, _) = gold();

        let mut bam = Reader::from_path(&"test/test.bam").expect("Error opening file.");
        let pileups = bam.pileup();
        for pileup in pileups.take(26) {
            let _pileup = pileup.expect("Expected successful pileup.");
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
        let mut bam = IndexedReader::from_path(&"test/test.bam").expect("Error opening file.");
        // read without fetch
        for pileup in bam.pileup() {
            pileup.unwrap();
        }
        // go back again
        let tid = bam.header().tid(b"CHROMOSOME_I").unwrap();
        bam.fetch((tid, 0, 5)).unwrap();
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
            .filter(|x| !x.is_empty() && x[0] != b'@')
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
        let mut bam = Reader::from_path(&Path::new("test/test.bam")).expect("Error opening file.");

        for (i, record) in bam.records().enumerate() {
            let rec = record.expect("Expected valid record");

            let cigar = rec.cigar();
            assert_eq!(*cigar, cigars[i]);
        }

        for (i, record) in bam.records().enumerate() {
            let mut rec = record.expect("Expected valid record");
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

        compare_inner_bam_cram_records(&cram_records, &bam_records);
    }

    #[test]
    fn test_write_cram() {
        // BAM file, records
        let bam_path = "./test/test_cram.bam";
        let ref_path = "./test/test_cram.fa";
        let mut bam_reader = Reader::from_path(bam_path).unwrap();
        let bam_records: Vec<Record> = bam_reader.records().map(|v| v.unwrap()).collect();

        // New CRAM file
        let tmp = tempfile::Builder::new()
            .prefix("rust-htslib")
            .tempdir()
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

            let mut cram_writer = Writer::from_path(&cram_path, &header, Format::Cram)
                .expect("Error opening CRAM file.");
            cram_writer.set_reference(ref_path).unwrap();

            // Write BAM records to CRAM file
            for rec in bam_records.iter() {
                cram_writer
                    .write(&rec)
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
            compare_inner_bam_cram_records(&cram_records, &bam_records);
        }

        tmp.close().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_compression_level_conversion() {
        // predefined compression levels
        assert_eq!(CompressionLevel::Uncompressed.convert().unwrap(), 0);
        assert_eq!(CompressionLevel::Fastest.convert().unwrap(), 1);
        assert_eq!(CompressionLevel::Maximum.convert().unwrap(), 9);

        // numeric compression levels
        for level in 0..=9 {
            assert_eq!(CompressionLevel::Level(level).convert().unwrap(), level);
        }
        // invalid levels
        assert!(CompressionLevel::Level(10).convert().is_err());
    }

    #[test]
    fn test_write_compression() {
        let tmp = tempfile::Builder::new()
            .prefix("rust-htslib")
            .tempdir()
            .expect("Cannot create temp dir");
        let input_bam_path = "test/test.bam";

        // test levels with decreasing compression factor
        let levels_to_test = vec![
            CompressionLevel::Maximum,
            CompressionLevel::Level(6),
            CompressionLevel::Fastest,
            CompressionLevel::Uncompressed,
        ];
        let file_sizes: Vec<_> = levels_to_test
            .iter()
            .map(|level| {
                let output_bam_path = tmp.path().join("test.bam");
                {
                    let mut reader = Reader::from_path(&input_bam_path).unwrap();
                    let header = Header::from_template(reader.header());
                    let mut writer =
                        Writer::from_path(&output_bam_path, &header, Format::Bam).unwrap();
                    writer.set_compression_level(*level).unwrap();
                    for record in reader.records() {
                        let r = record.unwrap();
                        writer.write(&r).unwrap();
                    }
                }
                fs::metadata(output_bam_path).unwrap().len()
            })
            .collect();

        // check that out BAM file sizes are in decreasing order, in line with the expected compression factor
        println!("testing compression leves: {:?}", levels_to_test);
        println!("got compressed sizes: {:?}", file_sizes);

        // libdeflate comes out with a slightly bigger file on Max compression
        // than on Level(6), so skip that check
        #[cfg(feature = "libdeflate")]
        assert!(file_sizes[1..].windows(2).all(|size| size[0] <= size[1]));

        #[cfg(not(feature = "libdeflate"))]
        assert!(file_sizes.windows(2).all(|size| size[0] <= size[1]));

        tmp.close().expect("Failed to delete temp dir");
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
        let bam_path = "./Cargo.toml";
        let bam_reader = Reader::from_path(bam_path);
        assert!(bam_reader.is_err());
    }

    #[test]
    fn test_sam_writer_example() {
        fn from_bam_with_filter<F>(bamfile: &str, samfile: &str, f: F) -> bool
        where
            F: Fn(&record::Record) -> Option<bool>,
        {
            let mut bam_reader = Reader::from_path(bamfile).unwrap(); // internal functions, just unwarp
            let header = header::Header::from_template(bam_reader.header());
            let mut sam_writer = Writer::from_path(samfile, &header, Format::Sam).unwrap();
            for record in bam_reader.records() {
                if record.is_err() {
                    return false;
                }
                let parsed = record.unwrap();
                match f(&parsed) {
                    None => return true,
                    Some(false) => {}
                    Some(true) => {
                        if sam_writer.write(&parsed).is_err() {
                            return false;
                        }
                    }
                }
            }
            true
        }
        use std::fs::File;
        use std::io::Read;
        let bamfile = "./test/bam2sam_test.bam";
        let samfile = "./test/bam2sam_out.sam";
        let expectedfile = "./test/bam2sam_expected.sam";
        let result = from_bam_with_filter(bamfile, samfile, |_| Some(true));
        assert!(result);
        let mut expected = Vec::new();
        let mut written = Vec::new();
        assert!(File::open(expectedfile)
            .unwrap()
            .read_to_end(&mut expected)
            .is_ok());
        assert!(File::open(samfile)
            .unwrap()
            .read_to_end(&mut written)
            .is_ok());
        assert_eq!(expected, written);
    }

    // #[cfg(feature = "curl")]
    // #[test]
    // fn test_http_connect() {
    //     let url: Url = Url::parse(
    //         "https://raw.githubusercontent.com/brainstorm/tiny-test-data/master/wgs/mt.bam",
    //     )
    //     .unwrap();
    //     let r = Reader::from_url(&url);
    //     println!("{:#?}", r);
    //     let r = r.unwrap();

    //     assert_eq!(r.header().target_names()[0], b"chr1");
    // }

    #[test]
    fn test_rc_records() {
        let (names, flags, seqs, quals, cigars) = gold();
        let mut bam = Reader::from_path(&Path::new("test/test.bam")).expect("Error opening file.");
        let del_len = [1, 1, 1, 1, 1, 100000];

        for (i, record) in bam.rc_records().enumerate() {
            //let rec = record.expect("Expected valid record");
            let rec = record.unwrap();
            println!("{}", str::from_utf8(rec.qname()).ok().unwrap());
            assert_eq!(rec.qname(), names[i]);
            assert_eq!(rec.flags(), flags[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);

            let cigar = rec.cigar();
            assert_eq!(*cigar, cigars[i]);

            let end_pos = cigar.end_pos();
            assert_eq!(end_pos, rec.pos() + 100 + del_len[i]);
            assert_eq!(
                cigar
                    .read_pos(end_pos as u32 - 10, false, false)
                    .unwrap()
                    .unwrap(),
                90
            );
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
    fn test_aux_arrays() {
        let bam_header = Header::new();
        let mut test_record = Record::from_sam(
            &mut HeaderView::from_header(&bam_header),
            "ali1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tFFFF".as_bytes(),
        )
        .unwrap();

        let array_i8: Vec<i8> = vec![std::i8::MIN, -1, 0, 1, std::i8::MAX];
        let array_u8: Vec<u8> = vec![std::u8::MIN, 0, 1, std::u8::MAX];
        let array_i16: Vec<i16> = vec![std::i16::MIN, -1, 0, 1, std::i16::MAX];
        let array_u16: Vec<u16> = vec![std::u16::MIN, 0, 1, std::u16::MAX];
        let array_i32: Vec<i32> = vec![std::i32::MIN, -1, 0, 1, std::i32::MAX];
        let array_u32: Vec<u32> = vec![std::u32::MIN, 0, 1, std::u32::MAX];
        let array_f32: Vec<f32> = vec![std::f32::MIN, 0.0, -0.0, 0.1, 0.99, std::f32::MAX];

        test_record
            .push_aux(b"XA", Aux::ArrayI8((&array_i8).into()))
            .unwrap();
        test_record
            .push_aux(b"XB", Aux::ArrayU8((&array_u8).into()))
            .unwrap();
        test_record
            .push_aux(b"XC", Aux::ArrayI16((&array_i16).into()))
            .unwrap();
        test_record
            .push_aux(b"XD", Aux::ArrayU16((&array_u16).into()))
            .unwrap();
        test_record
            .push_aux(b"XE", Aux::ArrayI32((&array_i32).into()))
            .unwrap();
        test_record
            .push_aux(b"XF", Aux::ArrayU32((&array_u32).into()))
            .unwrap();
        test_record
            .push_aux(b"XG", Aux::ArrayFloat((&array_f32).into()))
            .unwrap();

        {
            let tag = b"XA";
            if let Ok(Aux::ArrayI8(array)) = test_record.aux(tag) {
                // Retrieve aux array
                let aux_array_content = array.iter().collect::<Vec<_>>();
                assert_eq!(aux_array_content, array_i8);

                // Copy the stored aux array to another record
                {
                    let mut copy_test_record = test_record.clone();

                    // Pushing a field with an existing tag should fail
                    assert!(copy_test_record.push_aux(tag, Aux::I8(3)).is_err());

                    // Remove aux array from target record
                    copy_test_record.remove_aux(tag).unwrap();
                    assert!(copy_test_record.aux(tag).is_err());

                    // Copy array to target record
                    let src_aux = test_record.aux(tag).unwrap();
                    assert!(copy_test_record.push_aux(tag, src_aux).is_ok());
                    if let Ok(Aux::ArrayI8(array)) = copy_test_record.aux(tag) {
                        let aux_array_content_copied = array.iter().collect::<Vec<_>>();
                        assert_eq!(aux_array_content_copied, array_i8);
                    } else {
                        panic!("Aux tag not found");
                    }
                }
            } else {
                panic!("Aux tag not found");
            }
        }

        {
            let tag = b"XB";
            if let Ok(Aux::ArrayU8(array)) = test_record.aux(tag) {
                // Retrieve aux array
                let aux_array_content = array.iter().collect::<Vec<_>>();
                assert_eq!(aux_array_content, array_u8);

                // Copy the stored aux array to another record
                {
                    let mut copy_test_record = test_record.clone();

                    // Pushing a field with an existing tag should fail
                    assert!(copy_test_record.push_aux(tag, Aux::U8(3)).is_err());

                    // Remove aux array from target record
                    copy_test_record.remove_aux(tag).unwrap();
                    assert!(copy_test_record.aux(tag).is_err());

                    // Copy array to target record
                    let src_aux = test_record.aux(tag).unwrap();
                    assert!(copy_test_record.push_aux(tag, src_aux).is_ok());
                    if let Ok(Aux::ArrayU8(array)) = copy_test_record.aux(tag) {
                        let aux_array_content_copied = array.iter().collect::<Vec<_>>();
                        assert_eq!(aux_array_content_copied, array_u8);
                    } else {
                        panic!("Aux tag not found");
                    }
                }
            } else {
                panic!("Aux tag not found");
            }
        }

        {
            let tag = b"XC";
            if let Ok(Aux::ArrayI16(array)) = test_record.aux(tag) {
                // Retrieve aux array
                let aux_array_content = array.iter().collect::<Vec<_>>();
                assert_eq!(aux_array_content, array_i16);

                // Copy the stored aux array to another record
                {
                    let mut copy_test_record = test_record.clone();

                    // Pushing a field with an existing tag should fail
                    assert!(copy_test_record.push_aux(tag, Aux::I16(3)).is_err());

                    // Remove aux array from target record
                    copy_test_record.remove_aux(tag).unwrap();
                    assert!(copy_test_record.aux(tag).is_err());

                    // Copy array to target record
                    let src_aux = test_record.aux(tag).unwrap();
                    assert!(copy_test_record.push_aux(tag, src_aux).is_ok());
                    if let Ok(Aux::ArrayI16(array)) = copy_test_record.aux(tag) {
                        let aux_array_content_copied = array.iter().collect::<Vec<_>>();
                        assert_eq!(aux_array_content_copied, array_i16);
                    } else {
                        panic!("Aux tag not found");
                    }
                }
            } else {
                panic!("Aux tag not found");
            }
        }

        {
            let tag = b"XD";
            if let Ok(Aux::ArrayU16(array)) = test_record.aux(tag) {
                // Retrieve aux array
                let aux_array_content = array.iter().collect::<Vec<_>>();
                assert_eq!(aux_array_content, array_u16);

                // Copy the stored aux array to another record
                {
                    let mut copy_test_record = test_record.clone();

                    // Pushing a field with an existing tag should fail
                    assert!(copy_test_record.push_aux(tag, Aux::U16(3)).is_err());

                    // Remove aux array from target record
                    copy_test_record.remove_aux(tag).unwrap();
                    assert!(copy_test_record.aux(tag).is_err());

                    // Copy array to target record
                    let src_aux = test_record.aux(tag).unwrap();
                    assert!(copy_test_record.push_aux(tag, src_aux).is_ok());
                    if let Ok(Aux::ArrayU16(array)) = copy_test_record.aux(tag) {
                        let aux_array_content_copied = array.iter().collect::<Vec<_>>();
                        assert_eq!(aux_array_content_copied, array_u16);
                    } else {
                        panic!("Aux tag not found");
                    }
                }
            } else {
                panic!("Aux tag not found");
            }
        }

        {
            let tag = b"XE";
            if let Ok(Aux::ArrayI32(array)) = test_record.aux(tag) {
                // Retrieve aux array
                let aux_array_content = array.iter().collect::<Vec<_>>();
                assert_eq!(aux_array_content, array_i32);

                // Copy the stored aux array to another record
                {
                    let mut copy_test_record = test_record.clone();

                    // Pushing a field with an existing tag should fail
                    assert!(copy_test_record.push_aux(tag, Aux::I32(3)).is_err());

                    // Remove aux array from target record
                    copy_test_record.remove_aux(tag).unwrap();
                    assert!(copy_test_record.aux(tag).is_err());

                    // Copy array to target record
                    let src_aux = test_record.aux(tag).unwrap();
                    assert!(copy_test_record.push_aux(tag, src_aux).is_ok());
                    if let Ok(Aux::ArrayI32(array)) = copy_test_record.aux(tag) {
                        let aux_array_content_copied = array.iter().collect::<Vec<_>>();
                        assert_eq!(aux_array_content_copied, array_i32);
                    } else {
                        panic!("Aux tag not found");
                    }
                }
            } else {
                panic!("Aux tag not found");
            }
        }

        {
            let tag = b"XF";
            if let Ok(Aux::ArrayU32(array)) = test_record.aux(tag) {
                // Retrieve aux array
                let aux_array_content = array.iter().collect::<Vec<_>>();
                assert_eq!(aux_array_content, array_u32);

                // Copy the stored aux array to another record
                {
                    let mut copy_test_record = test_record.clone();

                    // Pushing a field with an existing tag should fail
                    assert!(copy_test_record.push_aux(tag, Aux::U32(3)).is_err());

                    // Remove aux array from target record
                    copy_test_record.remove_aux(tag).unwrap();
                    assert!(copy_test_record.aux(tag).is_err());

                    // Copy array to target record
                    let src_aux = test_record.aux(tag).unwrap();
                    assert!(copy_test_record.push_aux(tag, src_aux).is_ok());
                    if let Ok(Aux::ArrayU32(array)) = copy_test_record.aux(tag) {
                        let aux_array_content_copied = array.iter().collect::<Vec<_>>();
                        assert_eq!(aux_array_content_copied, array_u32);
                    } else {
                        panic!("Aux tag not found");
                    }
                }
            } else {
                panic!("Aux tag not found");
            }
        }

        {
            let tag = b"XG";
            if let Ok(Aux::ArrayFloat(array)) = test_record.aux(tag) {
                // Retrieve aux array
                let aux_array_content = array.iter().collect::<Vec<_>>();
                assert_eq!(aux_array_content, array_f32);

                // Copy the stored aux array to another record
                {
                    let mut copy_test_record = test_record.clone();

                    // Pushing a field with an existing tag should fail
                    assert!(copy_test_record.push_aux(tag, Aux::Float(3.0)).is_err());

                    // Remove aux array from target record
                    copy_test_record.remove_aux(tag).unwrap();
                    assert!(copy_test_record.aux(tag).is_err());

                    // Copy array to target record
                    let src_aux = test_record.aux(tag).unwrap();
                    assert!(copy_test_record.push_aux(tag, src_aux).is_ok());
                    if let Ok(Aux::ArrayFloat(array)) = copy_test_record.aux(tag) {
                        let aux_array_content_copied = array.iter().collect::<Vec<_>>();
                        assert_eq!(aux_array_content_copied, array_f32);
                    } else {
                        panic!("Aux tag not found");
                    }
                }
            } else {
                panic!("Aux tag not found");
            }
        }

        // Test via `Iterator` impl
        for item in test_record.aux_iter() {
            match item.unwrap() {
                (b"XA", Aux::ArrayI8(array)) => {
                    assert_eq!(&array.iter().collect::<Vec<_>>(), &array_i8);
                }
                (b"XB", Aux::ArrayU8(array)) => {
                    assert_eq!(&array.iter().collect::<Vec<_>>(), &array_u8);
                }
                (b"XC", Aux::ArrayI16(array)) => {
                    assert_eq!(&array.iter().collect::<Vec<_>>(), &array_i16);
                }
                (b"XD", Aux::ArrayU16(array)) => {
                    assert_eq!(&array.iter().collect::<Vec<_>>(), &array_u16);
                }
                (b"XE", Aux::ArrayI32(array)) => {
                    assert_eq!(&array.iter().collect::<Vec<_>>(), &array_i32);
                }
                (b"XF", Aux::ArrayU32(array)) => {
                    assert_eq!(&array.iter().collect::<Vec<_>>(), &array_u32);
                }
                (b"XG", Aux::ArrayFloat(array)) => {
                    assert_eq!(&array.iter().collect::<Vec<_>>(), &array_f32);
                }
                _ => {
                    panic!();
                }
            }
        }

        // Test via `PartialEq` impl
        assert_eq!(
            test_record.aux(b"XA").unwrap(),
            Aux::ArrayI8((&array_i8).into())
        );
        assert_eq!(
            test_record.aux(b"XB").unwrap(),
            Aux::ArrayU8((&array_u8).into())
        );
        assert_eq!(
            test_record.aux(b"XC").unwrap(),
            Aux::ArrayI16((&array_i16).into())
        );
        assert_eq!(
            test_record.aux(b"XD").unwrap(),
            Aux::ArrayU16((&array_u16).into())
        );
        assert_eq!(
            test_record.aux(b"XE").unwrap(),
            Aux::ArrayI32((&array_i32).into())
        );
        assert_eq!(
            test_record.aux(b"XF").unwrap(),
            Aux::ArrayU32((&array_u32).into())
        );
        assert_eq!(
            test_record.aux(b"XG").unwrap(),
            Aux::ArrayFloat((&array_f32).into())
        );
    }

    #[test]
    fn test_aux_scalars() {
        let bam_header = Header::new();
        let mut test_record = Record::from_sam(
            &mut HeaderView::from_header(&bam_header),
            "ali1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tFFFF".as_bytes(),
        )
        .unwrap();

        test_record.push_aux(b"XA", Aux::I8(i8::MIN)).unwrap();
        test_record.push_aux(b"XB", Aux::I8(i8::MAX)).unwrap();
        test_record.push_aux(b"XC", Aux::U8(u8::MIN)).unwrap();
        test_record.push_aux(b"XD", Aux::U8(u8::MAX)).unwrap();
        test_record.push_aux(b"XE", Aux::I16(i16::MIN)).unwrap();
        test_record.push_aux(b"XF", Aux::I16(i16::MAX)).unwrap();
        test_record.push_aux(b"XG", Aux::U16(u16::MIN)).unwrap();
        test_record.push_aux(b"XH", Aux::U16(u16::MAX)).unwrap();
        test_record.push_aux(b"XI", Aux::I32(i32::MIN)).unwrap();
        test_record.push_aux(b"XJ", Aux::I32(i32::MAX)).unwrap();
        test_record.push_aux(b"XK", Aux::U32(u32::MIN)).unwrap();
        test_record.push_aux(b"XL", Aux::U32(u32::MAX)).unwrap();
        test_record
            .push_aux(b"XM", Aux::Float(std::f32::consts::PI))
            .unwrap();
        test_record
            .push_aux(b"XN", Aux::Double(std::f64::consts::PI))
            .unwrap();
        test_record
            .push_aux(b"XO", Aux::String("Test str"))
            .unwrap();
        test_record.push_aux(b"XP", Aux::I8(0)).unwrap();

        let collected_aux_fields = test_record.aux_iter().collect::<Result<Vec<_>>>().unwrap();
        assert_eq!(
            collected_aux_fields,
            vec![
                (&b"XA"[..], Aux::I8(i8::MIN)),
                (&b"XB"[..], Aux::I8(i8::MAX)),
                (&b"XC"[..], Aux::U8(u8::MIN)),
                (&b"XD"[..], Aux::U8(u8::MAX)),
                (&b"XE"[..], Aux::I16(i16::MIN)),
                (&b"XF"[..], Aux::I16(i16::MAX)),
                (&b"XG"[..], Aux::U16(u16::MIN)),
                (&b"XH"[..], Aux::U16(u16::MAX)),
                (&b"XI"[..], Aux::I32(i32::MIN)),
                (&b"XJ"[..], Aux::I32(i32::MAX)),
                (&b"XK"[..], Aux::U32(u32::MIN)),
                (&b"XL"[..], Aux::U32(u32::MAX)),
                (&b"XM"[..], Aux::Float(std::f32::consts::PI)),
                (&b"XN"[..], Aux::Double(std::f64::consts::PI)),
                (&b"XO"[..], Aux::String("Test str")),
                (&b"XP"[..], Aux::I8(0)),
            ]
        );
    }

    #[test]
    fn test_aux_array_partial_eq() {
        use record::AuxArray;

        // Target types
        let one_data: Vec<i8> = vec![0, 1, 2, 3, 4, 5, 6];
        let one_aux_array = AuxArray::from(&one_data);

        let two_data: Vec<i8> = vec![0, 1, 2, 3, 4, 5];
        let two_aux_array = AuxArray::from(&two_data);

        assert_ne!(&one_data, &two_data);
        assert_ne!(&one_aux_array, &two_aux_array);

        let one_aux = Aux::ArrayI8(one_aux_array);
        let two_aux = Aux::ArrayI8(two_aux_array);
        assert_ne!(&one_aux, &two_aux);

        // Raw bytes
        let bam_header = Header::new();
        let mut test_record = Record::from_sam(
            &mut HeaderView::from_header(&bam_header),
            "ali1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tFFFF".as_bytes(),
        )
        .unwrap();

        test_record.push_aux(b"XA", one_aux).unwrap();
        test_record.push_aux(b"XB", two_aux).unwrap();

        // RawLeBytes == RawLeBytes
        assert_eq!(
            test_record.aux(b"XA").unwrap(),
            test_record.aux(b"XA").unwrap()
        );
        // RawLeBytes != RawLeBytes
        assert_ne!(
            test_record.aux(b"XA").unwrap(),
            test_record.aux(b"XB").unwrap()
        );

        // RawLeBytes == TargetType
        assert_eq!(
            test_record.aux(b"XA").unwrap(),
            Aux::ArrayI8((&one_data).into())
        );
        assert_eq!(
            test_record.aux(b"XB").unwrap(),
            Aux::ArrayI8((&two_data).into())
        );
        // RawLeBytes != TargetType
        assert_ne!(
            test_record.aux(b"XA").unwrap(),
            Aux::ArrayI8((&two_data).into())
        );
        assert_ne!(
            test_record.aux(b"XB").unwrap(),
            Aux::ArrayI8((&one_data).into())
        );
    }

    /// Test if both text and binary representations of a BAM header are in sync (#156)
    #[test]
    fn test_bam_header_sync() {
        let reader = Reader::from_path("test/test_issue_156_no_text.bam").unwrap();
        let header_hashmap = Header::from_template(reader.header()).to_hashmap();
        let header_refseqs = header_hashmap.get("SQ".into()).unwrap();
        assert_eq!(header_refseqs[0].get("SN").unwrap(), "ref_1",);
        assert_eq!(header_refseqs[0].get("LN").unwrap(), "10000000",);
    }

    #[test]
    fn test_idxstats_bam() {
        let mut reader = IndexedReader::from_path("test/test.bam").unwrap();
        let expected = vec![
            (0, 15072423, 6, 0),
            (1, 15279345, 0, 0),
            (2, 13783700, 0, 0),
            (3, 17493793, 0, 0),
            (4, 20924149, 0, 0),
            (-1, 0, 0, 0),
        ];
        let actual = reader.index_stats().unwrap();
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_number_mapped_and_unmapped_bam() {
        let reader = IndexedReader::from_path("test/test.bam").unwrap();
        let expected = (6, 0);
        let actual = reader.index().number_mapped_unmapped(0);
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_number_unmapped_global_bam() {
        let reader = IndexedReader::from_path("test/test_unmapped.bam").unwrap();
        let expected = 8;
        let actual = reader.index().number_unmapped();
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_idxstats_cram() {
        let mut reader = IndexedReader::from_path("test/test_cram.cram").unwrap();
        reader.set_reference("test/test_cram.fa").unwrap();
        let expected = vec![
            (0, 120, 2, 0),
            (1, 120, 2, 0),
            (2, 120, 2, 0),
            (-1, 0, 0, 0),
        ];
        let actual = reader.index_stats().unwrap();
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_slow_idxstats_cram() {
        let mut reader = IndexedReader::from_path("test/test_cram.cram").unwrap();
        reader.set_reference("test/test_cram.fa").unwrap();
        let expected = vec![
            (0, 120, 2, 0),
            (1, 120, 2, 0),
            (2, 120, 2, 0),
            (-1, 0, 0, 0),
        ];
        let actual = reader.index_stats().unwrap();
        assert_eq!(expected, actual);
    }

    // #[test]
    // fn test_number_mapped_and_unmapped_cram() {
    //     let mut reader = IndexedReader::from_path("test/test_cram.cram").unwrap();
    //     reader.set_reference("test/test_cram.fa").unwrap();
    //     let expected = (2, 0);
    //     let actual = reader.index().number_mapped_unmapped(0);
    //     assert_eq!(expected, actual);
    // }
    //
    // #[test]
    // fn test_number_unmapped_global_cram() {
    //     let mut reader = IndexedReader::from_path("test/test_unmapped.cram").unwrap();
    //     let expected = 8;
    //     let actual = reader.index().number_unmapped();
    //     assert_eq!(expected, actual);
    // }
}
