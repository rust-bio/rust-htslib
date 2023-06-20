// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::borrow::{Borrow, BorrowMut};
use std::f32;
use std::ffi;
use std::fmt;
use std::i32;
use std::marker::PhantomData;
use std::ops::Deref;
use std::os::raw::c_char;
use std::ptr;
use std::rc::Rc;
use std::slice;
use std::str;

use bio_types::genome;
use derive_new::new;
use ieee754::Ieee754;
use lazy_static::lazy_static;

use crate::bcf::header::{HeaderView, Id};
use crate::bcf::Error;
use crate::errors::Result;
use crate::htslib;

const MISSING_INTEGER: i32 = i32::MIN;
const VECTOR_END_INTEGER: i32 = i32::MIN + 1;

lazy_static! {
    static ref MISSING_FLOAT: f32 = Ieee754::from_bits(0x7F80_0001);
    static ref VECTOR_END_FLOAT: f32 = Ieee754::from_bits(0x7F80_0002);
}

/// Common methods for numeric INFO and FORMAT entries
pub trait Numeric {
    /// Return true if entry is a missing value
    fn is_missing(&self) -> bool;

    /// Return missing value for storage in BCF record.
    fn missing() -> Self;
}

impl Numeric for f32 {
    fn is_missing(&self) -> bool {
        self.bits() == MISSING_FLOAT.bits()
    }

    fn missing() -> f32 {
        *MISSING_FLOAT
    }
}

impl Numeric for i32 {
    fn is_missing(&self) -> bool {
        *self == MISSING_INTEGER
    }

    fn missing() -> i32 {
        MISSING_INTEGER
    }
}

trait NumericUtils {
    /// Return true if entry marks the end of the record.
    fn is_vector_end(&self) -> bool;
}

impl NumericUtils for f32 {
    fn is_vector_end(&self) -> bool {
        self.bits() == VECTOR_END_FLOAT.bits()
    }
}

impl NumericUtils for i32 {
    fn is_vector_end(&self) -> bool {
        *self == VECTOR_END_INTEGER
    }
}

/// A trait to allow for seamless use of bytes or integer identifiers for filters
pub trait FilterId {
    fn id_from_header(&self, header: &HeaderView) -> Result<Id>;
    fn is_pass(&self) -> bool;
}

impl FilterId for [u8] {
    fn id_from_header(&self, header: &HeaderView) -> Result<Id> {
        header.name_to_id(self)
    }
    fn is_pass(&self) -> bool {
        matches!(self, b"PASS" | b".")
    }
}

impl FilterId for Id {
    fn id_from_header(&self, _header: &HeaderView) -> Result<Id> {
        Ok(*self)
    }
    fn is_pass(&self) -> bool {
        *self == Id(0)
    }
}

/// A buffer for info or format data.
#[derive(Debug)]
pub struct Buffer {
    inner: *mut ::std::os::raw::c_void,
    len: i32,
}

impl Buffer {
    pub fn new() -> Self {
        Buffer {
            inner: ptr::null_mut(),
            len: 0,
        }
    }
}

impl Drop for Buffer {
    fn drop(&mut self) {
        unsafe {
            ::libc::free(self.inner as *mut ::libc::c_void);
        }
    }
}

#[derive(new, Debug)]
pub struct BufferBacked<'a, T: 'a + fmt::Debug, B: Borrow<Buffer> + 'a> {
    value: T,
    buffer: B,
    #[new(default)]
    phantom: PhantomData<&'a B>,
}

impl<'a, T: 'a + fmt::Debug, B: Borrow<Buffer> + 'a> Deref for BufferBacked<'a, T, B> {
    type Target = T;

    fn deref(&self) -> &T {
        &self.value
    }
}

impl<'a, T: 'a + fmt::Debug + fmt::Display, B: Borrow<Buffer> + 'a> fmt::Display
    for BufferBacked<'a, T, B>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.value, f)
    }
}

/// A VCF/BCF record.
/// New records can be created by the `empty_record` methods of [`bcf::Reader`](crate::bcf::Reader)
/// and [`bcf::Writer`](crate::bcf::Writer).
/// # Example
/// ```rust
/// use rust_htslib::bcf::{Format, Writer};
/// use rust_htslib::bcf::header::Header;
///
/// // Create minimal VCF header with a single sample
/// let mut header = Header::new();
/// header.push_sample("sample".as_bytes());
///
/// // Write uncompressed VCF to stdout with above header and get an empty record
/// let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
/// let mut record = vcf.empty_record();
/// ```
#[derive(Debug)]
pub struct Record {
    pub inner: *mut htslib::bcf1_t,
    header: Rc<HeaderView>,
}

impl Record {
    /// Construct record with reference to header `HeaderView`, for create-internal use.
    pub(crate) fn new(header: Rc<HeaderView>) -> Self {
        let inner = unsafe {
            let inner = htslib::bcf_init();
            // Always unpack record.
            htslib::bcf_unpack(inner, htslib::BCF_UN_ALL as i32);
            inner
        };
        Record { inner, header }
    }

    /// Force unpacking of internal record values.
    pub fn unpack(&mut self) {
        unsafe { htslib::bcf_unpack(self.inner, htslib::BCF_UN_ALL as i32) };
    }

    /// Return associated header.
    pub fn header(&self) -> &HeaderView {
        self.header.as_ref()
    }

    /// Set the record header.
    pub(crate) fn set_header(&mut self, header: Rc<HeaderView>) {
        self.header = header;
    }

    /// Return reference to the inner C struct.
    ///
    /// # Remarks
    ///
    /// Note that this function is only required as long as Rust-Htslib does not provide full
    /// access to all aspects of Htslib.
    pub fn inner(&self) -> &htslib::bcf1_t {
        unsafe { &*self.inner }
    }

    /// Return mutable reference to inner C struct.
    ///
    /// # Remarks
    ///
    /// Note that this function is only required as long as Rust-Htslib does not provide full
    /// access to all aspects of Htslib.
    pub fn inner_mut(&mut self) -> &mut htslib::bcf1_t {
        unsafe { &mut *self.inner }
    }

    /// Get the reference id of the record.
    ///
    /// To look up the contig name,
    /// use [`HeaderView::rid2name`](../header/struct.HeaderView.html#method.rid2name).
    ///
    /// # Returns
    ///
    /// - `Some(rid)` if the internal `rid` is set to a value that is not `-1`
    /// - `None` if the internal `rid` is set to `-1`
    pub fn rid(&self) -> Option<u32> {
        match self.inner().rid {
            -1 => None,
            rid => Some(rid as u32),
        }
    }

    /// Update the reference id of the record.
    ///
    /// To look up reference id for a contig name,
    /// use [`HeaderView::name2rid`](../header/struct.HeaderView.html#method.name2rid).
    ///
    /// # Example
    ///
    /// Example assumes we have a Record `record` from a VCF with a header containing region
    /// named `1`. See [module documentation](../index.html#example-writing) for how to set
    /// up VCF, header, and record.
    ///
    /// ```
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// # let mut header = Header::new();
    /// # let header_contig_line = r#"##contig=<ID=1,length=10>"#;
    /// # header.push_record(header_contig_line.as_bytes());
    /// # header.push_sample("test_sample".as_bytes());
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// let rid = record.header().name2rid(b"1").ok();
    /// record.set_rid(rid);
    /// assert_eq!(record.rid(), rid);
    /// let name = record.header().rid2name(record.rid().unwrap()).ok();
    /// assert_eq!(Some("1".as_bytes()), name);
    /// ```
    pub fn set_rid(&mut self, rid: Option<u32>) {
        match rid {
            Some(rid) => self.inner_mut().rid = rid as i32,
            None => self.inner_mut().rid = -1,
        }
    }

    /// Return **0-based** position
    pub fn pos(&self) -> i64 {
        self.inner().pos
    }

    /// Set **0-based** position
    pub fn set_pos(&mut self, pos: i64) {
        self.inner_mut().pos = pos;
    }

    /// Return the **0-based, exclusive** end position
    ///
    /// # Example
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Header, Writer};
    /// # use tempfile::NamedTempFile;
    /// # let tmp = NamedTempFile::new().unwrap();
    /// # let path = tmp.path();
    /// # let header = Header::new();
    /// # let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// let alleles: &[&[u8]] = &[b"AGG", b"TG"];
    /// record.set_alleles(alleles).expect("Failed to set alleles");
    /// record.set_pos(5);
    ///
    /// assert_eq!(record.end(), 8)
    /// ```
    pub fn end(&self) -> i64 {
        self.pos() + self.rlen()
    }

    /// Return the value of the ID column.
    ///
    /// When empty, returns `b".".to_vec()`.
    pub fn id(&self) -> Vec<u8> {
        if self.inner().d.id.is_null() {
            b".".to_vec()
        } else {
            let id = unsafe { ffi::CStr::from_ptr(self.inner().d.id) };
            id.to_bytes().to_vec()
        }
    }

    /// Update the ID string to the given value.
    pub fn set_id(&mut self, id: &[u8]) -> Result<()> {
        let c_str = ffi::CString::new(id).unwrap();
        if unsafe {
            htslib::bcf_update_id(
                self.header().inner,
                self.inner,
                c_str.as_ptr() as *mut c_char,
            )
        } == 0
        {
            Ok(())
        } else {
            Err(Error::BcfSetValues)
        }
    }

    /// Clear the ID column (set it to `"."`).
    pub fn clear_id(&mut self) -> Result<()> {
        let c_str = ffi::CString::new(&b"."[..]).unwrap();
        if unsafe {
            htslib::bcf_update_id(
                self.header().inner,
                self.inner,
                c_str.as_ptr() as *mut c_char,
            )
        } == 0
        {
            Ok(())
        } else {
            Err(Error::BcfSetValues)
        }
    }

    /// Add the ID string (the ID field is semicolon-separated), checking for duplicates.
    pub fn push_id(&mut self, id: &[u8]) -> Result<()> {
        let c_str = ffi::CString::new(id).unwrap();
        if unsafe {
            htslib::bcf_add_id(
                self.header().inner,
                self.inner,
                c_str.as_ptr() as *mut c_char,
            )
        } == 0
        {
            Ok(())
        } else {
            Err(Error::BcfSetValues)
        }
    }

    /// Return `Filters` iterator for enumerating all filters that have been set.
    ///
    /// A record having the `PASS` filter will return an empty `Filter` here.
    pub fn filters(&self) -> Filters<'_> {
        Filters::new(self)
    }

    /// Query whether the filter with the given ID has been set.
    ///
    /// This method can be used to check if a record passes filtering by using either `Id(0)`,
    /// `PASS` or `.`
    ///
    /// # Example
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Header, Writer};
    /// # use rust_htslib::bcf::header::Id;
    /// # use tempfile::NamedTempFile;
    /// # let tmp = tempfile::NamedTempFile::new().unwrap();
    /// # let path = tmp.path();
    /// let mut header = Header::new();
    /// header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
    /// # let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// assert!(record.has_filter("PASS".as_bytes()));
    /// assert!(record.has_filter(".".as_bytes()));
    /// assert!(record.has_filter(&Id(0)));
    ///
    /// record.push_filter("foo".as_bytes()).unwrap();
    /// assert!(record.has_filter("foo".as_bytes()));
    /// assert!(!record.has_filter("PASS".as_bytes()))
    /// ```
    pub fn has_filter<T: FilterId + ?Sized>(&self, flt_id: &T) -> bool {
        if flt_id.is_pass() && self.inner().d.n_flt == 0 {
            return true;
        }
        let id = match flt_id.id_from_header(self.header()) {
            Ok(i) => *i,
            Err(_) => return false,
        };
        for i in 0..(self.inner().d.n_flt as isize) {
            if unsafe { *self.inner().d.flt.offset(i) } == id as i32 {
                return true;
            }
        }
        false
    }

    /// Set the given filter IDs to the FILTER column.
    ///
    /// Setting an empty slice removes all filters and sets `PASS`.
    ///
    /// # Example
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Header, Writer};
    /// # use rust_htslib::bcf::header::Id;
    /// # use tempfile::NamedTempFile;
    /// # let tmp = tempfile::NamedTempFile::new().unwrap();
    /// # let path = tmp.path();
    /// let mut header = Header::new();
    /// header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
    /// header.push_record(br#"##FILTER=<ID=bar,Description="a horse walks into...">"#);
    /// # let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// let foo = record.header().name_to_id(b"foo").unwrap();
    /// let bar = record.header().name_to_id(b"bar").unwrap();
    /// assert!(record.has_filter("PASS".as_bytes()));
    /// let mut filters = vec![&foo, &bar];
    /// record.set_filters(&filters).unwrap();
    /// assert!(record.has_filter(&foo));
    /// assert!(record.has_filter(&bar));
    /// assert!(!record.has_filter("PASS".as_bytes()));
    /// filters.clear();
    /// record.set_filters(&filters).unwrap();
    /// assert!(record.has_filter("PASS".as_bytes()));
    /// assert!(!record.has_filter("foo".as_bytes()));
    /// // 'baz' isn't in the header
    /// assert!(record.set_filters(&["baz".as_bytes()]).is_err())
    /// ```
    ///
    /// # Errors
    /// If any of the filter IDs do not exist in the header, an [`Error::BcfUnknownID`] is returned.
    ///
    pub fn set_filters<T: FilterId + ?Sized>(&mut self, flt_ids: &[&T]) -> Result<()> {
        let mut ids: Vec<i32> = flt_ids
            .iter()
            .map(|id| id.id_from_header(self.header()).map(|id| *id as i32))
            .collect::<Result<Vec<i32>>>()?;
        unsafe {
            htslib::bcf_update_filter(
                self.header().inner,
                self.inner,
                ids.as_mut_ptr(),
                ids.len() as i32,
            );
        };
        Ok(())
    }

    /// Add the given filter to the FILTER column.
    ///
    /// If `flt_id` is `PASS` or `.` then all existing filters are removed first. Otherwise,
    /// any existing `PASS` filter is removed.
    ///
    /// # Example
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Header, Writer};
    /// # use tempfile::NamedTempFile;
    /// # let tmp = tempfile::NamedTempFile::new().unwrap();
    /// # let path = tmp.path();
    /// let mut header = Header::new();
    /// header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
    /// header.push_record(br#"##FILTER=<ID=bar,Description="dranks">"#);
    /// # let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// let foo = "foo".as_bytes();
    /// let bar = record.header().name_to_id(b"bar").unwrap();
    /// assert!(record.has_filter("PASS".as_bytes()));
    ///
    /// record.push_filter(foo).unwrap();
    /// record.push_filter(&bar).unwrap();
    /// assert!(record.has_filter(foo));
    /// assert!(record.has_filter(&bar));
    /// // filter must exist in the header
    /// assert!(record.push_filter("baz".as_bytes()).is_err())
    /// ```
    ///
    /// # Errors
    /// If the `flt_id` does not exist in the header, an [`Error::BcfUnknownID`] is returned.
    ///
    pub fn push_filter<T: FilterId + ?Sized>(&mut self, flt_id: &T) -> Result<()> {
        let id = flt_id.id_from_header(self.header())?;
        unsafe {
            htslib::bcf_add_filter(self.header().inner, self.inner, *id as i32);
        };
        Ok(())
    }

    /// Remove the given filter from the FILTER column.
    ///
    /// # Arguments
    ///
    /// - `flt_id` - The corresponding filter ID to remove.
    /// - `pass_on_empty` - Set to `PASS` when removing the last filter.
    ///
    /// # Example
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Header, Writer};
    /// # use tempfile::NamedTempFile;
    /// # let tmp = tempfile::NamedTempFile::new().unwrap();
    /// # let path = tmp.path();
    /// let mut header = Header::new();
    /// header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
    /// header.push_record(br#"##FILTER=<ID=bar,Description="a horse walks into...">"#);
    /// # let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// let foo = "foo".as_bytes();
    /// let bar = "bar".as_bytes();
    /// record.set_filters(&[foo, bar]).unwrap();
    /// assert!(record.has_filter(foo));
    /// assert!(record.has_filter(bar));
    ///
    /// record.remove_filter(foo, true).unwrap();
    /// assert!(!record.has_filter(foo));
    /// assert!(record.has_filter(bar));
    /// // 'baz' is not in the header
    /// assert!(record.remove_filter("baz".as_bytes(), true).is_err());
    ///
    /// record.remove_filter(bar, true).unwrap();
    /// assert!(!record.has_filter(bar));
    /// assert!(record.has_filter("PASS".as_bytes()));
    /// ```
    ///
    /// # Errors
    /// If the `flt_id` does not exist in the header, an [`Error::BcfUnknownID`] is returned.
    ///
    pub fn remove_filter<T: FilterId + ?Sized>(
        &mut self,
        flt_id: &T,
        pass_on_empty: bool,
    ) -> Result<()> {
        let id = flt_id.id_from_header(self.header())?;
        unsafe {
            htslib::bcf_remove_filter(
                self.header().inner,
                self.inner,
                *id as i32,
                pass_on_empty as i32,
            )
        };
        Ok(())
    }

    /// Get alleles strings.
    ///
    /// The first allele is the reference allele.
    pub fn alleles(&self) -> Vec<&[u8]> {
        unsafe { htslib::bcf_unpack(self.inner, htslib::BCF_UN_ALL as i32) };
        let n = self.inner().n_allele() as usize;
        let dec = self.inner().d;
        let alleles = unsafe { slice::from_raw_parts(dec.allele, n) };
        (0..n)
            .map(|i| unsafe { ffi::CStr::from_ptr(alleles[i]).to_bytes() })
            .collect()
    }

    /// Set alleles. The first allele is the reference allele.
    ///
    /// # Example
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// #
    /// # // Create minimal VCF header with a single sample
    /// # let mut header = Header::new();
    /// # header.push_sample("sample".as_bytes());
    /// #
    /// # // Write uncompressed VCF to stdout with above header and get an empty record
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// assert_eq!(record.allele_count(), 0);
    ///
    /// let alleles: &[&[u8]] = &[b"A", b"TG"];
    /// record.set_alleles(alleles).expect("Failed to set alleles");
    /// assert_eq!(record.allele_count(), 2)
    /// ```
    pub fn set_alleles(&mut self, alleles: &[&[u8]]) -> Result<()> {
        let cstrings: Vec<ffi::CString> = alleles
            .iter()
            .map(|vec| ffi::CString::new(*vec).unwrap())
            .collect();
        let mut ptrs: Vec<*const c_char> = cstrings
            .iter()
            .map(|cstr| cstr.as_ptr() as *const c_char)
            .collect();
        if unsafe {
            htslib::bcf_update_alleles(
                self.header().inner,
                self.inner,
                ptrs.as_mut_ptr(),
                alleles.len() as i32,
            )
        } == 0
        {
            Ok(())
        } else {
            Err(Error::BcfSetValues)
        }
    }

    /// Get variant quality.
    pub fn qual(&self) -> f32 {
        self.inner().qual
    }

    /// Set variant quality.
    pub fn set_qual(&mut self, qual: f32) {
        self.inner_mut().qual = qual;
    }

    pub fn info<'a>(&'a self, tag: &'a [u8]) -> Info<'a, Buffer> {
        self.info_shared_buffer(tag, Buffer::new())
    }

    /// Get the value of the given info tag.
    pub fn info_shared_buffer<'a, 'b, B: BorrowMut<Buffer> + Borrow<Buffer> + 'b>(
        &'a self,
        tag: &'a [u8],
        buffer: B,
    ) -> Info<'a, B> {
        Info {
            record: self,
            tag,
            buffer,
        }
    }

    /// Get the number of samples in the record.
    pub fn sample_count(&self) -> u32 {
        self.inner().n_sample()
    }

    /// Get the number of alleles, including reference allele.
    pub fn allele_count(&self) -> u32 {
        self.inner().n_allele()
    }

    /// Add/replace genotypes in FORMAT GT tag.
    ///
    /// # Arguments
    ///
    /// - `genotypes` - a flattened, two-dimensional array of GenotypeAllele,
    ///                 the first dimension contains one array for each sample.
    ///
    /// # Errors
    ///
    /// Returns error if GT tag is not present in header.
    ///
    /// # Example
    ///
    /// Example assumes we have a Record `record` from a VCF with a `GT` `FORMAT` tag.
    /// See [module documentation](../index.html#example-writing) for how to set up
    /// VCF, header, and record.
    ///
    /// ```
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// # use rust_htslib::bcf::record::GenotypeAllele;
    /// # let mut header = Header::new();
    /// # let header_contig_line = r#"##contig=<ID=1,length=10>"#;
    /// # header.push_record(header_contig_line.as_bytes());
    /// # let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
    /// # header.push_record(header_gt_line.as_bytes());
    /// # header.push_sample("test_sample".as_bytes());
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// let alleles = &[GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)];
    /// record.push_genotypes(alleles);
    /// assert_eq!("1/1", &format!("{}", record.genotypes().unwrap().get(0)));
    /// ```
    pub fn push_genotypes(&mut self, genotypes: &[GenotypeAllele]) -> Result<()> {
        let encoded: Vec<i32> = genotypes.iter().map(|gt| i32::from(*gt)).collect();
        self.push_format_integer(b"GT", &encoded)
    }

    /// Get genotypes as vector of one `Genotype` per sample.
    ///
    /// # Example
    /// Parsing genotype field (`GT` tag) from a VCF record:
    /// ```
    /// use crate::rust_htslib::bcf::{Reader, Read};
    /// let mut vcf = Reader::from_path(&"test/test_string.vcf").expect("Error opening file.");
    /// let expected = ["./1", "1|1", "0/1", "0|1", "1|.", "1/1"];
    /// for (rec, exp_gt) in vcf.records().zip(expected.iter()) {
    ///     let mut rec = rec.expect("Error reading record.");
    ///     let genotypes = rec.genotypes().expect("Error reading genotypes");
    ///     assert_eq!(&format!("{}", genotypes.get(0)), exp_gt);
    /// }
    /// ```
    pub fn genotypes(&self) -> Result<Genotypes<'_, Buffer>> {
        self.genotypes_shared_buffer(Buffer::new())
    }

    /// Get genotypes as vector of one `Genotype` per sample, using a given shared buffer
    /// to avoid unnecessary allocations.
    pub fn genotypes_shared_buffer<'a, B>(&self, buffer: B) -> Result<Genotypes<'a, B>>
    where
        B: BorrowMut<Buffer> + Borrow<Buffer> + 'a,
    {
        Ok(Genotypes {
            encoded: self.format_shared_buffer(b"GT", buffer).integer()?,
        })
    }

    /// Retrieve data for a `FORMAT` field
    ///
    /// # Example
    /// *Note: some boilerplate for the example is hidden for clarity. See [module documentation](../index.html#example-writing)
    /// for an example of the setup used here.*
    ///
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// #
    /// # // Create minimal VCF header with a single sample
    /// # let mut header = Header::new();
    /// header.push_sample(b"sample1").push_sample(b"sample2").push_record(br#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">"#);
    /// #
    /// # // Write uncompressed VCF to stdout with above header and get an empty record
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// record.push_format_integer(b"DP", &[20, 12]).expect("Failed to set DP format field");
    ///
    /// let read_depths = record.format(b"DP").integer().expect("Couldn't retrieve DP field");
    /// let sample1_depth = read_depths[0];
    /// assert_eq!(sample1_depth, &[20]);
    /// let sample2_depth = read_depths[1];
    /// assert_eq!(sample2_depth, &[12])
    /// ```
    ///
    /// # Errors
    /// **Attention:** the returned [`BufferBacked`] from [`integer()`](Format::integer)
    /// (`read_depths`), which holds the data, has to be kept in scope as long as the data is
    /// accessed. If parts of the data are accessed after the `BufferBacked` object is been
    /// dropped, you will access unallocated memory.
    pub fn format<'a>(&'a self, tag: &'a [u8]) -> Format<'a, Buffer> {
        self.format_shared_buffer(tag, Buffer::new())
    }

    /// Get the value of the given format tag for each sample.
    pub fn format_shared_buffer<'a, 'b, B: BorrowMut<Buffer> + Borrow<Buffer> + 'b>(
        &'a self,
        tag: &'a [u8],
        buffer: B,
    ) -> Format<'a, B> {
        Format::new(self, tag, buffer)
    }

    /// Add/replace an integer-typed FORMAT tag.
    ///
    /// # Arguments
    ///
    /// - `tag` - The tag's string.
    /// - `data` - a flattened, two-dimensional array, the first dimension contains one array
    ///            for each sample.
    ///
    /// # Errors
    ///
    /// Returns error if tag is not present in header.
    pub fn push_format_integer(&mut self, tag: &[u8], data: &[i32]) -> Result<()> {
        self.push_format(tag, data, htslib::BCF_HT_INT)
    }

    /// Add/replace a float-typed FORMAT tag.
    ///
    /// # Arguments
    ///
    /// - `tag` - The tag's string.
    /// - `data` - a flattened, two-dimensional array, the first dimension contains one array
    ///            for each sample.
    ///
    /// # Errors
    ///
    /// Returns error if tag is not present in header.
    ///
    /// # Example
    ///
    /// Example assumes we have a Record `record` from a VCF with an `AF` `FORMAT` tag.
    /// See [module documentation](../index.html#example-writing) for how to set up
    /// VCF, header, and record.
    ///
    /// ```
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// # use rust_htslib::bcf::record::GenotypeAllele;
    /// # let mut header = Header::new();
    /// # let header_contig_line = r#"##contig=<ID=1,length=10>"#;
    /// # header.push_record(header_contig_line.as_bytes());
    /// # let header_af_line = r#"##FORMAT=<ID=AF,Number=1,Type=Float,Description="Frequency">"#;
    /// # header.push_record(header_af_line.as_bytes());
    /// # header.push_sample("test_sample".as_bytes());
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// record.push_format_float(b"AF", &[0.5]);
    /// assert_eq!(0.5, record.format(b"AF").float().unwrap()[0][0]);
    /// ```
    pub fn push_format_float(&mut self, tag: &[u8], data: &[f32]) -> Result<()> {
        self.push_format(tag, data, htslib::BCF_HT_REAL)
    }

    /// Add/replace a single-char-typed FORMAT tag.
    ///
    /// # Arguments
    ///
    /// - `tag` - The tag's string.
    /// - `data` - a flattened, two-dimensional array, the first dimension contains one array
    ///            for each sample.
    ///
    /// # Errors
    ///
    /// Returns error if tag is not present in header.
    pub fn push_format_char(&mut self, tag: &[u8], data: &[u8]) -> Result<()> {
        self.push_format(tag, data, htslib::BCF_HT_STR)
    }

    /// Add a format tag. Data is a flattened two-dimensional array.
    /// The first dimension contains one array for each sample.
    fn push_format<T>(&mut self, tag: &[u8], data: &[T], ht: u32) -> Result<()> {
        let tag_c_str = ffi::CString::new(tag).unwrap();
        unsafe {
            if htslib::bcf_update_format(
                self.header().inner,
                self.inner,
                tag_c_str.as_ptr() as *mut c_char,
                data.as_ptr() as *const ::std::os::raw::c_void,
                data.len() as i32,
                ht as i32,
            ) == 0
            {
                Ok(())
            } else {
                Err(Error::BcfSetTag {
                    tag: str::from_utf8(tag).unwrap().to_owned(),
                })
            }
        }
    }

    // TODO: should we add convenience methods clear_format_*?

    /// Add a string-typed FORMAT tag. Note that genotypes are treated as a special case
    /// and cannot be added with this method. See instead [push_genotypes](#method.push_genotypes).
    ///
    /// # Arguments
    ///
    /// - `tag` - The tag's string.
    /// - `data` - a two-dimensional array, the first dimension contains one array
    ///            for each sample. Must be non-empty.
    ///
    /// # Errors
    ///
    /// Returns error if tag is not present in header.
    pub fn push_format_string<D: Borrow<[u8]>>(&mut self, tag: &[u8], data: &[D]) -> Result<()> {
        assert!(
            !data.is_empty(),
            "given string data must have at least 1 element"
        );
        let c_data = data
            .iter()
            .map(|s| ffi::CString::new(s.borrow()).unwrap())
            .collect::<Vec<ffi::CString>>();
        let c_ptrs = c_data
            .iter()
            .map(|s| s.as_ptr() as *mut i8)
            .collect::<Vec<*mut i8>>();
        let tag_c_str = ffi::CString::new(tag).unwrap();
        unsafe {
            if htslib::bcf_update_format_string(
                self.header().inner,
                self.inner,
                tag_c_str.as_ptr() as *mut c_char,
                c_ptrs.as_slice().as_ptr() as *mut *const c_char,
                data.len() as i32,
            ) == 0
            {
                Ok(())
            } else {
                Err(Error::BcfSetTag {
                    tag: str::from_utf8(tag).unwrap().to_owned(),
                })
            }
        }
    }

    /// Add/replace an integer-typed INFO entry.
    pub fn push_info_integer(&mut self, tag: &[u8], data: &[i32]) -> Result<()> {
        self.push_info(tag, data, htslib::BCF_HT_INT)
    }

    /// Remove the integer-typed INFO entry.
    pub fn clear_info_integer(&mut self, tag: &[u8]) -> Result<()> {
        self.push_info::<i32>(tag, &[], htslib::BCF_HT_INT)
    }

    /// Add/replace a float-typed INFO entry.
    pub fn push_info_float(&mut self, tag: &[u8], data: &[f32]) -> Result<()> {
        self.push_info(tag, data, htslib::BCF_HT_REAL)
    }

    /// Remove the float-typed INFO entry.
    pub fn clear_info_float(&mut self, tag: &[u8]) -> Result<()> {
        self.push_info::<u8>(tag, &[], htslib::BCF_HT_REAL)
    }

    /// Add/replace an INFO tag.
    ///
    /// # Arguments
    /// * `tag` - the tag to add/replace
    /// * `data` - the data to set
    /// * `ht` - the HTSLib type to use
    fn push_info<T>(&mut self, tag: &[u8], data: &[T], ht: u32) -> Result<()> {
        let tag_c_str = ffi::CString::new(tag).unwrap();
        unsafe {
            if htslib::bcf_update_info(
                self.header().inner,
                self.inner,
                tag_c_str.as_ptr() as *mut c_char,
                data.as_ptr() as *const ::std::os::raw::c_void,
                data.len() as i32,
                ht as i32,
            ) == 0
            {
                Ok(())
            } else {
                Err(Error::BcfSetTag {
                    tag: str::from_utf8(tag).unwrap().to_owned(),
                })
            }
        }
    }

    /// Set flag into the INFO column.
    pub fn push_info_flag(&mut self, tag: &[u8]) -> Result<()> {
        self.push_info_string_impl(tag, &[b""], htslib::BCF_HT_FLAG)
    }

    /// Remove the flag from the INFO column.
    pub fn clear_info_flag(&mut self, tag: &[u8]) -> Result<()> {
        self.push_info_string_impl(tag, &[], htslib::BCF_HT_FLAG)
    }

    /// Add/replace a string-typed INFO entry.
    pub fn push_info_string(&mut self, tag: &[u8], data: &[&[u8]]) -> Result<()> {
        self.push_info_string_impl(tag, data, htslib::BCF_HT_STR)
    }

    /// Remove the string field from the INFO column.
    pub fn clear_info_string(&mut self, tag: &[u8]) -> Result<()> {
        self.push_info_string_impl(tag, &[], htslib::BCF_HT_STR)
    }

    /// Add an string-valued INFO tag.
    fn push_info_string_impl(&mut self, tag: &[u8], data: &[&[u8]], ht: u32) -> Result<()> {
        let mut buf: Vec<u8> = Vec::new();
        for (i, &s) in data.iter().enumerate() {
            if i > 0 {
                buf.extend(b",");
            }
            buf.extend(s);
        }
        let c_str = ffi::CString::new(buf).unwrap();
        let len = if ht == htslib::BCF_HT_FLAG {
            data.len()
        } else {
            c_str.to_bytes().len()
        };
        let tag_c_str = ffi::CString::new(tag).unwrap();
        unsafe {
            if htslib::bcf_update_info(
                self.header().inner,
                self.inner,
                tag_c_str.as_ptr() as *mut c_char,
                c_str.as_ptr() as *const ::std::os::raw::c_void,
                len as i32,
                ht as i32,
            ) == 0
            {
                Ok(())
            } else {
                Err(Error::BcfSetTag {
                    tag: str::from_utf8(tag).unwrap().to_owned(),
                })
            }
        }
    }

    /// Remove unused alleles.
    pub fn trim_alleles(&mut self) -> Result<()> {
        match unsafe { htslib::bcf_trim_alleles(self.header().inner, self.inner) } {
            -1 => Err(Error::BcfRemoveAlleles),
            _ => Ok(()),
        }
    }

    pub fn remove_alleles(&mut self, remove: &[bool]) -> Result<()> {
        let rm_set = unsafe { htslib::kbs_init(remove.len() as u64) };

        for (i, &r) in remove.iter().enumerate() {
            if r {
                unsafe {
                    htslib::kbs_insert(rm_set, i as i32);
                }
            }
        }

        let ret = unsafe { htslib::bcf_remove_allele_set(self.header().inner, self.inner, rm_set) };

        unsafe {
            htslib::kbs_destroy(rm_set);
        }

        match ret {
            -1 => Err(Error::BcfRemoveAlleles),
            _ => Ok(()),
        }
    }

    /// Get the length of the reference allele. If the record has no reference allele, then the
    /// result will be `0`.
    ///
    /// # Example
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// #
    /// # // Create minimal VCF header with a single sample
    /// # let mut header = Header::new();
    /// # header.push_sample("sample".as_bytes());
    /// #
    /// # // Write uncompressed VCF to stdout with above header and get an empty record
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// # assert_eq!(record.rlen(), 0);
    /// let alleles: &[&[u8]] = &[b"AGG", b"TG"];
    /// record.set_alleles(alleles).expect("Failed to set alleles");
    /// assert_eq!(record.rlen(), 3)
    /// ```
    pub fn rlen(&self) -> i64 {
        self.inner().rlen
    }

    /// Clear all parts of the record. Useful if you plan to reuse a record object multiple times.
    ///
    /// # Example
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// #
    /// # // Create minimal VCF header with a single sample
    /// # let mut header = Header::new();
    /// # header.push_sample("sample".as_bytes());
    /// #
    /// # // Write uncompressed VCF to stdout with above header and get an empty record
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// let alleles: &[&[u8]] = &[b"AGG", b"TG"];
    /// record.set_alleles(alleles).expect("Failed to set alleles");
    /// record.set_pos(6);
    /// record.clear();
    /// assert_eq!(record.rlen(), 0);
    /// assert_eq!(record.pos(), 0)
    /// ```
    pub fn clear(&self) {
        unsafe { htslib::bcf_clear(self.inner) }
    }

    /// Provide short description of record for locating it in the BCF/VCF file.
    pub fn desc(&self) -> String {
        if let Some(rid) = self.rid() {
            if let Ok(contig) = self.header.rid2name(rid) {
                return format!("{}:{}", str::from_utf8(contig).unwrap(), self.pos());
            }
        }
        "".to_owned()
    }
}

impl Clone for Record {
    fn clone(&self) -> Self {
        let inner = unsafe {
            let inner = htslib::bcf_dup(self.inner);
            inner
        };
        Record {
            inner,
            header: self.header.clone(),
        }
    }
}

impl genome::AbstractLocus for Record {
    fn contig(&self) -> &str {
        str::from_utf8(
            self.header()
                .rid2name(self.rid().expect("rid not set"))
                .expect("unable to find rid in header"),
        )
        .expect("unable to interpret contig name as UTF-8")
    }

    fn pos(&self) -> u64 {
        self.pos() as u64
    }
}

/// Phased or unphased alleles, represented as indices.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GenotypeAllele {
    Unphased(i32),
    Phased(i32),
    UnphasedMissing,
    PhasedMissing,
}

impl GenotypeAllele {
    /// Decode given integer according to BCF standard.
    #[deprecated(
        since = "0.36.0",
        note = "Please use the conversion trait From<i32> for GenotypeAllele instead."
    )]
    pub fn from_encoded(encoded: i32) -> Self {
        match (encoded, encoded & 1) {
            (0, 0) => GenotypeAllele::UnphasedMissing,
            (1, 1) => GenotypeAllele::PhasedMissing,
            (e, 1) => GenotypeAllele::Phased((e >> 1) - 1),
            (e, 0) => GenotypeAllele::Unphased((e >> 1) - 1),
            _ => panic!("unexpected phasing type"),
        }
    }

    /// Get the index into the list of alleles.
    pub fn index(self) -> Option<u32> {
        match self {
            GenotypeAllele::Unphased(i) | GenotypeAllele::Phased(i) => Some(i as u32),
            GenotypeAllele::UnphasedMissing | GenotypeAllele::PhasedMissing => None,
        }
    }
}

impl fmt::Display for GenotypeAllele {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.index() {
            Some(a) => write!(f, "{}", a),
            None => write!(f, "."),
        }
    }
}

impl From<GenotypeAllele> for i32 {
    fn from(allele: GenotypeAllele) -> i32 {
        let (allele, phased) = match allele {
            GenotypeAllele::UnphasedMissing => (-1, 0),
            GenotypeAllele::PhasedMissing => (-1, 1),
            GenotypeAllele::Unphased(a) => (a, 0),
            GenotypeAllele::Phased(a) => (a, 1),
        };
        (allele + 1) << 1 | phased
    }
}

impl From<i32> for GenotypeAllele {
    fn from(encoded: i32) -> GenotypeAllele {
        match (encoded, encoded & 1) {
            (0, 0) => GenotypeAllele::UnphasedMissing,
            (1, 1) => GenotypeAllele::PhasedMissing,
            (e, 1) => GenotypeAllele::Phased((e >> 1) - 1),
            (e, 0) => GenotypeAllele::Unphased((e >> 1) - 1),
            _ => panic!("unexpected phasing type"),
        }
    }
}

custom_derive! {
    /// Genotype representation as a vector of `GenotypeAllele`.
    #[derive(NewtypeDeref, Debug, Clone, PartialEq, Eq, Hash)]
    pub struct Genotype(Vec<GenotypeAllele>);
}

impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let &Genotype(ref alleles) = self;
        write!(f, "{}", alleles[0])?;
        for a in &alleles[1..] {
            let sep = match a {
                GenotypeAllele::Phased(_) | GenotypeAllele::PhasedMissing => '|',
                GenotypeAllele::Unphased(_) | GenotypeAllele::UnphasedMissing => '/',
            };
            write!(f, "{}{}", sep, a)?;
        }
        Ok(())
    }
}

/// Lazy representation of genotypes, that does no computation until a particular genotype is queried.
#[derive(Debug)]
pub struct Genotypes<'a, B>
where
    B: Borrow<Buffer> + 'a,
{
    encoded: BufferBacked<'a, Vec<&'a [i32]>, B>,
}

impl<'a, B: Borrow<Buffer> + 'a> Genotypes<'a, B> {
    /// Get genotype of ith sample. So far, only supports diploid genotypes.
    ///
    /// Note that the result complies with the BCF spec. This means that the
    /// first allele will always be marked as `Unphased`. That is, if you have 1|1 in the VCF,
    /// this method will return `[Unphased(1), Phased(1)]`.
    pub fn get(&self, i: usize) -> Genotype {
        let igt = self.encoded[i];
        Genotype(igt.iter().map(|&e| GenotypeAllele::from(e)).collect())
    }
}

impl Drop for Record {
    fn drop(&mut self) {
        unsafe { htslib::bcf_destroy(self.inner) };
    }
}

unsafe impl Send for Record {}

unsafe impl Sync for Record {}

/// Info tag representation.
#[derive(Debug)]
pub struct Info<'a, B: BorrowMut<Buffer> + Borrow<Buffer>> {
    record: &'a Record,
    tag: &'a [u8],
    buffer: B,
}

impl<'a, 'b, B: BorrowMut<Buffer> + Borrow<Buffer> + 'b> Info<'a, B> {
    /// Short description of info tag.
    pub fn desc(&self) -> String {
        str::from_utf8(self.tag).unwrap().to_owned()
    }

    fn data(&mut self, data_type: u32) -> Result<Option<i32>> {
        let mut n: i32 = self.buffer.borrow().len;
        let c_str = ffi::CString::new(self.tag).unwrap();
        let ret = unsafe {
            htslib::bcf_get_info_values(
                self.record.header().inner,
                self.record.inner,
                c_str.as_ptr() as *mut c_char,
                &mut self.buffer.borrow_mut().inner,
                &mut n,
                data_type as i32,
            )
        };
        self.buffer.borrow_mut().len = n;

        match ret {
            -1 => Err(Error::BcfUndefinedTag { tag: self.desc() }),
            -2 => Err(Error::BcfUnexpectedType { tag: self.desc() }),
            -3 => Ok(None),
            ret => Ok(Some(ret)),
        }
    }

    /// Get integers from tag. `None` if tag not present in record.
    ///
    /// Import `bcf::record::Numeric` for missing value handling.
    ///
    /// **Attention:** the returned BufferBacked which holds the data has to be kept in scope
    /// as along as the data is accessed. If parts of the data are accessed while
    /// the BufferBacked object is already dropped, you will access unallocated
    /// memory.
    pub fn integer(mut self) -> Result<Option<BufferBacked<'b, &'b [i32], B>>> {
        self.data(htslib::BCF_HT_INT).map(|data| {
            data.map(|ret| {
                let values = unsafe {
                    slice::from_raw_parts(self.buffer.borrow().inner as *const i32, ret as usize)
                };
                BufferBacked::new(&values[..ret as usize], self.buffer)
            })
        })
    }

    /// Get floats from tag. `None` if tag not present in record.
    ///
    /// Import `bcf::record::Numeric` for missing value handling.
    ///
    /// **Attention:** the returned BufferBacked which holds the data has to be kept in scope
    /// as along as the data is accessed. If parts of the data are accessed while
    /// the BufferBacked object is already dropped, you will access unallocated
    /// memory.
    pub fn float(mut self) -> Result<Option<BufferBacked<'b, &'b [f32], B>>> {
        self.data(htslib::BCF_HT_REAL).map(|data| {
            data.map(|ret| {
                let values = unsafe {
                    slice::from_raw_parts(self.buffer.borrow().inner as *const f32, ret as usize)
                };
                BufferBacked::new(&values[..ret as usize], self.buffer)
            })
        })
    }

    /// Get flags from tag. `false` if not set.
    pub fn flag(&mut self) -> Result<bool> {
        self.data(htslib::BCF_HT_FLAG).map(|data| match data {
            Some(ret) => ret == 1,
            None => false,
        })
    }

    /// Get strings from tag. `None` if tag not present in record.
    ///
    /// **Attention:** the returned BufferBacked which holds the data has to be kept in scope
    /// as along as the data is accessed. If parts of the data are accessed while
    /// the BufferBacked object is already dropped, you will access unallocated
    /// memory.
    pub fn string(mut self) -> Result<Option<BufferBacked<'b, Vec<&'b [u8]>, B>>> {
        self.data(htslib::BCF_HT_STR).map(|data| {
            data.map(|ret| {
                BufferBacked::new(
                    unsafe {
                        slice::from_raw_parts(self.buffer.borrow().inner as *const u8, ret as usize)
                    }
                    .split(|c| *c == b',')
                    .map(|s| {
                        // stop at zero character
                        s.split(|c| *c == 0u8)
                            .next()
                            .expect("Bug: returned string should not be empty.")
                    })
                    .collect(),
                    self.buffer,
                )
            })
        })
    }
}

unsafe impl<'a, 'b, B: BorrowMut<Buffer> + Borrow<Buffer> + 'b> Send for Info<'a, B> {}

unsafe impl<'a, 'b, B: BorrowMut<Buffer> + Borrow<Buffer> + 'b> Sync for Info<'a, B> {}

fn trim_slice<T: PartialEq + NumericUtils>(s: &[T]) -> &[T] {
    s.split(|v| v.is_vector_end())
        .next()
        .expect("Bug: returned slice should not be empty.")
}

// Representation of per-sample data.
#[derive(Debug)]
pub struct Format<'a, B: BorrowMut<Buffer> + Borrow<Buffer>> {
    record: &'a Record,
    tag: &'a [u8],
    inner: *mut htslib::bcf_fmt_t,
    buffer: B,
}

impl<'a, 'b, B: BorrowMut<Buffer> + Borrow<Buffer> + 'b> Format<'a, B> {
    /// Create new format data in a given record.
    fn new(record: &'a Record, tag: &'a [u8], buffer: B) -> Format<'a, B> {
        let c_str = ffi::CString::new(tag).unwrap();
        let inner = unsafe {
            htslib::bcf_get_fmt(
                record.header().inner,
                record.inner,
                c_str.as_ptr() as *mut c_char,
            )
        };
        Format {
            record,
            tag,
            inner,
            buffer,
        }
    }

    /// Provide short description of format entry (just the tag name).
    pub fn desc(&self) -> String {
        str::from_utf8(self.tag).unwrap().to_owned()
    }

    pub fn inner(&self) -> &htslib::bcf_fmt_t {
        unsafe { &*self.inner }
    }

    pub fn inner_mut(&mut self) -> &mut htslib::bcf_fmt_t {
        unsafe { &mut *self.inner }
    }

    fn values_per_sample(&self) -> usize {
        self.inner().n as usize
    }

    /// Read and decode format data into a given type.
    fn data(&mut self, data_type: u32) -> Result<i32> {
        let mut n: i32 = self.buffer.borrow().len;
        let c_str = ffi::CString::new(self.tag).unwrap();
        let ret = unsafe {
            htslib::bcf_get_format_values(
                self.record.header().inner,
                self.record.inner,
                c_str.as_ptr() as *mut c_char,
                &mut self.buffer.borrow_mut().inner,
                &mut n,
                data_type as i32,
            )
        };
        self.buffer.borrow_mut().len = n;
        match ret {
            -1 => Err(Error::BcfUndefinedTag { tag: self.desc() }),
            -2 => Err(Error::BcfUnexpectedType { tag: self.desc() }),
            -3 => Err(Error::BcfMissingTag {
                tag: self.desc(),
                record: self.record.desc(),
            }),
            ret => Ok(ret),
        }
    }

    /// Get format data as integers.
    ///
    /// **Attention:** the returned BufferBacked which holds the data has to be kept in scope
    /// as long as the data is accessed. If parts of the data are accessed while
    /// the BufferBacked object is already dropped, you will access unallocated
    /// memory.
    pub fn integer(mut self) -> Result<BufferBacked<'b, Vec<&'b [i32]>, B>> {
        self.data(htslib::BCF_HT_INT).map(|ret| {
            BufferBacked::new(
                unsafe {
                    slice::from_raw_parts(
                        self.buffer.borrow_mut().inner as *const i32,
                        ret as usize,
                    )
                }
                .chunks(self.values_per_sample())
                .map(|s| trim_slice(s))
                .collect(),
                self.buffer,
            )
        })
    }

    /// Get format data as floats.
    ///
    /// **Attention:** the returned BufferBacked which holds the data has to be kept in scope
    /// as along as the data is accessed. If parts of the data are accessed while
    /// the BufferBacked object is already dropped, you will access unallocated
    /// memory.
    pub fn float(mut self) -> Result<BufferBacked<'b, Vec<&'b [f32]>, B>> {
        self.data(htslib::BCF_HT_REAL).map(|ret| {
            BufferBacked::new(
                unsafe {
                    slice::from_raw_parts(
                        self.buffer.borrow_mut().inner as *const f32,
                        ret as usize,
                    )
                }
                .chunks(self.values_per_sample())
                .map(|s| trim_slice(s))
                .collect(),
                self.buffer,
            )
        })
    }

    /// Get format data as byte slices. To obtain the values strings, use `std::str::from_utf8`.
    ///
    /// **Attention:** the returned BufferBacked which holds the data has to be kept in scope
    /// as along as the data is accessed. If parts of the data are accessed while
    /// the BufferBacked object is already dropped, you will access unallocated
    /// memory.
    pub fn string(mut self) -> Result<BufferBacked<'b, Vec<&'b [u8]>, B>> {
        self.data(htslib::BCF_HT_STR).map(|ret| {
            BufferBacked::new(
                unsafe {
                    slice::from_raw_parts(self.buffer.borrow_mut().inner as *const u8, ret as usize)
                }
                .chunks(self.values_per_sample())
                .map(|s| {
                    // stop at zero character
                    s.split(|c| *c == 0u8)
                        .next()
                        .expect("Bug: returned string should not be empty.")
                })
                .collect(),
                self.buffer,
            )
        })
    }
}

unsafe impl<'a, 'b, B: BorrowMut<Buffer> + Borrow<Buffer> + 'b> Send for Format<'a, B> {}

unsafe impl<'a, 'b, B: BorrowMut<Buffer> + Borrow<Buffer> + 'b> Sync for Format<'a, B> {}

#[derive(Debug)]
pub struct Filters<'a> {
    /// Reference to the `Record` to enumerate records for.
    record: &'a Record,
    /// Index of the next filter to return, if not at end.
    idx: i32,
}

impl<'a> Filters<'a> {
    pub fn new(record: &'a Record) -> Self {
        Filters { record, idx: 0 }
    }
}

impl<'a> Iterator for Filters<'a> {
    type Item = Id;

    fn next(&mut self) -> Option<Id> {
        if self.record.inner().d.n_flt <= self.idx {
            None
        } else {
            let i = self.idx as isize;
            self.idx += 1;
            Some(Id(unsafe { *self.record.inner().d.flt.offset(i) } as u32))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bcf::{Format, Header, Writer};
    use tempfile::NamedTempFile;

    #[test]
    fn test_missing_float() {
        let expected: u32 = 0x7F80_0001;
        assert_eq!(MISSING_FLOAT.bits(), expected);
    }

    #[test]
    fn test_vector_end_float() {
        let expected: u32 = 0x7F80_0002;
        assert_eq!(VECTOR_END_FLOAT.bits(), expected);
    }

    #[test]
    fn test_record_rlen() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        assert_eq!(record.rlen(), 0);
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        assert_eq!(record.rlen(), 3)
    }

    #[test]
    fn test_record_end() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record.set_pos(5);

        assert_eq!(record.end(), 8)
    }

    #[test]
    fn test_record_clear() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        header.push_sample("sample".as_bytes());
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record.set_pos(6);
        record.clear();

        assert_eq!(record.rlen(), 0);
        assert_eq!(record.sample_count(), 0);
        assert_eq!(record.pos(), 0)
    }

    #[test]
    fn test_record_clone() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record.set_pos(6);

        let mut cloned_record = record.clone();
        cloned_record.set_pos(5);

        assert_eq!(record.pos(), 6);
        assert_eq!(record.allele_count(), 2);
        assert_eq!(cloned_record.pos(), 5);
        assert_eq!(cloned_record.allele_count(), 2);
    }

    #[test]
    fn test_record_has_filter_pass_is_default() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let record = vcf.empty_record();

        assert!(record.has_filter("PASS".as_bytes()));
        assert!(record.has_filter(".".as_bytes()));
        assert!(record.has_filter(&Id(0)));
        assert!(!record.has_filter("foo".as_bytes()));
        assert!(!record.has_filter(&Id(2)));
    }

    #[test]
    fn test_record_has_filter_custom() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        record.push_filter("foo".as_bytes()).unwrap();

        assert!(record.has_filter("foo".as_bytes()));
        assert!(!record.has_filter("PASS".as_bytes()))
    }

    #[test]
    fn test_record_push_filter() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
        header.push_record(br#"##FILTER=<ID=bar,Description="dranks">"#);
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        assert!(record.has_filter("PASS".as_bytes()));
        record.push_filter("foo".as_bytes()).unwrap();
        let bar = record.header().name_to_id(b"bar").unwrap();
        record.push_filter(&bar).unwrap();
        assert!(record.has_filter("foo".as_bytes()));
        assert!(record.has_filter(&bar));
        assert!(!record.has_filter("PASS".as_bytes()));
        assert!(record.push_filter("baz".as_bytes()).is_err())
    }

    #[test]
    fn test_record_set_filters() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
        header.push_record(br#"##FILTER=<ID=bar,Description="a horse walks into...">"#);
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        assert!(record.has_filter("PASS".as_bytes()));
        record
            .set_filters(&["foo".as_bytes(), "bar".as_bytes()])
            .unwrap();
        assert!(record.has_filter("foo".as_bytes()));
        assert!(record.has_filter("bar".as_bytes()));
        assert!(!record.has_filter("PASS".as_bytes()));
        let filters: &[&Id] = &[];
        record.set_filters(filters).unwrap();
        assert!(record.has_filter("PASS".as_bytes()));
        assert!(!record.has_filter("foo".as_bytes()));
        assert!(record
            .set_filters(&["foo".as_bytes(), "baz".as_bytes()])
            .is_err())
    }

    #[test]
    fn test_record_remove_filter() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
        header.push_record(br#"##FILTER=<ID=bar,Description="a horse walks into...">"#);
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let foo = record.header().name_to_id(b"foo").unwrap();
        let bar = record.header().name_to_id(b"bar").unwrap();
        record.set_filters(&[&foo, &bar]).unwrap();
        assert!(record.has_filter(&foo));
        assert!(record.has_filter(&bar));
        record.remove_filter(&foo, true).unwrap();
        assert!(!record.has_filter(&foo));
        assert!(record.has_filter(&bar));
        assert!(record.remove_filter("baz".as_bytes(), true).is_err());
        record.remove_filter(&bar, true).unwrap();
        assert!(!record.has_filter(&bar));
        assert!(record.has_filter("PASS".as_bytes()));
    }
}
