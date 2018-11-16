// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::f32;
use std::ffi;
use std::fmt;
use std::i32;
use std::ptr;
use std::rc::Rc;
use std::slice;

use ieee754::Ieee754;
use itertools::Itertools;

use bcf::header::{HeaderView, Id};
use htslib;

const MISSING_INTEGER: i32 = i32::MIN;
const VECTOR_END_INTEGER: i32 = i32::MIN + 1;
lazy_static! {
    static ref MISSING_FLOAT: f32 = Ieee754::from_bits(0x7F800001);
    static ref VECTOR_END_FLOAT: f32 = Ieee754::from_bits(0x7F800002);
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

/// A BCF record.
/// New records can be created by the `empty_record` methods of `bcf::Reader` and `bcf::Writer`.
#[derive(Debug)]
pub struct Record {
    pub inner: *mut htslib::bcf1_t,
    header: Rc<HeaderView>,
    buffer: *mut ::std::os::raw::c_void,
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
        Record {
            inner: inner,
            header: header,
            buffer: ptr::null_mut(),
        }
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
    /// To look up the contig name, use `bcf::header::HeaderView::rid2name`.
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

    // Update the internal reference ID number.
    pub fn set_rid(&mut self, rid: &Option<u32>) {
        match rid {
            &Some(rid) => self.inner_mut().rid = rid as i32,
            &None => self.inner_mut().rid = -1,
        }
    }

    // Return 0-based position.
    pub fn pos(&self) -> u32 {
        self.inner().pos as u32
    }

    /// Set 0-based position.
    pub fn set_pos(&mut self, pos: i32) {
        self.inner_mut().pos = pos;
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
    pub fn set_id(&mut self, id: &[u8]) -> Result<(), IdWriteError> {
        if unsafe {
            htslib::bcf_update_id(
                self.header().inner,
                self.inner,
                ffi::CString::new(id).unwrap().as_ptr() as *mut i8,
            )
        } == 0
        {
            Ok(())
        } else {
            Err(IdWriteError::Some)
        }
    }

    /// Clear the ID column (set it to `"."`).
    pub fn clear_id(&mut self) -> Result<(), IdWriteError> {
        if unsafe {
            htslib::bcf_update_id(
                self.header().inner,
                self.inner,
                ffi::CString::new(".".as_bytes()).unwrap().as_ptr() as *mut i8,
            )
        } == 0
        {
            Ok(())
        } else {
            Err(IdWriteError::Some)
        }
    }

    /// Add the ID string (the ID field is semicolon-separated), checking for duplicates.
    pub fn push_id(&mut self, id: &[u8]) -> Result<(), IdWriteError> {
        if unsafe {
            htslib::bcf_add_id(
                self.header().inner,
                self.inner,
                ffi::CString::new(id).unwrap().as_ptr() as *mut i8,
            )
        } == 0
        {
            Ok(())
        } else {
            Err(IdWriteError::Some)
        }
    }

    /// Return `Filters` iterator for enumerating all filters that have been set.
    ///
    /// A record having the `PASS` filter will return an empty `Filter` here.
    pub fn filters(&self) -> Filters {
        Filters::new(self)
    }

    /// Query whether the filter with the given ID has been set.
    ///
    /// # Arguments
    ///
    /// - `flt_id` - The filter ID to query for.
    pub fn has_filter(&self, flt_id: &Id) -> bool {
        if **flt_id == 0 && self.inner().d.n_flt == 0 {
            return true;
        }
        for i in 0..(self.inner().d.n_flt as isize) {
            if unsafe { *self.inner().d.flt.offset(i) } == **flt_id as i32 {
                return true;
            }
        }
        false
    }

    /// Set the given filters IDs to the FILTER column.
    ///
    /// Setting an empty slice removes all filters.
    ///
    /// # Arguments
    ///
    /// - `flt_ids` - The identifiers of the filter values to set.
    pub fn set_filters(&mut self, flt_ids: &[Id]) {
        let mut flt_ids: Vec<i32> = flt_ids.iter().map(|x| **x as i32).collect();
        unsafe {
            htslib::bcf_update_filter(
                self.header().inner,
                self.inner,
                flt_ids.as_mut_ptr(),
                flt_ids.len() as i32,
            );
        }
    }

    /// Add the given filter to the FILTER column.
    ///
    /// If `val` corresponds to `"PASS"` then all existing filters are removed first. If other than
    /// `"PASS"`, then existing `"PASS"` is removed.
    ///
    /// # Arguments
    ///
    /// - `flt_id` - The corresponding filter ID value to add.
    pub fn push_filter(&mut self, flt_id: Id) {
        unsafe {
            htslib::bcf_add_filter(self.header().inner, self.inner, *flt_id as i32);
        }
    }

    /// Remove the given filter from the FILTER column.
    ///
    /// # Arguments
    ///
    /// - `val` - The corresponding filter ID to remove.
    /// - `pass_on_empty` - Set to "PASS" when removing the last value.
    pub fn remove_filter(&mut self, flt_id: Id, pass_on_empty: bool) {
        unsafe {
            htslib::bcf_remove_filter(
                self.header().inner,
                self.inner,
                *flt_id as i32,
                pass_on_empty as i32,
            );
        }
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

    /// Set alleles.
    pub fn set_alleles(&mut self, alleles: &[&[u8]]) -> Result<(), AlleleWriteError> {
        let cstrings: Vec<ffi::CString> = alleles
            .iter()
            .map(|vec| ffi::CString::new(*vec).unwrap())
            .collect();
        let mut ptrs: Vec<*const i8> = cstrings
            .iter()
            .map(|cstr| cstr.as_ptr() as *const i8)
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
            Err(AlleleWriteError::Some)
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

    /// Get the value of the given info tag.
    pub fn info<'a>(&'a mut self, tag: &'a [u8]) -> Info {
        Info {
            record: self,
            tag: tag,
        }
    }

    /// Get the number of samples.
    pub fn sample_count(&self) -> u32 {
        self.inner().n_sample()
    }

    /// Get the number of alleles, including reference allele.
    pub fn allele_count(&self) -> u32 {
        self.inner().n_allele()
    }

    // TODO fn push_genotypes(&mut self, Genotypes) {}?

    /// Get genotypes as vector of one `Genotype` per sample.
    pub fn genotypes(&mut self) -> Result<Genotypes, FormatReadError> {
        Ok(Genotypes {
            encoded: try!(self.format(b"GT").integer()),
        })
    }

    /// Get the value of the given format tag for each sample.
    pub fn format<'a>(&'a mut self, tag: &'a [u8]) -> Format {
        Format::new(self, tag)
    }

    /// Add an integer-typed FORMAT tag.
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
    pub fn push_format_integer(&mut self, tag: &[u8], data: &[i32]) -> Result<(), TagWriteError> {
        self.push_format(tag, data, htslib::BCF_HT_INT)
    }

    /// Add a float-typed FORMAT tag.
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
    pub fn push_format_float(&mut self, tag: &[u8], data: &[f32]) -> Result<(), TagWriteError> {
        self.push_format(tag, data, htslib::BCF_HT_REAL)
    }

    /// Add a char-typed FORMAT tag.
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
    pub fn push_format_char(&mut self, tag: &[u8], data: &[u8]) -> Result<(), TagWriteError> {
        self.push_format(tag, data, htslib::BCF_HT_STR)
    }

    /// Add a format tag. Data is a flattened two-dimensional array.
    /// The first dimension contains one array for each sample.
    fn push_format<T>(&mut self, tag: &[u8], data: &[T], ht: u32) -> Result<(), TagWriteError> {
        assert!(data.len() > 0);
        unsafe {
            if htslib::bcf_update_format(
                self.header().inner,
                self.inner,
                ffi::CString::new(tag).unwrap().as_ptr() as *mut i8,
                data.as_ptr() as *const ::std::os::raw::c_void,
                data.len() as i32,
                ht as i32,
            ) == 0
            {
                Ok(())
            } else {
                Err(TagWriteError::Some)
            }
        }
    }

    // TODO: should we add convenience methods clear_format_*?

    /// Add a string-typed FORMAT tag.
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
    pub fn push_format_string(&mut self, tag: &[u8], data: &[&[u8]]) -> Result<(), TagWriteError> {
        let c_data = data
            .iter()
            .map(|&s| ffi::CString::new(s).unwrap())
            .collect::<Vec<ffi::CString>>();
        let c_ptrs = c_data
            .iter()
            .map(|s| s.as_ptr() as *mut i8)
            .collect::<Vec<*mut i8>>();
        unsafe {
            if htslib::bcf_update_format_string(
                self.header().inner,
                self.inner,
                ffi::CString::new(tag).unwrap().as_ptr() as *mut i8,
                c_ptrs.as_slice().as_ptr() as *mut *const i8,
                data.len() as i32,
            ) == 0
            {
                Ok(())
            } else {
                Err(TagWriteError::Some)
            }
        }
    }

    /// Add an integer-typed INFO entry.
    pub fn push_info_integer(&mut self, tag: &[u8], data: &[i32]) -> Result<(), TagWriteError> {
        self.push_info(tag, data, htslib::BCF_HT_INT)
    }

    /// Remove the integer-typed INFO entry.
    pub fn clear_info_integer(&mut self, tag: &[u8]) -> Result<(), TagWriteError> {
        self.push_info::<i32>(tag, &[], htslib::BCF_HT_INT)
    }

    /// Add a float-typed INFO entry.
    pub fn push_info_float(&mut self, tag: &[u8], data: &[f32]) -> Result<(), TagWriteError> {
        self.push_info(tag, data, htslib::BCF_HT_REAL)
    }

    /// Remove the float-typed INFO entry.
    pub fn clear_info_float(&mut self, tag: &[u8]) -> Result<(), TagWriteError> {
        self.push_info::<u8>(tag, &[], htslib::BCF_HT_REAL)
    }

    /// Add a not INFO tag.
    fn push_info<T>(&mut self, tag: &[u8], data: &[T], ht: u32) -> Result<(), TagWriteError> {
        assert!(data.len() > 0);
        unsafe {
            if htslib::bcf_update_info(
                self.header().inner,
                self.inner,
                ffi::CString::new(tag).unwrap().as_ptr() as *mut i8,
                data.as_ptr() as *const ::std::os::raw::c_void,
                data.len() as i32,
                ht as i32,
            ) == 0
            {
                Ok(())
            } else {
                Err(TagWriteError::Some)
            }
        }
    }

    /// Set flag into the INFO column.
    pub fn push_info_flag(&mut self, tag: &[u8]) -> Result<(), TagWriteError> {
        self.push_info_string_impl(tag, &["".as_bytes()], htslib::BCF_HT_FLAG)
    }

    /// Remove the flag from the INFO column.
    pub fn clear_info_flag(&mut self, tag: &[u8]) -> Result<(), TagWriteError> {
        self.push_info_string_impl(tag, &[], htslib::BCF_HT_FLAG)
    }

    /// Add a string-typed INFO entry.
    pub fn push_info_string(&mut self, tag: &[u8], data: &[&[u8]]) -> Result<(), TagWriteError> {
        self.push_info_string_impl(tag, data, htslib::BCF_HT_STR)
    }

    /// Remove the string field from the INFO column.
    pub fn clear_info_string(&mut self, tag: &[u8]) -> Result<(), TagWriteError> {
        self.push_info_string_impl(tag, &[], htslib::BCF_HT_STR)
    }

    /// Add an string-valued INFO tag.
    fn push_info_string_impl(
        &mut self,
        tag: &[u8],
        data: &[&[u8]],
        ht: u32,
    ) -> Result<(), TagWriteError> {
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
        unsafe {
            if htslib::bcf_update_info(
                self.header().inner,
                self.inner,
                ffi::CString::new(tag).unwrap().as_ptr() as *mut i8,
                c_str.as_ptr() as *const ::std::os::raw::c_void,
                len as i32,
                ht as i32,
            ) == 0
            {
                Ok(())
            } else {
                Err(TagWriteError::Some)
            }
        }
    }

    /// Remove unused alleles.
    pub fn trim_alleles(&mut self) -> Result<(), RemoveAllelesError> {
        match unsafe { htslib::bcf_trim_alleles(self.header().inner, self.inner) } {
            -1 => Err(RemoveAllelesError::Some),
            _ => Ok(()),
        }
    }

    pub fn remove_alleles(&mut self, remove: &[bool]) -> Result<(), RemoveAllelesError> {
        let rm_set = unsafe { htslib::kbs_init(remove.len()) };

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
            -1 => Err(RemoveAllelesError::Some),
            _ => Ok(()),
        }
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
    pub fn index(&self) -> Option<u32> {
        match self {
            &GenotypeAllele::Unphased(i) => Some(i as u32),
            &GenotypeAllele::Phased(i) => Some(i as u32),
            &GenotypeAllele::UnphasedMissing => None,
            &GenotypeAllele::PhasedMissing => None,
        }
    }
}

impl fmt::Display for GenotypeAllele {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.index() {
            Some(a) => write!(f, "{}", a),
            None => write!(f, "."),
        }
    }
}

custom_derive! {
    /// Genotype representation as a vector of `GenotypeAllele`.
    #[derive(NewtypeDeref, Debug, Clone, PartialEq, Eq, Hash)]
    pub struct Genotype(Vec<GenotypeAllele>);
}

impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let &Genotype(ref alleles) = self;
        try!(write!(f, "{}", alleles[0]));
        for a in &alleles[1..] {
            let sep = match a {
                &GenotypeAllele::Phased(_) => '|',
                &GenotypeAllele::Unphased(_) => '/',
                &GenotypeAllele::UnphasedMissing => '/',
                &GenotypeAllele::PhasedMissing => '|',
            };
            try!(write!(f, "{}{}", sep, a));
        }
        Ok(())
    }
}

/// Lazy representation of genotypes, that does no computation until a particular genotype is queried.
#[derive(Debug, Clone)]
pub struct Genotypes<'a> {
    encoded: Vec<&'a [i32]>,
}

impl<'a> Genotypes<'a> {
    /// Get genotype of ith sample. So far, only supports diploid genotypes.
    ///
    /// Note that the result complies with the BCF spec. This means that the
    /// first allele will always be marked as `Unphased`. That is, if you have 1|1 in the VCF,
    /// this method will return `[Unphased(1), Phased(1)]`.
    pub fn get(&self, i: usize) -> Genotype {
        let igt = self.encoded[i];
        let gt = Genotype(
            igt.into_iter()
                .map(|&e| GenotypeAllele::from_encoded(e))
                .collect_vec(),
        );
        gt
    }
}

impl Drop for Record {
    fn drop(&mut self) {
        if !self.buffer.is_null() {
            unsafe { ::libc::free(self.buffer as *mut ::libc::c_void) };
        }
        unsafe { htslib::bcf_destroy(self.inner) };
    }
}

unsafe impl Send for Record {}
unsafe impl Sync for Record {}

/// Info tag representation.
#[derive(Debug)]
pub struct Info<'a> {
    record: &'a mut Record,
    tag: &'a [u8],
}

impl<'a> Info<'a> {
    fn data(&mut self, data_type: u32) -> Result<Option<(usize, i32)>, InfoReadError> {
        let mut n: i32 = 0;
        match unsafe {
            htslib::bcf_get_info_values(
                self.record.header().inner,
                self.record.inner,
                ffi::CString::new(self.tag).unwrap().as_ptr() as *mut i8,
                &mut self.record.buffer,
                &mut n,
                data_type as i32,
            )
        } {
            -1 => Err(InfoReadError::UndefinedTag),
            -2 => Err(InfoReadError::UnexpectedType),
            -3 => Ok(None),
            ret => Ok(Some((n as usize, ret))),
        }
    }

    /// Get integers from tag. `None` if tag not present in record.
    ///
    /// Import `bcf::record::Numeric` for missing value handling.
    pub fn integer(&mut self) -> Result<Option<&'a [i32]>, InfoReadError> {
        self.data(htslib::BCF_HT_INT).map(|data| {
            data.map(|(n, _)| {
                trim_slice(unsafe { slice::from_raw_parts(self.record.buffer as *const i32, n) })
            })
        })
    }

    /// Get floats from tag. `None` if tag not present in record.
    ///
    /// Import `bcf::record::Numeric` for missing value handling.
    pub fn float(&mut self) -> Result<Option<&'a [f32]>, InfoReadError> {
        self.data(htslib::BCF_HT_REAL).map(|data| {
            data.map(|(n, _)| {
                trim_slice(unsafe { slice::from_raw_parts(self.record.buffer as *const f32, n) })
            })
        })
    }

    /// Get flags from tag. `false` if not set.
    pub fn flag(&mut self) -> Result<bool, InfoReadError> {
        self.data(htslib::BCF_HT_FLAG).map(|data| match data {
            Some((_, ret)) => ret == 1,
            None => false,
        })
    }

    /// Get strings from tag. `None` if tag not present in record.
    pub fn string(&mut self) -> Result<Option<Vec<&'a [u8]>>, InfoReadError> {
        self.data(htslib::BCF_HT_STR).map(|data| {
            data.map(|(n, ret)| {
                unsafe { slice::from_raw_parts(self.record.buffer as *const u8, ret as usize) }
                    .chunks(n)
                    .map(|s| {
                        // stop at zero character
                        s.split(|c| *c == 0u8)
                            .next()
                            .expect("Bug: returned string should not be empty.")
                    }).collect()
            })
        })
    }
}

unsafe impl<'a> Send for Info<'a> {}
unsafe impl<'a> Sync for Info<'a> {}

fn trim_slice<T: PartialEq + NumericUtils>(s: &[T]) -> &[T] {
    s.split(|v| v.is_vector_end())
        .next()
        .expect("Bug: returned slice should not be empty.")
}

// Representation of per-sample data.
#[derive(Debug)]
pub struct Format<'a> {
    record: &'a mut Record,
    tag: &'a [u8],
    inner: *mut htslib::bcf_fmt_t,
}

impl<'a> Format<'a> {
    /// Create new format data in a given record.
    fn new(record: &'a mut Record, tag: &'a [u8]) -> Format<'a> {
        let inner = unsafe {
            htslib::bcf_get_fmt(
                record.header().inner,
                record.inner,
                ffi::CString::new(tag).unwrap().as_ptr() as *mut i8,
            )
        };
        Format {
            record: record,
            tag: tag,
            inner: inner,
        }
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
    fn data(&mut self, data_type: u32) -> Result<(usize, i32), FormatReadError> {
        let mut n: i32 = 0;
        match unsafe {
            htslib::bcf_get_format_values(
                self.record.header().inner,
                self.record.inner,
                ffi::CString::new(self.tag).unwrap().as_ptr() as *mut i8,
                &mut self.record.buffer,
                &mut n,
                data_type as i32,
            )
        } {
            -1 => Err(FormatReadError::UndefinedTag),
            -2 => Err(FormatReadError::UnexpectedType),
            -3 => Err(FormatReadError::MissingTag),
            ret => Ok((n as usize, ret)),
        }
    }

    /// Get format data as integers.
    pub fn integer(&mut self) -> Result<Vec<&'a [i32]>, FormatReadError> {
        self.data(htslib::BCF_HT_INT).map(|(n, _)| {
            unsafe { slice::from_raw_parts(self.record.buffer as *const i32, n) }
                .chunks(self.values_per_sample())
                .map(|s| trim_slice(s))
                .collect()
        })
    }

    /// Get format data as floats.
    pub fn float(&mut self) -> Result<Vec<&'a [f32]>, FormatReadError> {
        self.data(htslib::BCF_HT_REAL).map(|(n, _)| {
            unsafe { slice::from_raw_parts(self.record.buffer as *const f32, n) }
                .chunks(self.values_per_sample())
                .map(|s| trim_slice(s))
                .collect()
        })
    }

    /// Get format data as byte slices. To obtain the values strings, use `std::str::from_utf8`.
    pub fn string(&mut self) -> Result<Vec<&'a [u8]>, FormatReadError> {
        self.data(htslib::BCF_HT_STR).map(|(n, _)| {
            unsafe { slice::from_raw_parts(self.record.buffer as *const u8, n) }
                .chunks(self.values_per_sample())
                .map(|s| {
                    // stop at zero character
                    s.split(|c| *c == 0u8)
                        .next()
                        .expect("Bug: returned string should not be empty.")
                }).collect()
        })
    }
}

unsafe impl<'a> Send for Format<'a> {}
unsafe impl<'a> Sync for Format<'a> {}

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

quick_error! {
    #[derive(Debug, Clone)]
    pub enum IterFilterError {
        Some {
            description("problem enumerating FILTER entries")
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum InfoReadError {
        UndefinedTag {
            description("tag undefined in header")
        }
        UnexpectedType {
            description("tag type differs from header definition")
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum FormatReadError {
        UndefinedTag {
            description("tag undefined in header")
        }
        UnexpectedType {
            description("tag type differs from header definition")
        }
        MissingTag {
            description("tag missing from record")
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum TagWriteError {
        Some {
            description("error writing tag to record")
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum IdWriteError {
        Some {
            description("error writing ID to record")
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum AlleleWriteError {
        Some {
            description("error writing alleles to record")
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum FilterWriteError {
        Some {
            description("error writing filters to record")
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum RemoveAllelesError {
        Some {
            description("error trimming alleles")
        }
    }
}
