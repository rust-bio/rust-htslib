// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::ptr;
use std::slice;
use std::ffi;
use std::i32;
use std::f32;
use std::fmt;

use ieee754::Ieee754;
use itertools::Itertools;

use htslib;

const MISSING_INTEGER: i32 = i32::MIN;
const VECTOR_END_INTEGER: i32 = i32::MIN + 1;
lazy_static!{
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
pub struct Record {
    pub inner: *mut htslib::vcf::bcf1_t,
    pub header: *mut htslib::vcf::bcf_hdr_t,
    buffer: *mut ::libc::c_void,
}


impl Record {
    pub fn new() -> Self {
        let inner = unsafe { htslib::vcf::bcf_init() };
        Record { inner: inner, header: ptr::null_mut(), buffer: ptr::null_mut() }
    }

    pub fn inner(&self) -> &htslib::vcf::bcf1_t {
        unsafe { &*self.inner }
    }

    pub fn inner_mut(&mut self) -> &mut htslib::vcf::bcf1_t {
        unsafe { &mut *self.inner }
    }

    pub fn rid(&self) -> Option<u32> {
        match self.inner().rid {
            -1  => None,
            rid => Some(rid as u32)
        }
    }

    // 0-based position.
    pub fn pos(&self) -> u32 {
        self.inner().pos as u32
    }


    /// Set 0-based position.
    pub fn set_pos(&mut self, pos: i32) {
        self.inner_mut().pos = pos;
    }

    /// Get alleles. The first allele is the reference allele.
    pub fn alleles(&self) -> Vec<&[u8]> {
        unsafe { htslib::vcf::bcf_unpack(self.inner, htslib::vcf::BCF_UN_STR) };
        let n = self.inner().n_allele as usize;
        let dec = self.inner().d;
        let alleles = unsafe { slice::from_raw_parts(dec.allele, n) };
        (0..n).map(|i| unsafe { ffi::CStr::from_ptr(alleles[i]).to_bytes() }).collect()
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
        Info { record: self, tag: tag }
    }

    /// Get the number of samples.
    pub fn sample_count(&self) -> u32 {
        self.inner().n_fmt_n_sample >> 8
    }

    /// Get the number of alleles, including reference allele.
    pub fn allele_count(&self) -> u16 {
        self.inner().n_allele
    }

    /// Get genotypes as vector of one `Genotype` per sample.
    pub fn genotypes(&mut self) -> Result<Genotypes, FormatReadError> {
        Ok(Genotypes {
            encoded: try!(self.format(b"GT").integer())
        })
    }

    /// Get the value of the given format tag for each sample.
    pub fn format<'a>(&'a mut self, tag: &'a [u8]) -> Format {
        Format::new(self, tag)
    }

    /// Add an integer format tag. Data is a flattened two-dimensional array.
    /// The first dimension contains one array for each sample.
    /// Returns error if tag is not present in header.
    pub fn push_format_integer(&mut self, tag: &[u8], data: &[i32]) -> Result<(), TagWriteError> {
        self.push_format(tag, data, htslib::vcf::BCF_HT_INT)
    }

    /// Add a float format tag. Data is a flattened two-dimensional array.
    /// The first dimension contains one array for each sample.
    /// Returns error if tag is not present in header.
    pub fn push_format_float(&mut self, tag: &[u8], data: &[f32]) -> Result<(), TagWriteError> {
        self.push_format(tag, data, htslib::vcf::BCF_HT_REAL)
    }

    /// Add a format tag. Data is a flattened two-dimensional array.
    /// The first dimension contains one array for each sample.
    fn push_format<T>(&mut self, tag: &[u8], data: &[T], ht: i32) -> Result<(), TagWriteError> {
        assert!(data.len() > 0);
        unsafe {
            if htslib::vcf::bcf_update_format(
                self.header,
                self.inner,
                ffi::CString::new(tag).unwrap().as_ptr() as *mut i8,
                data.as_ptr() as *const ::libc::c_void,
                data.len() as i32,
                ht
            ) == 0 {
                Ok(())
            }
            else {
                Err(TagWriteError::Some)
            }
        }
    }

    /// Add an integer info tag.
    pub fn push_info_integer(&mut self, tag: &[u8], data: &[i32]) -> Result<(), TagWriteError> {
        self.push_info(tag, data, htslib::vcf::BCF_HT_INT)
    }

    /// Add a float info tag.
    pub fn push_info_float(&mut self, tag: &[u8], data: &[f32]) -> Result<(), TagWriteError> {
        self.push_info(tag, data, htslib::vcf::BCF_HT_REAL)
    }

    /// Add an info tag.
    pub fn push_info<T>(&mut self, tag: &[u8], data: &[T], ht: i32) -> Result<(), TagWriteError> {
        assert!(data.len() > 0);
        unsafe {
            if htslib::vcf::bcf_update_info(
                self.header,
                self.inner,
                ffi::CString::new(tag).unwrap().as_ptr() as *mut i8,
                data.as_ptr() as *const ::libc::c_void,
                data.len() as i32,
                ht
            ) == 0 {
                Ok(())
            }
            else {
                Err(TagWriteError::Some)
            }
        }
    }

    /// Remove unused alleles.
    pub fn trim_alleles(&mut self) -> Result<(), TrimAllelesError> {
        match unsafe { htslib::vcfutils::bcf_trim_alleles(self.header, self.inner) } {
            -1 => Err(TrimAllelesError::Some),
            _  => Ok(())
        }
    }
}


/// Phased or unphased alleles, represented as indices.
#[derive(Debug)]
pub enum GenotypeAllele {
    Unphased(i32),
    Phased(i32),
    UnphasedMissing,
    PhasedMissing
}


impl GenotypeAllele {
    /// Decode given integer according to BCF standard.
    pub fn from_encoded(encoded: i32) -> Self {
        match (encoded, encoded & 1) {
            (0, 0) => GenotypeAllele::UnphasedMissing,
            (1, 1) => GenotypeAllele::PhasedMissing,
            (e, 1) => GenotypeAllele::Phased((e >> 1) - 1),
            (e, 0) => GenotypeAllele::Unphased((e >> 1) - 1),
            _ => panic!("unexpected phasing type")
        }
    }

    /// Get the index into the list of alleles.
    pub fn index(&self) -> Option<u32> {
        match self {
            &GenotypeAllele::Unphased(i) => Some(i as u32),
            &GenotypeAllele::Phased(i) => Some(i as u32),
            &GenotypeAllele::UnphasedMissing => None,
            &GenotypeAllele::PhasedMissing => None
        }
    }
}


impl fmt::Display for GenotypeAllele {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.index() {
            Some(a) => write!(f, "{}", a),
            None => write!(f, ".")
        }
    }
}


custom_derive! {
    /// Genotype representation as a vector of `GenotypeAllele`.
    #[derive(NewtypeDeref, Debug)]
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
                &GenotypeAllele::PhasedMissing => '|'
            };
            try!(write!(f, "{}{}", sep, a));
        }
        Ok(())
    }
}


/// Lazy representation of genotypes, that does no computation until a particular genotype is queried.
#[derive(Debug)]
pub struct Genotypes<'a> {
    encoded: Vec<&'a [i32]>
}


impl<'a> Genotypes<'a> {

    /// Get genotype of ith sample. So far, only supports diploid genotypes.
    ///
    /// Note that the result complies with the BCF spec. This means that the
    /// first allele will always be marked as `Unphased`. That is, if you have 1|1 in the VCF,
    /// this method will return `[Unphased(1), Phased(1)]`.
    pub fn get(&self, i: usize) -> Genotype {
        let igt = self.encoded[i];
        let gt = Genotype(igt.into_iter().map(|&e| GenotypeAllele::from_encoded(e)).collect_vec());
        gt
    }
}


impl Drop for Record {
    fn drop(&mut self) {
        if !self.buffer.is_null() {
            unsafe { ::libc::free(self.buffer) };
        }
        unsafe { htslib::vcf::bcf_destroy(self.inner) };
    }
}


unsafe impl Send for Record {}
unsafe impl Sync for Record {}


pub struct Info<'a> {
    record: &'a mut Record,
    tag: &'a [u8],
}


impl<'a> Info<'a> {
    fn data(&mut self, data_type: i32) -> Result<Option<(usize, i32)>, InfoReadError> {
        let mut n: i32 = 0;
        match unsafe {
            htslib::vcf::bcf_get_info_values(
                self.record.header,
                self.record.inner,
                ffi::CString::new(self.tag).unwrap().as_ptr() as *mut i8,
                &mut self.record.buffer,
                &mut n,
                data_type
            )
        } {
            -1 => Err(InfoReadError::UndefinedTag),
            -2 => Err(InfoReadError::UnexpectedType),
            -3 => Ok(None),
            ret  => Ok(Some((n as usize, ret))),
        }
    }

    /// Get integers from tag. `None` if tag not present in record.
    /// Import `bcf::record::Numeric` for missing value handling.
    pub fn integer(&mut self) -> Result<Option<&'a [i32]>, InfoReadError> {
        self.data(htslib::vcf::BCF_HT_INT).map(|data| data.map(|(n, _)| {
            trim_slice(
                unsafe { slice::from_raw_parts(self.record.buffer as *const i32, n) }
            )
        }))
    }

    /// Get mutable integers from tag. `None` if tag not present in record.
    /// Import `bcf::record::Numeric` for missing value handling.
    pub fn integer_mut(&mut self) -> Result<Option<&'a mut [i32]>, InfoReadError> {
        self.data(htslib::vcf::BCF_HT_INT).map(|data| data.map(|(n, _)| {
            unsafe { slice::from_raw_parts_mut(self.record.buffer as *mut i32, n) }
        }))
    }

    /// Get floats from tag. `None` if tag not present in record.
    /// Import `bcf::record::Numeric` for missing value handling.
    pub fn float(&mut self) -> Result<Option<&'a [f32]>, InfoReadError> {
        self.data(htslib::vcf::BCF_HT_REAL).map(|data| data.map(|(n, _)| {
            trim_slice(
                unsafe { slice::from_raw_parts(self.record.buffer as *const f32, n) }
            )
        }))
    }

    /// Get mutable floats from tag. `None` if tag not present in record.
    /// Import `bcf::record::Numeric` for missing value handling.
    pub fn float_mut(&mut self) -> Result<Option<&'a mut [f32]>, InfoReadError> {
        self.data(htslib::vcf::BCF_HT_REAL).map(|data| data.map(|(n, _)| {
            unsafe { slice::from_raw_parts_mut(self.record.buffer as *mut f32, n) }
        }))
    }

    pub fn flag(&mut self) -> Result<bool, InfoReadError> {
        self.data(htslib::vcf::BCF_HT_FLAG).map(|data| {
            match data {
                Some((_, ret)) => ret == 1,
                None => false
            }
        })
    }

    /// Get strings from tag. `None` if tag not present in record.
    pub fn string(&mut self) -> Result<Option<Vec<&'a [u8]>>, InfoReadError> {
        self.data(htslib::vcf::BCF_HT_STR).map(|data| data.map(|(n, ret)| {
            unsafe {
                slice::from_raw_parts(self.record.buffer as *const u8, ret as usize)
            }.chunks(n).map(|s| {
                // stop at zero character
                s.split(|c| *c == 0u8).next().expect("Bug: returned string should not be empty.")
            }).collect()
        }))
    }

    /// Get mutable strings from tag. `None` if tag not present in record.
    pub fn string_mut(&mut self) -> Result<Option<Vec<&'a mut [u8]>>, InfoReadError> {
        self.data(htslib::vcf::BCF_HT_STR).map(|data| data.map(|(n, ret)| {
            unsafe {
                slice::from_raw_parts_mut(self.record.buffer as *mut u8, ret as usize)
            }.chunks_mut(n).collect()
        }))
    }
}


unsafe impl<'a> Send for Info<'a> {}
unsafe impl<'a> Sync for Info<'a> {}


fn trim_slice<T: PartialEq + NumericUtils>(s: &[T]) -> &[T] {
    s.split(|v| v.is_vector_end()).next().expect("Bug: returned slice should not be empty.")
}


// TODO implement format.
pub struct Format<'a> {
    record: &'a mut Record,
    tag: &'a [u8],
    inner: *mut htslib::vcf::bcf_fmt_t,
}


impl<'a> Format<'a> {
    /// Create new format data in a given record.
    fn new(record: &'a mut Record, tag: &'a [u8]) -> Format<'a> {
        let inner = unsafe { htslib::vcf::bcf_get_fmt(
            record.header,
            record.inner,
            ffi::CString::new(tag).unwrap().as_ptr() as *mut i8
        ) };
        Format { record: record, tag: tag, inner: inner }
    }

    pub fn inner(&self) -> &htslib::vcf::bcf_fmt_t {
        unsafe { &*self.inner }
    }

    pub fn inner_mut(&mut self) -> &mut htslib::vcf::bcf_fmt_t {
        unsafe { &mut *self.inner }
    }

    fn values_per_sample(&self) -> usize {
        self.inner().n as usize
    }

    /// Read and decode format data into a given type.
    fn data(&mut self, data_type: i32) -> Result<(usize, i32), FormatReadError> {
        let mut n: i32 = 0;
        match unsafe {
            htslib::vcf::bcf_get_format_values(
                self.record.header,
                self.record.inner,
                ffi::CString::new(self.tag).unwrap().as_ptr() as *mut i8,
                &mut self.record.buffer,
                &mut n,
                data_type
            )
        } {
            -1 => Err(FormatReadError::UndefinedTag),
            -2 => Err(FormatReadError::UnexpectedType),
            -3 => Err(FormatReadError::MissingTag),
            ret  => Ok((n as usize, ret)),
        }
    }

    /// Get format data as integers.
    pub fn integer(&mut self) -> Result<Vec<&'a [i32]>, FormatReadError> {
        self.data(htslib::vcf::BCF_HT_INT).map(|(n, _)| {
            unsafe {
                slice::from_raw_parts(self.record.buffer as *const i32, n)
            }.chunks(self.values_per_sample()).map(|s| trim_slice(s)).collect()
        })
    }

    /// Get format data as mutable integers.
    pub fn integer_mut(&mut self) -> Result<Vec<&'a mut [i32]>, FormatReadError> {
        self.data(htslib::vcf::BCF_HT_INT).map(|(n, _)| {
            unsafe {
                slice::from_raw_parts_mut(self.record.buffer as *mut i32, n)
            }.chunks_mut(self.values_per_sample()).collect()
        })
    }

    /// Get format data as floats.
    pub fn float(&mut self) -> Result<Vec<&'a [f32]>, FormatReadError> {
        self.data(htslib::vcf::BCF_HT_REAL).map(|(n, _)| {
            unsafe {
                slice::from_raw_parts(self.record.buffer as *const f32, n)
            }.chunks(self.values_per_sample()).map(|s| trim_slice(s)).collect()
        })
    }

    pub fn float_mut(&mut self) -> Result<Vec<&'a mut [f32]>, FormatReadError> {
        self.data(htslib::vcf::BCF_HT_REAL).map(|(n, _)| {
            unsafe {
                slice::from_raw_parts_mut(self.record.buffer as *mut f32, n)
            }.chunks_mut(self.values_per_sample()).collect()
        })
    }

    pub fn string(&mut self) -> Result<Vec<&'a [u8]>, FormatReadError> {
        self.data(htslib::vcf::BCF_HT_STR).map(|(n, _)| {
            unsafe {
                slice::from_raw_parts(self.record.buffer as *const u8, n)
            }.chunks(self.values_per_sample()).map(|s| {
                // stop at zero character
                s.split(|c| *c == 0u8).next().expect("Bug: returned string should not be empty.")
            }).collect()
        })
    }

    pub fn string_mut(&mut self) -> Result<Vec<&'a mut [u8]>, FormatReadError> {
        self.data(htslib::vcf::BCF_HT_STR).map(|(n, _)| {
            unsafe {
                slice::from_raw_parts_mut(self.record.buffer as *mut u8, n)
            }.chunks_mut(self.values_per_sample()).collect()
        })
    }
}


unsafe impl<'a> Send for Format<'a> {}
unsafe impl<'a> Sync for Format<'a> {}


quick_error! {
    #[derive(Debug)]
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
    #[derive(Debug)]
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
    #[derive(Debug)]
    pub enum TagWriteError {
        Some {
            description("error writing tag to record")
        }
    }
}


quick_error! {
    #[derive(Debug)]
    pub enum TrimAllelesError {
        Some {
            description("error trimming alleles")
        }
    }
}
