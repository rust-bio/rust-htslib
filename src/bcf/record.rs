// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::ptr;
use std::slice;
use std::ffi;
use std::i32;
use std::f32;

use htslib;

pub const MISSING_INTEGER: i32 = i32::MIN;
pub const MISSING_FLOAT: f32 = f32::NAN;
pub const VECTOR_END_INTEGER: i32 = i32::MIN + 1;
pub const VECTOR_END_FLOAT: f32 = 2.139095e+09;

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


    pub fn set_pos(&mut self, pos: i32) {
        self.inner_mut().pos = pos;
    }

    pub fn alleles(&self) -> Vec<&[u8]> {
        unsafe { htslib::vcf::bcf_unpack(self.inner, htslib::vcf::BCF_UN_STR) };
        let n = self.inner().n_allele as usize;
        let dec = self.inner().d;
        let alleles = unsafe { slice::from_raw_parts(dec.allele, n) };
        (0..n).map(|i| unsafe { ffi::CStr::from_ptr(alleles[i]).to_bytes() }).collect()
    }

    pub fn qual(&self) -> f32 {
        self.inner().qual
    }

    pub fn set_qual(&mut self, qual: f32) {
        self.inner_mut().qual = qual;
    }

    /// Get the value of the given info tag.
    pub fn info<'a>(&'a mut self, tag: &'a [u8]) -> Info {
        Info { record: self, tag: tag }
    }

    pub fn sample_count(&self) -> u32 {
        self.inner().n_fmt_n_sample >> 8
    }

    pub fn allele_count(&self) -> u16 {
        self.inner().n_allele
    }

    /// Get the value of the given format tag for each sample.
    pub fn format<'a>(&'a mut self, tag: &'a [u8]) -> Format {
        Format::new(self, tag)
    }

    /// Add an integer format tag. Data is a flattened two-dimensional array.
    /// The first dimension contains one array for each sample.
    /// Returns error if tag is not present in header.
    pub fn push_format_integer(&mut self, tag: &[u8], data: &[i32]) -> Result<(), ()> {
        self.push_format(tag, data, htslib::vcf::BCF_HT_INT)
    }

    /// Add a float format tag. Data is a flattened two-dimensional array.
    /// The first dimension contains one array for each sample.
    /// Returns error if tag is not present in header.
    pub fn push_format_float(&mut self, tag: &[u8], data: &[f32]) -> Result<(), ()> {
        self.push_format(tag, data, htslib::vcf::BCF_HT_REAL)
    }

    /// Add a format tag. Data is a flattened two-dimensional array.
    /// The first dimension contains one array for each sample.
    fn push_format<T>(&mut self, tag: &[u8], data: &[T], ht: i32) -> Result<(), ()> {
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
                Err(())
            }
        }
    }

    /// Add a info tag.
    pub fn push_info<T>(&mut self, tag: &[u8], data: &[T], ht: i32) -> Result<(), ()> {
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
                Err(())
            }
        }
    }

    /// Remove unused alleles.
    pub fn trim_alleles(&mut self) -> Result<(), ()> {
        match unsafe { htslib::vcfutils::bcf_trim_alleles(self.header, self.inner) } {
            -1 => Err(()),
            _  => Ok(())
        }
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


pub struct Info<'a> {
    record: &'a mut Record,
    tag: &'a [u8],
}


impl<'a> Info<'a> {
    fn data(&mut self, data_type: i32) -> Result<(usize, i32), TagError> {
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
            -1 => Err(TagError::UndefinedTag),
            -2 => Err(TagError::UnexpectedType),
            -3 => Err(TagError::MissingTag),
            ret  => Ok((n as usize, ret)),
        }
    }

    pub fn integer(&mut self) -> Result<&[i32], TagError> {
        self.data(htslib::vcf::BCF_HT_INT).map(|(n, _)| {
            trim_slice(
                unsafe { slice::from_raw_parts(self.record.buffer as *const i32, n) },
                VECTOR_END_INTEGER
            )
        })
    }

    pub fn integer_mut(&mut self) -> Result<&mut [i32], TagError> {
        self.data(htslib::vcf::BCF_HT_INT).map(|(n, _)| {
            unsafe { slice::from_raw_parts_mut(self.record.buffer as *mut i32, n) }
        })
    }

    pub fn float(&mut self) -> Result<&[f32], TagError> {
        self.data(htslib::vcf::BCF_HT_REAL).map(|(n, _)| {
            trim_slice(
                unsafe { slice::from_raw_parts(self.record.buffer as *const f32, n) },
                VECTOR_END_FLOAT
            )
        })
    }

    pub fn float_mut(&mut self) -> Result<&mut [f32], TagError> {
        self.data(htslib::vcf::BCF_HT_REAL).map(|(n, _)| {
            unsafe { slice::from_raw_parts_mut(self.record.buffer as *mut f32, n) }
        })
    }

    pub fn flag(&mut self) -> Result<bool, TagError> {
        self.data(htslib::vcf::BCF_HT_FLAG).map(|(_, ret)| {
            ret == 1
        })
    }

    pub fn string(&mut self) -> Result<Vec<&[u8]>, TagError> {
        self.data(htslib::vcf::BCF_HT_STR).map(|(n, ret)| {
            unsafe {
                slice::from_raw_parts(self.record.buffer as *const u8, ret as usize)
            }.chunks(n).map(|s| {
                // stop at zero character
                s.split(|c| *c == 0u8).next().expect("Bug: returned string should not be empty.")
            }).collect()
        })
    }

    pub fn string_mut(&mut self) -> Result<Vec<&mut [u8]>, TagError> {
        self.data(htslib::vcf::BCF_HT_STR).map(|(n, ret)| {
            unsafe {
                slice::from_raw_parts_mut(self.record.buffer as *mut u8, ret as usize)
            }.chunks_mut(n).collect()
        })
    }
}


fn trim_slice<T: PartialEq>(s: &[T], end_value: T) -> &[T] {
    s.split(|v| *v == end_value).next().expect("Bug: returned slice should not be empty.")
}


// TODO implement format.
pub struct Format<'a> {
    record: &'a mut Record,
    tag: &'a [u8],
    inner: *mut htslib::vcf::bcf_fmt_t,
}


impl<'a> Format<'a> {
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

    fn data(&mut self, data_type: i32) -> Result<(usize, i32), TagError> {
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
            -1 => Err(TagError::UndefinedTag),
            -2 => Err(TagError::UnexpectedType),
            -3 => Err(TagError::MissingTag),
            ret  => Ok((n as usize, ret)),
        }
    }

    pub fn integer(&mut self) -> Result<Vec<&[i32]>, TagError> {
        self.data(htslib::vcf::BCF_HT_INT).map(|(n, _)| {
            unsafe {
                slice::from_raw_parts(self.record.buffer as *const i32, n)
            }.chunks(self.values_per_sample()).map(|s| trim_slice(s, VECTOR_END_INTEGER)).collect()
        })
    }

    pub fn integer_mut(&mut self) -> Result<Vec<&mut [i32]>, TagError> {
        self.data(htslib::vcf::BCF_HT_INT).map(|(n, _)| {
            unsafe {
                slice::from_raw_parts_mut(self.record.buffer as *mut i32, n)
            }.chunks_mut(self.values_per_sample()).collect()
        })
    }

    pub fn float(&mut self) -> Result<Vec<&[f32]>, TagError> {
        self.data(htslib::vcf::BCF_HT_REAL).map(|(n, _)| {
            unsafe {
                slice::from_raw_parts(self.record.buffer as *const f32, n)
            }.chunks(self.values_per_sample()).map(|s| trim_slice(s, VECTOR_END_FLOAT)).collect()
        })
    }

    pub fn float_mut(&mut self) -> Result<Vec<&mut [f32]>, TagError> {
        self.data(htslib::vcf::BCF_HT_REAL).map(|(n, _)| {
            unsafe {
                slice::from_raw_parts_mut(self.record.buffer as *mut f32, n)
            }.chunks_mut(self.values_per_sample()).collect()
        })
    }

    pub fn string(&mut self) -> Result<Vec<&[u8]>, TagError> {
        self.data(htslib::vcf::BCF_HT_STR).map(|(n, _)| {
            unsafe {
                slice::from_raw_parts(self.record.buffer as *const u8, n)
            }.chunks(self.values_per_sample()).map(|s| {
                // stop at zero character
                s.split(|c| *c == 0u8).next().expect("Bug: returned string should not be empty.")
            }).collect()
        })
    }

    pub fn string_mut(&mut self) -> Result<Vec<&mut [u8]>, TagError> {
        self.data(htslib::vcf::BCF_HT_STR).map(|(n, _)| {
            unsafe {
                slice::from_raw_parts_mut(self.record.buffer as *mut u8, n)
            }.chunks_mut(self.values_per_sample()).collect()
        })
    }
}


#[derive(Debug)]
pub enum TagError {
    UndefinedTag,
    UnexpectedType,
    MissingTag,
}
