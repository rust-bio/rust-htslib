// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::ptr;
use std::slice;
use std::ffi;

use htslib;


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

    #[inline]
    fn inner(&self) -> htslib::vcf::bcf1_t {
        unsafe { *self.inner }
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

    /// Get the value of the given info tag.
    pub fn info<'a>(&'a mut self, tag: &'a [u8]) -> Info {
        Info { record: self, tag: tag }
    }

    pub fn sample_count(&self) -> u32 {
        self.inner().n_fmt_n_sample & 0xffffff
    }

    pub fn allele_count(&self) -> u16 {
        self.inner().n_allele
    }

    /// Get the value of the given format tag for each sample.
    pub fn format<'a>(&'a mut self, tag: &'a [u8]) -> Format {
        Format::new(self, tag)
    }

    pub fn set_qual(&mut self, qual: f32) {
        self.inner().qual = qual;
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

    pub fn integer(&mut self) -> Result<&mut [i32], TagError> {
        self.data(htslib::vcf::BCF_HT_INT).map(|(n, _)| {
            unsafe { slice::from_raw_parts_mut(self.record.buffer as *mut i32, n) }
        })
    }

    pub fn float(&mut self) -> Result<&mut [f32], TagError> {
        self.data(htslib::vcf::BCF_HT_REAL).map(|(n, _)| {
            unsafe { slice::from_raw_parts_mut(self.record.buffer as *mut f32, n) }
        })
    }

    pub fn flag(&mut self) -> Result<bool, TagError> {
        self.data(htslib::vcf::BCF_HT_FLAG).map(|(_, ret)| {
            ret == 1
        })
    }

    pub fn string(&mut self) -> Result<&mut [u8], TagError> {
        self.data(htslib::vcf::BCF_HT_STR).map(|(_, ret)| {
            unsafe { slice::from_raw_parts_mut(self.record.buffer as *mut u8, ret as usize) }
        })
    }
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

    fn values_per_sample(&self) -> usize {
        unsafe { (*self.inner).n as usize }
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

    pub fn integer(&mut self) -> Result<Vec<&mut [i32]>, TagError> {
        self.data(htslib::vcf::BCF_HT_INT).map(|(n, _)| {
            unsafe {
                slice::from_raw_parts_mut(self.record.buffer as *mut i32, n)
            }.chunks_mut(self.values_per_sample()).collect()
        })
    }

    pub fn float(&mut self) -> Result<Vec<&mut [f32]>, TagError> {
        self.data(htslib::vcf::BCF_HT_REAL).map(|(n, _)| {
            unsafe {
                slice::from_raw_parts_mut(self.record.buffer as *mut f32, n)
            }.chunks_mut(self.values_per_sample()).collect()
        })
    }
}


pub enum TagError {
    UndefinedTag,
    UnexpectedType,
    MissingTag,
}
