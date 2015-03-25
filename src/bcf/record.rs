// Copyright 2014 Christopher Schröder, Johannes Köster.
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
}


impl Record {
    pub fn new() -> Self {
        let inner = unsafe { htslib::vcf::bcf_init() };
        Record { inner: inner, header: ptr::null_mut() }
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

    pub fn qual(&self) -> f32 {
        self.inner().qual
    }

    pub fn info<'a>(&'a self, tag: &'a [u8]) -> Info {
        Info { record: self, tag: tag }
    }
}


impl Drop for Record {
    fn drop(&mut self) {
        unsafe { htslib::vcf::bcf_destroy(self.inner) };
    }
}


pub struct Info<'a> {
    record: &'a Record,
    tag: &'a [u8],
}


impl<'a> Info<'a> {
    fn data(&self, data_type: i32) -> Result<(*mut ::libc::c_void, usize, i32), InfoError> {
        let mut data = ptr::null_mut();
        let mut n: i32 = 0;
        match unsafe {
            htslib::vcf::bcf_get_info_values(
                self.record.header,
                self.record.inner,
                ffi::CString::new(self.tag).unwrap().as_ptr() as *mut i8,
                &mut data,
                &mut n,
                data_type
            )
        } {
            -1 => Err(InfoError::UndefinedTag),
            -2 => Err(InfoError::UnexpectedType),
            -3 => Err(InfoError::MissingTag),
            ret  => Ok((data, n as usize, ret)),
        }
    }

    pub fn integer(&self) -> Result<&[i32], InfoError> {
        self.data(htslib::vcf::BCF_HT_INT).map(|(data, n, _)| {
            unsafe { slice::from_raw_parts(data as *mut i32, n) }
        })
    }

    pub fn float(&self) -> Result<&[f32], InfoError> {
        self.data(htslib::vcf::BCF_HT_REAL).map(|(data, n, _)| {
            unsafe { slice::from_raw_parts(data as *mut f32, n) }
        })
    }

    pub fn flag(&self) -> Result<bool, InfoError> {
        self.data(htslib::vcf::BCF_HT_FLAG).map(|(_, _, ret)| {
            ret == 1
        })
    }

    pub fn string(&self) -> Result<&[u8], InfoError> {
        self.data(htslib::vcf::BCF_HT_STR).map(|(data, _, ret)| {
            unsafe { slice::from_raw_parts(data as *mut u8, ret as usize) }
        })
    }
}


pub enum InfoError {
    UndefinedTag,
    UnexpectedType,
    MissingTag,
}
