// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::ptr;
use std::slice;
use std::ffi;
use std::vec;

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

    pub fn qual(&self) -> f32 {
        self.inner().qual
    }

    pub fn info<'a>(&'a mut self, tag: &'a [u8]) -> Info {
        Info { record: self, tag: tag }
    }

    pub fn format<'a>(&'a mut self, tag: &'a [u8]) -> Format {
        Format { record: self, tag: tag }
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


pub struct Info<'a> {
    record: &'a mut Record,
    tag: &'a [u8],
}


impl<'a> Info<'a> {
    fn data(&mut self, data_type: i32) -> Result<(usize, i32), InfoError> {
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
            -1 => Err(InfoError::UndefinedTag),
            -2 => Err(InfoError::UnexpectedType),
            -3 => Err(InfoError::MissingTag),
            ret  => Ok((n as usize, ret)),
        }
    }

    pub fn integer(&mut self) -> Result<&[i32], InfoError> {
        self.data(htslib::vcf::BCF_HT_INT).map(|(n, _)| {
            unsafe { slice::from_raw_parts(self.record.buffer as *mut i32, n) }
        })
    }

    pub fn float(&mut self) -> Result<&[f32], InfoError> {
        self.data(htslib::vcf::BCF_HT_REAL).map(|(n, _)| {
            unsafe { slice::from_raw_parts(self.record.buffer as *mut f32, n) }
        })
    }

    pub fn flag(&mut self) -> Result<bool, InfoError> {
        self.data(htslib::vcf::BCF_HT_FLAG).map(|(_, ret)| {
            ret == 1
        })
    }

    pub fn string(&mut self) -> Result<&[u8], InfoError> {
        self.data(htslib::vcf::BCF_HT_STR).map(|(_, ret)| {
            unsafe { slice::from_raw_parts(self.record.buffer as *mut u8, ret as usize) }
        })
    }
}


// TODO implement format.
pub struct Format<'a> {
    record: &'a Record,
    tag: &'a [u8],
}


pub enum InfoError {
    UndefinedTag,
    UnexpectedType,
    MissingTag,
}
