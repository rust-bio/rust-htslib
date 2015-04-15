// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


use std::slice;
use std::ffi;

use htslib;


/// A BCF header.
pub struct Header {
    pub inner: *mut htslib::vcf::bcf_hdr_t,
}


impl Header {
    /// Create a new header.
    pub fn new() -> Self {
        Header { inner: unsafe { htslib::vcf::bcf_hdr_init(ffi::CString::new(&b"w"[..]).unwrap().as_ptr()) } }
    }

    pub fn with_template(header: &HeaderView) -> Self {
        Header { inner: unsafe { htslib::vcf::bcf_hdr_dup(header.inner) } }
    }

    pub fn push_sample(&mut self, sample: &[u8]) -> &mut Self {
        unsafe { htslib::vcf::bcf_hdr_add_sample(self.inner, ffi::CString::new(sample).unwrap().as_ptr()) };
        self
    }

    /// Add a record to the header.
    pub fn push_record(&mut self, record: &[u8]) -> &mut Self {
        unsafe { htslib::vcf::bcf_hdr_append(self.inner, ffi::CString::new(record).unwrap().as_ptr()) };
        self
    }
}


impl Drop for Header {
    fn drop(&mut self) {
        unsafe { htslib::vcf::bcf_hdr_destroy(self.inner) };
    }
}



pub struct HeaderView {
    pub inner: *mut htslib::vcf::bcf_hdr_t,
}


impl HeaderView {
    pub fn new(inner: *mut htslib::vcf::bcf_hdr_t) -> Self {
        HeaderView { inner: inner }
    }

    #[inline]
    fn inner(&self) -> htslib::vcf::bcf_hdr_t {
        unsafe { (*self.inner) }
    }

    pub fn sample_count(&self) -> u32 {
        self.inner().n[htslib::vcf::BCF_DT_SAMPLE as usize] as u32
    }

    pub fn samples(&self) -> Vec<&[u8]> {
        let names = unsafe { slice::from_raw_parts(self.inner().samples, self.sample_count() as usize) };
        names.iter().map(|name| unsafe { ffi::CStr::from_ptr(*name).to_bytes() }).collect()
    }
}
