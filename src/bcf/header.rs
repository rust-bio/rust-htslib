// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


use std::slice;
use std::ffi;
use std::str;

use htslib;


pub type SampleSubset = Vec<i32>;


/// A BCF header.
pub struct Header {
    pub inner: *mut htslib::vcf::bcf_hdr_t,
    pub subset: Option<SampleSubset>,
}


impl Header {
    /// Create a new header.
    pub fn new() -> Self {
        Header {
            inner: unsafe { htslib::vcf::bcf_hdr_init(ffi::CString::new(&b"w"[..]).unwrap().as_ptr()) },
            subset: None,
        }
    }

    pub fn with_template(header: &HeaderView) -> Self {
        Header { inner: unsafe { htslib::vcf::bcf_hdr_dup(header.inner) }, subset: None }
    }

    pub fn subset_template(header: &HeaderView, samples: &[&[u8]]) -> Result<Self, SubsetError> {
        let mut imap = vec![0; samples.len()];
        let names: Vec<_> = samples.iter().map(|&s| ffi::CString::new(s).unwrap()).collect();
        let name_pointers: Vec<_> = names.iter().map(|s| s.as_ptr() as *mut i8).collect();
        let inner = unsafe {
            htslib::vcf::bcf_hdr_subset(header.inner, samples.len() as i32, name_pointers.as_ptr(), imap.as_mut_ptr() as *mut i32)
        };
        if inner.is_null() {
            Err(SubsetError::DuplicateSampleName)
        }
        else {
            Ok(Header { inner: inner, subset: Some(imap) })
        }
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

    pub fn remove_info(&mut self, tag: &[u8]) -> &mut Self {
        unsafe {
            htslib::vcf::bcf_hdr_remove(self.inner, htslib::vcf::BCF_HL_INFO, tag.as_ptr() as *const i8);
        }
        self
    }

    pub fn remove_format(&mut self, tag: &[u8]) -> &mut Self {
        unsafe {
            htslib::vcf::bcf_hdr_remove(self.inner, htslib::vcf::BCF_HL_FMT, tag.as_ptr() as *const i8);
        }
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

    pub fn rid2name(&self, rid: u32) -> &[u8] {
        unsafe {
            let dict = self.inner().id[htslib::vcf::BCF_DT_CTG as usize];
            let ptr = (*dict.offset(rid as isize)).key;
            ffi::CStr::from_ptr(ptr).to_bytes()
        }
    }

    pub fn name2rid(&self, name: &[u8]) -> Result<u32, RidError> {
        unsafe {
            match htslib::vcf::bcf_hdr_id2int(
                self.inner,
                htslib::vcf::BCF_DT_CTG,
                ffi::CString::new(name).unwrap().as_ptr() as *mut i8
            ) {
                -1 => Err(RidError::UnknownSequence(str::from_utf8(name).unwrap().to_owned())),
                i  => Ok(i as u32)
            }
        }
    }

    pub fn info_type(&self, tag: &[u8]) -> Result<(TagType, TagLength), TagTypeError> {
        self.tag_type(tag, htslib::vcf::BCF_HL_INFO)
    }

    pub fn format_type(&self, tag: &[u8]) -> Result<(TagType, TagLength), TagTypeError> {
        self.tag_type(tag, htslib::vcf::BCF_HL_FMT)
    }

    fn tag_type(&self, tag: &[u8], hdr_type: ::libc::c_int) -> Result<(TagType, TagLength), TagTypeError> {
        let (_type, length) = unsafe {
            let id = htslib::vcf::bcf_hdr_id2int(
                self.inner,
                htslib::vcf::BCF_DT_ID,
                ffi::CString::new(tag).unwrap().as_ptr() as *mut i8
            );
            if id < 0 {
                return Err(TagTypeError::UndefinedTag(str::from_utf8(tag).unwrap().to_owned()));
            }
            let n = (*self.inner).n[htslib::vcf::BCF_DT_ID as usize] as usize;
            let entry = slice::from_raw_parts((*self.inner).id[htslib::vcf::BCF_DT_ID as usize], n);
            let d = (*entry[id as usize].val).info[hdr_type as usize];
            (d >> 4 & 0xf, d >> 8 & 0xf)
        };
        let _type = match _type as ::libc::c_int {
            htslib::vcf::BCF_HT_FLAG => TagType::Flag,
            htslib::vcf::BCF_HT_INT => TagType::Integer,
            htslib::vcf::BCF_HT_REAL => TagType::Float,
            htslib::vcf::BCF_HT_STR => TagType::String,
            _ => return Err(TagTypeError::UnexpectedTagType)
        };
        let length = match length as ::libc::c_int {
            htslib::vcf::BCF_VL_FIXED => TagLength::Fixed,
            htslib::vcf::BCF_VL_VAR => TagLength::Variable,
            htslib::vcf::BCF_VL_A => TagLength::AltAlleles,
            htslib::vcf::BCF_VL_R => TagLength::Alleles,
            htslib::vcf::BCF_VL_G => TagLength::Genotypes,
            _ => return Err(TagTypeError::UnexpectedTagType)
        };

        Ok((_type, length))
    }
}


impl Clone for HeaderView {
    fn clone(&self) -> Self {
        HeaderView {
            inner: unsafe { htslib::vcf::bcf_hdr_dup(self.inner) }
        }
    }
}


impl Drop for HeaderView {
    fn drop(&mut self) {
        unsafe {
            htslib::vcf::bcf_hdr_destroy(self.inner);
        }
    }
}


pub enum TagType {
    Flag,
    Integer,
    Float,
    String
}


pub enum TagLength {
    Fixed,
    AltAlleles,
    Alleles,
    Genotypes,
    Variable
}


quick_error! {
    #[derive(Debug)]
    pub enum RidError {
        UnknownSequence(name: String) {
            description("unknown sequence")
            display("sequence {} not found in header", name)
        }
    }
}


quick_error! {
    #[derive(Debug)]
    pub enum SubsetError {
        DuplicateSampleName {
            description("duplicate sample name when subsetting header")
        }
    }
}


quick_error! {
    #[derive(Debug)]
    pub enum TagTypeError {
        UnexpectedTagType {
            description("unexpected tag type in header")
        }
        UndefinedTag(name: String) {
            description("undefined tag")
            display("tag {} is undefined in header", name)
        }
    }
}
