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
#[derive(Debug)]
pub struct Header {
    pub inner: *mut htslib::bcf_hdr_t,
    pub subset: Option<SampleSubset>,
}


impl Header {
    /// Create a new header.
    pub fn new() -> Self {
        Header {
            inner: unsafe { htslib::bcf_hdr_init(ffi::CString::new(&b"w"[..]).unwrap().as_ptr()) },
            subset: None,
        }
    }

    pub fn with_template(header: &HeaderView) -> Self {
        Header { inner: unsafe { htslib::bcf_hdr_dup(header.inner) }, subset: None }
    }

    pub fn subset_template(header: &HeaderView, samples: &[&[u8]]) -> Result<Self, SubsetError> {
        let mut imap = vec![0; samples.len()];
        let names: Vec<_> = samples.iter().map(|&s| ffi::CString::new(s).unwrap()).collect();
        let name_pointers: Vec<_> = names.iter().map(|s| s.as_ptr() as *mut i8).collect();
        let inner = unsafe {
            htslib::bcf_hdr_subset(
                header.inner,
                samples.len() as i32,
                name_pointers.as_ptr() as *const *const i8,
                imap.as_mut_ptr() as *mut i32
            )
        };
        if inner.is_null() {
            Err(SubsetError::DuplicateSampleName)
        }
        else {
            Ok(Header { inner: inner, subset: Some(imap) })
        }
    }

    pub fn push_sample(&mut self, sample: &[u8]) -> &mut Self {
        unsafe { htslib::bcf_hdr_add_sample(self.inner, ffi::CString::new(sample).unwrap().as_ptr()) };
        self
    }

    /// Add a record to the header.
    pub fn push_record(&mut self, record: &[u8]) -> &mut Self {
        unsafe { htslib::bcf_hdr_append(self.inner, ffi::CString::new(record).unwrap().as_ptr()) };
        self
    }

    pub fn remove_info(&mut self, tag: &[u8]) -> &mut Self {
        unsafe {
            htslib::bcf_hdr_remove(self.inner, htslib::BCF_HL_INFO as i32, tag.as_ptr() as *const i8);
        }
        self
    }

    pub fn remove_format(&mut self, tag: &[u8]) -> &mut Self {
        unsafe {
            htslib::bcf_hdr_remove(
                self.inner,
                htslib::BCF_HL_FMT as i32,
                tag.as_ptr() as *const i8
            );
        }
        self
    }
}


impl Drop for Header {
    fn drop(&mut self) {
        unsafe { htslib::bcf_hdr_destroy(self.inner) };
    }
}


#[derive(Debug)]
pub struct HeaderView {
    pub inner: *mut htslib::bcf_hdr_t,
}


impl HeaderView {
    pub fn new(inner: *mut htslib::bcf_hdr_t) -> Self {
        HeaderView { inner: inner }
    }

    #[inline]
    fn inner(&self) -> htslib::bcf_hdr_t {
        unsafe { (*self.inner) }
    }

    pub fn sample_count(&self) -> u32 {
        self.inner().n[htslib::BCF_DT_SAMPLE as usize] as u32
    }

    pub fn samples(&self) -> Vec<&[u8]> {
        let names = unsafe { slice::from_raw_parts(self.inner().samples, self.sample_count() as usize) };
        names.iter().map(|name| unsafe { ffi::CStr::from_ptr(*name).to_bytes() }).collect()
    }

    pub fn rid2name(&self, rid: u32) -> &[u8] {
        unsafe {
            let dict = self.inner().id[htslib::BCF_DT_CTG as usize];
            let ptr = (*dict.offset(rid as isize)).key;
            ffi::CStr::from_ptr(ptr).to_bytes()
        }
    }

    pub fn name2rid(&self, name: &[u8]) -> Result<u32, RidError> {
        unsafe {
            match htslib::bcf_hdr_id2int(
                self.inner,
                htslib::BCF_DT_CTG as i32,
                ffi::CString::new(name).unwrap().as_ptr() as *mut i8
            ) {
                -1 => Err(RidError::UnknownSequence(str::from_utf8(name).unwrap().to_owned())),
                i  => Ok(i as u32)
            }
        }
    }

    pub fn info_type(&self, tag: &[u8]) -> Result<(TagType, TagLength), TagTypeError> {
        self.tag_type(tag, htslib::BCF_HL_INFO)
    }

    pub fn format_type(&self, tag: &[u8]) -> Result<(TagType, TagLength), TagTypeError> {
        self.tag_type(tag, htslib::BCF_HL_FMT)
    }

    fn tag_type(&self, tag: &[u8], hdr_type: ::libc::c_uint) -> Result<(TagType, TagLength), TagTypeError> {
        let (_type, length) = unsafe {
            let id = htslib::bcf_hdr_id2int(
                self.inner,
                htslib::BCF_DT_ID as i32,
                ffi::CString::new(tag).unwrap().as_ptr() as *mut i8
            );
            if id < 0 {
                return Err(TagTypeError::UndefinedTag(str::from_utf8(tag).unwrap().to_owned()));
            }
            let n = (*self.inner).n[htslib::BCF_DT_ID as usize] as usize;
            let entry = slice::from_raw_parts((*self.inner).id[htslib::BCF_DT_ID as usize], n);
            let d = (*entry[id as usize].val).info[hdr_type as usize];
            (d >> 4 & 0xf, d >> 8 & 0xf)
        };
        let _type = match _type as ::libc::c_uint {
            htslib::BCF_HT_FLAG => TagType::Flag,
            htslib::BCF_HT_INT => TagType::Integer,
            htslib::BCF_HT_REAL => TagType::Float,
            htslib::BCF_HT_STR => TagType::String,
            _ => return Err(TagTypeError::UnexpectedTagType)
        };
        let length = match length as ::libc::c_uint {
            htslib::BCF_VL_FIXED => TagLength::Fixed,
            htslib::BCF_VL_VAR => TagLength::Variable,
            htslib::BCF_VL_A => TagLength::AltAlleles,
            htslib::BCF_VL_R => TagLength::Alleles,
            htslib::BCF_VL_G => TagLength::Genotypes,
            _ => return Err(TagTypeError::UnexpectedTagType)
        };

        Ok((_type, length))
    }

    /// Convert string ID (e.g., for a `FILTER` value) to its numeric identifier.
    pub fn id2int(&self, id: &[u8]) -> Result<u32, IdError> {
        unsafe {
            match htslib::bcf_hdr_id2int(
                self.inner,
                htslib::BCF_DT_ID as i32,
                ffi::CString::new(id).unwrap().as_ptr() as *const i8
            ) {
                -1 => Err(IdError::UnknownID(str::from_utf8(id).unwrap().to_owned())),
                i => Ok(i as u32),
            }
        }
    }

    /// Convert integer representing an identifier (e.g., a `FILTER` value) to its string
    /// name.bam
    pub fn int2id(&self, id: i32) -> Vec<u8> {
        let key = unsafe { ffi::CStr::from_ptr(
            (*(*self.inner).id[htslib::BCF_DT_ID as usize].offset(id as isize)).key) };
        key.to_bytes().to_vec()
    }

    /// Convert string sample name to its numeric identifier.
    pub fn sample2int(&self, id: &[u8]) -> Result<u32, SampleError> {
        unsafe {
            match htslib::bcf_hdr_id2int(
                self.inner,
                htslib::BCF_DT_SAMPLE as i32,
                ffi::CString::new(id).unwrap().as_ptr() as *const i8
            ) {
                -1 => Err(SampleError::UnknownSample(str::from_utf8(id).unwrap().to_owned())),
                i => Ok(i as u32),
            }
        }
    }

    /// Convert integer representing an contig to its name.
    pub fn int2sample(&self, id: i32) -> Vec<u8> {
        let key = unsafe { ffi::CStr::from_ptr(
            (*(*self.inner).id[htslib::BCF_DT_SAMPLE as usize].offset(id as isize)).key) };
        key.to_bytes().to_vec()
    }
}


impl Clone for HeaderView {
    fn clone(&self) -> Self {
        HeaderView {
            inner: unsafe { htslib::bcf_hdr_dup(self.inner) }
        }
    }
}


impl Drop for HeaderView {
    fn drop(&mut self) {
        unsafe {
            htslib::bcf_hdr_destroy(self.inner);
        }
    }
}


#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum TagType {
    Flag,
    Integer,
    Float,
    String
}


#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum TagLength {
    Fixed,
    AltAlleles,
    Alleles,
    Genotypes,
    Variable
}


quick_error! {
    #[derive(Debug, Clone)]
    pub enum RidError {
        UnknownSequence(name: String) {
            description("unknown sequence")
            display("sequence {} not found in header", name)
        }
    }
}


quick_error! {
    #[derive(Debug, Clone)]
    pub enum IdError {
        UnknownID(name: String) {
            description("unknown ID")
            display("ID {} not found in header", name)
        }
    }
}


quick_error! {
    #[derive(Debug, Clone)]
    pub enum SampleError {
        UnknownSample(name: String) {
            description("unknown sample")
            display("sample {} not found in header", name)
        }
    }
}


quick_error! {
    #[derive(Debug, Clone)]
    pub enum SubsetError {
        DuplicateSampleName {
            description("duplicate sample name when subsetting header")
        }
    }
}


quick_error! {
    #[derive(Debug, Clone)]
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
