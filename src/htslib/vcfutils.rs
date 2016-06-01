#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]
#![allow(non_snake_case)]


use htslib::vcf::{bcf1_t, bcf_hdr_t};


extern "C" {
    pub fn bcf_trim_alleles(hdr: *const bcf_hdr_t, line: *mut bcf1_t) -> ::libc::c_int;
}
