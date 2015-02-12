use htslib;
use std::ffi::c_str_to_bytes;
use std::slice::from_raw_parts;
use std::ffi::CString;


pub struct Record<'a> {
    b: &'a htslib::bam1_t,
    data: &'a [u8],
}

impl<'a> Record<'a>{
    pub fn new(b: &'a htslib::bam1_t) -> Record {
        Record { b: b, data: unsafe { from_raw_parts((*b).data, b.l_data as usize) } }
    }

    pub fn qname(&'a self) -> &'a [u8] {
        self.data[1..10].as_slice()
    }
}
