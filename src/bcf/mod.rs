
use std::ffi;
use std::path;
use std::ffi::AsOsStr;
use std::os::unix::prelude::OsStrExt;


pub mod record;

use htslib;


pub struct Reader {
    inner: *mut htslib::vcf::htsFile,
    header: *mut htslib::vcf::bcf_hdr_t,
}

impl Reader {
   pub fn new<P: path::AsPath>(path: &P) -> Self {
        let htsfile = bcf_open(path, b"r");
        let header = unsafe { htslib::vcf::bcf_hdr_read(htsfile) };
        Reader { inner: htsfile, header: header }
    }

    pub fn read(&self, record: &mut record::Record) -> Result<(), ReadError> {
        match unsafe { htslib::vcf::bcf_read(self.inner, self.header, record.inner) } {
            -1 => Err(ReadError::NoMoreRecord),
            0  => Ok(()),
            _  => Err(ReadError::Invalid),
        }
    }

    pub fn records(&self) -> Records {
        Records { reader: self }
    }
}


pub struct Records<'a> {
    reader: &'a Reader,
}


impl<'a> Iterator for Records<'a> {
    type Item = Result<record::Record, ReadError>;

    fn next(&mut self) -> Option<Result<record::Record, ReadError>> {
        let mut record = record::Record::new();
        match self.reader.read(&mut record) {
            Err(ReadError::NoMoreRecord) => None,
            Err(e)                       => Some(Err(e)),
            Ok(())                       => Some(Ok(record)),
        }
    }
}


/// Wrapper for opening a BCF file.
fn bcf_open<P: path::AsPath>(path: &P, mode: &[u8]) -> *mut htslib::vcf::htsFile {
    unsafe {
        htslib::vcf::hts_open(
            path.as_path().as_os_str().to_cstring().unwrap().as_ptr(),
            ffi::CString::new(mode).unwrap().as_ptr()
        )
    }
}


pub enum ReadError {
    Truncated,
    Invalid,
    NoMoreRecord,
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read() {
        let bcf = Reader::new(&"test.bcf");
        for (i, rec) in bcf.records().enumerate() {
            let record = rec.ok().expect("Error reading record.");
            assert_eq!(record.rid().expect("Error reading rid."), 0);
            assert_eq!(record.pos(), 10021 + i as u32);
            assert_eq!(record.qual(), 0f32);
        }
    }
}
