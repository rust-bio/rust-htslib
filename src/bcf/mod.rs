
use std::ffi;
use std::convert::AsRef;
use std::path::Path;
use std::fmt;
use std::error::Error;

pub mod record;
pub mod header;

use htslib;
use bcf::header::{HeaderView, SampleSubset};
use utils;

pub use bcf::header::Header;
pub use bcf::record::Record;


pub struct Reader {
    inner: *mut htslib::vcf::htsFile,
    pub header: HeaderView,
}


unsafe impl Send for Reader {}


impl Reader {
   pub fn new<P: AsRef<Path>>(path: &P) -> Result<Self, BCFError> {
        let htsfile = try!(bcf_open(path, b"r"));
        let header = unsafe { htslib::vcf::bcf_hdr_read(htsfile) };
        Ok(Reader { inner: htsfile, header: HeaderView::new(header) })
    }

    pub fn read(&self, record: &mut record::Record) -> Result<(), ReadError> {
        match unsafe { htslib::vcf::bcf_read(self.inner, self.header.inner, record.inner) } {
            0  => {
                record.header = self.header.inner;
                Ok(())
            },
            -1 => Err(ReadError::NoMoreRecord),
            _  => Err(ReadError::Invalid),
        }
    }

    pub fn records(&self) -> Records {
        Records { reader: self }
    }
}


impl Drop for Reader {
    fn drop(&mut self) {
        unsafe {
            htslib::vcf::bcf_hdr_destroy(self.header.inner);
            htslib::vcf::hts_close(self.inner);
        }
    }
}


pub struct Writer {
    inner: *mut htslib::vcf::htsFile,
    pub header: HeaderView,
    subset: Option<SampleSubset>,
}


unsafe impl Send for Writer {}


impl Writer {
    pub fn new<P: AsRef<Path>>(path: &P, header: &Header, uncompressed: bool, vcf: bool) -> Result<Self, BCFError> {
        let mode: &[u8] = match (uncompressed, vcf) {
            (true, true)   => b"w",
            (false, true)  => b"wz",
            (true, false)  => b"wu",
            (false, false) => b"wb",
        };

        let htsfile = try!(bcf_open(path, mode));
        unsafe { htslib::vcf::bcf_hdr_write(htsfile, header.inner) };
        Ok(Writer {
            inner: htsfile,
            header: HeaderView::new(unsafe { htslib::vcf::bcf_hdr_dup(header.inner) }),
            subset: header.subset.clone()
        })
    }

    /// Translate record to header of this writer.
    pub fn translate(&mut self, record: &mut record::Record) {
        unsafe {
            htslib::vcf::bcf_translate(self.header.inner, record.header, record.inner);
        }
        record.header = self.header.inner;
    }

    /// Subset samples of record to match header of this writer.
    pub fn subset(&mut self, record: &mut record::Record) {
        match self.subset {
            Some(ref mut subset) => unsafe {
                htslib::vcf::bcf_subset(self.header.inner, record.inner, subset.len() as i32, subset.as_mut_ptr());
            },
            None         => ()
        }
    }

    pub fn write(&mut self, record: &record::Record) -> Result<(), ()> {
        if unsafe { htslib::vcf::bcf_write(self.inner, self.header.inner, record.inner) } == -1 {
            Err(())
        }
        else {
            Ok(())
        }
    }
}


impl Drop for Writer {
    fn drop(&mut self) {
        unsafe {
            htslib::vcf::bcf_hdr_destroy(self.header.inner);
            htslib::vcf::hts_close(self.inner);
        }
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


#[derive(Debug)]
pub enum BCFError {
    InvalidPath
}


impl fmt::Display for BCFError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.description().fmt(f)
    }
}


impl Error for BCFError {
    fn description(&self) -> &str {
        match self {
            &BCFError::InvalidPath => "invalid path"
        }
    }
}


/// Wrapper for opening a BCF file.
fn bcf_open<P: AsRef<Path>>(path: &P, mode: &[u8]) -> Result<*mut htslib::vcf::htsFile, BCFError> {
    if let Some(p) = utils::path_to_cstring(path) {
        Ok(unsafe {
            htslib::vcf::hts_open(
                p.as_ptr(),
                ffi::CString::new(mode).unwrap().as_ptr()
            )
        })
    }
    else {
        Err(BCFError::InvalidPath)
    }
}


#[derive(Debug)]
pub enum ReadError {
    Invalid,
    NoMoreRecord,
}

impl ReadError {
    /// Returns true if no record has been read because the end of the file was reached.
    pub fn is_eof(&self) -> bool {
        match self {
            &ReadError::NoMoreRecord => true,
            _ => false
        }
    }
}


impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.description().fmt(f)
    }
}


impl Error for ReadError {
    fn description(&self) -> &str {
        match self {
            &ReadError::Invalid => "invalid record",
            &ReadError::NoMoreRecord => "no more record"
        }
    }
}


#[cfg(test)]
mod tests {
    extern crate tempdir;
    use super::*;
    use std::path::Path;

    fn _test_read<P: AsRef<Path>>(path: &P) {
        let bcf = Reader::new(path).ok().expect("Error opening file.");
        assert_eq!(bcf.header.samples(), [b"NA12878.subsample-0.25-0"]);

        for (i, rec) in bcf.records().enumerate() {
            let mut record = rec.ok().expect("Error reading record.");
            assert_eq!(record.sample_count(), 1);

            assert_eq!(record.rid().expect("Error reading rid."), 0);
            assert_eq!(record.pos(), 10021 + i as u32);
            assert_eq!(record.qual(), 0f32);
            assert_eq!(record.info(b"MQ0F").float().ok().expect("Error reading info."), [1.0]);
            if i == 59 {
                assert_eq!(record.info(b"SGB").float().ok().expect("Error reading info."), [-0.379885]);
            }
            // the artificial "not observed" allele is present in each record.
            assert_eq!(record.alleles().iter().last().unwrap(), b"<X>");

            let mut fmt = record.format(b"PL");
            let pl = fmt.integer().ok().expect("Error reading format.");
            assert_eq!(pl.len(), 1);
            if i == 59 {
                assert_eq!(pl[0].len(), 6);
            }
            else {
                assert_eq!(pl[0].len(), 3);
            }
        }
    }

    #[test]
    fn test_read() {
        _test_read(&"test/test.bcf");
    }

    #[test]
    fn test_write() {
        let bcf = Reader::new(&"test/test_multi.bcf").ok().expect("Error opening file.");
        let tmp = tempdir::TempDir::new("rust-htslib").ok().expect("Cannot create temp dir");
        let bcfpath = tmp.path().join("test.bcf");
        println!("{:?}", bcfpath);
        {
            let header = Header::subset_template(&bcf.header, &[b"NA12878.subsample-0.25-0"]).ok().expect("Error subsetting samples.");
            let mut writer = Writer::new(&bcfpath, &header, false, false).ok().expect("Error opening file.");
            for rec in bcf.records() {
                let mut record = rec.ok().expect("Error reading record.");
                writer.translate(&mut record);
                writer.subset(&mut record);
                record.trim_alleles().ok().expect("Error trimming alleles.");
                writer.write(&record).ok().expect("Error writing record");
            }
        }
        {
            _test_read(&bcfpath);
        }
        tmp.close().ok().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_strings() {
        let vcf = Reader::new(&"test/test_string.vcf").ok().expect("Error opening file.");
        let fs1 = [&b"LongString1"[..], &b"LongString2"[..], &b"."[..], &b"LongString4"[..], &b"evenlength"[..], &b"ss6"[..]];
        for (i, rec) in vcf.records().enumerate() {
            println!("record {}", i);
            let mut record = rec.ok().expect("Error reading record.");
            assert_eq!(record.info(b"S1").string().ok().expect("Error reading string.")[0], format!("string{}", i + 1).as_bytes());
            println!("{}", String::from_utf8_lossy(record.format(b"FS1").string().ok().expect("Error reading string.")[0]));
            assert_eq!(record.format(b"FS1").string().ok().expect("Error reading string.")[0], fs1[i]);
        }
    }

    #[test]
    fn test_missing() {
        let vcf = Reader::new(&"test/test_missing.vcf").ok().expect("Error opening file.");
        let fn4 = [&[record::MISSING_INTEGER, record::MISSING_INTEGER, record::MISSING_INTEGER, record::MISSING_INTEGER][..], &[record::MISSING_INTEGER][..]];
        let f1 = [false, true];
        for (i, rec) in vcf.records().enumerate() {
            let mut record = rec.ok().expect("Error reading record.");
            assert_eq!(record.info(b"F1").float().ok().expect("Error reading float.")[0].is_nan(), f1[i]);
            assert_eq!(record.format(b"FN4").integer().ok().expect("Error reading integer.")[1], fn4[i]);
            println!("{:?}", record.format(b"FF4").float().ok().expect("Error reading float.")[1]);
            assert!(record.format(b"FF4").float().ok().expect("Error reading float.")[1].iter().all(|v| v.is_nan()));
        }
    }
}
