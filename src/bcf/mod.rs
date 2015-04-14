
use std::ffi;
use std::ffi::AsOsStr;
use std::convert::AsRef;
use std::path::Path;
use std::slice;

pub mod record;

use htslib;


pub struct HeaderView {
    inner: *mut htslib::vcf::bcf_hdr_t,
}


impl HeaderView {
    fn new(inner: *mut htslib::vcf::bcf_hdr_t) -> Self {
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


pub struct Reader {
    inner: *mut htslib::vcf::htsFile,
    header: HeaderView,
}

impl Reader {
   pub fn new<P: AsRef<Path>>(path: &P) -> Self {
        let htsfile = bcf_open(path, b"r");
        let header = unsafe { htslib::vcf::bcf_hdr_read(htsfile) };
        Reader { inner: htsfile, header: HeaderView::new(header) }
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


pub struct Writer<'a> {
    inner: *mut htslib::vcf::htsFile,
    header: &'a HeaderView,
}


impl<'a> Writer<'a> {
    pub fn with_template<P: AsRef<Path>>(template: &'a Reader, path: &P) -> Self {
        let htsfile = bcf_open(path, b"w");
        unsafe { htslib::vcf::bcf_hdr_write(htsfile, template.header.inner) };
        Writer { inner: htsfile, header: &template.header }
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


impl<'a> Drop for Writer<'a> {
    fn drop(&mut self) {
        unsafe {
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


/// Wrapper for opening a BCF file.
fn bcf_open<P: AsRef<Path>>(path: &P, mode: &[u8]) -> *mut htslib::vcf::htsFile {
    unsafe {
        htslib::vcf::hts_open(
            path.as_ref().as_os_str().to_cstring().unwrap().as_ptr(),
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
    extern crate tempdir;
    use super::*;
    use std::path::Path;

    fn _test_read<P: AsRef<Path>>(path: &P) {
        let bcf = Reader::new(path);
        assert_eq!(bcf.header.samples(), [b"NA12878.subsample-0.25-0"]);

        for (i, rec) in bcf.records().enumerate() {
            let mut record = rec.ok().expect("Error reading record.");

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
        _test_read(&"test.bcf");
    }

    #[test]
    fn test_write() {
        let bcf = Reader::new(&"test.bcf");
        let tmp = tempdir::TempDir::new("rust-htslib").ok().expect("Cannot create temp dir");
        let bcfpath = tmp.path().join("test.bcf");
        println!("{:?}", bcfpath);
        {
            let mut writer = Writer::with_template(&bcf, &bcfpath);
            for record in bcf.records() {
                writer.write(&record.ok().expect("Error reading record.")).ok().expect("Error writing record");
            }
        }
        {
            _test_read(&bcfpath);
        }
        tmp.close().ok().expect("Failed to delete temp dir");
    }
}
