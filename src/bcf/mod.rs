// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::ffi;
use std::path::Path;
use std::rc::Rc;

use url::Url;

pub mod record;
pub mod header;
pub mod buffer;

use htslib;
use bcf::header::{Id, HeaderView, SampleSubset};

pub use bcf::header::Header;
pub use bcf::record::Record;
pub use bcf::buffer::RecordBuffer;


/// Redefinition of corresponding `#define` in `vcf.h.`.
#[allow(non_upper_case_globals)]
pub const GT_MISSING: i32 = 0;


/// A VCF/BCF reader.
#[derive(Debug)]
pub struct Reader {
    inner: *mut htslib::htsFile,
    header: Rc<HeaderView>,
}


unsafe impl Send for Reader {}


/// Implementation for `Reader::set_threads()` and `Writer::set_threads`.
pub fn set_threads(hts_file: *mut htslib::htsFile, n_threads: usize)
        -> Result<(), ThreadingError> {
    assert!(n_threads > 0, "n_threads must be > 0");

    let r = unsafe { htslib::hts_set_threads(hts_file, n_threads as i32) };
    if r != 0 {
        Err(ThreadingError::Some)
    } else {
        Ok(())
    }
}


impl Reader {
    /// Create a new reader from a given path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, BCFPathError> {
        match path.as_ref().to_str() {
            Some(p) if path.as_ref().exists() => {
                Ok(try!(Self::new(p.as_bytes())))
            },
            _ => {
                Err(BCFPathError::InvalidPath)
            }
        }
    }

    /// Create a new reader from a given URL.
    pub fn from_url(url: &Url) -> Result<Self, BCFError> {
        Self::new(url.as_str().as_bytes())
    }

    /// Create a new reader from standard input.
    pub fn from_stdin() -> Result<Self, BCFError> {
        Self::new(b"-")
    }

    fn new(path: &[u8]) -> Result<Self, BCFError> {
        let htsfile = try!(bcf_open(path, b"r"));
        let header = unsafe { htslib::bcf_hdr_read(htsfile) };
        Ok(Reader { inner: htsfile, header: Rc::new(HeaderView::new(header)) })
    }

    /// Get header of the read VCF/BCF file.
    pub fn header(&self) -> &HeaderView {
        &self.header
    }

    /// Create empty record for this reader.
    /// The record can be reused multiple times.
    pub fn empty_record(&self) -> Record {
        record::Record::new(self.header.clone())
    }

    /// Read the next record.
    ///
    /// # Arguments
    /// * record - an empty record, that can be created with `bcf::Reader::empty_record`.
    pub fn read(&mut self, record: &mut record::Record) -> Result<(), ReadError> {
        match unsafe { htslib::bcf_read(self.inner, self.header.inner, record.inner) } {
            0  => {
                record.set_header(self.header.clone());
                Ok(())
            },
            -1 => Err(ReadError::NoMoreRecord),
            _  => Err(ReadError::Invalid),
        }
    }

    /// Return an iterator over all records of the VCF/BCF file.
    pub fn records(&mut self) -> Records {
        Records { reader: self }
    }

    /// Activate multi-threaded BCF read support in htslib. This should permit faster
    /// reading of large BCF files.
    ///
    /// # Arguments
    ///
    /// * `n_threads` - number of extra background reader threads to use, must be `> 0`.
    pub fn set_threads(&mut self, n_threads: usize) -> Result<(), ThreadingError> {
        set_threads(self.inner, n_threads)
    }
}


impl Drop for Reader {
    fn drop(&mut self) {
        unsafe {
            htslib::hts_close(self.inner);
        }
    }
}


/// A VCF/BCF writer.
#[derive(Debug)]
pub struct Writer {
    inner: *mut htslib::htsFile,
    header: Rc<HeaderView>,
    subset: Option<SampleSubset>,
}


unsafe impl Send for Writer {}


impl Writer {
    /// Create a new writer that writes to the given path.
    ///
    /// # Arguments
    ///
    /// * `path` - the path.
    /// * `header` - header definition to use
    /// * `uncompressed` - disable compression
    /// * `vcf` - write VCF instead of BCF
    pub fn from_path<P: AsRef<Path>>(path: P, header: &Header, uncompressed: bool, vcf: bool) -> Result<Self, BCFPathError> {
        if let Some(p) = path.as_ref().to_str() {
            Ok(try!(Self::new(p.as_bytes(), header, uncompressed, vcf)))
        } else {
            Err(BCFPathError::InvalidPath)
        }
    }

    pub fn from_url(url: &Url, header: &Header, uncompressed: bool, vcf: bool) -> Result<Self, BCFError> {
        Self::new(url.as_str().as_bytes(), header, uncompressed, vcf)
    }

    pub fn from_stdout(header: &Header, uncompressed: bool, vcf: bool) -> Result<Self, BCFError> {
        Self::new(b"-", header, uncompressed, vcf)
    }

    fn new(path: &[u8], header: &Header, uncompressed: bool, vcf: bool) -> Result<Self, BCFError> {
        let mode: &[u8] = match (uncompressed, vcf) {
            (true, true)   => b"w",
            (false, true)  => b"wz",
            (true, false)  => b"wu",
            (false, false) => b"wb",
        };

        let htsfile = try!(bcf_open(path, mode));
        unsafe { htslib::bcf_hdr_write(htsfile, header.inner) };
        Ok(Writer {
            inner: htsfile,
            header: Rc::new(HeaderView::new(unsafe { htslib::bcf_hdr_dup(header.inner) })),
            subset: header.subset.clone()
        })
    }

    pub fn header(&self) -> &HeaderView {
        &self.header
    }

    /// Create empty record for writing to this writer.
    /// The record can be reused multiple times.
    pub fn empty_record(&self) -> Record {
        record::Record::new(self.header.clone())
    }

    /// Translate record to header of this writer.
    pub fn translate(&mut self, record: &mut record::Record) {
        unsafe {
            htslib::bcf_translate(self.header.inner, record.header().inner, record.inner);
        }
        record.set_header(self.header.clone());
    }

    /// Subset samples of record to match header of this writer.
    pub fn subset(&mut self, record: &mut record::Record) {
        match self.subset {
            Some(ref mut subset) => unsafe {
                htslib::bcf_subset(self.header.inner, record.inner, subset.len() as i32, subset.as_mut_ptr());
            },
            None         => ()
        }
    }

    pub fn write(&mut self, record: &record::Record) -> Result<(), WriteError> {
        if unsafe { htslib::bcf_write(self.inner, self.header.inner, record.inner) } == -1 {
            Err(WriteError::Some)
        }
        else {
            Ok(())
        }
    }

    /// Activate multi-threaded BCF write support in htslib. This should permit faster
    /// writing of large BCF files.
    ///
    /// # Arguments
    ///
    /// * `n_threads` - number of extra background writer threads to use, must be `> 0`.
    pub fn set_threads(&mut self, n_threads: usize) -> Result<(), ThreadingError> {
        set_threads(self.inner, n_threads)
    }
}


impl Drop for Writer {
    fn drop(&mut self) {
        unsafe {
            htslib::hts_close(self.inner);
        }
    }
}


#[derive(Debug)]
pub struct Records<'a> {
    reader: &'a mut Reader,
}


impl<'a> Iterator for Records<'a> {
    type Item = Result<record::Record, ReadError>;

    fn next(&mut self) -> Option<Result<record::Record, ReadError>> {
        let mut record = self.reader.empty_record();
        match self.reader.read(&mut record) {
            Err(ReadError::NoMoreRecord) => None,
            Err(e)                       => Some(Err(e)),
            Ok(())                       => Some(Ok(record)),
        }
    }
}


quick_error! {
    #[derive(Debug, Clone)]
    pub enum BCFError {
        Some {
            description("error reading BCF/VCF file")
        }
    }
}


quick_error! {
    #[derive(Debug, Clone)]
    pub enum BCFPathError {
        InvalidPath {
            description("invalid path")
        }
        BCFError(err: BCFError) {
            from()
        }
    }
}


quick_error! {
    #[derive(Debug, Clone)]
    pub enum ThreadingError {
        Some {
            description("error setting threads for multi-threaded I/O")
        }
    }
}


/// Wrapper for opening a BCF file.
fn bcf_open(path: &[u8], mode: &[u8]) -> Result<*mut htslib::htsFile, BCFError> {
    let p = ffi::CString::new(path).unwrap();
    let ret = unsafe {
        htslib::hts_open(
            p.as_ptr(),
            ffi::CString::new(mode).unwrap().as_ptr()
        )
    };
    if ret.is_null() {
        Err(BCFError::Some)
    } else {
        Ok(ret)
    }
}


quick_error! {
    #[derive(Debug, Clone)]
    pub enum ReadError {
        Invalid {
            description("invalid record")
        }
        NoMoreRecord {
            description("no more record")
        }
    }
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


quick_error! {
    #[derive(Debug, Clone)]
    pub enum WriteError {
        Some {
            description("failed to write record")
        }
    }
}


#[cfg(test)]
mod tests {
    extern crate tempdir;
    use super::*;
    use std::path::Path;
    use bcf::record::Numeric;

    fn _test_read<P: AsRef<Path>>(path: &P) {
        let mut bcf = Reader::from_path(path).ok().expect("Error opening file.");
        assert_eq!(bcf.header.samples(), [b"NA12878.subsample-0.25-0"]);

        for (i, rec) in bcf.records().enumerate() {
            let mut record = rec.ok().expect("Error reading record.");
            assert_eq!(record.sample_count(), 1);

            assert_eq!(record.rid().expect("Error reading rid."), 0);
            assert_eq!(record.pos(), 10021 + i as u32);
            assert_eq!(record.qual(), 0f32);
            assert_eq!(record.info(b"MQ0F").float().ok().expect("Error reading info.").expect("Missing tag"), [1.0]);
            if i == 59 {
                assert_eq!(record.info(b"SGB").float().ok().expect("Error reading info.").expect("Missing tag"), [-0.379885]);
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
    fn test_reader_set_threads() {
        let path = &"test/test.bcf";
        let mut bcf = Reader::from_path(path).ok().expect("Error opening file.");
        bcf.set_threads(2).unwrap();
    }

    #[test]
    fn test_writer_set_threads() {
        let path = &"test/test.bcf";
        let tmp = tempdir::TempDir::new("rust-htslib").ok().expect("Cannot create temp dir");
        let bcfpath = tmp.path().join("test.bcf");
        let bcf = Reader::from_path(path).ok().expect("Error opening file.");
        let header = Header::subset_template(&bcf.header, &[b"NA12878.subsample-0.25-0"]).ok().expect("Error subsetting samples.");
        let mut writer = Writer::from_path(&bcfpath, &header, false, false).ok().expect("Error opening file.");
        writer.set_threads(2).unwrap();
    }

    #[test]
    fn test_write() {
        let mut bcf = Reader::from_path(&"test/test_multi.bcf").ok().expect("Error opening file.");
        let tmp = tempdir::TempDir::new("rust-htslib").ok().expect("Cannot create temp dir");
        let bcfpath = tmp.path().join("test.bcf");
        println!("{:?}", bcfpath);
        {
            let header = Header::subset_template(&bcf.header, &[b"NA12878.subsample-0.25-0"]).ok().expect("Error subsetting samples.");
            let mut writer = Writer::from_path(&bcfpath, &header, false, false).ok().expect("Error opening file.");
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
        let mut vcf = Reader::from_path(&"test/test_string.vcf").ok().expect("Error opening file.");
        let fs1 = [&b"LongString1"[..], &b"LongString2"[..], &b"."[..], &b"LongString4"[..], &b"evenlength"[..], &b"ss6"[..]];
        for (i, rec) in vcf.records().enumerate() {
            println!("record {}", i);
            let mut record = rec.ok().expect("Error reading record.");
            assert_eq!(record.info(b"S1").string().ok().expect("Error reading string.").expect("Missing tag")[0], format!("string{}", i + 1).as_bytes());
            println!("{}", String::from_utf8_lossy(record.format(b"FS1").string().ok().expect("Error reading string.")[0]));
            assert_eq!(record.format(b"FS1").string().ok().expect("Error reading string.")[0], fs1[i]);
        }
    }

    #[test]
    fn test_missing() {
        let mut vcf = Reader::from_path(&"test/test_missing.vcf").ok().expect("Error opening file.");
        let fn4 = [&[i32::missing(), i32::missing(), i32::missing(), i32::missing()][..], &[i32::missing()][..]];
        let f1 = [false, true];
        for (i, rec) in vcf.records().enumerate() {
            let mut record = rec.ok().expect("Error reading record.");
            assert_eq!(record.info(b"F1").float().ok().expect("Error reading float.").expect("Missing tag")[0].is_nan(), f1[i]);
            assert_eq!(record.format(b"FN4").integer().ok().expect("Error reading integer.")[1], fn4[i]);
            assert!(record.format(b"FF4").float().ok().expect("Error reading float.")[1].iter().all(|&v| v.is_missing()));
        }
    }

    #[test]
    fn test_genotypes() {
        let mut vcf = Reader::from_path(&"test/test_string.vcf").ok().expect("Error opening file.");
        let expected = ["./1", "1|1", "0/1", "0|1", "1|.", "1/1"];
        for (rec, exp_gt) in vcf.records().zip(expected.into_iter()) {
            let mut rec = rec.ok().expect("Error reading record.");
            let genotypes = rec.genotypes().expect("Error reading genotypes");
            assert_eq!(&format!("{}", genotypes.get(0)), exp_gt);
        }
    }

    #[test]
    fn test_header_ids() {
        let vcf = Reader::from_path(&"test/test_string.vcf").ok().expect("Error opening file.");
        let header = &vcf.header();
        use bcf::header::Id;

        assert_eq!(header.id_to_name(Id(4)), b"GT");
        assert_eq!(header.name_to_id(b"GT").unwrap(), Id(4));
        assert!(header.name_to_id(b"XX").is_err());
    }

    #[test]
    fn test_header_samples() {
        let vcf = Reader::from_path(&"test/test_string.vcf").ok().expect("Error opening file.");
        let header = &vcf.header();

        assert_eq!(header.id_to_sample(Id(0)), b"one");
        assert_eq!(header.id_to_sample(Id(1)), b"two");
        assert_eq!(header.sample_to_id(b"one").unwrap(), Id(0));
        assert_eq!(header.sample_to_id(b"two").unwrap(), Id(1));
        assert!(header.sample_to_id(b"three").is_err());
    }
}
