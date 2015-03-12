// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


pub mod record;
pub mod header;

use std::ffi;
use std::path;
use std::ffi::AsOsStr;
use std::os::unix::prelude::OsStrExt;

use htslib;


/// A trait for a BAM reader with a read method.
pub trait BAMRead {
    /// Read next BAM record into given record.
    /// Use this method in combination with a single allocated record to avoid the reallocations
    /// occurring with the iterator.
    ///
    /// # Arguments
    ///
    /// * `record` - the record to be filled
    fn read(&self, record: &mut record::Record) -> Result<(), ReadError>;
}


/// A BAM reader.
pub struct BAMReader {
    f: *mut htslib::Struct_BGZF,
    header: *mut htslib::bam_hdr_t,
}


impl BAMReader {
    /// Create a new BAMReader.
    ///
    /// # Arguments
    ///
    /// * `path` - the path. Use "-" for stdin.
    pub fn new<P: path::AsPath>(path: P) -> Self {
        let f = bgzf_open(path, b"r");
        let header = unsafe { htslib::bam_hdr_read(f) };
        BAMReader { f : f, header : header }
    }

    /// Iterator over the records of the BAM file.
    pub fn records(self) -> Records {
        Records { bam: self }
    }
}


impl BAMRead for BAMReader {
    fn read(&self, record: &mut record::Record) -> Result<(), ReadError> {
        match unsafe { htslib::bam_read1(self.f, &mut record.inner) } {
            -1 => Err(ReadError::EOF),
            -2 => Err(ReadError::Truncated),
            -4 => Err(ReadError::Invalid),
            _  => Ok(())
        }
    }
}


impl Drop for BAMReader {
    fn drop(&mut self) {
        unsafe {
            htslib::bam_hdr_destroy(self.header);
            htslib::bgzf_close(self.f);
        }
    }
}


/// A BAM writer.
pub struct BAMWriter {
    f: *mut htslib::Struct_BGZF,
    header: *mut htslib::bam_hdr_t,
}


impl BAMWriter {
    /// Create a new BAM file.
    ///
    /// # Arguments
    ///
    /// * `path` - the path. Use "-" for stdin.
    /// * `header` - header definition to use
    pub fn new<P: path::AsPath>(path: &P, header: &header::Header) -> Self {
        let f = bgzf_open(path, b"w");
        let header_record = unsafe {
            let header_string = header.to_bytes();
            htslib::sam_hdr_parse(
                header_string.len() as i32,
                ffi::CString::new(header_string).unwrap().as_ptr()
            )
        };
        unsafe { htslib::bam_hdr_write(f, header_record); }
        BAMWriter { f: f, header: header_record }
    }

    /// Create a new BAM file from template.
    ///
    /// # Arguments
    ///
    /// * `path` - the path. Use "-" for stdin.
    /// * `template` - the template BAM. Use "-" for stdin.
    pub fn with_template<P: path::AsPath, T: path::AsPath>(path: &P, template: &T) -> Self {
        let t = bgzf_open(template, b"r");
        let header = unsafe { htslib::bam_hdr_read(t) };
        unsafe { htslib::bgzf_close(t); }
        let f = bgzf_open(path, b"w");
        unsafe { htslib::bam_hdr_write(f, header); }
        BAMWriter { f: f, header: header }
    }

    /// Write record to BAM.
    ///
    /// # Arguments
    ///
    /// * `record` - the record to write
    pub fn write(&mut self, record: &record::Record) -> Result<(), ()> {
        if unsafe { htslib::bam_write1(self.f, &record.inner) } == -1 {
            Err(())
        }
        else {
            Ok(())
        }
    }
}


impl Drop for BAMWriter {
    fn drop(&mut self) {
        unsafe {
            htslib::bam_hdr_destroy(self.header);
            htslib::bgzf_close(self.f);
        }
    }
}


/// Iterator over the records of a BAM.
pub struct Records {
    bam: BAMReader
}


impl Iterator for Records {
    type Item = Result<record::Record, ReadError>;

    fn next(&mut self) -> Option<Result<record::Record, ReadError>> {
        let mut record = record::Record::new();
        match self.bam.read(&mut record) {
            Err(ReadError::EOF) => None,
            Ok(())   => Some(Ok(record)),
            Err(err) => Some(Err(err))
        }
    }
}


pub enum ReadError {
    Truncated,
    Invalid,
    EOF,
}


/// Wrapper for opening a BAM file.
fn bgzf_open<P: path::AsPath>(path: &P, mode: &[u8]) -> *mut htslib::Struct_BGZF {
    unsafe {
        htslib::bgzf_open(
            path.as_path().as_os_str().to_cstring().unwrap().as_ptr(),
            ffi::CString::new(mode).unwrap().as_ptr()
        )
    }
}


#[cfg(test)]
mod tests {
    use super::record::*;
    use super::*;
    use std::str;

    #[test]
    fn test_read() {
        let names = [b"I", b"II.14978392", b"III", b"IV", b"V", b"VI"];
        let flags = [16u16, 16u16, 16u16, 16u16, 16u16, 2048u16];
        let seqs = [
            b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA",
            b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA",
            b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA",
            b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA",
            b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA",
            b"ACTAAGCCTAAGCCTAAGCCTAAGCCAATTATCGATTTCTGAAAAAATTATCGAATTTTCTAGAAATTTTGCAAATTTTTTCATAAAATTATCGATTTTA",
        ];
        let cigars = [
            [Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)],
            [Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)],
            [Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)],
            [Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)],
            [Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)],
            [Cigar::Match(27), Cigar::Del(100000), Cigar::Match(73)],
        ];

        let bam = BAMReader::new("test.bam");

        for (i, record) in bam.records().enumerate() {
            let rec = record.ok().expect("Expected valid record");
            println!("{}", str::from_utf8(rec.qname()).ok().unwrap());
            assert_eq!(rec.qname(), names[i]);
            assert_eq!(rec.flags(), flags[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);
            assert_eq!(rec.cigar(), cigars[i]);
        }
    }

    #[test]
    fn test_set_record() {
        let qname = b"I";
        let cigar = [Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)];
        let seq =  b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGC";
        let qual = b"JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ";

        let mut rec = record::Record::new();
        rec.set_reverse();
        rec.set(qname, &cigar, seq, qual);
        rec.push_aux(b"NM", &Aux::Integer(15));
        assert_eq!(rec.qname(), qname);
        assert_eq!(rec.cigar(), cigar);
        assert_eq!(rec.seq().as_bytes(), seq);
        assert_eq!(rec.qual(), qual);
        assert!(rec.is_reverse());
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));
    }
}
