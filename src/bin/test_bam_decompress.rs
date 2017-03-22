extern crate rust_htslib;
use self::rust_htslib::sam::SAMWriter;
use self::rust_htslib::bam::{Reader, Read};

use std::io::Write;

/// Read bam file. For each record apply f to it, and write to sam file if f returned Some(true), skip record if Some(false) if None then terminate iteration
///
/// # Arguments
///
/// * `bamfile` - the bam file to read from
/// * `samfile` - the sam file to write
/// * `f` - the predicate to apply
use self::rust_htslib::bam::header;
use self::rust_htslib::bam::record;
pub fn from_bam_with_filter<'a, 'b, F>(bamfile:&'a str, samfile:&'b str, f:F) -> bool where F:Fn(&record::Record) -> Option<bool> {
    let bam_reader = if bamfile != "-" {
        Reader::from_path(bamfile).unwrap()
    } else {
        Reader::from_stdin().unwrap()
    };
    let header = header::Header::from_template(bam_reader.header());
    let mut sam_writer = if samfile != "-" {
            SAMWriter::from_path(samfile, &header).unwrap()
        } else {
            SAMWriter::from_stdout(&header).unwrap()
        };
    for record in bam_reader.records() {
        if record.is_err() {
            return false;
        } 
        let parsed = record.unwrap();
        match f(&parsed) {
            None => return true,
            Some(false) => {},
            Some(true) => if let Err(_) = sam_writer.write(&parsed) {
                return false;
            }
        }
    }
    true
}

pub fn main(){
  let example_bam = "./test/bam2sam_test.bam";
  let result = from_bam_with_filter(example_bam, &"-", |_|{Some(true)});	
  assert!(writeln!(&mut std::io::stderr(), "{:?}", result).is_ok());
}
