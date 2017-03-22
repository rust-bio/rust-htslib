extern crate rust_htslib;
use self::rust_htslib::sam::SAMWriter;

use std::io::Write;


pub fn main(){
  let example_bam = "./test/bam2sam_test.bam";
  let result = SAMWriter::from_bam_with_filter(example_bam, &"-", |_|{Some(true)});	
  assert!(writeln!(&mut std::io::stderr(), "{:?}", result).is_ok());
}
