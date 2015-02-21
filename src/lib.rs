#![feature(libc)]
extern crate libc;
pub mod htslib;
pub mod bam;

#[cfg(test)]
mod tests {
    use htslib;
    use bam;
    use std::str;

    #[test]
    fn test_record() {
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

        let samfile = bam::Samfile::new(b"test.bam");

        for (((record, &name), &flag), &seq) in samfile.records()
                                                       .zip(names.iter())
                                                       .zip(flags.iter())
                                                       .zip(seqs.iter()) {
            let rec = record.ok().expect("Expected valid record");
            println!("{}", str::from_utf8(rec.qname()).ok().unwrap());
            //println!("{}", str::from_utf8(rec.seq()).ok().unwrap());
            assert_eq!(rec.qname(), name);
            assert_eq!(rec.flag(), flag);
            assert_eq!(rec.seq(), seq);
        }


//        println!("{}", String::from_utf8_lossy(record.aux(b"MD").ok().unwrap().string()).as_slice());
//        println!("{}", record.aux(b"SM").ok().unwrap().integer());

    }
}
