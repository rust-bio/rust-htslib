#![feature(libc)]
extern crate libc;
pub mod htslib;
pub mod bam;

#[cfg(test)]
mod tests {
    use htslib;
    use bam;

    #[test]
    fn test_record() {
        let names = [b"I", b"II.14978392", b"III", b"IV", b"V", b"VI"];
        let flags = [16u16, 16u16, 16u16, 16u16, 16u16, 2048u16];
        let f = bam::Samfile::new(b"test.bam");
        let records: Vec<bam::Record> = f.take(6).collect();
        assert!(records.len() == 6);

        for ((record, &name), &flag) in records.iter().zip(names.iter()).zip(flags.iter()) {
            assert_eq!(record.qname(), name);
            assert_eq!(record.flag(), flag);
            println!("{:?}", String::from_utf8_lossy(record.qname()));
            println!("{:?}", record.flag());
        }


//        println!("{}", String::from_utf8_lossy(record.aux(b"MD").ok().unwrap().string()).as_slice());
//        println!("{}", record.aux(b"SM").ok().unwrap().integer());

    }
}
