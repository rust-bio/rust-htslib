use std::slice;
use std::ffi::CStr;
use std::mem;
use std::iter;
use libc;

use htslib;

#[derive(Debug)]
pub enum Aux<'a> {
    Integer(i32),
    String(&'a [u8]),
    Float(f64),
    Char(u8),
}


impl<'a> Aux<'a> {
    pub fn string(&self) -> &'a [u8] {
        match *self {
            Aux::String(x) => x,
            _ => panic!("not a string"),
        }
    }

    pub fn float(&self) -> f64 {
        match *self {
            Aux::Float(x) => x,
            _ => panic!("not a float"),
        }
    }

    pub fn integer(&self) -> i32 {
        match *self {
            Aux::Integer(x) => x,
            _ => panic!("not an integer"),
        }
    }

    pub fn char(&self) -> u8 {
        match *self {
            Aux::Char(x) => x,
            _ => panic!("not a character"),
        }
    }
}


pub struct Record {
    b: htslib::bam1_t,
}


impl Record {
    pub fn new() -> Self {
        let mut b = unsafe { *htslib::bam_init1() };
        b.m_data = 0;
        Record { b: b }
    }

    fn data(&self) -> &[u8] {
        unsafe { slice::from_raw_parts(self.b.data, self.b.l_data as usize) }
    }

    pub fn tid(&self) -> i32 {
        self.b.core.tid
    }

    pub fn set_tid(&mut self, tid: i32) {
        self.b.core.tid = tid;
    }

    pub fn pos(&self) -> i32 {
        self.b.core.pos
    }

    pub fn set_pos(&mut self, pos: i32) {
        self.b.core.pos = pos;
    }

    pub fn bin(&self) -> u16 {
        self.b.core.bin
    }

    pub fn set_bin(&mut self, bin: u16) {
        self.b.core.bin = bin;
    }

    pub fn mapq(&self) -> u8 {
        self.b.core.qual
    }

    pub fn set_mapq(&mut self, mapq: u8) {
        self.b.core.qual = mapq;
    }

    pub fn flags(&self) -> u16 {
        self.b.core.flag
    }

    pub fn set_flags(&mut self, flags: u16) {
        self.b.core.flag = flags;
    }

    pub fn mtid(&self) -> i32 {
        self.b.core.mtid
    }

    pub fn set_mtid(&mut self, mtid: i32) {
        self.b.core.mtid = mtid;
    }

    pub fn mpos(&self) -> i32 {
        self.b.core.mpos
    }

    pub fn set_mpos(&mut self, mpos: i32) {
        self.b.core.mpos = mpos;
    }

    pub fn insert_size(&self) -> i32 {
        self.b.core.isize
    }

    pub fn set_insert_size(&mut self, insert_size: i32) {
        self.b.core.isize = insert_size;
    }

    fn qname_len(&self) -> usize {
        self.b.core.l_qname as usize
    }

    pub fn qname(&self) -> &[u8] {
        &self.data()[..self.qname_len()-1] // -1 ignores the termination symbol
    }

    pub fn set(&mut self, qname: &[u8], cigar: &[Cigar], seq: &[u8], qual: &[u8]) {
        self.b.l_data = (qname.len() + 1 + cigar.len() * 4 + seq.len() / 2 + qual.len()) as i32;

        if self.b.m_data < self.b.l_data {
            
            self.b.m_data = self.b.l_data;
            self.b.m_data += 32 - self.b.m_data % 32;
            unsafe {
                self.b.data = libc::funcs::c95::stdlib::realloc(
                    self.b.data as *mut libc::types::common::c95::c_void, self.b.m_data as u64
                ) as *mut u8;
            }
        }

        let mut data = unsafe { slice::from_raw_parts_mut(self.b.data, self.b.l_data as usize) };
        // qname
        slice::bytes::copy_memory(data, qname);
        data[qname.len()] = b'\0';
        let mut i = qname.len() + 1;
        self.b.core.l_qname = i as u8;

        // cigar
        {
            let mut cigar_data = unsafe {
                 slice::from_raw_parts_mut(data[i..].as_ptr() as *mut u32, cigar.len())
            };
            for (i, c) in cigar.iter().enumerate() {
                cigar_data[i] = c.encode();
            }
            self.b.core.n_cigar = cigar.len() as u16;
            i += cigar.len() * 4;
        }

        // seq
        {
            for j in iter::range_step(0, seq.len(), 2) {
                data[i + j / 2] = ENCODE_BASE[seq[j] as usize] << 4 | ENCODE_BASE[seq[j + 1] as usize];
            }
            self.b.core.l_qseq = seq.len() as i32;
            i += (seq.len() + 1) / 2;
        }

        // qual
        slice::bytes::copy_memory(&mut data[i..], qual);
    }

    fn cigar_len(&self) -> usize {
        self.b.core.n_cigar as usize
    }

    pub fn cigar(&self) -> Vec<Cigar> {
        let raw = unsafe { slice::from_raw_parts(self.data()[self.qname_len()..].as_ptr() as *const u32, self.cigar_len()) };
        raw.iter().map(|&c| {
            let len = c >> 4;
            match c & 0b1111 {
                0 => Cigar::Match(len),
                1 => Cigar::Ins(len),
                2 => Cigar::Del(len),
                3 => Cigar::RefSkip(len),
                4 => Cigar::SoftClip(len),
                5 => Cigar::HardClip(len),
                6 => Cigar::Pad(len),
                7 => Cigar::Equal(len),
                8 => Cigar::Diff(len),
                9 => Cigar::Back(len),
                _ => panic!("Unexpected cigar type"),
            }
        }).collect()
    }

    fn seq_len(&self) -> usize {
        self.b.core.l_qseq as usize
    }

    pub fn seq(&self) -> Seq {
        Seq { 
            encoded: &self.data()
                        [self.qname_len() + self.cigar_len()*4..]
                        [..(self.seq_len() + 1) / 2]
        }
    }

    pub fn qual(&self) -> &[u8] {
        &self.data()[self.qname_len() + self.cigar_len()*4 + (self.seq_len()+1)/2..][..self.seq_len()]
    }

    pub fn aux(&self, name: &[u8]) -> Option<Aux> {
        let aux = unsafe { htslib::bam_aux_get(&self.b, name.as_ptr() as *mut i8 ) };

        unsafe {
            match *aux {
                b'c'|b'C'|b's'|b'S'|b'i'|b'I' => Some(Aux::Integer(htslib::bam_aux2i(aux))),
                b'f'|b'd' => Some(Aux::Float(htslib::bam_aux2f(aux))),
                b'A' => Some(Aux::Char(htslib::bam_aux2A(aux) as u8)),
                b'Z'|b'H' => {
                    let f = aux.offset(1) as *const i8;
                    let x = CStr::from_ptr(f).to_bytes();
                    Some(Aux::String(mem::copy_lifetime(self, x)))
                },
                _ => None,
            }
        }
    }

    pub fn is_paired(&self) -> bool {
        (self.flags() & 1u16) > 0
    }

    pub fn set_is_paired(&mut self) {
        self.b.core.flag | 1u16;
    }

    pub fn is_proper_pair(&self) -> bool{
        (self.flags() & 2u16) > 0
    }

    pub fn is_unmapped(&self) -> bool{
        (self.flags() & 4u16) > 0
    }

    pub fn mate_is_unmapped(&self) -> bool{
        (self.flags() & 8u16) > 0
    }

    pub fn is_reverse(&self) -> bool{
        (self.flags() & 16u16) > 0
    }

    pub fn mate_is_reverse(&self) -> bool{
        (self.flags() & 32u16) > 0
    }

    pub fn is_first_in_pair(&self) -> bool{
        (self.flags() & 64u16) > 0
    }

    pub fn is_second_in_pair(&self) -> bool{
        (self.flags() & 128u16) > 0
    }

    pub fn is_secondary(&self) -> bool{
        (self.flags() & 256u16) > 0
    }

    pub fn is_quality_check_failed(&self) -> bool{
        (self.flags() & 512u16) > 0
    }

    pub fn is_duplicate(&self) -> bool{
        (self.flags() & 1024u16) > 0
    }

    pub fn is_supplementary(&self) -> bool{
        (self.flags() & 2048u16) > 0
    }
}


impl Drop for Record {
    fn drop(&mut self) {
        unsafe { libc::funcs::c95::stdlib::free(self.b.data as *mut libc::types::common::c95::c_void) };
    }
}


pub struct Bamfile {
    f: *mut htslib::Struct_BGZF,
    header: *mut htslib::bam_hdr_t,
}


impl Bamfile{
     pub fn new(filename: &[u8]) -> Self {
        let f = unsafe { htslib::bgzf_open(filename.as_ptr() as *const i8, b"r\0".as_ptr() as *const i8) };
        let header = unsafe { htslib::bam_hdr_read(f) };
        Bamfile { f : f, header : header }
    }

    pub fn read(&self, record: &mut Record) -> Result<(), Error> {
        match unsafe { htslib::bam_read1(self.f, &mut record.b as *mut htslib::bam1_t) } {
            -1 => Err(Error::EOF),
            -2 => Err(Error::Truncated),
            -4 => Err(Error::Invalid),
            _  => Ok(())
        }
    }

    pub fn records(self) -> Records {
        Records { bam: self }
    }
}


impl Drop for Bamfile {
    fn drop(&mut self) {
        unsafe {
            htslib::bam_hdr_destroy(self.header);
            htslib::bgzf_close(self.f);
        }
    }
}


pub struct Records {
    bam: Bamfile
}


impl Iterator for Records {
    type Item = Result<Record, Error>;

    fn next(&mut self) -> Option<Result<Record, Error>> {
        let mut record = Record::new();
        match self.bam.read(&mut record) {
            Err(Error::EOF) => None,
            Ok(())   => Some(Ok(record)),
            Err(err) => Some(Err(err))
        }
    }
}


static DECODE_BASE: &'static [u8] = b"=ACMGRSVTWYHKDBN";
static ENCODE_BASE: [u8; 256] = [
15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0,15,15,
15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
15,15, 5, 6, 8,15, 7, 9, 15,10,15,15, 15,15,15,15,
15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
15,15, 5, 6, 8,15, 7, 9, 15,10,15,15, 15,15,15,15,
15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
];


pub struct Seq<'a> {
    pub encoded: &'a [u8]
}


impl<'a> Seq<'a> {
    #[inline]
    pub fn encoded_base(&self, i: usize) -> u8 {
        (self.encoded[i / 2] >> ((! i & 1) << 2)) & 0b1111
    }

    #[inline]
    pub fn base(&self, i: usize) -> u8 {
        DECODE_BASE[self.encoded_base(i) as usize]
    }

    pub fn as_bytes(&self) -> Vec<u8> {
        (0..self.len()).map(|i| self.base(i)).collect()
    }

    pub fn len(&self) -> usize {
        self.encoded.len() * 2
    }
}


#[derive(Debug)]
#[derive(PartialEq)]
pub enum Cigar {
    Match(u32),  // M
    Ins(u32),  // I
    Del(u32),  // D
    RefSkip(u32),  // N
    SoftClip(u32),  // S
    HardClip(u32),  // H
    Pad(u32),  // P
    Equal(u32),  // =
    Diff(u32),  // X
    Back(u32)  // B
}


impl Cigar {
    fn encode(&self) -> u32 {
        match *self {
            Cigar::Match(len)    => len << 4 | 0,
            Cigar::Ins(len)      => len << 4 | 1,
            Cigar::Del(len)      => len << 4 | 2,
            Cigar::RefSkip(len)  => len << 4 | 3,
            Cigar::SoftClip(len) => len << 4 | 4,
            Cigar::HardClip(len) => len << 4 | 5,
            Cigar::Pad(len)      => len << 4 | 6,
            Cigar::Equal(len)    => len << 4 | 7,
            Cigar::Diff(len)     => len << 4 | 8,
            Cigar::Back(len)     => len << 4 | 9,
        }
    }
}


pub enum Error {
    Truncated,
    Invalid,
    EOF,
}


#[cfg(test)]
mod tests {
    use super::*;
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
        let cigars = [
            [Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)],
            [Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)],
            [Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)],
            [Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)],
            [Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)],
            [Cigar::Match(27), Cigar::Del(100000), Cigar::Match(73)],
        ];

        let bam = Bamfile::new(b"test.bam");

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

        let mut rec = Record::new();
        rec.set(qname, &cigar, seq, qual);
        assert_eq!(rec.qname(), qname);
        assert_eq!(rec.cigar(), cigar);
        assert_eq!(rec.seq().as_bytes(), seq);
        //assert_eq!(rec.qual(), qual);
    }
}
