// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


use std::slice;
use std::ffi;
use std::ops;

use itertools::Itertools;

use htslib;
use utils;


/// A macro creating methods for flag access.
macro_rules! flag {
    ($get:ident, $set:ident, $bit:expr) => (
        pub fn $get(&self) -> bool {
            self.inner().core.flag & $bit != 0
        }

        pub fn $set(&mut self) {
            self.inner_mut().core.flag |= $bit;
        }
    )
}


/// A BAM record.
pub struct Record {
    pub inner: *mut htslib::bam1_t,
    own: bool,
}


unsafe impl Send for Record {}


impl Record {
    /// Create an empty BAM record.
    pub fn new() -> Self {
        let inner = unsafe { htslib::bam_init1() };
        let mut record = Record { inner: inner, own: true };
        record.inner_mut().m_data = 0;
        record
    }

    pub fn from_inner(inner: *mut htslib::bam1_t) -> Self {
        Record { inner: inner, own: false }
    }

    fn data(&self) -> &[u8] {
        unsafe { slice::from_raw_parts(self.inner().data, self.inner().l_data as usize) }
    }

    #[inline]
    pub fn inner_mut(&mut self) -> &mut htslib::bam1_t {
        unsafe { &mut *self.inner }
    }

    #[inline]
    pub fn inner(&self) -> &htslib::bam1_t {
        unsafe { &*self.inner }
    }

    /// Get target id.
    pub fn tid(&self) -> i32 {
        self.inner().core.tid
    }

    /// Set target id.
    pub fn set_tid(&mut self, tid: i32) {
        self.inner_mut().core.tid = tid;
    }

    /// Get position (0-based).
    pub fn pos(&self) -> i32 {
        self.inner().core.pos
    }

    /// Set position (0-based).
    pub fn set_pos(&mut self, pos: i32) {
        self.inner_mut().core.pos = pos;
    }

    pub fn bin(&self) -> u16 {
        self.inner().core.bin
    }

    pub fn set_bin(&mut self, bin: u16) {
        self.inner_mut().core.bin = bin;
    }

    /// Get MAPQ.
    pub fn mapq(&self) -> u8 {
        self.inner().core.qual
    }

    /// Set MAPQ.
    pub fn set_mapq(&mut self, mapq: u8) {
        self.inner_mut().core.qual = mapq;
    }

    /// Get raw flags.
    pub fn flags(&self) -> u16 {
        self.inner().core.flag
    }

    /// Set raw flags.
    pub fn set_flags(&mut self, flags: u16) {
        self.inner_mut().core.flag = flags;
    }

    /// Unset all flags.
    pub fn unset_flags(&mut self) {
        self.inner_mut().core.flag = 0;
    }

    /// Get target id of mate.
    pub fn mtid(&self) -> i32 {
        self.inner().core.mtid
    }

    /// Set target id of mate.
    pub fn set_mtid(&mut self, mtid: i32) {
        self.inner_mut().core.mtid = mtid;
    }

    /// Get mate position.
    pub fn mpos(&self) -> i32 {
        self.inner().core.mpos
    }

    /// Set mate position.
    pub fn set_mpos(&mut self, mpos: i32) {
        self.inner_mut().core.mpos = mpos;
    }

    /// Get insert size.
    pub fn insert_size(&self) -> i32 {
        self.inner().core.isize
    }

    /// Set insert size.
    pub fn set_insert_size(&mut self, insert_size: i32) {
        self.inner_mut().core.isize = insert_size;
    }

    fn qname_len(&self) -> usize {
        self.inner().core.l_qname as usize
    }

    /// Get qname (read name).
    pub fn qname(&self) -> &[u8] {
        &self.data()[..self.qname_len()-1] // -1 ignores the termination symbol
    }

    /// Set variable length data (qname, cigar, seq, qual).
    pub fn set(&mut self, qname: &[u8], cigar: &[Cigar], seq: &[u8], qual: &[u8]) {
        self.inner_mut().l_data = (qname.len() + 1 + cigar.len() * 4 + seq.len() / 2 + qual.len()) as i32;

        if self.inner().m_data < self.inner().l_data {

            self.inner_mut().m_data = self.inner().l_data;
            self.inner_mut().m_data += 32 - self.inner().m_data % 32;
            unsafe {
                self.inner_mut().data = ::libc::realloc(
                    self.inner().data as *mut ::libc::c_void, self.inner().m_data as usize
                ) as *mut u8;
            }
        }

        let mut data = unsafe { slice::from_raw_parts_mut((*self.inner).data, self.inner().l_data as usize) };
        // qname
        utils::copy_memory(qname, data);
        data[qname.len()] = b'\0';
        let mut i = qname.len() + 1;
        self.inner_mut().core.l_qname = i as u8;

        // cigar
        {
            let mut cigar_data = unsafe {
                 slice::from_raw_parts_mut(data[i..].as_ptr() as *mut u32, cigar.len())
            };
            for (i, c) in cigar.iter().enumerate() {
                cigar_data[i] = c.encode();
            }
            self.inner_mut().core.n_cigar = cigar.len() as u16;
            i += cigar.len() * 4;
        }

        // seq
        {
            for j in (0..seq.len()).step(2) {
                data[i + j / 2] = ENCODE_BASE[seq[j] as usize] << 4 | ENCODE_BASE[seq[j + 1] as usize];
            }
            self.inner_mut().core.l_qseq = seq.len() as i32;
            i += (seq.len() + 1) / 2;
        }

        // qual
        utils::copy_memory(qual, &mut data[i..]);
    }

    fn cigar_len(&self) -> usize {
        self.inner().core.n_cigar as usize
    }

    /// Get cigar sequence.
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
        self.inner().core.l_qseq as usize
    }

    /// Get read sequence.
    pub fn seq(&self) -> Seq {
        Seq {
            encoded: &self.data()
                        [self.qname_len() + self.cigar_len()*4..]
                        [..(self.seq_len() + 1) / 2],
            len: self.seq_len()
        }
    }

    /// Get base qualities.
    pub fn qual(&self) -> &[u8] {
        &self.data()[self.qname_len() + self.cigar_len()*4 + (self.seq_len()+1)/2..][..self.seq_len()]
    }

    /// Get auxiliary data (tags).
    pub fn aux(&self, tag: &[u8]) -> Option<Aux> {
        let aux = unsafe { htslib::bam_aux_get(self.inner, ffi::CString::new(tag).unwrap().as_ptr() as *mut i8 ) };

        unsafe {
            if aux.is_null() {
                return None;
            }
            match *aux {
                b'c'|b'C'|b's'|b'S'|b'i'|b'I' => Some(Aux::Integer(htslib::bam_aux2i(aux))),
                b'f'|b'd' => Some(Aux::Float(htslib::bam_aux2f(aux))),
                b'A' => Some(Aux::Char(htslib::bam_aux2A(aux) as u8)),
                b'Z'|b'H' => {
                    let f = aux.offset(1) as *const i8;
                    let x = ffi::CStr::from_ptr(f).to_bytes();
                    Some(Aux::String(x))
                },
                _ => None,
            }
        }
    }

    /// Add auxiliary data.
    pub fn push_aux(&mut self, tag: &[u8], value: &Aux) {
        let ctag = tag.as_ptr() as *mut i8;
        unsafe {
            match *value {
                Aux::Integer(v) => htslib::bam_aux_append(self.inner, ctag, b'i' as i8, 4, [v].as_mut_ptr() as *mut u8),
                Aux::Float(v) => htslib::bam_aux_append(self.inner, ctag, b'f' as i8, 4, [v].as_mut_ptr() as *mut u8),
                Aux::Char(v) => htslib::bam_aux_append(self.inner, ctag, b'A' as i8, 1, [v].as_mut_ptr() as *mut u8),
                Aux::String(v) => htslib::bam_aux_append(
                    self.inner,
                    ctag,
                    b'Z' as i8,
                    (v.len() + 1) as i32,
                    ffi::CString::new(v).unwrap().as_ptr() as *mut u8
                ),
            }
        }
    }

    flag!(is_paired, set_paired, 1u16);
    flag!(is_proper_pair, set_proper_pair, 2u16);
    flag!(is_unmapped, set_unmapped, 4u16);
    flag!(is_mate_unmapped, set_mate_unmapped, 8u16);
    flag!(is_reverse, set_reverse, 16u16);
    flag!(is_mate_reverse, set_mate_reverse, 32u16);
    flag!(is_first_in_pair, set_first_in_pair, 64u16);
    flag!(is_second_in_pair, set_second_in_pair, 128u16);
    flag!(is_secondary, set_secondary, 256u16);
    flag!(is_quality_check_failed, set_quality_check_failed, 512u16);
    flag!(is_duplicate, set_duplicate, 1024u16);
    flag!(is_supplementary, set_supplementary, 2048u16);
}


impl Drop for Record {
    fn drop(&mut self) {
        if self.own {
            unsafe { htslib::bam_destroy1(self.inner) };
        }
    }
}


/// Auxiliary record data.
#[derive(Debug)]
#[derive(PartialEq)]
pub enum Aux<'a> {
    Integer(i32),
    String(&'a [u8]),
    Float(f64),
    Char(u8),
}


impl<'a> Aux<'a> {
    /// Get string from aux data (panics if not a string).
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
    pub encoded: &'a [u8],
    len: usize
}


impl<'a> Seq<'a> {
    #[inline]
    pub fn encoded_base(&self, i: usize) -> u8 {
        (self.encoded[i / 2] >> ((! i & 1) << 2)) & 0b1111
    }

    pub fn as_bytes(&self) -> Vec<u8> {
        (0..self.len()).map(|i| self[i]).collect()
    }

    pub fn len(&self) -> usize {
        self.len
    }
}


impl<'a> ops::Index<usize> for Seq<'a> {
    type Output = u8;

    fn index(&self, index: usize) -> &u8 {
        &DECODE_BASE[self.encoded_base(index) as usize]
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
