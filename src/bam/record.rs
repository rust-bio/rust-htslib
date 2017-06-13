// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


use std::slice;
use std::ffi;
use std::ops;
use std::fmt;
use std::error::Error;

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


quick_error! {
    #[derive(Debug)]
    pub enum CigarError {
        UnsupportedOperation(msg: String) {
            description("Unsupported CIGAR operation")
            display(x) -> ("{}: {}", x.description(), msg)
        }
        UnexpectedOperation(msg: String) {
            description("CIGAR operation not allowed at this point")
            display(x) -> ("{}: {}", x.description(), msg)
        }
    }
}


/// A BAM record.
pub struct Record {
    pub inner: *mut htslib::bam1_t,
    own: bool
}


unsafe impl Send for Record {}
unsafe impl Sync for Record {}


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

    /// Get qname (read name). Complexity: O(1).
    pub fn qname(&self) -> &[u8] {
        &self.data()[..self.qname_len()-1] // -1 ignores the termination symbol
    }

    /// Set variable length data (qname, cigar, seq, qual).
    pub fn set(&mut self, qname: &[u8], cigar: &CigarString, seq: &[u8], qual: &[u8]) {
        self.inner_mut().l_data = (qname.len() + 1 + cigar.len() * 4 + ((seq.len() as f32 / 2.0).ceil() as usize) + qual.len()) as i32;

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
                data[i + j / 2] = ENCODE_BASE[seq[j] as usize] << 4 | (if j + 1 < seq.len() { ENCODE_BASE[seq[j + 1] as usize] } else { 0 });
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

    fn raw_cigar(&self) -> &[u32] {
        unsafe { slice::from_raw_parts(self.data()[self.qname_len()..].as_ptr() as *const u32, self.cigar_len()) }
    }

    /// Get cigar string. Complexity: O(k) with k being the length of the cigar string.
    pub fn cigar(&self) -> CigarStringView {
        let raw = self.raw_cigar();
        CigarStringView {
            inner: CigarString(raw.iter().map(|&c| {
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
            }).collect()),
            pos: self.pos()
        }
    }

    fn seq_len(&self) -> usize {
        self.inner().core.l_qseq as usize
    }

    /// Get read sequence. Complexity: O(1).
    pub fn seq(&self) -> Seq {
        Seq {
            encoded: &self.data()
                        [self.qname_len() + self.cigar_len()*4..]
                        [..(self.seq_len() + 1) / 2],
            len: self.seq_len()
        }
    }

    /// Get base qualities (PHRED-scaled probability that base is wrong).
    /// This does not entail any offsets, hence the qualities can be used directly without
    /// e.g. subtracting 33. Complexity: O(1).
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

    // Delete auxiliary tag.
    pub fn remove_aux(&self, tag: &[u8]) -> bool {
        let aux = unsafe { htslib::bam_aux_get(self.inner, ffi::CString::new(tag).unwrap().as_ptr() as *mut i8 ) };
        unsafe {
            if aux.is_null() {
                false
            } else {
                htslib::bam_aux_del(self.inner, aux);
                true
            }
        }
    }

    flag!(is_paired, set_paired, 1u16);
    flag!(is_proper_pair, set_proper_pair, 2u16);
    flag!(is_unmapped, set_unmapped, 4u16);
    flag!(is_mate_unmapped, set_mate_unmapped, 8u16);
    flag!(is_reverse, set_reverse, 16u16);
    flag!(is_mate_reverse, set_mate_reverse, 32u16);
    flag!(is_first_in_template, set_first_in_template, 64u16);
    flag!(is_last_in_template, set_last_in_template, 128u16);
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


unsafe impl<'a> Send for Aux<'a> {}
unsafe impl<'a> Sync for Aux<'a> {}


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


/// The sequence of a record.
pub struct Seq<'a> {
    pub encoded: &'a [u8],
    len: usize
}


impl<'a> Seq<'a> {
    /// Return encoded base. Complexity: O(1).
    #[inline]
    pub fn encoded_base(&self, i: usize) -> u8 {
        (self.encoded[i / 2] >> ((! i & 1) << 2)) & 0b1111
    }

    /// Return decoded sequence. Complexity: O(m) with m being the read length.
    pub fn as_bytes(&self) -> Vec<u8> {
        (0..self.len()).map(|i| self[i]).collect()
    }

    /// Return length (in bases) of the sequence.
    pub fn len(&self) -> usize {
        self.len
    }
}


impl<'a> ops::Index<usize> for Seq<'a> {
    type Output = u8;

    /// Return decoded base at given position within read. Complexity: O(1).
    fn index(&self, index: usize) -> &u8 {
        &DECODE_BASE[self.encoded_base(index) as usize]
    }
}


unsafe impl<'a> Send for Seq<'a> {}
unsafe impl<'a> Sync for Seq<'a> {}


#[derive(PartialEq, Eq, Debug, Clone)]
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

    /// Return the length the CIGAR.
    pub fn len(&self) -> u32 {
        match *self {
            Cigar::Match(len)    => len,
            Cigar::Ins(len)      => len,
            Cigar::Del(len)      => len,
            Cigar::RefSkip(len)  => len,
            Cigar::SoftClip(len) => len,
            Cigar::HardClip(len) => len,
            Cigar::Pad(len)      => len,
            Cigar::Equal(len)    => len,
            Cigar::Diff(len)     => len,
            Cigar::Back(len)     => len
        }
    }

    /// Return the character representing the CIGAR.
    pub fn char(&self) -> char {
        match *self {
            Cigar::Match(_)    => 'M',
            Cigar::Ins(_)      => 'I',
            Cigar::Del(_)      => 'D',
            Cigar::RefSkip(_)  => 'N',
            Cigar::SoftClip(_) => 'S',
            Cigar::HardClip(_) => 'H',
            Cigar::Pad(_)      => 'P',
            Cigar::Equal(_)    => '=',
            Cigar::Diff(_)     => 'X',
            Cigar::Back(_)     => 'B'
        }
    }
}


impl fmt::Display for Cigar {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        fmt.write_fmt(format_args!("{}{}", self.len(), self.char()))
    }
}


unsafe impl Send for Cigar {}
unsafe impl Sync for Cigar {}


custom_derive! {
    /// A CIGAR string. This type wraps around a `Vec<Cigar>`.
    ///
    /// # Example
    ///
    /// ```
    /// use rust_htslib::bam::record::{Cigar, CigarString};
    ///
    /// let cigar = CigarString(vec![Cigar::Match(100), Cigar::SoftClip(10)]);
    ///
    /// // access by index
    /// assert_eq!(cigar[0], Cigar::Match(100));
    /// // format into classical string representation
    /// assert_eq!(format!("{}", cigar), "100M10S");
    /// // iterate
    /// for op in &cigar {
    ///    println!("{}", op);
    /// }
    /// ```
    #[derive(NewtypeDeref,
             NewtypeIndex(usize),
             NewtypeIndexMut(usize),
             NewtypeFrom,
             PartialEq,
             Eq,
             NewtypeDebug,
             Clone
    )]
    pub struct CigarString(pub Vec<Cigar>);
}


impl<'a> CigarString {
    pub fn iter(&'a self) -> ::std::slice::Iter<'a, Cigar> {
        self.into_iter()
    }
}


impl<'a> IntoIterator for &'a CigarString {
    type Item = &'a Cigar;
    type IntoIter = ::std::slice::Iter<'a, Cigar>;

    fn into_iter(self) -> Self::IntoIter {
        (&(self.0)).into_iter()
    }
}


impl fmt::Display for CigarString {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        for op in self {
            fmt.write_fmt(format_args!("{}{}", op.len(), op.char()))?;
        }
        Ok(())
    }
}


#[derive(Eq, PartialEq, Clone, Debug)]
pub struct CigarStringView {
    inner: CigarString,
    pos: i32
}


impl CigarStringView {
    /// Get end position of alignment.
    pub fn end_pos(&self) -> i32 {
        let mut pos = self.pos;
        for c in self {
            match c {
                &Cigar::Match(l) | &Cigar::RefSkip(l) | &Cigar::Del(l) |
                &Cigar::Equal(l) | &Cigar::Diff(l) => pos += l as i32,
                &Cigar::Back(l) => pos -= l as i32,
                // these don't add to end_pos on reference
                &Cigar::Ins(_) | &Cigar::SoftClip(_) | &Cigar::HardClip(_) | &Cigar::Pad(_) => ()
            }
        }
        pos
    }

    /// For a given position in the reference, get corresponding position within read.
    /// If reference position is outside of the read alignment, return None.
    ///
    /// # Arguments
    ///
    /// * `ref_pos` - the reference position
    /// * `include_softclips` - if true, softclips will be considered as matches or mismatches
    /// * `include_dels` - if true, positions within deletions will be considered (start of deletion will be returned)
    ///
    pub fn read_pos(
        &self,
        ref_pos: u32,
        include_softclips: bool,
        include_dels: bool
    ) -> Result<Option<u32>, CigarError> {
        let mut rpos = self.pos as u32; // reference position
        let mut qpos = 0u32; // position within read
        let mut j = 0; // index into cigar operation vector

        // find first cigar operation referring to qpos = 0 (and thus bases in record.seq()),
        // because all augmentations of qpos and rpos before that are invalid
        for (i, c) in self.iter().enumerate() {
            match c {
                &Cigar::Match(_) |
                &Cigar::Diff(_)  |
                &Cigar::Equal(_) |
                // this is unexpected, but bwa + GATK indel realignment can produce insertions
                // before matching positions
                &Cigar::Ins(_)   |
                &Cigar::Pad(_)   => { //TODO: double check if we do not want Pad to be an UnsupportedOperation, e.g. with the following:
//                    return Err( CigarError::UnsupportedOperation("Pad operation definition in SAMv1 spec is unclear.".to_string()) );
                    j = i;
                    break;
                },
                &Cigar::SoftClip(l) => {
                    j = i;
                    if include_softclips {
                        // Alignment starts with softclip and we want to include it in the
                        // projection of the reference position. However, the POS field does not
                        // include the softclip. Hence we have to subtract its length.
                        rpos = rpos.saturating_sub(l);
                    }
                    break;
                },
                &Cigar::Del(_) => {
                    return Err(CigarError::UnexpectedOperation(
                        "'deletion' (D) found before any operation describing read sequence".to_owned()
                    ));
                },
                &Cigar::Back(_) => {
                    return Err(CigarError::UnsupportedOperation(
                        "'back' (B) operation is deprecated according to htslib/bam_plcmd.c and is not in SAMv1 spec".to_owned()
                    ));
                },
                &Cigar::RefSkip(_) => {
                    return Err(CigarError::UnexpectedOperation(
                        "'reference skip' (N) found before any operation describing read sequence".to_owned()
                    ));
                },
                &Cigar::HardClip(_) => () // just skip HardClips, they're not in the read sequence any more
            }
        }

        let contains_ref_pos = |cigar_op_start: u32, cigar_op_length: u32| {
            cigar_op_start <= ref_pos && cigar_op_start + cigar_op_length > ref_pos
        };

        while rpos <= ref_pos && j < self.len() {
            match &self[j] {
                // potential SNV evidence
                &Cigar::Match(l) | &Cigar::Diff(l) | &Cigar::Equal(l)
                if contains_ref_pos(rpos, l) => {
                    // difference between desired position and first position of current cigar
                    // operation
                    qpos += ref_pos - rpos;
                    return Ok(Some(qpos));
                },
                &Cigar::SoftClip(l) if include_softclips && contains_ref_pos(rpos, l) => {
                    qpos += ref_pos - rpos;
                    return Ok(Some(qpos));
                },
                &Cigar::Del(l) if include_dels && contains_ref_pos(rpos, l) => {
                    // qpos shall resemble the start of the deletion
                    return Ok(Some(qpos));
                },
                // for others, just increase pos and qpos as needed
                &Cigar::Match(l)   |
                &Cigar::Diff(l)    |
                &Cigar::Equal(l)   => {
                    rpos += l;
                    qpos += l;
                    j += 1;
                },
                &Cigar::SoftClip(l) => {
                    qpos += l;
                    j += 1;
                    if include_softclips {
                        rpos += l;
                    }
                },
                &Cigar::Ins(l)  => {
                    qpos += l;
                    j += 1;
                },
                &Cigar::RefSkip(l) |
                &Cigar::Del(l) |
                &Cigar::Pad(l) => { //TODO: double check if we do not want Pad to be an UnsupportedOperation, e.g. with the following:
//                    return Err( CigarError::UnsupportedOperation("Pad operation definition in SAMv1 spec is unclear.".to_string()) );
                    rpos += l;
                    j += 1;
                },
                &Cigar::HardClip(_) => return Ok(None),
                &Cigar::Back(_) => {
                    return Err(CigarError::UnsupportedOperation(
                        "'back' (B) operation is deprecated according to htslib/bam_plcmd.c and is not in SAMv1 spec".to_owned()
                    ));
                }
            }
        }

        Ok(None)
    }
}


impl ops::Deref for CigarStringView {
    type Target = CigarString;

    fn deref(&self) -> &CigarString {
        &self.inner
    }
}


impl ops::Index<usize> for CigarStringView {
    type Output = Cigar;

    fn index(&self, index: usize) -> &Cigar {
        self.inner.index(index)
    }
}


impl ops::IndexMut<usize> for CigarStringView {

    fn index_mut(&mut self, index: usize) -> &mut Cigar {
        self.inner.index_mut(index)
    }
}


impl<'a> CigarStringView {
    pub fn iter(&'a self) -> ::std::slice::Iter<'a, Cigar> {
        self.inner.into_iter()
    }
}


impl<'a> IntoIterator for &'a CigarStringView {
    type Item = &'a Cigar;
    type IntoIter = ::std::slice::Iter<'a, Cigar>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.into_iter()
    }
}


impl fmt::Display for CigarStringView {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        self.inner.fmt(fmt)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_string() {
        let cigar = CigarString(vec![Cigar::Match(100), Cigar::SoftClip(10)]);

        assert_eq!(cigar[0], Cigar::Match(100));
        assert_eq!(format!("{}", cigar), "100M10S");
        for op in &cigar {
            println!("{}", op);
        }
    }
}
