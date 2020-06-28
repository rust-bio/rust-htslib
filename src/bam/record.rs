// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::convert::TryFrom;
use std::ffi;
use std::fmt;
use std::mem::{size_of, MaybeUninit};
use std::ops;
use std::rc::Rc;
use std::slice;
use std::str;
use std::str::FromStr;
use std::u32;

use lazy_static::lazy_static;
use regex::Regex;

use crate::bam::errors::Result;
use crate::bam::Error;
use crate::bam::HeaderView;
use crate::htslib;
use crate::utils;

use bio_types::alignment::{Alignment, AlignmentMode, AlignmentOperation};
use bio_types::genome;
use bio_types::sequence::SequenceRead;
use bio_types::strand::ReqStrand;

/// A macro creating methods for flag access.
macro_rules! flag {
    ($get:ident, $set:ident, $unset:ident, $bit:expr) => {
        pub fn $get(&self) -> bool {
            self.inner().core.flag & $bit != 0
        }

        pub fn $set(&mut self) {
            self.inner_mut().core.flag |= $bit;
        }

        pub fn $unset(&mut self) {
            self.inner_mut().core.flag &= !$bit;
        }
    };
}

/// A BAM record.
pub struct Record {
    pub inner: htslib::bam1_t,
    own: bool,
    cigar: Option<CigarStringView>,
    header: Option<Rc<HeaderView>>,
}

unsafe impl Send for Record {}
unsafe impl Sync for Record {}

impl Clone for Record {
    fn clone(&self) -> Self {
        let mut copy = Record::new();
        unsafe { htslib::bam_copy1(copy.inner_ptr_mut(), self.inner_ptr()) };
        copy
    }
}

impl PartialEq for Record {
    fn eq(&self, other: &Record) -> bool {
        self.tid() == other.tid()
            && self.pos() == other.pos()
            && self.bin() == other.bin()
            && self.mapq() == other.mapq()
            && self.flags() == other.flags()
            && self.mtid() == other.mtid()
            && self.mpos() == other.mpos()
            && self.insert_size() == other.insert_size()
            && self.data() == other.data()
            && self.inner().core.l_extranul == other.inner().core.l_extranul
    }
}

impl Eq for Record {}

impl fmt::Debug for Record {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        fmt.write_fmt(format_args!(
            "Record(tid: {}, pos: {})",
            self.tid(),
            self.pos()
        ))
    }
}

impl Default for Record {
    fn default() -> Self {
        Self::new()
    }
}

#[inline]
fn extranul_from_qname(qname: &[u8]) -> usize {
    let qlen = qname.len() + 1;
    if qlen % 4 != 0 {
        4 - qlen % 4
    } else {
        0
    }
}

impl Record {
    /// Create an empty BAM record.
    pub fn new() -> Self {
        Record {
            inner: unsafe { MaybeUninit::zeroed().assume_init() },
            own: true,
            cigar: None,
            header: None,
        }
    }

    pub fn from_inner(from: *mut htslib::bam1_t) -> Self {
        Record {
            inner: {
                #[allow(clippy::uninit_assumed_init)]
                let mut inner = unsafe { MaybeUninit::uninit().assume_init() };
                unsafe {
                    ::libc::memcpy(
                        &mut inner as *mut htslib::bam1_t as *mut ::libc::c_void,
                        from as *const ::libc::c_void,
                        size_of::<htslib::bam1_t>(),
                    );
                }
                inner
            },
            own: false,
            cigar: None,
            header: None,
        }
    }

    // Create a BAM record from a line SAM text. SAM slice need not be 0-terminated.
    pub fn from_sam(header_view: &HeaderView, sam: &[u8]) -> Result<Record> {
        let mut record = Self::new();

        let mut sam_copy = Vec::with_capacity(sam.len() + 1);
        sam_copy.extend(sam);
        sam_copy.push(0);

        let mut sam_string = htslib::kstring_t {
            s: sam_copy.as_ptr() as *mut i8,
            l: sam_copy.len() as u64,
            m: sam_copy.len() as u64,
        };

        let succ = unsafe {
            htslib::sam_parse1(
                &mut sam_string,
                header_view.inner_ptr() as *mut htslib::bam_hdr_t,
                record.inner_ptr_mut(),
            )
        };

        if succ == 0 {
            Ok(record)
        } else {
            Err(Error::ParseSAM {
                rec: str::from_utf8(&sam_copy).unwrap().to_owned(),
            })
        }
    }

    pub fn set_header(&mut self, header: Rc<HeaderView>) {
        self.header = Some(header);
    }

    pub(super) fn data(&self) -> &[u8] {
        unsafe { slice::from_raw_parts(self.inner().data, self.inner().l_data as usize) }
    }

    #[inline]
    pub fn inner_mut(&mut self) -> &mut htslib::bam1_t {
        &mut self.inner
    }

    #[inline]
    pub(super) fn inner_ptr_mut(&mut self) -> *mut htslib::bam1_t {
        &mut self.inner as *mut htslib::bam1_t
    }

    #[inline]
    pub fn inner(&self) -> &htslib::bam1_t {
        &self.inner
    }

    #[inline]
    pub(super) fn inner_ptr(&self) -> *const htslib::bam1_t {
        &self.inner as *const htslib::bam1_t
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
    pub fn pos(&self) -> i64 {
        self.inner().core.pos
    }

    /// Set position (0-based).
    pub fn set_pos(&mut self, pos: i64) {
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

    /// Get strand information from record flags.
    pub fn strand(&mut self) -> ReqStrand {
        let reverse = self.flags() & 0x10 != 0;
        if reverse {
            ReqStrand::Reverse
        } else {
            ReqStrand::Forward
        }
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
    pub fn mpos(&self) -> i64 {
        self.inner().core.mpos
    }

    /// Set mate position.
    pub fn set_mpos(&mut self, mpos: i64) {
        self.inner_mut().core.mpos = mpos;
    }

    /// Get insert size.
    pub fn insert_size(&self) -> i64 {
        self.inner().core.isize
    }

    /// Set insert size.
    pub fn set_insert_size(&mut self, insert_size: i64) {
        self.inner_mut().core.isize = insert_size;
    }

    fn qname_capacity(&self) -> usize {
        self.inner().core.l_qname as usize
    }

    fn qname_len(&self) -> usize {
        // discount all trailing zeros (the default one and extra nulls)
        self.qname_capacity() - 1 - self.inner().core.l_extranul as usize
    }

    /// Get qname (read name). Complexity: O(1).
    pub fn qname(&self) -> &[u8] {
        &self.data()[..self.qname_len()]
    }

    /// Set the variable length data buffer
    pub fn set_data(&mut self, new_data: &[u8]) {
        self.cigar = None;

        self.inner_mut().l_data = new_data.len() as i32;
        if (self.inner().m_data as i32) < self.inner().l_data {
            // Verbosity due to lexical borrowing
            let l_data = self.inner().l_data;
            self.realloc_var_data(l_data as usize);
        }

        // Copy new data into buffer
        let data =
            unsafe { slice::from_raw_parts_mut(self.inner.data, self.inner().l_data as usize) };
        utils::copy_memory(new_data, data);
    }

    /// Set variable length data (qname, cigar, seq, qual).
    /// The aux data is left unchanged.
    /// `qual` is Phred-scaled quality values, without any offset.
    /// NOTE: seq.len() must equal qual.len() or this method
    /// will panic. If you don't have quality values use
    /// `let quals = vec![ 255 as u8; seq.len()];` as a placeholder that will
    /// be recognized as missing QVs by `samtools`.
    pub fn set(&mut self, qname: &[u8], cigar: Option<&CigarString>, seq: &[u8], qual: &[u8]) {
        assert!(qname.len() < 255);
        assert_eq!(seq.len(), qual.len(), "seq.len() must equal qual.len()");

        self.cigar = None;

        let cigar_width = if let Some(cigar_string) = cigar {
            cigar_string.len()
        } else {
            0
        } * 4;
        let q_len = qname.len() + 1;
        let extranul = extranul_from_qname(qname);

        let orig_aux_offset = self.qname_capacity()
            + 4 * self.cigar_len()
            + (self.seq_len() + 1) / 2
            + self.seq_len();
        let new_aux_offset = q_len + extranul + cigar_width + (seq.len() + 1) / 2 + qual.len();
        assert!(orig_aux_offset <= self.inner.l_data as usize);
        let aux_len = self.inner.l_data as usize - orig_aux_offset;
        self.inner_mut().l_data = (new_aux_offset + aux_len) as i32;
        if (self.inner().m_data as i32) < self.inner().l_data {
            // Verbosity due to lexical borrowing
            let l_data = self.inner().l_data;
            self.realloc_var_data(l_data as usize);
        }

        // Copy the aux data.
        if aux_len > 0 && orig_aux_offset != new_aux_offset {
            let data =
                unsafe { slice::from_raw_parts_mut(self.inner.data, self.inner().m_data as usize) };
            data.copy_within(orig_aux_offset..orig_aux_offset + aux_len, new_aux_offset);
        }

        let data =
            unsafe { slice::from_raw_parts_mut(self.inner.data, self.inner().l_data as usize) };

        // qname
        utils::copy_memory(qname, data);
        for i in 0..=extranul {
            data[qname.len() + i] = b'\0';
        }
        let mut i = q_len + extranul;
        self.inner_mut().core.l_qname = i as u16;
        self.inner_mut().core.l_extranul = extranul as u8;

        // cigar
        if let Some(cigar_string) = cigar {
            let cigar_data = unsafe {
                //cigar is always aligned to 4 bytes (see extranul above) - so this is safe
                #[allow(clippy::cast_ptr_alignment)]
                slice::from_raw_parts_mut(data[i..].as_ptr() as *mut u32, cigar_string.len())
            };
            for (i, c) in cigar_string.iter().enumerate() {
                cigar_data[i] = c.encode();
            }
            self.inner_mut().core.n_cigar = cigar_string.len() as u32;
            i += cigar_string.len() * 4;
        } else {
            self.inner_mut().core.n_cigar = 0;
        };

        // seq
        {
            for j in (0..seq.len()).step_by(2) {
                data[i + j / 2] = ENCODE_BASE[seq[j] as usize] << 4
                    | (if j + 1 < seq.len() {
                        ENCODE_BASE[seq[j + 1] as usize]
                    } else {
                        0
                    });
            }
            self.inner_mut().core.l_qseq = seq.len() as i32;
            i += (seq.len() + 1) / 2;
        }

        // qual
        utils::copy_memory(qual, &mut data[i..]);
    }

    /// Replace current qname with a new one.
    pub fn set_qname(&mut self, new_qname: &[u8]) {
        // 251 + 1NUL is the max 32-bit aligned value that fits in u8
        assert!(new_qname.len() < 252);

        let old_q_len = self.qname_capacity();
        // We're going to add a terminal NUL
        let extranul = extranul_from_qname(new_qname);
        let new_q_len = new_qname.len() + 1 + extranul;

        // Length of data after qname
        let other_len = self.inner_mut().l_data - old_q_len as i32;

        if new_q_len < old_q_len && self.inner().l_data > (old_q_len as i32) {
            self.inner_mut().l_data -= (old_q_len - new_q_len) as i32;
        } else if new_q_len > old_q_len {
            self.inner_mut().l_data += (new_q_len - old_q_len) as i32;

            // Reallocate if necessary
            if (self.inner().m_data as i32) < self.inner().l_data {
                // Verbosity due to lexical borrowing
                let l_data = self.inner().l_data;
                self.realloc_var_data(l_data as usize);
            }
        }

        if new_q_len != old_q_len {
            // Move other data to new location
            unsafe {
                let data = slice::from_raw_parts_mut(self.inner.data, self.inner().l_data as usize);

                ::libc::memmove(
                    data.as_mut_ptr().add(new_q_len) as *mut ::libc::c_void,
                    data.as_mut_ptr().add(old_q_len) as *mut ::libc::c_void,
                    other_len as usize,
                );
            }
        }

        // Copy qname data
        let data =
            unsafe { slice::from_raw_parts_mut(self.inner.data, self.inner().l_data as usize) };
        utils::copy_memory(new_qname, data);
        for i in 0..=extranul {
            data[new_q_len - i - 1] = b'\0';
        }
        self.inner_mut().core.l_qname = new_q_len as u16;
        self.inner_mut().core.l_extranul = extranul as u8;
    }

    fn realloc_var_data(&mut self, new_len: usize) {
        // pad request
        let new_len = new_len as u32;
        let new_request = new_len + 32 - (new_len % 32);

        let ptr = unsafe {
            ::libc::realloc(
                self.inner().data as *mut ::libc::c_void,
                new_request as usize,
            ) as *mut u8
        };

        if ptr.is_null() {
            panic!("ran out of memory in rust_htslib trying to realloc");
        }

        // don't update m_data until we know we have
        // a successful allocation.
        self.inner_mut().m_data = new_request;
        self.inner_mut().data = ptr;

        // we now own inner.data
        self.own = true;
    }

    pub fn cigar_len(&self) -> usize {
        self.inner().core.n_cigar as usize
    }

    /// Get reference to raw cigar string representation (as stored in BAM file).
    /// Usually, the method `Record::cigar` should be used instead.
    pub fn raw_cigar(&self) -> &[u32] {
        //cigar is always aligned to 4 bytes - so this is safe
        #[allow(clippy::cast_ptr_alignment)]
        unsafe {
            slice::from_raw_parts(
                self.data()[self.qname_capacity()..].as_ptr() as *const u32,
                self.cigar_len(),
            )
        }
    }

    /// Return unpacked cigar string. This will create a fresh copy the Cigar data.
    pub fn cigar(&self) -> CigarStringView {
        match self.cigar {
            Some(ref c) => c.clone(),
            None => self.unpack_cigar(),
        }
    }

    // Return unpacked cigar string. This returns None unless you have first called `bam::Record::cache_cigar`.
    pub fn cigar_cached(&self) -> Option<&CigarStringView> {
        self.cigar.as_ref()
    }

    /// Decode the cigar string and cache it inside the `Record`
    pub fn cache_cigar(&mut self) {
        self.cigar = Some(self.unpack_cigar())
    }

    /// Unpack cigar string. Complexity: O(k) with k being the length of the cigar string.
    fn unpack_cigar(&self) -> CigarStringView {
        CigarString(
            self.raw_cigar()
                .iter()
                .map(|&c| {
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
                        _ => panic!("Unexpected cigar operation"),
                    }
                })
                .collect(),
        )
        .into_view(self.pos())
    }

    fn seq_len(&self) -> usize {
        self.inner().core.l_qseq as usize
    }

    fn seq_data(&self) -> &[u8] {
        let offset = self.qname_capacity() + self.cigar_len() * 4;
        &self.data()[offset..][..(self.seq_len() + 1) / 2]
    }

    /// Get read sequence. Complexity: O(1).
    pub fn seq(&self) -> Seq<'_> {
        Seq {
            encoded: self.seq_data(),
            len: self.seq_len(),
        }
    }

    /// Get base qualities (PHRED-scaled probability that base is wrong).
    /// This does not entail any offsets, hence the qualities can be used directly without
    /// e.g. subtracting 33. Complexity: O(1).
    pub fn qual(&self) -> &[u8] {
        &self.data()[self.qname_capacity() + self.cigar_len() * 4 + (self.seq_len() + 1) / 2..]
            [..self.seq_len()]
    }

    /// Get auxiliary data (tags).
    pub fn aux(&self, tag: &[u8]) -> Option<Aux<'_>> {
        let c_str = ffi::CString::new(tag).unwrap();
        let aux = unsafe {
            htslib::bam_aux_get(
                &self.inner as *const htslib::bam1_t,
                c_str.as_ptr() as *mut i8,
            )
        };

        unsafe {
            if aux.is_null() {
                return None;
            }
            match *aux {
                b'c' | b'C' | b's' | b'S' | b'i' | b'I' => {
                    Some(Aux::Integer(htslib::bam_aux2i(aux) as i64))
                }
                b'f' | b'd' => Some(Aux::Float(htslib::bam_aux2f(aux))),
                b'A' => Some(Aux::Char(htslib::bam_aux2A(aux) as u8)),
                b'Z' | b'H' => {
                    let f = aux.offset(1) as *const i8;
                    let x = ffi::CStr::from_ptr(f).to_bytes();
                    Some(Aux::String(x))
                }
                _ => None,
            }
        }
    }

    /// Add auxiliary data.
    pub fn push_aux(&mut self, tag: &[u8], value: &Aux<'_>) {
        let ctag = tag.as_ptr() as *mut i8;
        let ret = unsafe {
            match *value {
                Aux::Integer(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'i' as i8,
                    4,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::Float(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'f' as i8,
                    4,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::Char(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'A' as i8,
                    1,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::String(v) => {
                    let c_str = ffi::CString::new(v).unwrap();
                    htslib::bam_aux_append(
                        self.inner_ptr_mut(),
                        ctag,
                        b'Z' as i8,
                        (v.len() + 1) as i32,
                        c_str.as_ptr() as *mut u8,
                    )
                }
            }
        };

        if ret < 0 {
            panic!("htslib ran out of memory in push_aux");
        }
    }

    // Delete auxiliary tag.
    pub fn remove_aux(&mut self, tag: &[u8]) -> bool {
        let c_str = ffi::CString::new(tag).unwrap();
        let aux = unsafe {
            htslib::bam_aux_get(
                &self.inner as *const htslib::bam1_t,
                c_str.as_ptr() as *mut i8,
            )
        };
        unsafe {
            if aux.is_null() {
                false
            } else {
                htslib::bam_aux_del(self.inner_ptr_mut(), aux);
                true
            }
        }
    }

    flag!(is_paired, set_paired, unset_paired, 1u16);
    flag!(is_proper_pair, set_proper_pair, unset_proper_pair, 2u16);
    flag!(is_unmapped, set_unmapped, unset_unmapped, 4u16);
    flag!(
        is_mate_unmapped,
        set_mate_unmapped,
        unset_mate_unmapped,
        8u16
    );
    flag!(is_reverse, set_reverse, unset_reverse, 16u16);
    flag!(is_mate_reverse, set_mate_reverse, unset_mate_reverse, 32u16);
    flag!(
        is_first_in_template,
        set_first_in_template,
        unset_first_in_template,
        64u16
    );
    flag!(
        is_last_in_template,
        set_last_in_template,
        unset_last_in_template,
        128u16
    );
    flag!(is_secondary, set_secondary, unset_secondary, 256u16);
    flag!(
        is_quality_check_failed,
        set_quality_check_failed,
        unset_quality_check_failed,
        512u16
    );
    flag!(is_duplicate, set_duplicate, unset_duplicate, 1024u16);
    flag!(
        is_supplementary,
        set_supplementary,
        unset_supplementary,
        2048u16
    );
}

impl Drop for Record {
    fn drop(&mut self) {
        if self.own {
            unsafe { ::libc::free(self.inner.data as *mut ::libc::c_void) }
        }
    }
}

impl SequenceRead for Record {
    fn name(&self) -> &[u8] {
        self.qname()
    }

    fn base(&self, i: usize) -> u8 {
        decode_base(encoded_base(self.seq_data(), i))
    }

    fn base_qual(&self, i: usize) -> u8 {
        self.qual()[i]
    }

    fn len(&self) -> usize {
        self.seq_len()
    }
}

impl genome::AbstractInterval for Record {
    /// Return contig name. Panics if record does not know its header (which happens if it has not been read from a file).
    fn contig(&self) -> &str {
        let tid = self.tid();
        if tid < 0 {
            panic!("invalid tid, must be at least zero");
        }
        str::from_utf8(
            self.header
                .as_ref()
                .expect(
                    "header must be set (this is the case if the record has been read from a file)",
                )
                .tid2name(tid as u32),
        )
        .expect("unable to interpret contig name as UTF-8")
    }

    /// Return genomic range covered by alignment. Panics if `Record::cache_cigar()` has not been called first or `Record::pos()` is less than zero.
    fn range(&self) -> ops::Range<genome::Position> {
        let end_pos = self
            .cigar_cached()
            .expect("cigar has not been cached yet, call cache_cigar() first")
            .end_pos() as u64;

        if self.pos() < 0 {
            panic!("invalid position, must be positive")
        }

        self.pos() as u64..end_pos
    }
}

/// Auxiliary record data.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Aux<'a> {
    Integer(i64),
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

    pub fn integer(&self) -> i64 {
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

static DECODE_BASE: &[u8] = b"=ACMGRSVTWYHKDBN";
static ENCODE_BASE: [u8; 256] = [
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    1, 2, 4, 8, 15, 15, 15, 15, 15, 15, 15, 15, 15, 0, 15, 15, 15, 1, 14, 2, 13, 15, 15, 4, 11, 15,
    15, 12, 15, 3, 15, 15, 15, 15, 5, 6, 8, 15, 7, 9, 15, 10, 15, 15, 15, 15, 15, 15, 15, 1, 14, 2,
    13, 15, 15, 4, 11, 15, 15, 12, 15, 3, 15, 15, 15, 15, 5, 6, 8, 15, 7, 9, 15, 10, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
];

fn encoded_base(encoded_seq: &[u8], i: usize) -> u8 {
    (encoded_seq[i / 2] >> ((!i & 1) << 2)) & 0b1111
}

fn decode_base(base: u8) -> u8 {
    DECODE_BASE[base as usize]
}

/// The sequence of a record.
#[derive(Debug, Copy, Clone)]
pub struct Seq<'a> {
    pub encoded: &'a [u8],
    len: usize,
}

impl<'a> Seq<'a> {
    /// Return encoded base. Complexity: O(1).
    #[inline]
    pub fn encoded_base(&self, i: usize) -> u8 {
        encoded_base(self.encoded, i)
    }

    /// Return decoded sequence. Complexity: O(m) with m being the read length.
    pub fn as_bytes(&self) -> Vec<u8> {
        (0..self.len()).map(|i| self[i]).collect()
    }

    /// Return length (in bases) of the sequence.
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
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

#[derive(PartialEq, Eq, Debug, Clone, Copy, Hash)]
pub enum Cigar {
    Match(u32),    // M
    Ins(u32),      // I
    Del(u32),      // D
    RefSkip(u32),  // N
    SoftClip(u32), // S
    HardClip(u32), // H
    Pad(u32),      // P
    Equal(u32),    // =
    Diff(u32),     // X
}

impl Cigar {
    fn encode(self) -> u32 {
        match self {
            Cigar::Match(len) => len << 4, // | 0,
            Cigar::Ins(len) => len << 4 | 1,
            Cigar::Del(len) => len << 4 | 2,
            Cigar::RefSkip(len) => len << 4 | 3,
            Cigar::SoftClip(len) => len << 4 | 4,
            Cigar::HardClip(len) => len << 4 | 5,
            Cigar::Pad(len) => len << 4 | 6,
            Cigar::Equal(len) => len << 4 | 7,
            Cigar::Diff(len) => len << 4 | 8,
        }
    }

    /// Return the length of the CIGAR.
    pub fn len(self) -> u32 {
        match self {
            Cigar::Match(len) => len,
            Cigar::Ins(len) => len,
            Cigar::Del(len) => len,
            Cigar::RefSkip(len) => len,
            Cigar::SoftClip(len) => len,
            Cigar::HardClip(len) => len,
            Cigar::Pad(len) => len,
            Cigar::Equal(len) => len,
            Cigar::Diff(len) => len,
        }
    }

    pub fn is_empty(self) -> bool {
        self.len() == 0
    }

    /// Return the character representing the CIGAR.
    pub fn char(self) -> char {
        match self {
            Cigar::Match(_) => 'M',
            Cigar::Ins(_) => 'I',
            Cigar::Del(_) => 'D',
            Cigar::RefSkip(_) => 'N',
            Cigar::SoftClip(_) => 'S',
            Cigar::HardClip(_) => 'H',
            Cigar::Pad(_) => 'P',
            Cigar::Equal(_) => '=',
            Cigar::Diff(_) => 'X',
        }
    }
}

impl fmt::Display for Cigar {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
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
             Clone,
             Hash
    )]
    pub struct CigarString(pub Vec<Cigar>);
}

impl CigarString {
    /// Create a `CigarStringView` from this CigarString at position `pos`
    pub fn into_view(self, pos: i64) -> CigarStringView {
        CigarStringView::new(self, pos)
    }

    /// Calculate the bam cigar from the alignment struct. x is the target string
    /// and y is the reference. `hard_clip` controls how unaligned read bases are encoded in the
    /// cigar string. Set to true to use the hard clip (`H`) code, or false to use soft clip
    /// (`S`) code. See the [SAM spec](https://samtools.github.io/hts-specs/SAMv1.pdf) for more details.
    pub fn from_alignment(alignment: &Alignment, hard_clip: bool) -> Self {
        match alignment.mode {
            AlignmentMode::Global => {
                panic!(" Bam cigar fn not supported for Global Alignment mode")
            }
            AlignmentMode::Local => panic!(" Bam cigar fn not supported for Local Alignment mode"),
            _ => {}
        }

        let mut cigar = Vec::new();
        if alignment.operations.is_empty() {
            return CigarString(cigar);
        }

        let add_op = |op: AlignmentOperation, length: u32, cigar: &mut Vec<Cigar>| match op {
            AlignmentOperation::Del => cigar.push(Cigar::Del(length)),
            AlignmentOperation::Ins => cigar.push(Cigar::Ins(length)),
            AlignmentOperation::Subst => cigar.push(Cigar::Diff(length)),
            AlignmentOperation::Match => cigar.push(Cigar::Equal(length)),
            _ => {}
        };

        if alignment.xstart > 0 {
            cigar.push(if hard_clip {
                Cigar::HardClip(alignment.xstart as u32)
            } else {
                Cigar::SoftClip(alignment.xstart as u32)
            });
        }

        let mut last = alignment.operations[0];
        let mut k = 1u32;
        for &op in alignment.operations[1..].iter() {
            if op == last {
                k += 1;
            } else {
                add_op(last, k, &mut cigar);
                k = 1;
            }
            last = op;
        }
        add_op(last, k, &mut cigar);
        if alignment.xlen > alignment.xend {
            cigar.push(if hard_clip {
                Cigar::HardClip((alignment.xlen - alignment.xend) as u32)
            } else {
                Cigar::SoftClip((alignment.xlen - alignment.xend) as u32)
            });
        }

        CigarString(cigar)
    }
}

impl TryFrom<&[u8]> for CigarString {
    type Error = Error;

    /// Create a CigarString from given bytes.
    fn try_from(text: &[u8]) -> Result<Self> {
        Self::try_from(str::from_utf8(text).map_err(|_| Error::ParseCigar {
            msg: "unable to parse as UTF8".to_owned(),
        })?)
    }
}

impl TryFrom<&str> for CigarString {
    type Error = Error;

    /// Create a CigarString from given str.
    fn try_from(text: &str) -> Result<Self> {
        lazy_static! {
            // regex for a cigar string operation
            static ref OP_RE: Regex = Regex::new("^(?P<n>[0-9]+)(?P<op>[MIDNSHP=X])").unwrap();
        }
        let mut inner = Vec::new();
        let mut i = 0;
        while i < text.len() {
            if let Some(caps) = OP_RE.captures(&text[i..]) {
                let n = &caps["n"];
                let op = &caps["op"];
                i += n.len() + op.len();
                let n = u32::from_str(n).map_err(|_| Error::ParseCigar {
                    msg: "expected integer".to_owned(),
                })?;
                inner.push(match op {
                    "M" => Cigar::Match(n),
                    "I" => Cigar::Ins(n),
                    "D" => Cigar::Del(n),
                    "N" => Cigar::RefSkip(n),
                    "H" => Cigar::HardClip(n),
                    "S" => Cigar::SoftClip(n),
                    "P" => Cigar::Pad(n),
                    "=" => Cigar::Equal(n),
                    "X" => Cigar::Diff(n),
                    op => {
                        return Err(Error::ParseCigar {
                            msg: format!("operation {} not expected", op),
                        });
                    }
                });
            } else {
                return Err(Error::ParseCigar {
                    msg: "expected cigar operation [0-9]+[MIDNSHP=X]".to_owned(),
                });
            }
        }

        Ok(CigarString(inner))
    }
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
        (&(self.0)).iter()
    }
}

impl fmt::Display for CigarString {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        for op in self {
            fmt.write_fmt(format_args!("{}{}", op.len(), op.char()))?;
        }
        Ok(())
    }
}

#[derive(Eq, PartialEq, Clone, Debug)]
pub struct CigarStringView {
    inner: CigarString,
    pos: i64,
}

impl CigarStringView {
    /// Construct a new CigarStringView from a CigarString at a position
    pub fn new(c: CigarString, pos: i64) -> CigarStringView {
        CigarStringView { inner: c, pos }
    }

    /// Get (exclusive) end position of alignment.
    pub fn end_pos(&self) -> i64 {
        let mut pos = self.pos;
        for c in self {
            match c {
                Cigar::Match(l)
                | Cigar::RefSkip(l)
                | Cigar::Del(l)
                | Cigar::Equal(l)
                | Cigar::Diff(l) => pos += *l as i64,
                // these don't add to end_pos on reference
                Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => (),
            }
        }
        pos
    }

    /// Get number of bases softclipped at the beginning of the alignment.
    pub fn leading_softclips(&self) -> i64 {
        self.first().map_or(0, |cigar| {
            if let Cigar::SoftClip(s) = cigar {
                *s as i64
            } else {
                0
            }
        })
    }

    /// Get number of bases softclipped at the end of the alignment.
    pub fn trailing_softclips(&self) -> i64 {
        self.last().map_or(0, |cigar| {
            if let Cigar::SoftClip(s) = cigar {
                *s as i64
            } else {
                0
            }
        })
    }

    /// Get number of bases hardclipped at the beginning of the alignment.
    pub fn leading_hardclips(&self) -> i64 {
        self.first().map_or(0, |cigar| {
            if let Cigar::HardClip(s) = cigar {
                *s as i64
            } else {
                0
            }
        })
    }

    /// Get number of bases hardclipped at the end of the alignment.
    pub fn trailing_hardclips(&self) -> i64 {
        self.last().map_or(0, |cigar| {
            if let Cigar::HardClip(s) = cigar {
                *s as i64
            } else {
                0
            }
        })
    }

    /// For a given position in the reference, get corresponding position within read.
    /// If reference position is outside of the read alignment, return None.
    ///
    /// # Arguments
    ///
    /// * `ref_pos` - the reference position
    /// * `include_softclips` - if true, softclips will be considered as matches or mismatches
    /// * `include_dels` - if true, positions within deletions will be considered (first reference matching read position after deletion will be returned)
    ///
    pub fn read_pos(
        &self,
        ref_pos: u32,
        include_softclips: bool,
        include_dels: bool,
    ) -> Result<Option<u32>> {
        let mut rpos = self.pos as u32; // reference position
        let mut qpos = 0u32; // position within read
        let mut j = 0; // index into cigar operation vector

        // find first cigar operation referring to qpos = 0 (and thus bases in record.seq()),
        // because all augmentations of qpos and rpos before that are invalid
        for (i, c) in self.iter().enumerate() {
            match c {
                Cigar::Match(_) |
                Cigar::Diff(_)  |
                Cigar::Equal(_) |
                // this is unexpected, but bwa + GATK indel realignment can produce insertions
                // before matching positions
                Cigar::Ins(_) => {
                    j = i;
                    break;
                },
                Cigar::SoftClip(l) => {
                    j = i;
                    if include_softclips {
                        // Alignment starts with softclip and we want to include it in the
                        // projection of the reference position. However, the POS field does not
                        // include the softclip. Hence we have to subtract its length.
                        rpos = rpos.saturating_sub(*l);
                    }
                    break;
                },
                Cigar::Del(_) => {
                    return Err(Error::UnexpectedCigarOperation {
                        msg: "'deletion' (D) found before any operation describing read sequence".to_owned()
                    });
                },
                Cigar::RefSkip(_) => {
                    return Err(Error::UnexpectedCigarOperation {
                        msg: "'reference skip' (N) found before any operation describing read sequence".to_owned()
                    });
                },
                Cigar::HardClip(_) if i > 0 && i < self.len()-1 => {
                    return Err(Error::UnexpectedCigarOperation{
                        msg: "'hard clip' (H) found in between operations, contradicting SAMv1 spec that hard clips can only be at the ends of reads".to_owned()
                    });
                },
                // if we have reached the end of the CigarString with only pads and hard clips, we have no read position matching the variant
                Cigar::Pad(_) | Cigar::HardClip(_) if i == self.len()-1 => return Ok(None),
                // skip leading HardClips and Pads, as they consume neither read sequence nor reference sequence
                Cigar::Pad(_) | Cigar::HardClip(_) => ()
            }
        }

        let contains_ref_pos = |cigar_op_start: u32, cigar_op_length: u32| {
            cigar_op_start <= ref_pos && cigar_op_start + cigar_op_length > ref_pos
        };

        while rpos <= ref_pos && j < self.len() {
            match self[j] {
                // potential SNV evidence
                Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) if contains_ref_pos(rpos, l) => {
                    // difference between desired position and first position of current cigar
                    // operation
                    qpos += ref_pos - rpos;
                    return Ok(Some(qpos));
                }
                Cigar::SoftClip(l) if include_softclips && contains_ref_pos(rpos, l) => {
                    qpos += ref_pos - rpos;
                    return Ok(Some(qpos));
                }
                Cigar::Del(l) if include_dels && contains_ref_pos(rpos, l) => {
                    // qpos shall resemble the start of the deletion
                    return Ok(Some(qpos));
                }
                // for others, just increase pos and qpos as needed
                Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) => {
                    rpos += l;
                    qpos += l;
                    j += 1;
                }
                Cigar::SoftClip(l) => {
                    qpos += l;
                    j += 1;
                    if include_softclips {
                        rpos += l;
                    }
                }
                Cigar::Ins(l) => {
                    qpos += l;
                    j += 1;
                }
                Cigar::RefSkip(l) | Cigar::Del(l) => {
                    rpos += l;
                    j += 1;
                }
                Cigar::Pad(_) => {
                    j += 1;
                }
                Cigar::HardClip(_) if j < self.len() - 1 => {
                    return Err(Error::UnexpectedCigarOperation{
                        msg: "'hard clip' (H) found in between operations, contradicting SAMv1 spec that hard clips can only be at the ends of reads".to_owned()
                    });
                }
                Cigar::HardClip(_) => return Ok(None),
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
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
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

    #[test]
    fn test_cigar_read_pos() {
        let vpos = 5; // variant position

        // Ignore leading HardClip
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c01: 7H                 M  M
        // qpos:                  00 01
        let c01 = CigarString(vec![Cigar::HardClip(7), Cigar::Match(2)]).into_view(4);
        assert_eq!(c01.read_pos(vpos, false, false).unwrap(), Some(1));

        // Skip leading SoftClip or use as pre-POS matches
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c02: 5H2S         M  M  M  M  M  M
        // qpos:  00        02 03 04 05 06 07
        // c02: 5H     S  S  M  M  M  M  M  M
        // qpos:      00 01 02 03 04 05 06 07
        let c02 = CigarString(vec![Cigar::SoftClip(2), Cigar::Match(6)]).into_view(2);
        assert_eq!(c02.read_pos(vpos, false, false).unwrap(), Some(5));
        assert_eq!(c02.read_pos(vpos, true, false).unwrap(), Some(5));

        // Skip leading SoftClip returning None for unmatched reference positiong or use as
        // pre-POS matches
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c03:  3S                      M  M
        // qpos: 00                     03 04
        // c03:                 S  S  S  M  M
        // qpos:               00 01 02 03 04
        let c03 = CigarString(vec![Cigar::SoftClip(3), Cigar::Match(6)]).into_view(6);
        assert_eq!(c03.read_pos(vpos, false, false).unwrap(), None);
        assert_eq!(c03.read_pos(vpos, true, false).unwrap(), Some(2));

        // Skip leading Insertion before variant position
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c04:  3I                X  X  X
        // qpos: 00               03 04 05
        let c04 = CigarString(vec![Cigar::Ins(3), Cigar::Diff(3)]).into_view(4);
        assert_eq!(c04.read_pos(vpos, true, false).unwrap(), Some(4));

        // Matches and deletion before variant position
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c05:        =  =  D  D  X  =  =
        // qpos:      00 01       02 03 04 05
        let c05 = CigarString(vec![
            Cigar::Equal(2),
            Cigar::Del(2),
            Cigar::Diff(1),
            Cigar::Equal(2),
        ])
        .into_view(0);
        assert_eq!(c05.read_pos(vpos, true, false).unwrap(), Some(3));

        // single nucleotide Deletion covering variant position
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c06:                 =  =  D  X  X
        // qpos:               00 01    02 03
        let c06 = CigarString(vec![Cigar::Equal(2), Cigar::Del(1), Cigar::Diff(2)]).into_view(3);
        assert_eq!(c06.read_pos(vpos, false, true).unwrap(), Some(2));
        assert_eq!(c06.read_pos(vpos, false, false).unwrap(), None);

        // three nucleotide Deletion covering variant position
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c07:              =  =  D  D  D  M  M
        // qpos:            00 01          02 03
        let c07 = CigarString(vec![Cigar::Equal(2), Cigar::Del(3), Cigar::Match(2)]).into_view(2);
        assert_eq!(c07.read_pos(vpos, false, true).unwrap(), Some(2));
        assert_eq!(c07.read_pos(vpos, false, false).unwrap(), None);

        // three nucleotide RefSkip covering variant position
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c08:              =  X  N  N  N  M  M
        // qpos:            00 01          02 03
        let c08 = CigarString(vec![
            Cigar::Equal(1),
            Cigar::Diff(1),
            Cigar::RefSkip(3),
            Cigar::Match(2),
        ])
        .into_view(2);
        assert_eq!(c08.read_pos(vpos, false, true).unwrap(), None);
        assert_eq!(c08.read_pos(vpos, false, false).unwrap(), None);

        // internal hard clip before variant pos
        // ref:       00 01 02 03    04 05 06 07 08 09 10 11 12 13 14 15
        // var:                          V
        // c09: 3H           =  = 3H  =  =
        // qpos:            00 01    02 03
        let c09 = CigarString(vec![
            Cigar::HardClip(3),
            Cigar::Equal(2),
            Cigar::HardClip(3),
            Cigar::Equal(2),
        ])
        .into_view(2);
        assert_eq!(c09.read_pos(vpos, false, true).is_err(), true);

        // Deletion right before variant position
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c10:           M  M  D  D  M  M
        // qpos:         00 01       02 03
        let c10 = CigarString(vec![Cigar::Match(2), Cigar::Del(2), Cigar::Match(2)]).into_view(1);
        assert_eq!(c10.read_pos(vpos, false, false).unwrap(), Some(2));

        // Insertion right before variant position
        // ref:       00 01 02 03 04    05 06 07 08 09 10 11 12 13 14 15
        // var:                          V
        // c11:                 M  M 3I  M
        // qpos:               00 01 02 05 06
        let c11 = CigarString(vec![Cigar::Match(2), Cigar::Ins(3), Cigar::Match(2)]).into_view(3);
        assert_eq!(c11.read_pos(vpos, false, false).unwrap(), Some(5));

        // Insertion right after variant position
        // ref:       00 01 02 03 04 05    06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c12:                 M  M  M 2I  =
        // qpos:               00 01 02 03 05
        let c12 = CigarString(vec![Cigar::Match(3), Cigar::Ins(2), Cigar::Equal(1)]).into_view(3);
        assert_eq!(c12.read_pos(vpos, false, false).unwrap(), Some(2));

        // Deletion right after variant position
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c13:                 M  M  M  D  =
        // qpos:               00 01 02    03
        let c13 = CigarString(vec![Cigar::Match(3), Cigar::Del(1), Cigar::Equal(1)]).into_view(3);
        assert_eq!(c13.read_pos(vpos, false, false).unwrap(), Some(2));

        // A messy and complicated example, including a Pad operation
        let vpos2 = 15;
        // ref:       00    01 02    03 04 05    06 07 08 09 10 11 12 13 14 15
        // var:                                                           V
        // c14: 5H3S   = 2P  M  X 3I  M  M  D 2I  =  =  N  N  N  M  M  M  =  =  5S2H
        // qpos:  00  03    04 05 06 09 10    11 13 14          15 16 17 18 19
        let c14 = CigarString(vec![
            Cigar::HardClip(5),
            Cigar::SoftClip(3),
            Cigar::Equal(1),
            Cigar::Pad(2),
            Cigar::Match(1),
            Cigar::Diff(1),
            Cigar::Ins(3),
            Cigar::Match(2),
            Cigar::Del(1),
            Cigar::Ins(2),
            Cigar::Equal(2),
            Cigar::RefSkip(3),
            Cigar::Match(3),
            Cigar::Equal(2),
            Cigar::SoftClip(5),
            Cigar::HardClip(2),
        ])
        .into_view(0);
        assert_eq!(c14.read_pos(vpos2, false, false).unwrap(), Some(19));

        // HardClip after Pad
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c15: 5P1H            =  =  =
        // qpos:               00 01 02
        let c15 =
            CigarString(vec![Cigar::Pad(5), Cigar::HardClip(1), Cigar::Equal(3)]).into_view(3);
        assert_eq!(c15.read_pos(vpos, false, false).is_err(), true);

        // only HardClip and Pad operations
        // c16: 7H5P2H
        let c16 =
            CigarString(vec![Cigar::HardClip(7), Cigar::Pad(5), Cigar::HardClip(2)]).into_view(3);
        assert_eq!(c16.read_pos(vpos, false, false).unwrap(), None);
    }

    #[test]
    fn test_clone() {
        let mut rec = Record::new();
        rec.set_pos(300);
        rec.set_qname(b"read1");
        let clone = rec.clone();
        assert_eq!(rec, clone);
    }

    #[test]
    fn test_flags() {
        let mut rec = Record::new();

        rec.set_paired();
        assert_eq!(rec.is_paired(), true);

        rec.set_supplementary();
        assert_eq!(rec.is_supplementary(), true);
        assert_eq!(rec.is_supplementary(), true);

        rec.unset_paired();
        assert_eq!(rec.is_paired(), false);
        assert_eq!(rec.is_supplementary(), true);

        rec.unset_supplementary();
        assert_eq!(rec.is_paired(), false);
        assert_eq!(rec.is_supplementary(), false);
    }

    #[test]
    fn test_cigar_parse() {
        let cigar = "1S20M1D2I3X1=2H";
        let parsed = CigarString::try_from(cigar).unwrap();
        assert_eq!(parsed.to_string(), cigar);
    }
}

#[cfg(test)]
mod alignment_cigar_tests {
    use super::*;
    use bio_types::alignment::AlignmentOperation::{Del, Ins, Match, Subst, Xclip, Yclip};
    use bio_types::alignment::{Alignment, AlignmentMode};

    #[test]
    fn test_cigar() {
        let alignment = Alignment {
            score: 5,
            xstart: 3,
            ystart: 0,
            xend: 9,
            yend: 10,
            ylen: 10,
            xlen: 10,
            operations: vec![Match, Match, Match, Subst, Ins, Ins, Del, Del],
            mode: AlignmentMode::Semiglobal,
        };
        assert_eq!(alignment.cigar(false), "3S3=1X2I2D1S");
        assert_eq!(
            CigarString::from_alignment(&alignment, false).0,
            vec![
                Cigar::SoftClip(3),
                Cigar::Equal(3),
                Cigar::Diff(1),
                Cigar::Ins(2),
                Cigar::Del(2),
                Cigar::SoftClip(1),
            ]
        );

        let alignment = Alignment {
            score: 5,
            xstart: 0,
            ystart: 5,
            xend: 4,
            yend: 10,
            ylen: 10,
            xlen: 5,
            operations: vec![Yclip(5), Match, Subst, Subst, Ins, Del, Del, Xclip(1)],
            mode: AlignmentMode::Custom,
        };
        assert_eq!(alignment.cigar(false), "1=2X1I2D1S");
        assert_eq!(alignment.cigar(true), "1=2X1I2D1H");
        assert_eq!(
            CigarString::from_alignment(&alignment, false).0,
            vec![
                Cigar::Equal(1),
                Cigar::Diff(2),
                Cigar::Ins(1),
                Cigar::Del(2),
                Cigar::SoftClip(1),
            ]
        );
        assert_eq!(
            CigarString::from_alignment(&alignment, true).0,
            vec![
                Cigar::Equal(1),
                Cigar::Diff(2),
                Cigar::Ins(1),
                Cigar::Del(2),
                Cigar::HardClip(1),
            ]
        );

        let alignment = Alignment {
            score: 5,
            xstart: 0,
            ystart: 5,
            xend: 3,
            yend: 8,
            ylen: 10,
            xlen: 3,
            operations: vec![Yclip(5), Subst, Match, Subst, Yclip(2)],
            mode: AlignmentMode::Custom,
        };
        assert_eq!(alignment.cigar(false), "1X1=1X");
        assert_eq!(
            CigarString::from_alignment(&alignment, false).0,
            vec![Cigar::Diff(1), Cigar::Equal(1), Cigar::Diff(1)]
        );

        let alignment = Alignment {
            score: 5,
            xstart: 0,
            ystart: 5,
            xend: 3,
            yend: 8,
            ylen: 10,
            xlen: 3,
            operations: vec![Subst, Match, Subst],
            mode: AlignmentMode::Semiglobal,
        };
        assert_eq!(alignment.cigar(false), "1X1=1X");
        assert_eq!(
            CigarString::from_alignment(&alignment, false).0,
            vec![Cigar::Diff(1), Cigar::Equal(1), Cigar::Diff(1)]
        );
    }
}
