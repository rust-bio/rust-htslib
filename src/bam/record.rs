// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::convert::TryFrom;
use std::convert::TryInto;
use std::ffi;
use std::fmt;
use std::marker::PhantomData;
use std::mem::{size_of, MaybeUninit};
use std::ops;
use std::os::raw::c_char;
use std::rc::Rc;
use std::slice;
use std::str;

use byteorder::{LittleEndian, ReadBytesExt};

use crate::bam::Error;
use crate::bam::HeaderView;
use crate::errors::Result;
use crate::htslib;
use crate::utils;
#[cfg(feature = "serde_feature")]
use serde::{self, Deserialize, Serialize};

use bio_types::alignment::{Alignment, AlignmentMode, AlignmentOperation};
use bio_types::genome;
use bio_types::sequence::SequenceRead;
use bio_types::sequence::SequenceReadPairOrientation;
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
        let mut record = Record {
            inner: unsafe { MaybeUninit::zeroed().assume_init() },
            own: true,
            cigar: None,
            header: None,
        };
        // The read/query name needs to be set as empty to properly initialize
        // the record
        record.set_qname(b"");
        // Developer note: these are needed so the returned record is properly
        // initialized as unmapped.
        record.set_unmapped();
        record.set_tid(-1);
        record.set_pos(-1);
        record.set_mpos(-1);
        record.set_mtid(-1);
        record
    }

    pub fn from_inner(from: *mut htslib::bam1_t) -> Self {
        Record {
            inner: {
                #[allow(clippy::uninit_assumed_init, invalid_value)]
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
            s: sam_copy.as_ptr() as *mut c_char,
            l: sam_copy.len(),
            m: sam_copy.len(),
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
            Err(Error::BamParseSAM {
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
    pub fn strand(&self) -> ReqStrand {
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
        self.inner().core.isize_
    }

    /// Set insert size.
    pub fn set_insert_size(&mut self, insert_size: i64) {
        self.inner_mut().core.isize_ = insert_size;
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
            + self.seq_len().div_ceil(2)
            + self.seq_len();
        let new_aux_offset = q_len + extranul + cigar_width + seq.len().div_ceil(2) + qual.len();
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
                data[i + j / 2] = (ENCODE_BASE[seq[j] as usize] << 4)
                    | (if j + 1 < seq.len() {
                        ENCODE_BASE[seq[j + 1] as usize]
                    } else {
                        0
                    });
            }
            self.inner_mut().core.l_qseq = seq.len() as i32;
            i += seq.len().div_ceil(2);
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

    /// Replace current cigar with a new one.
    pub fn set_cigar(&mut self, new_cigar: Option<&CigarString>) {
        self.cigar = None;

        let qname_data_len = self.qname_capacity();
        let old_cigar_data_len = self.cigar_len() * 4;

        // Length of data after cigar
        let other_data_len = self.inner_mut().l_data - (qname_data_len + old_cigar_data_len) as i32;

        let new_cigar_len = match new_cigar {
            Some(x) => x.len(),
            None => 0,
        };
        let new_cigar_data_len = new_cigar_len * 4;

        if new_cigar_data_len < old_cigar_data_len {
            self.inner_mut().l_data -= (old_cigar_data_len - new_cigar_data_len) as i32;
        } else if new_cigar_data_len > old_cigar_data_len {
            self.inner_mut().l_data += (new_cigar_data_len - old_cigar_data_len) as i32;

            // Reallocate if necessary
            if (self.inner().m_data as i32) < self.inner().l_data {
                // Verbosity due to lexical borrowing
                let l_data = self.inner().l_data;
                self.realloc_var_data(l_data as usize);
            }
        }

        if new_cigar_data_len != old_cigar_data_len {
            // Move other data to new location
            unsafe {
                ::libc::memmove(
                    self.inner.data.add(qname_data_len + new_cigar_data_len) as *mut ::libc::c_void,
                    self.inner.data.add(qname_data_len + old_cigar_data_len) as *mut ::libc::c_void,
                    other_data_len as usize,
                );
            }
        }

        // Copy cigar data
        if let Some(cigar_string) = new_cigar {
            let cigar_data = unsafe {
                #[allow(clippy::cast_ptr_alignment)]
                slice::from_raw_parts_mut(
                    self.inner.data.add(qname_data_len) as *mut u32,
                    cigar_string.len(),
                )
            };
            for (i, c) in cigar_string.iter().enumerate() {
                cigar_data[i] = c.encode();
            }
        }
        self.inner_mut().core.n_cigar = new_cigar_len as u32;
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

    pub fn seq_len(&self) -> usize {
        self.inner().core.l_qseq as usize
    }

    fn seq_data(&self) -> &[u8] {
        let offset = self.qname_capacity() + self.cigar_len() * 4;
        &self.data()[offset..][..self.seq_len().div_ceil(2)]
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
        &self.data()[self.qname_capacity() + self.cigar_len() * 4 + self.seq_len().div_ceil(2)..]
            [..self.seq_len()]
    }

    /// Look up an auxiliary field by its tag.
    ///
    /// Only the first two bytes of a given tag are used for the look-up of a field.
    /// See [`Aux`] for more details.
    pub fn aux(&self, tag: &[u8]) -> Result<Aux<'_>> {
        let c_str = ffi::CString::new(tag).map_err(|_| Error::BamAuxStringError)?;
        let aux = unsafe {
            htslib::bam_aux_get(
                &self.inner as *const htslib::bam1_t,
                c_str.as_ptr() as *mut c_char,
            )
        };
        unsafe { Self::read_aux_field(aux).map(|(aux_field, _length)| aux_field) }
    }

    unsafe fn read_aux_field<'a>(aux: *const u8) -> Result<(Aux<'a>, usize)> {
        const TAG_LEN: isize = 2;
        // Used for skipping type identifier
        const TYPE_ID_LEN: isize = 1;

        if aux.is_null() {
            return Err(Error::BamAuxTagNotFound);
        }

        let (data, type_size) = match *aux {
            b'A' => {
                let type_size = size_of::<u8>();
                (Aux::Char(*aux.offset(TYPE_ID_LEN)), type_size)
            }
            b'c' => {
                let type_size = size_of::<i8>();
                (Aux::I8(*aux.offset(TYPE_ID_LEN).cast::<i8>()), type_size)
            }
            b'C' => {
                let type_size = size_of::<u8>();
                (Aux::U8(*aux.offset(TYPE_ID_LEN)), type_size)
            }
            b's' => {
                let type_size = size_of::<i16>();
                (
                    Aux::I16(
                        slice::from_raw_parts(aux.offset(TYPE_ID_LEN), type_size)
                            .read_i16::<LittleEndian>()
                            .map_err(|_| Error::BamAuxParsingError)?,
                    ),
                    type_size,
                )
            }
            b'S' => {
                let type_size = size_of::<u16>();
                (
                    Aux::U16(
                        slice::from_raw_parts(aux.offset(TYPE_ID_LEN), type_size)
                            .read_u16::<LittleEndian>()
                            .map_err(|_| Error::BamAuxParsingError)?,
                    ),
                    type_size,
                )
            }
            b'i' => {
                let type_size = size_of::<i32>();
                (
                    Aux::I32(
                        slice::from_raw_parts(aux.offset(TYPE_ID_LEN), type_size)
                            .read_i32::<LittleEndian>()
                            .map_err(|_| Error::BamAuxParsingError)?,
                    ),
                    type_size,
                )
            }
            b'I' => {
                let type_size = size_of::<u32>();
                (
                    Aux::U32(
                        slice::from_raw_parts(aux.offset(TYPE_ID_LEN), type_size)
                            .read_u32::<LittleEndian>()
                            .map_err(|_| Error::BamAuxParsingError)?,
                    ),
                    type_size,
                )
            }
            b'f' => {
                let type_size = size_of::<f32>();
                (
                    Aux::Float(
                        slice::from_raw_parts(aux.offset(TYPE_ID_LEN), type_size)
                            .read_f32::<LittleEndian>()
                            .map_err(|_| Error::BamAuxParsingError)?,
                    ),
                    type_size,
                )
            }
            b'd' => {
                let type_size = size_of::<f64>();
                (
                    Aux::Double(
                        slice::from_raw_parts(aux.offset(TYPE_ID_LEN), type_size)
                            .read_f64::<LittleEndian>()
                            .map_err(|_| Error::BamAuxParsingError)?,
                    ),
                    type_size,
                )
            }
            b'Z' | b'H' => {
                let c_str = ffi::CStr::from_ptr(aux.offset(TYPE_ID_LEN).cast::<c_char>());
                let rust_str = c_str.to_str().map_err(|_| Error::BamAuxParsingError)?;
                (Aux::String(rust_str), c_str.to_bytes_with_nul().len())
            }
            b'B' => {
                const ARRAY_INNER_TYPE_LEN: isize = 1;
                const ARRAY_COUNT_LEN: isize = 4;

                // Used for skipping metadata
                let array_data_offset = TYPE_ID_LEN + ARRAY_INNER_TYPE_LEN + ARRAY_COUNT_LEN;

                let length =
                    slice::from_raw_parts(aux.offset(TYPE_ID_LEN + ARRAY_INNER_TYPE_LEN), 4)
                        .read_u32::<LittleEndian>()
                        .map_err(|_| Error::BamAuxParsingError)? as usize;

                // Return tuples of an `Aux` enum and the length of data + metadata in bytes
                let (array_data, array_size) = match *aux.offset(TYPE_ID_LEN) {
                    b'c' => (
                        Aux::ArrayI8(AuxArray::<'a, i8>::from_bytes(slice::from_raw_parts(
                            aux.offset(array_data_offset),
                            length,
                        ))),
                        length,
                    ),
                    b'C' => (
                        Aux::ArrayU8(AuxArray::<'a, u8>::from_bytes(slice::from_raw_parts(
                            aux.offset(array_data_offset),
                            length,
                        ))),
                        length,
                    ),
                    b's' => (
                        Aux::ArrayI16(AuxArray::<'a, i16>::from_bytes(slice::from_raw_parts(
                            aux.offset(array_data_offset),
                            length * size_of::<i16>(),
                        ))),
                        length * std::mem::size_of::<i16>(),
                    ),
                    b'S' => (
                        Aux::ArrayU16(AuxArray::<'a, u16>::from_bytes(slice::from_raw_parts(
                            aux.offset(array_data_offset),
                            length * size_of::<u16>(),
                        ))),
                        length * std::mem::size_of::<u16>(),
                    ),
                    b'i' => (
                        Aux::ArrayI32(AuxArray::<'a, i32>::from_bytes(slice::from_raw_parts(
                            aux.offset(array_data_offset),
                            length * size_of::<i32>(),
                        ))),
                        length * std::mem::size_of::<i32>(),
                    ),
                    b'I' => (
                        Aux::ArrayU32(AuxArray::<'a, u32>::from_bytes(slice::from_raw_parts(
                            aux.offset(array_data_offset),
                            length * size_of::<u32>(),
                        ))),
                        length * std::mem::size_of::<u32>(),
                    ),
                    b'f' => (
                        Aux::ArrayFloat(AuxArray::<f32>::from_bytes(slice::from_raw_parts(
                            aux.offset(array_data_offset),
                            length * size_of::<f32>(),
                        ))),
                        length * std::mem::size_of::<f32>(),
                    ),
                    _ => {
                        return Err(Error::BamAuxUnknownType);
                    }
                };
                (
                    array_data,
                    // Offset: array-specific metadata + array size
                    ARRAY_INNER_TYPE_LEN as usize + ARRAY_COUNT_LEN as usize + array_size,
                )
            }
            _ => {
                return Err(Error::BamAuxUnknownType);
            }
        };

        // Offset: metadata + type size
        Ok((data, TAG_LEN as usize + TYPE_ID_LEN as usize + type_size))
    }

    /// Returns an iterator over the auxiliary fields of the record.
    ///
    /// When an error occurs, the `Err` variant will be returned
    /// and the iterator will not be able to advance anymore.
    pub fn aux_iter(&self) -> AuxIter<'_> {
        AuxIter {
            // In order to get to the aux data section of a `bam::Record`
            // we need to skip fields in front of it
            aux: &self.data()[
                // NUL terminated read name:
                self.qname_capacity()
                // CIGAR (uint32_t):
                + self.cigar_len() * std::mem::size_of::<u32>()
                // Read sequence (4-bit encoded):
                + self.seq_len().div_ceil(2)
                // Base qualities (char):
                + self.seq_len()..],
        }
    }

    /// Add auxiliary data.
    pub fn push_aux(&mut self, tag: &[u8], value: Aux<'_>) -> Result<()> {
        // Don't allow pushing aux data when the given tag is already present in the record.
        // `htslib` seems to allow this (for non-array values), which can lead to problems
        // since retrieving aux fields consumes &[u8; 2] and yields one field only.
        if self.aux(tag).is_ok() {
            return Err(Error::BamAuxTagAlreadyPresent);
        }

        let ctag = tag.as_ptr() as *mut c_char;
        let ret = unsafe {
            match value {
                Aux::Char(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'A' as c_char,
                    size_of::<u8>() as i32,
                    [v].as_mut_ptr(),
                ),
                Aux::I8(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'c' as c_char,
                    size_of::<i8>() as i32,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::U8(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'C' as c_char,
                    size_of::<u8>() as i32,
                    [v].as_mut_ptr(),
                ),
                Aux::I16(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b's' as c_char,
                    size_of::<i16>() as i32,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::U16(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'S' as c_char,
                    size_of::<u16>() as i32,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::I32(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'i' as c_char,
                    size_of::<i32>() as i32,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::U32(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'I' as c_char,
                    size_of::<u32>() as i32,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::Float(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'f' as c_char,
                    size_of::<f32>() as i32,
                    [v].as_mut_ptr() as *mut u8,
                ),
                // Not part of specs but implemented in `htslib`:
                Aux::Double(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'd' as c_char,
                    size_of::<f64>() as i32,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::String(v) => {
                    let c_str = ffi::CString::new(v).map_err(|_| Error::BamAuxStringError)?;
                    htslib::bam_aux_append(
                        self.inner_ptr_mut(),
                        ctag,
                        b'Z' as c_char,
                        (v.len() + 1) as i32,
                        c_str.as_ptr() as *mut u8,
                    )
                }
                Aux::HexByteArray(v) => {
                    let c_str = ffi::CString::new(v).map_err(|_| Error::BamAuxStringError)?;
                    htslib::bam_aux_append(
                        self.inner_ptr_mut(),
                        ctag,
                        b'H' as c_char,
                        (v.len() + 1) as i32,
                        c_str.as_ptr() as *mut u8,
                    )
                }
                // Not sure it's safe to cast an immutable slice to a mutable pointer in the following branches
                Aux::ArrayI8(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'c',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'c',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayU8(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'C',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'C',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayI16(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b's',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b's',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayU16(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'S',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'S',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayI32(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'i',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'i',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayU32(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'I',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'I',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayFloat(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'f',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'f',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
            }
        };

        if ret < 0 {
            Err(Error::BamAux)
        } else {
            Ok(())
        }
    }

    // Delete auxiliary tag.
    pub fn remove_aux(&mut self, tag: &[u8]) -> Result<()> {
        let c_str = ffi::CString::new(tag).map_err(|_| Error::BamAuxStringError)?;
        let aux = unsafe {
            htslib::bam_aux_get(
                &self.inner as *const htslib::bam1_t,
                c_str.as_ptr() as *mut c_char,
            )
        };
        unsafe {
            if aux.is_null() {
                Err(Error::BamAuxTagNotFound)
            } else {
                htslib::bam_aux_del(self.inner_ptr_mut(), aux);
                Ok(())
            }
        }
    }

    /// Access the base modifications associated with this Record through the MM tag.
    /// Example:
    /// ```
    ///    use rust_htslib::bam::{Read, Reader, Record};
    ///    let mut bam = Reader::from_path("test/base_mods/MM-orient.sam").unwrap();
    ///    let mut mod_count = 0;
    ///    for r in bam.records() {
    ///        let record = r.unwrap();
    ///        if let Ok(mods) = record.basemods_iter() {
    ///            // print metadata for the modifications present in this record
    ///            for mod_code in mods.recorded() {
    ///                if let Ok(mod_metadata) = mods.query_type(*mod_code) {
    ///                    println!("mod found with code {}/{} flags: [{} {} {}]",
    ///                              mod_code, *mod_code as u8 as char,
    ///                              mod_metadata.strand, mod_metadata.implicit, mod_metadata.canonical as u8 as char);
    ///                }
    ///            }
    ///
    ///            // iterate over the modifications in this record
    ///            // the modifications are returned as a tuple with the
    ///            // position within SEQ and an hts_base_mod struct
    ///            for res in mods {
    ///                if let Ok( (position, m) ) = res {
    ///                    println!("{} {},{}", position, m.modified_base as u8 as char, m.qual);
    ///                    mod_count += 1;
    ///                }
    ///            }
    ///        };
    ///    }
    ///    assert_eq!(mod_count, 14);
    /// ```
    pub fn basemods_iter(&self) -> Result<BaseModificationsIter<'_>> {
        BaseModificationsIter::new(self)
    }

    /// An iterator that returns all of the modifications for each position as a vector.
    /// This is useful for the case where multiple possible modifications can be annotated
    /// at a single position (for example a C could be 5-mC or 5-hmC)
    pub fn basemods_position_iter(&self) -> Result<BaseModificationsPositionIter<'_>> {
        BaseModificationsPositionIter::new(self)
    }

    /// Infer read pair orientation from record. Returns `SequenceReadPairOrientation::None` if record
    /// is not paired, mates are not mapping to the same contig, or mates start at the
    /// same position.
    pub fn read_pair_orientation(&self) -> SequenceReadPairOrientation {
        if self.is_paired()
            && !self.is_unmapped()
            && !self.is_mate_unmapped()
            && self.tid() == self.mtid()
        {
            if self.pos() == self.mpos() {
                // both reads start at the same position, we cannot decide on the orientation.
                return SequenceReadPairOrientation::None;
            }

            let (pos_1, pos_2, fwd_1, fwd_2) = if self.is_first_in_template() {
                (
                    self.pos(),
                    self.mpos(),
                    !self.is_reverse(),
                    !self.is_mate_reverse(),
                )
            } else {
                (
                    self.mpos(),
                    self.pos(),
                    !self.is_mate_reverse(),
                    !self.is_reverse(),
                )
            };

            if pos_1 < pos_2 {
                match (fwd_1, fwd_2) {
                    (true, true) => SequenceReadPairOrientation::F1F2,
                    (true, false) => SequenceReadPairOrientation::F1R2,
                    (false, true) => SequenceReadPairOrientation::R1F2,
                    (false, false) => SequenceReadPairOrientation::R1R2,
                }
            } else {
                match (fwd_2, fwd_1) {
                    (true, true) => SequenceReadPairOrientation::F2F1,
                    (true, false) => SequenceReadPairOrientation::F2R1,
                    (false, true) => SequenceReadPairOrientation::R2F1,
                    (false, false) => SequenceReadPairOrientation::R2R1,
                }
            }
        } else {
            SequenceReadPairOrientation::None
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
        *decode_base_unchecked(encoded_base(self.seq_data(), i))
    }

    fn base_qual(&self, i: usize) -> u8 {
        self.qual()[i]
    }

    fn len(&self) -> usize {
        self.seq_len()
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
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

/// Auxiliary record data
///
/// The specification allows a wide range of types to be stored as an auxiliary data field of a BAM record.
///
/// Please note that the [`Aux::Double`] variant is _not_ part of the specification, but it is supported by `htslib`.
///
/// # Examples
///
/// ```
/// use rust_htslib::{
///     bam,
///     bam::record::{Aux, AuxArray},
///     errors::Error,
/// };
///
/// //Set up BAM record
/// let bam_header = bam::Header::new();
/// let mut record = bam::Record::from_sam(
///     &mut bam::HeaderView::from_header(&bam_header),
///     "ali1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tFFFF".as_bytes(),
/// )
/// .unwrap();
///
/// // Add an integer field
/// let aux_integer_field = Aux::I32(1234);
/// record.push_aux(b"XI", aux_integer_field).unwrap();
///
/// match record.aux(b"XI") {
///     Ok(value) => {
///         // Typically, callers expect an aux field to be of a certain type.
///         // If that's not the case, the value can be `match`ed exhaustively.
///         if let Aux::I32(v) = value {
///             assert_eq!(v, 1234);
///         }
///     }
///     Err(e) => {
///         panic!("Error reading aux field: {}", e);
///     }
/// }
///
/// // Add an array field
/// let array_like_data = vec![0.4, 0.3, 0.2, 0.1];
/// let slice_of_data = &array_like_data;
/// let aux_array: AuxArray<f32> = slice_of_data.into();
/// let aux_array_field = Aux::ArrayFloat(aux_array);
/// record.push_aux(b"XA", aux_array_field).unwrap();
///
/// if let Ok(Aux::ArrayFloat(array)) = record.aux(b"XA") {
///     let read_array = array.iter().collect::<Vec<_>>();
///     assert_eq!(read_array, array_like_data);
/// } else {
///     panic!("Could not read array data");
/// }
/// ```
#[derive(Debug, PartialEq)]
pub enum Aux<'a> {
    Char(u8),
    I8(i8),
    U8(u8),
    I16(i16),
    U16(u16),
    I32(i32),
    U32(u32),
    Float(f32),
    Double(f64), // Not part of specs but implemented in `htslib`
    String(&'a str),
    HexByteArray(&'a str),
    ArrayI8(AuxArray<'a, i8>),
    ArrayU8(AuxArray<'a, u8>),
    ArrayI16(AuxArray<'a, i16>),
    ArrayU16(AuxArray<'a, u16>),
    ArrayI32(AuxArray<'a, i32>),
    ArrayU32(AuxArray<'a, u32>),
    ArrayFloat(AuxArray<'a, f32>),
}

unsafe impl Send for Aux<'_> {}
unsafe impl Sync for Aux<'_> {}

/// Types that can be used in aux arrays.
pub trait AuxArrayElement: Copy {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self>;
}

impl AuxArrayElement for i8 {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        std::io::Cursor::new(bytes).read_i8().ok()
    }
}
impl AuxArrayElement for u8 {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        std::io::Cursor::new(bytes).read_u8().ok()
    }
}
impl AuxArrayElement for i16 {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        std::io::Cursor::new(bytes).read_i16::<LittleEndian>().ok()
    }
}
impl AuxArrayElement for u16 {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        std::io::Cursor::new(bytes).read_u16::<LittleEndian>().ok()
    }
}
impl AuxArrayElement for i32 {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        std::io::Cursor::new(bytes).read_i32::<LittleEndian>().ok()
    }
}
impl AuxArrayElement for u32 {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        std::io::Cursor::new(bytes).read_u32::<LittleEndian>().ok()
    }
}
impl AuxArrayElement for f32 {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        std::io::Cursor::new(bytes).read_f32::<LittleEndian>().ok()
    }
}

/// Provides access to aux arrays.
///
/// Provides methods to either retrieve single elements or an iterator over the
/// array.
///
/// This type is used for wrapping both, array data that was read from a
/// BAM record and slices of data that are going to be stored in one.
///
/// In order to be able to add an `AuxArray` field to a BAM record, `AuxArray`s
/// can be constructed via the `From` trait which is implemented for all
/// supported types (see [`AuxArrayElement`] for a list).
///
/// # Examples
///
/// ```
/// use rust_htslib::{
///     bam,
///     bam::record::{Aux, AuxArray},
/// };
///
/// //Set up BAM record
/// let bam_header = bam::Header::new();
/// let mut record = bam::Record::from_sam(
///     &mut bam::HeaderView::from_header(&bam_header),
///     "ali1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tFFFF".as_bytes(),
/// ).unwrap();
///
/// let data = vec![0.4, 0.3, 0.2, 0.1];
/// let slice_of_data = &data;
/// let aux_array: AuxArray<f32> = slice_of_data.into();
/// let aux_field = Aux::ArrayFloat(aux_array);
/// record.push_aux(b"XA", aux_field);
///
/// if let Ok(Aux::ArrayFloat(array)) = record.aux(b"XA") {
///     // Retrieve the second element from the array
///     assert_eq!(array.get(1).unwrap(), 0.3);
///     // Iterate over the array and collect it into a `Vec`
///     let read_array = array.iter().collect::<Vec<_>>();
///     assert_eq!(read_array, data);
/// } else {
///     panic!("Could not read array data");
/// }
/// ```
#[derive(Debug)]
pub enum AuxArray<'a, T> {
    TargetType(AuxArrayTargetType<'a, T>),
    RawLeBytes(AuxArrayRawLeBytes<'a, T>),
}

impl<T> PartialEq<AuxArray<'_, T>> for AuxArray<'_, T>
where
    T: AuxArrayElement + PartialEq,
{
    fn eq(&self, other: &AuxArray<'_, T>) -> bool {
        use AuxArray::*;
        match (self, other) {
            (TargetType(v), TargetType(v_other)) => v == v_other,
            (RawLeBytes(v), RawLeBytes(v_other)) => v == v_other,
            (TargetType(_), RawLeBytes(_)) => self.iter().eq(other.iter()),
            (RawLeBytes(_), TargetType(_)) => self.iter().eq(other.iter()),
        }
    }
}

/// Create AuxArrays from slices of allowed target types.
impl<'a, I, T> From<&'a T> for AuxArray<'a, I>
where
    I: AuxArrayElement,
    T: AsRef<[I]> + ?Sized,
{
    fn from(src: &'a T) -> Self {
        AuxArray::TargetType(AuxArrayTargetType {
            slice: src.as_ref(),
        })
    }
}

impl<'a, T> AuxArray<'a, T>
where
    T: AuxArrayElement,
{
    /// Returns the element at a position or None if out of bounds.
    pub fn get(&self, index: usize) -> Option<T> {
        match self {
            AuxArray::TargetType(v) => v.get(index),
            AuxArray::RawLeBytes(v) => v.get(index),
        }
    }

    /// Returns the number of elements in the array.
    pub fn len(&self) -> usize {
        match self {
            AuxArray::TargetType(a) => a.len(),
            AuxArray::RawLeBytes(a) => a.len(),
        }
    }

    /// Returns true if the array contains no elements.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns an iterator over the array.
    pub fn iter(&self) -> AuxArrayIter<'_, T> {
        AuxArrayIter {
            index: 0,
            array: self,
        }
    }

    /// Create AuxArrays from raw byte slices borrowed from `bam::Record`.
    fn from_bytes(bytes: &'a [u8]) -> Self {
        Self::RawLeBytes(AuxArrayRawLeBytes {
            slice: bytes,
            phantom_data: PhantomData,
        })
    }
}

/// Encapsulates slice of target type.
#[doc(hidden)]
#[derive(Debug, PartialEq)]
pub struct AuxArrayTargetType<'a, T> {
    slice: &'a [T],
}

impl<T> AuxArrayTargetType<'_, T>
where
    T: AuxArrayElement,
{
    fn get(&self, index: usize) -> Option<T> {
        self.slice.get(index).copied()
    }

    fn len(&self) -> usize {
        self.slice.len()
    }
}

/// Encapsulates slice of raw bytes to prevent it from being accidentally accessed.
#[doc(hidden)]
#[derive(Debug, PartialEq)]
pub struct AuxArrayRawLeBytes<'a, T> {
    slice: &'a [u8],
    phantom_data: PhantomData<T>,
}

impl<T> AuxArrayRawLeBytes<'_, T>
where
    T: AuxArrayElement,
{
    fn get(&self, index: usize) -> Option<T> {
        let type_size = std::mem::size_of::<T>();
        if index * type_size + type_size > self.slice.len() {
            return None;
        }
        T::from_le_bytes(&self.slice[index * type_size..][..type_size])
    }

    fn len(&self) -> usize {
        self.slice.len() / std::mem::size_of::<T>()
    }
}

/// Aux array iterator
///
/// This struct is created by the [`AuxArray::iter`] method.
pub struct AuxArrayIter<'a, T> {
    index: usize,
    array: &'a AuxArray<'a, T>,
}

impl<T> Iterator for AuxArrayIter<'_, T>
where
    T: AuxArrayElement,
{
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        let value = self.array.get(self.index);
        self.index += 1;
        value
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let array_length = self.array.len() - self.index;
        (array_length, Some(array_length))
    }
}

/// Auxiliary data iterator
///
/// This struct is created by the [`Record::aux_iter`] method.
///
/// This iterator returns `Result`s that wrap tuples containing
/// a slice which represents the two-byte tag (`&[u8; 2]`) as
/// well as an `Aux` enum that wraps the associated value.
///
/// When an error occurs, the `Err` variant will be returned
/// and the iterator will not be able to advance anymore.
pub struct AuxIter<'a> {
    aux: &'a [u8],
}

impl<'a> Iterator for AuxIter<'a> {
    type Item = Result<(&'a [u8], Aux<'a>)>;

    fn next(&mut self) -> Option<Self::Item> {
        // We're finished
        if self.aux.is_empty() {
            return None;
        }
        // Incomplete aux data
        if (1..=3).contains(&self.aux.len()) {
            // In the case of an error, we can not safely advance in the aux data, so we terminate the Iteration
            self.aux = &[];
            return Some(Err(Error::BamAuxParsingError));
        }
        let tag = &self.aux[..2];
        Some(unsafe {
            let data_ptr = self.aux[2..].as_ptr();
            Record::read_aux_field(data_ptr)
                .map(|(aux, offset)| {
                    self.aux = &self.aux[offset..];
                    (tag, aux)
                })
                .inspect_err(|_e| {
                    // In the case of an error, we can not safely advance in the aux data, so we terminate the Iteration
                    self.aux = &[];
                })
        })
    }
}

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

#[inline]
fn encoded_base(encoded_seq: &[u8], i: usize) -> u8 {
    (encoded_seq[i / 2] >> ((!i & 1) << 2)) & 0b1111
}

#[inline]
unsafe fn encoded_base_unchecked(encoded_seq: &[u8], i: usize) -> u8 {
    (encoded_seq.get_unchecked(i / 2) >> ((!i & 1) << 2)) & 0b1111
}

#[inline]
fn decode_base_unchecked(base: u8) -> &'static u8 {
    unsafe { DECODE_BASE.get_unchecked(base as usize) }
}

/// The sequence of a record.
#[derive(Debug, Copy, Clone)]
pub struct Seq<'a> {
    pub encoded: &'a [u8],
    len: usize,
}

impl Seq<'_> {
    /// Return encoded base. Complexity: O(1).
    #[inline]
    pub fn encoded_base(&self, i: usize) -> u8 {
        encoded_base(self.encoded, i)
    }

    /// Return encoded base. Complexity: O(1).
    ///
    /// # Safety
    ///
    /// TODO
    #[inline]
    pub unsafe fn encoded_base_unchecked(&self, i: usize) -> u8 {
        encoded_base_unchecked(self.encoded, i)
    }

    /// Obtain decoded base without performing bounds checking.
    /// Use index based access seq()[i], for checked, safe access.
    /// Complexity: O(1).
    ///
    /// # Safety
    ///
    /// TODO
    #[inline]
    pub unsafe fn decoded_base_unchecked(&self, i: usize) -> u8 {
        *decode_base_unchecked(self.encoded_base_unchecked(i))
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

impl ops::Index<usize> for Seq<'_> {
    type Output = u8;

    /// Return decoded base at given position within read. Complexity: O(1).
    fn index(&self, index: usize) -> &u8 {
        decode_base_unchecked(self.encoded_base(index))
    }
}

unsafe impl Send for Seq<'_> {}
unsafe impl Sync for Seq<'_> {}

#[cfg_attr(feature = "serde_feature", derive(Serialize, Deserialize))]
#[derive(PartialEq, PartialOrd, Eq, Debug, Clone, Copy, Hash)]
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
            Cigar::Ins(len) => (len << 4) | 1,
            Cigar::Del(len) => (len << 4) | 2,
            Cigar::RefSkip(len) => (len << 4) | 3,
            Cigar::SoftClip(len) => (len << 4) | 4,
            Cigar::HardClip(len) => (len << 4) | 5,
            Cigar::Pad(len) => (len << 4) | 6,
            Cigar::Equal(len) => (len << 4) | 7,
            Cigar::Diff(len) => (len << 4) | 8,
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
    #[cfg_attr(feature = "serde_feature", derive(Serialize, Deserialize))]
    #[derive(NewtypeDeref,
            NewtypeDerefMut,
             NewtypeIndex(usize),
             NewtypeIndexMut(usize),
             NewtypeFrom,
             PartialEq,
             PartialOrd,
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

    /// Create a CigarString from given &[u8].
    /// # Example
    /// ```
    /// use rust_htslib::bam::record::*;
    /// use rust_htslib::bam::record::CigarString;
    /// use rust_htslib::bam::record::Cigar::*;
    /// use std::convert::TryFrom;
    ///
    /// let cigar_str = "2H10M5X3=2H".as_bytes();
    /// let cigar = CigarString::try_from(cigar_str)
    ///     .expect("Unable to parse cigar string.");
    /// let expected_cigar = CigarString(vec![
    ///     HardClip(2),
    ///     Match(10),
    ///     Diff(5),
    ///     Equal(3),
    ///     HardClip(2),
    /// ]);
    /// assert_eq!(cigar, expected_cigar);
    /// ```
    fn try_from(bytes: &[u8]) -> Result<Self> {
        let mut inner = Vec::new();
        let mut i = 0;
        let text_len = bytes.len();
        while i < text_len {
            let mut j = i;
            while j < text_len && bytes[j].is_ascii_digit() {
                j += 1;
            }
            // check that length is provided
            if i == j {
                return Err(Error::BamParseCigar {
                    msg: "Expected length before cigar operation [0-9]+[MIDNSHP=X]".to_owned(),
                });
            }
            // get the length of the operation
            let s = str::from_utf8(&bytes[i..j]).map_err(|_| Error::BamParseCigar {
                msg: format!("Invalid utf-8 bytes '{:?}'.", &bytes[i..j]),
            })?;
            let n = s.parse().map_err(|_| Error::BamParseCigar {
                msg: format!("Unable to parse &str '{:?}' to u32.", s),
            })?;
            // get the operation
            let op = &bytes[j];
            inner.push(match op {
                b'M' => Cigar::Match(n),
                b'I' => Cigar::Ins(n),
                b'D' => Cigar::Del(n),
                b'N' => Cigar::RefSkip(n),
                b'H' => {
                    if i == 0 || j + 1 == text_len {
                        Cigar::HardClip(n)
                    } else {
                        return Err(Error::BamParseCigar {
                            msg: "Hard clipping ('H') is only valid at the start or end of a cigar."
                                .to_owned(),
                        });
                    }
                }
                b'S' => {
                    if i == 0
                        || j + 1 == text_len
                        || bytes[i-1] == b'H'
                        || bytes[j+1..].iter().all(|c| c.is_ascii_digit() || *c == b'H') {
                        Cigar::SoftClip(n)
                    } else {
                        return Err(Error::BamParseCigar {
                        msg: "Soft clips ('S') can only have hard clips ('H') between them and the end of the CIGAR string."
                            .to_owned(),
                        });
                    }
                },
                b'P' => Cigar::Pad(n),
                b'=' => Cigar::Equal(n),
                b'X' => Cigar::Diff(n),
                op => {
                    return Err(Error::BamParseCigar {
                        msg: format!("Expected cigar operation [MIDNSHP=X] but got [{}]", op),
                    })
                }
            });
            i = j + 1;
        }
        Ok(CigarString(inner))
    }
}

impl TryFrom<&str> for CigarString {
    type Error = Error;

    /// Create a CigarString from given &str.
    /// # Example
    /// ```
    /// use rust_htslib::bam::record::*;
    /// use rust_htslib::bam::record::CigarString;
    /// use rust_htslib::bam::record::Cigar::*;
    /// use std::convert::TryFrom;
    ///
    /// let cigar_str = "2H10M5X3=2H";
    /// let cigar = CigarString::try_from(cigar_str)
    ///     .expect("Unable to parse cigar string.");
    /// let expected_cigar = CigarString(vec![
    ///     HardClip(2),
    ///     Match(10),
    ///     Diff(5),
    ///     Equal(3),
    ///     HardClip(2),
    /// ]);
    /// assert_eq!(cigar, expected_cigar);
    /// ```
    fn try_from(text: &str) -> Result<Self> {
        let bytes = text.as_bytes();
        if text.chars().count() != bytes.len() {
            return Err(Error::BamParseCigar {
                msg: "CIGAR string contained non-ASCII characters, which are not valid. Valid are [0-9MIDNSHP=X].".to_owned(),
            });
        }
        CigarString::try_from(bytes)
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
        self.0.iter()
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

// Get number of leading/trailing softclips if a CigarString taking hardclips into account
fn calc_softclips<'a>(mut cigar: impl DoubleEndedIterator<Item = &'a Cigar>) -> i64 {
    match (cigar.next(), cigar.next()) {
        (Some(Cigar::HardClip(_)), Some(Cigar::SoftClip(s))) | (Some(Cigar::SoftClip(s)), _) => {
            *s as i64
        }
        _ => 0,
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

    /// Get the start position of the alignment (0-based).
    pub fn pos(&self) -> i64 {
        self.pos
    }

    /// Get number of bases softclipped at the beginning of the alignment.
    pub fn leading_softclips(&self) -> i64 {
        calc_softclips(self.iter())
    }

    /// Get number of bases softclipped at the end of the alignment.
    pub fn trailing_softclips(&self) -> i64 {
        calc_softclips(self.iter().rev())
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
                Cigar::Del(l) => {
                    // METHOD: leading deletions can happen in case of trimmed reads where
                    // a primer has been removed AFTER read mapping.
                    // Example: 24M8I8D18M9S before trimming, 32H8D18M9S after trimming
                    // with fgbio. While leading deletions should be impossible with
                    // normal read mapping, they make perfect sense with primer trimming
                    // because the mapper still had the evidence to decide in favor of
                    // the deletion via the primer sequence.
                    rpos += l;
                },
                Cigar::RefSkip(_) => {
                    return Err(Error::BamUnexpectedCigarOperation {
                        msg: "'reference skip' (N) found before any operation describing read sequence".to_owned()
                    });
                },
                Cigar::HardClip(_) if i > 0 && i < self.len()-1 => {
                    return Err(Error::BamUnexpectedCigarOperation{
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
                    return Err(Error::BamUnexpectedCigarOperation{
                        msg: "'hard clip' (H) found in between operations, contradicting SAMv1 spec that hard clips can only be at the ends of reads".to_owned()
                    });
                }
                Cigar::HardClip(_) => return Ok(None),
            }
        }

        Ok(None)
    }

    /// transfer ownership of the Cigar out of the CigarView
    pub fn take(self) -> CigarString {
        self.inner
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

pub struct BaseModificationMetadata {
    pub strand: i32,
    pub implicit: i32,
    pub canonical: u8,
}

/// struct containing the internal state required to access
/// the base modifications for a bam::Record
pub struct BaseModificationState<'a> {
    record: &'a Record,
    state: *mut htslib::hts_base_mod_state,
    buffer: Vec<htslib::hts_base_mod>,
    buffer_pos: i32,
}

impl BaseModificationState<'_> {
    /// Initialize a new BaseModification struct from a bam::Record
    /// This function allocates memory for the state structure
    /// and initializes the iterator to the start of the modification
    /// records.
    fn new(r: &Record) -> Result<BaseModificationState<'_>> {
        let mut bm = unsafe {
            BaseModificationState {
                record: r,
                state: hts_sys::hts_base_mod_state_alloc(),
                buffer: Vec::new(),
                buffer_pos: -1,
            }
        };

        if bm.state.is_null() {
            panic!("Unable to allocate memory for hts_base_mod_state");
        }

        // parse the MM tag to initialize the state
        unsafe {
            let ret = hts_sys::bam_parse_basemod(bm.record.inner_ptr(), bm.state);
            if ret != 0 {
                return Err(Error::BamBaseModificationTagNotFound);
            }
        }

        let types = bm.recorded();
        bm.buffer.reserve(types.len());
        Ok(bm)
    }

    pub fn buffer_next_mods(&mut self) -> Result<usize> {
        unsafe {
            let ret = hts_sys::bam_next_basemod(
                self.record.inner_ptr(),
                self.state,
                self.buffer.as_mut_ptr(),
                self.buffer.capacity() as i32,
                &mut self.buffer_pos,
            );

            if ret < 0 {
                return Err(Error::BamBaseModificationIterationFailed);
            }

            // the htslib API won't write more than buffer.capacity() mods to the output array but it will
            // return the actual number of modifications found. We return an error to the caller
            // in the case where there was insufficient storage to return all mods.
            if ret as usize > self.buffer.capacity() {
                return Err(Error::BamBaseModificationTooManyMods);
            }

            // we read the modifications directly into the vector, which does
            // not update the length so needs to be manually set
            self.buffer.set_len(ret as usize);

            Ok(ret as usize)
        }
    }

    /// Return an array containing the modification codes listed for this record.
    /// Positive values are ascii character codes (eg m), negative values are chEBI codes.
    pub fn recorded<'a>(&self) -> &'a [i32] {
        unsafe {
            let mut n: i32 = 0;
            let data_ptr: *const i32 = hts_sys::bam_mods_recorded(self.state, &mut n);

            // htslib should not return a null pointer, even when there are no base mods
            if data_ptr.is_null() {
                panic!("Unable to obtain pointer to base modifications");
            }
            assert!(n >= 0);
            slice::from_raw_parts(data_ptr, n as usize)
        }
    }

    /// Return metadata for the specified character code indicating the strand
    /// the base modification was called on, whether the tag uses implicit mode
    /// and the ascii code for the canonical base.
    /// If there are multiple modifications with the same code this will return the data
    /// for the first mod.  See https://github.com/samtools/htslib/issues/1635
    pub fn query_type(&self, code: i32) -> Result<BaseModificationMetadata> {
        unsafe {
            let mut strand: i32 = 0;
            let mut implicit: i32 = 0;
            // This may be i8 or u8 in hts_sys.
            let mut canonical: c_char = 0;

            let ret = hts_sys::bam_mods_query_type(
                self.state,
                code,
                &mut strand,
                &mut implicit,
                &mut canonical,
            );
            if ret == -1 {
                Err(Error::BamBaseModificationTypeNotFound)
            } else {
                Ok(BaseModificationMetadata {
                    strand,
                    implicit,
                    canonical: canonical.try_into().unwrap(),
                })
            }
        }
    }
}

impl Drop for BaseModificationState<'_> {
    fn drop<'a>(&mut self) {
        unsafe {
            hts_sys::hts_base_mod_state_free(self.state);
        }
    }
}

/// Iterator over the base modifications that returns
/// a vector for all of the mods at each position
pub struct BaseModificationsPositionIter<'a> {
    mod_state: BaseModificationState<'a>,
}

impl BaseModificationsPositionIter<'_> {
    fn new(r: &Record) -> Result<BaseModificationsPositionIter<'_>> {
        let state = BaseModificationState::new(r)?;
        Ok(BaseModificationsPositionIter { mod_state: state })
    }

    pub fn recorded<'a>(&self) -> &'a [i32] {
        self.mod_state.recorded()
    }

    pub fn query_type(&self, code: i32) -> Result<BaseModificationMetadata> {
        self.mod_state.query_type(code)
    }
}

impl Iterator for BaseModificationsPositionIter<'_> {
    type Item = Result<(i32, Vec<hts_sys::hts_base_mod>)>;

    fn next(&mut self) -> Option<Self::Item> {
        let ret = self.mod_state.buffer_next_mods();

        // Three possible things happened in buffer_next_mods:
        // 1. the htslib API call was successful but there are no more mods
        // 2. ths htslib API call was successful and we read some mods
        // 3. the htslib API call failed, we propogate the error wrapped in an option
        match ret {
            Ok(num_mods) => {
                if num_mods == 0 {
                    None
                } else {
                    let data = (self.mod_state.buffer_pos, self.mod_state.buffer.clone());
                    Some(Ok(data))
                }
            }
            Err(e) => Some(Err(e)),
        }
    }
}

/// Iterator over the base modifications that returns
/// the next modification found, one by one
pub struct BaseModificationsIter<'a> {
    mod_state: BaseModificationState<'a>,
    buffer_idx: usize,
}

impl BaseModificationsIter<'_> {
    fn new(r: &Record) -> Result<BaseModificationsIter<'_>> {
        let state = BaseModificationState::new(r)?;
        Ok(BaseModificationsIter {
            mod_state: state,
            buffer_idx: 0,
        })
    }

    pub fn recorded<'a>(&self) -> &'a [i32] {
        self.mod_state.recorded()
    }

    pub fn query_type(&self, code: i32) -> Result<BaseModificationMetadata> {
        self.mod_state.query_type(code)
    }
}

impl Iterator for BaseModificationsIter<'_> {
    type Item = Result<(i32, hts_sys::hts_base_mod)>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.buffer_idx == self.mod_state.buffer.len() {
            // need to use the internal state to read the next
            // set of modifications into the buffer
            let ret = self.mod_state.buffer_next_mods();

            match ret {
                Ok(num_mods) => {
                    if num_mods == 0 {
                        // done iterating
                        return None;
                    } else {
                        // we read some mods, reset the position in the buffer then fall through
                        self.buffer_idx = 0;
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        // if we got here when there are mods buffered that we haven't emitted yet
        assert!(self.buffer_idx < self.mod_state.buffer.len());
        let data = (
            self.mod_state.buffer_pos,
            self.mod_state.buffer[self.buffer_idx],
        );
        self.buffer_idx += 1;
        Some(Ok(data))
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
    fn test_cigar_string_view_pos() {
        let cigar = CigarString(vec![Cigar::Match(100), Cigar::SoftClip(10)]).into_view(5);
        assert_eq!(cigar.pos(), 5);
    }

    #[test]
    fn test_cigar_string_leading_softclips() {
        let cigar = CigarString(vec![Cigar::SoftClip(10), Cigar::Match(100)]).into_view(0);
        assert_eq!(cigar.leading_softclips(), 10);
        let cigar2 = CigarString(vec![
            Cigar::HardClip(5),
            Cigar::SoftClip(10),
            Cigar::Match(100),
        ])
        .into_view(0);
        assert_eq!(cigar2.leading_softclips(), 10);
    }

    #[test]
    fn test_cigar_string_trailing_softclips() {
        let cigar = CigarString(vec![Cigar::Match(100), Cigar::SoftClip(10)]).into_view(0);
        assert_eq!(cigar.trailing_softclips(), 10);
        let cigar2 = CigarString(vec![
            Cigar::Match(100),
            Cigar::SoftClip(10),
            Cigar::HardClip(5),
        ])
        .into_view(0);
        assert_eq!(cigar2.trailing_softclips(), 10);
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
    use crate::bam::{Read, Reader};
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

    #[test]
    fn test_read_orientation_f1r2() {
        let mut bam = Reader::from_path("test/test_paired.sam").unwrap();

        for res in bam.records() {
            let record = res.unwrap();
            assert_eq!(
                record.read_pair_orientation(),
                SequenceReadPairOrientation::F1R2
            );
        }
    }

    #[test]
    fn test_read_orientation_f2r1() {
        let mut bam = Reader::from_path("test/test_nonstandard_orientation.sam").unwrap();

        for res in bam.records() {
            let record = res.unwrap();
            assert_eq!(
                record.read_pair_orientation(),
                SequenceReadPairOrientation::F2R1
            );
        }
    }

    #[test]
    fn test_read_orientation_supplementary() {
        let mut bam = Reader::from_path("test/test_orientation_supplementary.sam").unwrap();

        for res in bam.records() {
            let record = res.unwrap();
            assert_eq!(
                record.read_pair_orientation(),
                SequenceReadPairOrientation::F2R1
            );
        }
    }

    #[test]
    pub fn test_cigar_parsing_non_ascii_error() {
        let cigar_str = "43ጷ";
        let expected_error = Err(Error::BamParseCigar {
                msg: "CIGAR string contained non-ASCII characters, which are not valid. Valid are [0-9MIDNSHP=X].".to_owned(),
            });

        let result = CigarString::try_from(cigar_str);
        assert_eq!(expected_error, result);
    }

    #[test]
    pub fn test_cigar_parsing() {
        // parsing test cases
        let cigar_strs = [
            "1H10M4D100I300N1102=10P25X11S", // test every cigar opt
            "100M",                          // test a single op
            "",                              // test empty input
            "1H1=1H",                        // test simple hardclip
            "1S1=1S",                        // test simple softclip
            "11H11S11=11S11H",               // test complex softclip
            "10H",
            "10S",
        ];
        // expected results
        let cigars = [
            CigarString(vec![
                Cigar::HardClip(1),
                Cigar::Match(10),
                Cigar::Del(4),
                Cigar::Ins(100),
                Cigar::RefSkip(300),
                Cigar::Equal(1102),
                Cigar::Pad(10),
                Cigar::Diff(25),
                Cigar::SoftClip(11),
            ]),
            CigarString(vec![Cigar::Match(100)]),
            CigarString(vec![]),
            CigarString(vec![
                Cigar::HardClip(1),
                Cigar::Equal(1),
                Cigar::HardClip(1),
            ]),
            CigarString(vec![
                Cigar::SoftClip(1),
                Cigar::Equal(1),
                Cigar::SoftClip(1),
            ]),
            CigarString(vec![
                Cigar::HardClip(11),
                Cigar::SoftClip(11),
                Cigar::Equal(11),
                Cigar::SoftClip(11),
                Cigar::HardClip(11),
            ]),
            CigarString(vec![Cigar::HardClip(10)]),
            CigarString(vec![Cigar::SoftClip(10)]),
        ];
        // compare
        for (&cigar_str, truth) in cigar_strs.iter().zip(cigars.iter()) {
            let cigar_parse = CigarString::try_from(cigar_str)
                .unwrap_or_else(|_| panic!("Unable to parse cigar: {}", cigar_str));
            assert_eq!(&cigar_parse, truth);
        }
    }
}

#[cfg(test)]
mod basemod_tests {
    use crate::bam::{Read, Reader};

    #[test]
    pub fn test_count_recorded() {
        let mut bam = Reader::from_path("test/base_mods/MM-double.sam").unwrap();

        for r in bam.records() {
            let record = r.unwrap();
            if let Ok(mods) = record.basemods_iter() {
                let n = mods.recorded().len();
                assert_eq!(n, 3);
            };
        }
    }

    #[test]
    pub fn test_query_type() {
        let mut bam = Reader::from_path("test/base_mods/MM-orient.sam").unwrap();

        let mut n_fwd = 0;
        let mut n_rev = 0;

        for r in bam.records() {
            let record = r.unwrap();
            if let Ok(mods) = record.basemods_iter() {
                for mod_code in mods.recorded() {
                    if let Ok(mod_metadata) = mods.query_type(*mod_code) {
                        if mod_metadata.strand == 0 {
                            n_fwd += 1;
                        }
                        if mod_metadata.strand == 1 {
                            n_rev += 1;
                        }
                    }
                }
            };
        }
        assert_eq!(n_fwd, 2);
        assert_eq!(n_rev, 2);
    }

    #[test]
    pub fn test_mod_iter() {
        let mut bam = Reader::from_path("test/base_mods/MM-double.sam").unwrap();
        let expected_positions = [1, 7, 12, 13, 13, 22, 30, 31];
        let mut i = 0;

        for r in bam.records() {
            let record = r.unwrap();
            for res in record.basemods_iter().unwrap().flatten() {
                let (position, _m) = res;
                assert_eq!(position, expected_positions[i]);
                i += 1;
            }
        }
    }

    #[test]
    pub fn test_position_iter() {
        let mut bam = Reader::from_path("test/base_mods/MM-double.sam").unwrap();
        let expected_positions = [1, 7, 12, 13, 22, 30, 31];
        let expected_counts = [1, 1, 1, 2, 1, 1, 1];
        let mut i = 0;

        for r in bam.records() {
            let record = r.unwrap();
            for res in record.basemods_position_iter().unwrap().flatten() {
                let (position, elements) = res;
                assert_eq!(position, expected_positions[i]);
                assert_eq!(elements.len(), expected_counts[i]);
                i += 1;
            }
        }
    }
}
