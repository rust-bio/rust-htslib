use std::fmt;
use std::fmt::Write;
use std::str;
use std::vec::IntoIter;

use bio_types::alignment::{Alignment, AlignmentMode, AlignmentOperation};

use bam;
use bam::record::Cigar;

/// Alignment position based on information in the CIGAR and MD aux
/// field for a BAM record. The `CigarMDPos` entries for a BAM record
/// represent all information from these fields -- together with the
/// sequence itself, these can fully reconstruct the read-vs-reference
/// alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CigarMDPos {
    Match {
        read_seq_pos: u32,
        ref_pos: u32,
    },
    Mismatch {
        ref_nt: u8,
        read_seq_pos: u32,
        ref_pos: u32,
    },
    Insert {
        read_seq_pos: u32,
        ref_pos_next: u32,
    },
    Delete {
        ref_nt: u8,
        read_seq_pos_next: u32,
        ref_pos: u32,
    },
    SoftClip {
        read_seq_pos: u32,
    },
}

#[allow(unused_variables)]
impl CigarMDPos {
    /// Zero-based offset in the read sequence as reported in the BAM
    /// alignment. _Note_ that these positions will run from the last
    /// to the first base in the original input sequence for
    /// reverse-strand alignments.
    ///
    /// `None is returned for read deletions.
    pub fn read_seq_pos(&self) -> Option<u32> {
        match self {
            CigarMDPos::Match {
                read_seq_pos,
                ref_pos,
            } => Some(*read_seq_pos),
            CigarMDPos::Mismatch {
                ref_nt,
                read_seq_pos,
                ref_pos,
            } => Some(*read_seq_pos),
            CigarMDPos::Insert {
                read_seq_pos,
                ref_pos_next,
            } => Some(*read_seq_pos),
            CigarMDPos::Delete {
                ref_nt,
                read_seq_pos_next,
                ref_pos,
            } => None,
            CigarMDPos::SoftClip { read_seq_pos } => Some(*read_seq_pos),
        }
    }

    /// Zero-based offset in the read sequence as reported in the BAM
    /// alignment. For read deletions, the offset of the next aligned
    /// read sequence nucleotide is reported. _Note_ that these
    /// positions will run from the last to the first base in the
    /// original input sequence for reverse-strand alignments.
    pub fn read_seq_pos_or_next(&self) -> u32 {
        match self {
            CigarMDPos::Match {
                read_seq_pos,
                ref_pos,
            } => *read_seq_pos,
            CigarMDPos::Mismatch {
                ref_nt,
                read_seq_pos,
                ref_pos,
            } => *read_seq_pos,
            CigarMDPos::Insert {
                read_seq_pos,
                ref_pos_next,
            } => *read_seq_pos,
            CigarMDPos::Delete {
                ref_nt,
                read_seq_pos_next,
                ref_pos,
            } => *read_seq_pos_next,
            CigarMDPos::SoftClip { read_seq_pos } => *read_seq_pos,
        }
    }

    /// Zero-based offset in the reference sequence.
    ///
    /// `None` is returned for read insertions and
    /// soft clipped positions.
    pub fn ref_pos(&self) -> Option<u32> {
        match self {
            CigarMDPos::Match {
                read_seq_pos,
                ref_pos,
            } => Some(*ref_pos),
            CigarMDPos::Mismatch {
                ref_nt,
                read_seq_pos,
                ref_pos,
            } => Some(*ref_pos),
            CigarMDPos::Insert {
                read_seq_pos,
                ref_pos_next,
            } => None,
            CigarMDPos::Delete {
                ref_nt,
                read_seq_pos_next,
                ref_pos,
            } => Some(*ref_pos),
            CigarMDPos::SoftClip { read_seq_pos } => None,
        }
    }

    /// Zero-based offset in the reference sequence. For read
    /// insertions, the offset of the next aligned reference sequence
    /// nucleotide is reported.
    ///
    /// `None` is returned for soft clipped positions.
    pub fn ref_pos_or_next(&self) -> Option<u32> {
        match self {
            CigarMDPos::Match {
                read_seq_pos,
                ref_pos,
            } => Some(*ref_pos),
            CigarMDPos::Mismatch {
                ref_nt,
                read_seq_pos,
                ref_pos,
            } => Some(*ref_pos),
            CigarMDPos::Insert {
                read_seq_pos,
                ref_pos_next,
            } => Some(*ref_pos_next),
            CigarMDPos::Delete {
                ref_nt,
                read_seq_pos_next,
                ref_pos,
            } => Some(*ref_pos),
            CigarMDPos::SoftClip { read_seq_pos } => None,
        }
    }

    /// Read nucleotide sequence as reported in the BAM
    /// alignment. _Note_ that this is the complement of the read
    /// sequence itself for reverse-strand alignments.
    ///
    /// `None` is returned for read deletions.
    ///
    /// # Arguments
    ///
    /// `record` is the associated BAM alignment
    pub fn read_nt(&self, record: &bam::record::Record) -> Option<u8> {
        self.read_seq_pos().map(|pos| record.seq()[pos as usize])
    }

    /// Read nucleotide quality score as reported in the BAM
    /// alignment.
    ///
    /// `None` is returned at read deletion positions.
    ///
    /// # Arguments
    ///
    /// `record` is the associated BAM alignment
    pub fn read_qual(&self, record: &bam::record::Record) -> Option<u8> {
        self.read_seq_pos().map(|pos| record.qual()[pos as usize])
    }

    /// Reference nucleotide sequence as reported in the BAM alignment.
    ///
    /// `None` is returned for read insertions and soft clipped
    /// positions.
    ///
    /// # Arguments
    ///
    /// `record` is the associated BAM alignment
    pub fn ref_nt(&self, record: &bam::record::Record) -> Option<u8> {
        match self {
            CigarMDPos::Match {
                read_seq_pos,
                ref_pos,
            } => self.read_nt(record),
            CigarMDPos::Mismatch {
                ref_nt,
                read_seq_pos,
                ref_pos,
            } => Some(*ref_nt),
            CigarMDPos::Insert {
                read_seq_pos,
                ref_pos_next,
            } => None,
            CigarMDPos::Delete {
                ref_nt,
                read_seq_pos_next,
                ref_pos,
            } => Some(*ref_nt),
            CigarMDPos::SoftClip { read_seq_pos } => None,
        }
    }

    /// Zero-based offset for this position in the original input read
    /// sequence. Position 0 is the first nucleotide from the input
    /// read sequence for either forward or reverse strand alignments.
    ///
    /// `None` is returned for read deletions.
    pub fn read_pos_on_read(&self, record: &bam::record::Record) -> Option<u32> {
        self.read_seq_pos().map(|seq_pos| {
            if record.is_reverse() {
                record.qual().len() as u32 - (1 + seq_pos)
            } else {
                seq_pos
            }
        })
    }

    /// Nucleotide in the original input read sequence, taking into
    /// account the fact that reverse-strand alignments report the
    /// reverse-complement of the input read sequence.
    ///
    /// `None` is returned for read deletions.
    pub fn read_nt_on_read(&self, record: &bam::record::Record) -> Option<u8> {
        self.read_nt(record).map(|nt| {
            if record.is_reverse() {
                fast_compl(nt)
            } else {
                nt
            }
        })
    }

    /// Nucleotide for the reference sequence that is matched against
    /// `read_onread_nt()`. _Note_ that this is the complement of the
    /// nucleotide in the reference sequence (as per `ref_nt`) for
    /// reverse-strand alignments.
    ///
    /// `None` is returned for read insertions and soft-clipped
    /// positions.
    pub fn ref_nt_on_read(&self, record: &bam::record::Record) -> Option<u8> {
        self.ref_nt(record).map(|nt| {
            if record.is_reverse() {
                fast_compl(nt)
            } else {
                nt
            }
        })
    }

    /// Character for read sequence line in pretty-printed alignment
    /// format. This character is the read sequence base, in
    /// lower-case for soft-clipped positions, or a `-` for read
    /// deletions.
    pub fn read_line_char(&self, record: &bam::record::Record) -> char {
        if self.ref_pos_or_next().is_none() {
            self.read_nt(record)
                .map_or('-', |nt| nt.to_ascii_lowercase() as char)
        } else {
            self.read_nt(record).map_or('-', |nt| nt as char)
        }
    }

    /// Character for middle match line in pretty-printed alignment
    /// format. This is a vertical bar at match positions and a space
    /// otherwise.
    pub fn match_line_char(&self, record: &bam::record::Record) -> char {
        let read = self.read_nt(record);
        if read.is_some() && read == self.ref_nt(record) {
            '|'
        } else {
            ' '
        }
    }

    /// Character for reference sequence line in pretty-printed
    /// alignment format. This character is the reference sequence
    /// base when present, a `-` for read insertions, and a space for
    /// soft-clipped positions.
    pub fn ref_line_char(&self, record: &bam::record::Record) -> char {
        self.ref_nt(record).map_or_else(
            || {
                if self.ref_pos_or_next().is_none() {
                    ' '
                } else {
                    '-'
                }
            },
            |nt| nt as char,
        )
    }
}

/// Iterator over the `CigarMDPos` positions represented by a BAM
/// record.
///
/// Note that these positions are generated in reference sequence
/// order. For a reverse-strand alignment, they run from the last to
/// the first sequenced base.
#[derive(Debug)]
pub struct CigarMDIter<I, J> {
    md_iter: I,
    md_curr: Option<MatchDesc>,
    cigar_iter: J,
    cigar_curr: Option<Cigar>,
    ref_pos_curr: u32,
    read_pos_curr: u32,
}

impl CigarMDIter<IntoIter<MatchDesc>, IntoIter<Cigar>> {
    /// Create a new iterator for a BAM record.
    ///
    /// # Arguments
    ///
    /// * `record` is the BAM record whose alignment will be extracted
    ///
    /// # Errors An error variant is returned when `record` has no MD
    /// aux field or when there's an error extracting this field.
    pub fn new_from_record(record: &bam::record::Record) -> Result<Self, MDAlignError> {
        Ok(Self::new(
            MDString::new_from_record(record)?,
            (*record.cigar()).to_vec(),
            record.pos() as u32,
        ))
    }
}

impl<I: Iterator<Item = MatchDesc>, J: Iterator<Item = Cigar>> CigarMDIter<I, J> {
    /// Create a new iterator that consumes `MatchDesc` and `Cigar`
    /// entries from iterators in order to yield `CigarMDPos` entries
    /// describing the alignment.
    ///
    /// # Arguments
    ///
    /// `mdstring` can be converted into an iterator over `MatchDesc`
    /// entries
    ///
    /// `cigarstring` can be converted into an iterator over `Cigar`
    /// entires
    ///
    /// `startpos` is the starting position of the alignment on the
    /// reference sequence
    pub fn new<S, T>(mdstring: S, cigarstring: T, startpos: u32) -> Self
    where
        S: IntoIterator<IntoIter = I, Item = MatchDesc>,
        T: IntoIterator<IntoIter = J, Item = Cigar>,
    {
        let mut md_iter = mdstring.into_iter();
        let md_curr = md_iter.next();

        let mut cigar_iter = cigarstring.into_iter();
        let cigar_curr = cigar_iter.next();

        CigarMDIter {
            md_iter: md_iter,
            md_curr: md_curr,
            cigar_iter: cigar_iter,
            cigar_curr: cigar_curr,
            ref_pos_curr: startpos,
            read_pos_curr: 0,
        }
    }

    // Utility function that yields the next CigarMDPos.
    // Requires the cigar stack is non-empty
    // Requires that 0-length matches and non-yielding cigar entries
    //   are cleared from the tops of those respective stacks.
    fn next_with_some(&mut self) -> Result<CigarMDPos, MDAlignError> {
        let next_cigar;

        let res = match self.cigar_curr.as_mut().unwrap() {
            Cigar::Match(ref mut ciglen) => {
                let next_md;
                let mmm = match self.md_curr
                    .as_mut()
                    .ok_or_else(|| MDAlignError::MDvsCIGAR)?
                {
                    MatchDesc::Matches(ref mut mdlen) => {
                        let pos = CigarMDPos::Match {
                            read_seq_pos: self.read_pos_curr,
                            ref_pos: self.ref_pos_curr,
                        };
                        self.read_pos_curr += 1;
                        self.ref_pos_curr += 1;
                        *mdlen -= 1;
                        next_md = *mdlen == 0;
                        *ciglen -= 1;
                        next_cigar = *ciglen == 0;
                        pos
                    }
                    MatchDesc::Mismatch(ref_nt) => {
                        let pos = CigarMDPos::Mismatch {
                            ref_nt: *ref_nt,
                            read_seq_pos: self.read_pos_curr,
                            ref_pos: self.ref_pos_curr,
                        };
                        self.read_pos_curr += 1;
                        self.ref_pos_curr += 1;
                        next_md = true;
                        *ciglen -= 1;
                        next_cigar = *ciglen == 0;
                        pos
                    }
                    _ => {
                        return Err(MDAlignError::MDvsCIGAR);
                    }
                };
                if next_md {
                    self.md_curr = self.md_iter.next();
                }
                mmm
            }
            Cigar::Equal(ref mut ciglen) => {
                let next_md;
                let m = match self.md_curr
                    .as_mut()
                    .ok_or_else(|| MDAlignError::MDvsCIGAR)?
                {
                    MatchDesc::Matches(ref mut mdlen) => {
                        let pos = CigarMDPos::Match {
                            read_seq_pos: self.read_pos_curr,
                            ref_pos: self.ref_pos_curr,
                        };
                        self.read_pos_curr += 1;
                        self.ref_pos_curr += 1;

                        *mdlen -= 1;
                        next_md = *mdlen == 0;
                        *ciglen -= 1;
                        next_cigar = *ciglen == 0;
                        pos
                    }
                    _ => return Err(MDAlignError::MDvsCIGAR),
                };
                if next_md {
                    self.md_curr = self.md_iter.next();
                }
                m
            }
            Cigar::Diff(ref mut ciglen) => {
                let next_md;
                let mm = match self.md_curr
                    .as_mut()
                    .ok_or_else(|| MDAlignError::MDvsCIGAR)?
                {
                    MatchDesc::Mismatch(ref_nt) => {
                        let pos = CigarMDPos::Mismatch {
                            ref_nt: *ref_nt,
                            read_seq_pos: self.read_pos_curr,
                            ref_pos: self.ref_pos_curr,
                        };
                        self.read_pos_curr += 1;
                        self.ref_pos_curr += 1;
                        next_md = true;
                        *ciglen -= 1;
                        next_cigar = *ciglen == 0;
                        pos
                    }
                    _ => return Err(MDAlignError::MDvsCIGAR),
                };
                if next_md {
                    self.md_curr = self.md_iter.next();
                }
                mm
            }
            Cigar::Ins(ref mut len) => {
                let mut pos = CigarMDPos::Insert {
                    read_seq_pos: self.read_pos_curr,
                    ref_pos_next: self.ref_pos_curr,
                };
                self.read_pos_curr += 1;
                *len -= 1;
                next_cigar = *len == 0;
                pos
            }
            Cigar::Del(ref mut len) => {
                let next_md;

                let del = match self.md_curr.as_mut() {
                    Some(MatchDesc::Deletion(ref mut ref_nts)) => {
                        if ref_nts.len() > 0 {
                            let ref_nt = ref_nts.remove(0);
                            next_md = ref_nts.is_empty();
                            let pos = CigarMDPos::Delete {
                                ref_nt: ref_nt,
                                ref_pos: self.ref_pos_curr,
                                read_seq_pos_next: self.read_pos_curr,
                            };
                            self.ref_pos_curr += 1;
                            pos
                        } else {
                            return Err(MDAlignError::MDvsCIGAR);
                        }
                    }
                    _ => return Err(MDAlignError::MDvsCIGAR),
                };
                *len -= 1;
                next_cigar = *len == 0;
                if next_md {
                    self.md_curr = self.md_iter.next();
                }
                del
            }
            Cigar::SoftClip(ref mut len) => {
                let pos = CigarMDPos::SoftClip {
                    read_seq_pos: self.read_pos_curr,
                };
                self.read_pos_curr += 1;
                *len -= 1;
                next_cigar = *len == 0;
                pos
            }
            Cigar::RefSkip(_) => panic!("RefSkip in next_with_some"),
            Cigar::HardClip(_) => panic!("HardClip in next_with_some"),
            Cigar::Pad(_) => panic!("Pad in next_with_some"),
        };

        if next_cigar {
            self.cigar_curr = self.cigar_iter.next();
        }

        Ok(res)
    }

    /// Collect the alignment positions from this record into an
    /// `Alignment`. This is always a semi-global alignment that
    /// includes the soft clipping from the read as `Xclip` operations.
    pub fn collect_into_alignment(self) -> Result<Alignment, MDAlignError> {
        let mds_result: Result<Vec<CigarMDPos>, MDAlignError> = self.collect();
        let mds = mds_result?;
        let mut iter = mds.as_slice().into_iter();

        let mut ops = Vec::new();

        let mut curr = iter.next();

        let mut left_clip = 0;
        while curr.ok_or_else(|| MDAlignError::EmptyAlign)?
            .ref_pos_or_next()
            .is_none()
        {
            left_clip += 1;
            curr = iter.next();
        }
        ops.push(AlignmentOperation::Xclip(left_clip));

        let (xstart, ystart) = if let Some(md) = curr {
            (md.read_seq_pos_or_next(), md.ref_pos_or_next().unwrap())
        } else {
            panic!("No CigarMDPos remaining after left clip")
        };

        let mut xend = xstart;
        let mut yend = ystart;

        while let Some(md) = curr {
            match md {
                CigarMDPos::Match {
                    read_seq_pos,
                    ref_pos,
                } => {
                    ops.push(AlignmentOperation::Match);
                    xend = *read_seq_pos;
                    yend = *ref_pos;
                }
                CigarMDPos::Mismatch {
                    ref_nt: _,
                    read_seq_pos,
                    ref_pos,
                } => {
                    ops.push(AlignmentOperation::Subst);
                    xend = *read_seq_pos;
                    yend = *ref_pos;
                }
                CigarMDPos::Insert {
                    read_seq_pos,
                    ref_pos_next: _,
                } => {
                    ops.push(AlignmentOperation::Ins);
                    xend = *read_seq_pos;
                }
                CigarMDPos::Delete {
                    ref_nt: _,
                    read_seq_pos_next: _,
                    ref_pos,
                } => {
                    ops.push(AlignmentOperation::Del);
                    yend = *ref_pos;
                }
                CigarMDPos::SoftClip { read_seq_pos: _ } => break,
            }
            curr = iter.next();
        }

        let mut right_clip = 0;
        let mut xlen = xend;
        while let Some(md) = curr {
            match md {
                CigarMDPos::SoftClip { read_seq_pos } => {
                    xlen = *read_seq_pos;
                }
                _ => {
                    return Err(MDAlignError::BadMD);
                }
            };
            right_clip += 1;
            curr = iter.next();
        }
        ops.push(AlignmentOperation::Xclip(right_clip));

        Ok(Alignment {
            score: 0,
            xstart: xstart as usize,
            ystart: ystart as usize,
            xend: xend as usize,
            yend: yend as usize,
            ylen: (1 + yend - ystart) as usize,
            xlen: xlen as usize,
            operations: ops,
            mode: AlignmentMode::Semiglobal,
        })
    }
}

impl<I: Iterator<Item = MatchDesc>, J: Iterator<Item = Cigar>> Iterator for CigarMDIter<I, J> {
    type Item = Result<CigarMDPos, MDAlignError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Process non-yielding cigar entries
        loop {
            let handled = match self.cigar_curr {
                Some(Cigar::RefSkip(len)) => {
                    self.ref_pos_curr += len;
                    true
                }
                Some(Cigar::HardClip(_len)) => true,
                Some(Cigar::Pad(_len)) => true,
                _ => false,
            };
            if handled {
                self.cigar_curr = self.cigar_iter.next();
            } else {
                break;
            }
        }

        // Consume 0-length match
        while self.md_curr == Some(MatchDesc::Matches(0)) {
            self.md_curr = self.md_iter.next();
        }

        if self.cigar_curr.is_none() {
            None
        } else {
            Some(
                self.next_with_some()
                    .map_err(|e| ::std::convert::From::from(e)),
            )
        }
    }
}

/// Abstract representation of entries in an MD aux field
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MatchDesc {
    Matches(u32),
    Mismatch(u8),
    Deletion(Vec<u8>),
}

impl MatchDesc {
    /// Create a new, matching MD entry.
    ///
    /// # Arguments
    ///
    /// * `len` is the length of the match
    pub fn new_matches(len: u32) -> Self {
        MatchDesc::Matches(len)
    }

    /// Create a new mismatched MD entry.
    ///
    /// # Arguments
    ///
    /// * `refnt` is the reference base at the mismatched position
    ///
    /// # Errors
    ///
    /// If `refnt` is not an upper-case ASCII character, an error
    /// variant is returned.
    pub fn new_mismatch(refnt: u8) -> Result<Self, MDAlignError> {
        if refnt.is_ascii_uppercase() {
            Ok(MatchDesc::Mismatch(refnt))
        } else {
            Err(MDAlignError::BadMD)
        }
    }

    /// Create a new read-deletion MD entry
    ///
    /// # Arguments
    ///
    /// * `refnts` are the reference bases at the deletion position
    ///
    /// # Errors
    ///
    /// If `refnts` is empty or contains invalid (non upper-case
    /// ASCII) characters, an error variant is returned.
    pub fn new_deletion(refnts: &[u8]) -> Result<Self, MDAlignError> {
        if refnts.len() > 0 && refnts.iter().all(|&c| c.is_ascii_uppercase()) {
            Ok(MatchDesc::Deletion(refnts.to_vec()))
        } else {
            Err(MDAlignError::BadMD)
        }
    }

    /// True if this `MatchDesc` represents a run of matches
    pub fn is_matches(&self) -> bool {
        match self {
            MatchDesc::Matches(_) => true,
            _ => false,
        }
    }

    /// True if this `MatchDesc` represents a mismatch
    pub fn is_mismatch(&self) -> bool {
        match self {
            MatchDesc::Mismatch(_) => true,
            _ => false,
        }
    }

    /// True if this `MatchDesc` represents a read deletion relative to the reference
    pub fn is_deletion(&self) -> bool {
        match self {
            MatchDesc::Deletion(_) => true,
            _ => false,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MDString(pub Vec<MatchDesc>);

impl MDString {
    /// Create an `MDString` corresponding to the MD aux field on a record.
    ///
    /// # Arguments
    ///
    /// * `record` is the source of the MD aux field to parse
    ///
    /// # Errors
    ///
    /// If `record` has no MD aux field, or if the value of the aux field is malformed,
    /// an error variant will be returned.
    pub fn new_from_record(record: &bam::record::Record) -> Result<Self, MDAlignError> {
        Self::new(record.aux_md().ok_or_else(|| MDAlignError::NoMD)?)
    }

    /// Create an `MDString` by parseing a bytestring as an MD aux field description.
    ///
    /// # Arguments
    ///
    /// * `mdstring` is a bytestring MD aux field
    ///
    /// # Errors
    ///
    /// If the `mdstring` is not a well-formed MD aux field, an error variant is returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use rust_htslib::bam::md_align::{MatchDesc,MDString,MDAlignError};
    /// # fn try_main() -> Result<(), MDAlignError> {
    /// assert_eq!(MDString::new(b"10A5^AC6")?,
    ///            MDString(vec![MatchDesc::Matches(10),
    ///             MatchDesc::Mismatch(b'A'),
    ///             MatchDesc::Matches(5),
    ///             MatchDesc::Deletion(b"AC".to_vec()),
    ///             MatchDesc::Matches(6)]));
    /// # Ok( () )
    /// # }
    /// # fn main() { try_main().unwrap(); }
    /// ```
    // This can fail, to it can't be From<&[u8]>. Consider TryFrom<&[u8]> when stabilized
    pub fn new(mdstring: &[u8]) -> Result<Self, MDAlignError> {
        let mut mdvec = Vec::new();

        let mut mdrest = mdstring;

        while !mdrest.is_empty() {
            let ch = mdrest.first().unwrap();
            if ch.is_ascii_digit() {
                let endpos = mdrest
                    .iter()
                    .position(|&c| !c.is_ascii_digit())
                    .unwrap_or(mdrest.len());
                let numstr = str::from_utf8(&mdrest[0..endpos])?;
                mdrest = &mdrest[endpos..];
                mdvec.push(MatchDesc::Matches(numstr.parse()?))
            } else if *ch == b'^' {
                let endpos = mdrest
                    .iter()
                    .skip(1)
                    .position(|&c| !c.is_ascii_uppercase())
                    .unwrap_or(mdrest.len() - 1) + 1;
                mdvec.push(MatchDesc::Deletion(mdrest[1..endpos].to_vec()));
                mdrest = &mdrest[endpos..];
            } else if ch.is_ascii_uppercase() {
                mdrest = &mdrest[1..];
                mdvec.push(MatchDesc::Mismatch(*ch));
            } else {
                return Err(MDAlignError::BadMD);
            }
        }

        Ok(MDString(mdvec))
    }
}

impl fmt::Display for MDString {
    /// Convert a vector of `MatchDesc` entries into an MD aux field
    /// string. The string will be partly normalized by merging
    /// adjacent matches (but not deletions), which guarantees
    /// unambiguous parsing. Zero-length matches will be added between
    /// a mismatch and a subsequent mismatch, deletion, or
    /// end-of-string. This ensures that the resulting MD field string
    /// matches the regexp
    ///
    /// `[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*`
    ///
    /// # Arguments
    ///
    /// * `mds` are a vector of `MatchDesc` entries    
    ///
    /// # Examples
    ///
    /// ```
    /// use rust_htslib::bam::md_align::{MatchDesc,MDString,MDAlignError};
    /// # fn try_main() -> Result<(), MDAlignError> {
    /// let mdvec = vec![MatchDesc::Matches(10),
    ///                  MatchDesc::Mismatch(b'A'),
    ///                  MatchDesc::Matches(5),
    ///                  MatchDesc::Deletion(b"AC".to_vec()),
    ///                  MatchDesc::Matches(6)];
    /// let mdstr = MDString(mdvec);
    /// assert_eq!(mdstr.to_string(), "10A5^AC6");
    /// let unnorm_vec = vec![MatchDesc::Matches(10),
    ///                       MatchDesc::Mismatch(b'A'),
    ///                       MatchDesc::Mismatch(b'G'),
    ///                       MatchDesc::Matches(2),
    ///                       MatchDesc::Matches(2),
    ///                       MatchDesc::Deletion(b"AC".to_vec()),
    ///                       MatchDesc::Mismatch(b'T'),
    ///                       MatchDesc::Matches(5)];
    /// let unnorm_str = MDString(unnorm_vec);
    /// assert_eq!(unnorm_str.to_string(), "10A0G4^AC0T5");
    /// # Ok( () )
    /// # }
    /// # fn main() { try_main().unwrap(); }
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        let mut match_accum = 0;

        let mut md_iter = self.0.iter().peekable();
        while let Some(md) = md_iter.next() {
            match md {
                MatchDesc::Matches(mlen) => {
                    if md_iter.peek().map_or(false, |next| next.is_matches()) {
                        match_accum += mlen;
                    } else {
                        write!(f, "{}", match_accum + mlen)?;
                        match_accum = 0;
                    }
                }
                MatchDesc::Mismatch(refnt) => {
                    if md_iter.peek().map_or(true, |next| !next.is_matches()) {
                        write!(f, "{}0", *refnt as char)?;
                    } else {
                        write!(f, "{}", *refnt as char)?;
                    }
                }
                MatchDesc::Deletion(refnts) => {
                    if md_iter.peek().map_or(true, |next| !next.is_matches()) {
                        write!(f, "^{}0", str::from_utf8(refnts).unwrap())?;
                    } else {
                        write!(f, "^{}", str::from_utf8(refnts).unwrap())?;
                    }
                }
            };
        }

        Ok(())
    }
}

impl<'a> IntoIterator for &'a MDString {
    type Item = &'a MatchDesc;
    type IntoIter = ::std::slice::Iter<'a, MatchDesc>;

    fn into_iter(self) -> Self::IntoIter {
        (&self.0).into_iter()
    }
}

impl IntoIterator for MDString {
    type Item = MatchDesc;
    type IntoIter = ::std::vec::IntoIter<MatchDesc>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

quick_error! {
    #[derive(Debug,Clone)]
    pub enum MDAlignError {
        NoMD {
            description("no MD aux field")
        }
        BadMD {
            description("bad MD value")
        }
        MDvsCIGAR {
            description("MD inconsistent with CIGAR")
        }
        BadSeqLen {
            description("Sequence/quality length inconsistent with MD/CIGAR")
        }
        EmptyAlign {
            description("Alignment has no positions")
        }
        ParseInt(err: ::std::num::ParseIntError) {
            from()
        }
        Utf8(err: ::std::str::Utf8Error) {
            from()
        }
    }
}

fn fast_compl(nt: u8) -> u8 {
    unsafe { *COMPL.get_unchecked((nt as usize) & 0x7f) }
}

const NT_NONE: u8 = b'*';
#[cfg_attr(rustfmt, rustfmt_skip)]
static COMPL: [u8; 128] =
    [ NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE,
      NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE,
// 0x10
      NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE,
      NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, 
// 0x20
      NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE,
      NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, 
// 0x30
      NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE,
      NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, 
// 0x40
      NT_NONE, b'T',    NT_NONE, b'G',    NT_NONE, NT_NONE, NT_NONE, b'C',
      NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, b'N',    NT_NONE, 
// 0x50
      NT_NONE, NT_NONE, NT_NONE, NT_NONE, b'A',    NT_NONE, NT_NONE, NT_NONE,
      NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, 
// 0x60
      NT_NONE, b'T',    NT_NONE, b'G',    NT_NONE, NT_NONE, NT_NONE, b'C',
      NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, b'N',    NT_NONE, 
// 0x70
      NT_NONE, NT_NONE, NT_NONE, NT_NONE, b'A',    NT_NONE, NT_NONE, NT_NONE,
      NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE, NT_NONE];

#[cfg(test)]
mod tests {
    use super::super::*;
    use super::*;
    use std::path::Path;

    const N_TEST: usize = 14;

    static TEST_READ_NAMES: [&[u8]; N_TEST] = [
        b"K00364:89:HWTF2BBXX:2:1102:3376:38293",
        b"K00364:89:HWTF2BBXX:2:1102:3092:38293",
        b"K00364:89:HWTF2BBXX:2:1102:3356:38293",
        b"K00364:89:HWTF2BBXX:2:1102:3640:38293",
        b"K00364:89:HWTF2BBXX:2:1102:20994:38293",
        b"K00364:89:HWTF2BBXX:2:1102:18852:38310",
        b"K00364:89:HWTF2BBXX:2:1102:7395:38363",
        b"K00364:89:HWTF2BBXX:2:1102:12560:38310",
        b"K00364:89:HWTF2BBXX:2:1102:13017:38328",
        b"K00364:89:HWTF2BBXX:2:1102:31324:38293",
        b"K00364:89:HWTF2BBXX:2:1102:27630:38293",
        b"K00364:89:HWTF2BBXX:2:1102:15240:38310",
        b"K00364:89:HWTF2BBXX:2:1102:3062:38310",
        b"K00364:89:HWTF2BBXX:2:1102:30371:38293",
    ];

    static TEST_CIGAR: [&str; N_TEST] = [
        "151M",
        "151M",
        "151M",
        "151M",
        "149M2S",
        "149M2S",
        "2S149M",
        "2S149M",
        "44M2D107M",
        "133M2D11M7S",
        "47M1I103M",
        "103M1I47M",
        "53M99N98M",
        "48M1998N103M",
    ];

    static TEST_MD: [&[u8]; N_TEST] = [
        b"151",
        b"151",
        b"41T87T21",
        b"137A13",
        b"96T52",
        b"93A0A54",
        b"127T21",
        b"149",
        b"44^GT107",
        b"8A68A26A15A12^TG1T9",
        b"76T73",
        b"150",
        b"151",
        b"75A75",
    ];

    static TEST_ALIGN_CIGAR: [&str; N_TEST] = [
        "151=",
        "151=",
        "41=1X87=1X21=",
        "137=1X13=",
        "96=1X52=2S",
        "93=2X54=2S",
        "2S127=1X21=",
        "2S149=",
        "44=2D107=",
        "8=1X68=1X26=1X15=1X12=2D1=1X9=7S",
        "47=1I29=1X73=",
        "103=1I47=",
        "151=",
        "75=1X75=",
    ];

    static TEST_READ_LINE: [&str; N_TEST] =
        ["TATCTCTGCTATGGTCGATGCTGCTAAGGATCAAGTTGCCCTAAGAATGCCAGAATTGATTCCAGTCTTGTCTGAAACCATGTGGGACACCAAGAAGGAAGTCAAGGCTGCTGCTACTGCCGCCATGACCAAGGCTACCGAAACTGTTGAC",
         "CATAGCAGTGGTTTCAGAAATGTTGGCAGACATTCTGTGGAAAACAGTGAAGTCACCGTTACCCAAGGTGTGGTGCAACAACAATTGCTTAGCTTGAGCAGAGATGGATGGGACACCAACAACGTGCAAAACACCGACGTGTTCAGCGTAA",
         "CGTTGACTATGATCAAGTACCAGACGAAATTCAATCAAGAACTGCGTTATGCTTAATGCATAATATGAAGAGAATTATGGACTCTAGCTTTAAATTAGATGAATTGACTGTACAAGATGATTTAATATTCGATGGTTCCTTGATAACTGAG",
         "CTGCAAAGCTATATAATGCCACGGCCTTTGGAGAAGACGTGGATGGGGATGTTGGAGCAGTGAATCTGCCCTGACTAGTTTGCGGTGTTGAAGCGGACGAGATTCTAGATTTAGAAAATCTGTTCGACAAGTCATCGGCTTCACGGTCTTT",
         "ACGACACAGGAAGACATGTCATTGTTACAACGGGTGTGGGGCAACATCAAATGTGGGCTGCTCAACACTGGACATGGAGAAATCCACATACTTTCACCACATCAGGTGGTTTAGGTACGATGGGTTACGGTCTCCCTGCCGCCATCGGTtc",
         "TGCGGAACAGTTTTATTATTGGTACCACCCGTACTGGATATTGGTACGTTTGTATGATTAGTCTCATTGTCACTGTACGAGTCTGAGTGTCTGGGATCTTTAGATTTACTGGCGTGCGACGACTCATGTGTGTTAGATTGGGACATGGGgg",
         "tcCCGAACCAGCTTTACAAAATAAGAATTGTATAGAAGCGACAGTAATGCAGTCCAAGGAACGCCCTAATGACAAGATCATACCAACCAAAACCGAAAAAAACGATTTTGGAATAGGCACTCAATGGTTCGAACGCAAACAAATATCAAGA",
         "tcTCGGTTGGTGAAATCAAGTACTTAAAAGACCTAACAAACTCTACAATTTGCGCCTTTCTTTCCCCTAATAGAACATGTTAGATTTATAAGCGAAAACGGAAAAGAAAACTACAGCCCAGCAATAAGAGTGCAAATCAAAGTGTGCATGA",
         "TTCCTCTACTTCTTACGCTATCAAGAAGAAGGATGAATTGGAAC--GTTGCCAAGTCTAACCGTTAAGAAGCTAAAAAAAGTGAAAGATTTTCAATATTACATATGTTTTTCTTATTATTACTTTTATCTTCTTTGTACACTCTATTGTTTAA",
         "GCCCAAATGAAATGTTTTCCTTTACAGTACCATTCATTATCCATGGAACTTGTGAAACATAAGCAACAGAACCATGAGCGGTGGCGAAACCTTTAACCCTGAATGGATCACCTAACATGCGTGACAATAGAGC--TATTACCACTGagatcgg",
         "GCTGGTGTTAAGGCTGAATGATTCGATCAAAAGCGAAGTATTTTTTTCTTCCTAGTTTTTATGTAAAAATATATTCACTTTTTTCTCTACGCAGTTGCTGCAGTGTTATAAGAAAAATAGAATGAAATAAGATAAAGCTAGACGAATTGTG",
         "CGTCAGCATCGTGGTAGTCAAAAAGTTCATCTGCACCGTACTCTTTCAACAATTTTTCATGTTTACGAGAAGCAACGACGATGATCTTGCTGAAACCGTTTAGTTTTTTTTGCCAATTGAATAAGCATCTGGCCAACAGCAGTGGCACCAC",
         "ACAAACCATGAATCCTTTGTAGGATCAATAAACGATCCTGAAGATGGACCTGCTGTTGCACATACTGTTTATTTGGCTGCCTTGGTATACCTCGTGTTTTTCGTATTCTGTGGGTTCCAAGTTTACCTAGCCAGAAGAAAACCTTCGATCG",
         "CGACATATGGAGATACTTTATTTCCTTTTCTTAATTATTAACGTATACCTATAAATTAACAAAGTATCTAAACAAGATACATAAGTGTACTCAAACTGAGTAGAATCGTCGATTAAACTTCCTTCTCCTTTTAAAAATTAAAAACAGCAAA"];

    static TEST_REF_LINE: [&str; N_TEST] =
        ["TATCTCTGCTATGGTCGATGCTGCTAAGGATCAAGTTGCCCTAAGAATGCCAGAATTGATTCCAGTCTTGTCTGAAACCATGTGGGACACCAAGAAGGAAGTCAAGGCTGCTGCTACTGCCGCCATGACCAAGGCTACCGAAACTGTTGAC",
         "CATAGCAGTGGTTTCAGAAATGTTGGCAGACATTCTGTGGAAAACAGTGAAGTCACCGTTACCCAAGGTGTGGTGCAACAACAATTGCTTAGCTTGAGCAGAGATGGATGGGACACCAACAACGTGCAAAACACCGACGTGTTCAGCGTAA",
         "CGTTGACTATGATCAAGTACCAGACGAAATTCAATCAAGAATTGCGTTATGCTTAATGCATAATATGAAGAGAATTATGGACTCTAGCTTTAAATTAGATGAATTGACTGTACAAGATGATTTAATATTTGATGGTTCCTTGATAACTGAG",
         "CTGCAAAGCTATATAATGCCACGGCCTTTGGAGAAGACGTGGATGGGGATGTTGGAGCAGTGAATCTGCCCTGACTAGTTTGCGGTGTTGAAGCGGACGAGATTCTAGATTTAGAAAATCTGTTCGACAAGTCATCGACTTCACGGTCTTT",
         "ACGACACAGGAAGACATGTCATTGTTACAACGGGTGTGGGGCAACATCAAATGTGGGCTGCTCAACACTGGACATGGAGAAATCCACATACTTTCATCACATCAGGTGGTTTAGGTACGATGGGTTACGGTCTCCCTGCCGCCATCGGT  ",
         "TGCGGAACAGTTTTATTATTGGTACCACCCGTACTGGATATTGGTACGTTTGTATGATTAGTCTCATTGTCACTGTACGAGTCTGAGTGTCTGAAATCTTTAGATTTACTGGCGTGCGACGACTCATGTGTGTTAGATTGGGACATGGG  ",
         "  CCGAACCAGCTTTACAAAATAAGAATTGTATAGAAGCGACAGTAATGCAGTCCAAGGAACGCCCTAATGACAAGATCATACCAACCAAAACCGAAAAAAACGATTTTGGAATAGGCACTCAATGGTTTGAACGCAAACAAATATCAAGA",
         "  TCGGTTGGTGAAATCAAGTACTTAAAAGACCTAACAAACTCTACAATTTGCGCCTTTCTTTCCCCTAATAGAACATGTTAGATTTATAAGCGAAAACGGAAAAGAAAACTACAGCCCAGCAATAAGAGTGCAAATCAAAGTGTGCATGA",
         "TTCCTCTACTTCTTACGCTATCAAGAAGAAGGATGAATTGGAACGTGTTGCCAAGTCTAACCGTTAAGAAGCTAAAAAAAGTGAAAGATTTTCAATATTACATATGTTTTTCTTATTATTACTTTTATCTTCTTTGTACACTCTATTGTTTAA",
         "GCCCAAATAAAATGTTTTCCTTTACAGTACCATTCATTATCCATGGAACTTGTGAAACATAAGCAACAGAACCATGAACGGTGGCGAAACCTTTAACCCTGAATAGATCACCTAACATGCATGACAATAGAGCTGTTTTACCACTG       ",
         "GCTGGTGTTAAGGCTGAATGATTCGATCAAAAGCGAAGTATTTTTTT-TTCCTAGTTTTTATGTAAAAATATATTCATTTTTTTCTCTACGCAGTTGCTGCAGTGTTATAAGAAAAATAGAATGAAATAAGATAAAGCTAGACGAATTGTG",
         "CGTCAGCATCGTGGTAGTCAAAAAGTTCATCTGCACCGTACTCTTTCAACAATTTTTCATGTTTACGAGAAGCAACGACGATGATCTTGCTGAAACCGTTTAG-TTTTTTTGCCAATTGAATAAGCATCTGGCCAACAGCAGTGGCACCAC",
         "ACAAACCATGAATCCTTTGTAGGATCAATAAACGATCCTGAAGATGGACCTGCTGTTGCACATACTGTTTATTTGGCTGCCTTGGTATACCTCGTGTTTTTCGTATTCTGTGGGTTCCAAGTTTACCTAGCCAGAAGAAAACCTTCGATCG",
         "CGACATATGGAGATACTTTATTTCCTTTTCTTAATTATTAACGTATACCTATAAATTAACAAAGTATCTAAACAAAATACATAAGTGTACTCAAACTGAGTAGAATCGTCGATTAAACTTCCTTCTCCTTTTAAAAATTAAAAACAGCAAA"];

    #[test]
    fn test_md_align() {
        let mut bam = Reader::from_path(&Path::new("test/test_md_align.bam"))
            .ok()
            .expect("Error opening file.");

        for (i, record) in bam.records().enumerate() {
            let mut rec = record.ok().expect("Error reading BAM record");
            assert_eq!(rec.qname(), TEST_READ_NAMES[i]);

            assert_eq!(rec.aux_md().unwrap(), TEST_MD[i]);
            assert_eq!(rec.cigar().to_string(), TEST_CIGAR[i]);

            let cigar_md_iter = CigarMDIter::new_from_record(&rec)
                .ok()
                .expect("Creating CigarMDIter");
            let cigar_md_res: Result<Vec<CigarMDPos>, MDAlignError> = cigar_md_iter.collect();
            let cigar_md_vec = cigar_md_res.ok().expect("Unable to create Vec<CigarMdPos>");
            let read_line: String = cigar_md_vec
                .iter()
                .map(|rap| rap.read_line_char(&rec))
                .collect();
            assert_eq!(read_line, TEST_READ_LINE[i]);
            let ref_line: String = raps.iter().map(|rap| rap.ref_line_char(&rec)).collect();
            assert_eq!(ref_line, TEST_REF_LINE[i]);

            let md_iter = CigarMDIter::new_from_record(&rec)
                .ok()
                .expect("Creating CigarMDIter");
            let alignment = md_iter
                .collect_into_alignment()
                .ok()
                .expect("collecting into alignment");
            let a_cigar = alignment.cigar(false);
            assert_eq!(a_cigar, TEST_ALIGN_CIGAR[i]);
        }
    }
}
