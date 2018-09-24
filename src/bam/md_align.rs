use std::fmt::Write;
use std::str;

use bio_types::alignment::{Alignment, AlignmentMode, AlignmentOperation};

use bam;
use bam::record::{Cigar, CigarString};

/// Alignment position on a BAM record, reconstructed using read
/// sequence information along with the Cigar string and the MD Aux
/// field.
pub struct MDAlignPos<'a> {
    record: &'a bam::record::Record,
    pos: CigarMDPos,
}

impl<'a> MDAlignPos<'a> {
    fn new(record: &'a bam::record::Record, pos: CigarMDPos) -> Self {
        MDAlignPos {
            record: record,
            pos: pos,
        }
    }

    fn read_nt_at(&self, at: u32) -> u8 {
        self.record.seq()[at as usize]
    }

    fn qual_at(&self, at: u32) -> u8 {
        self.record.qual()[at as usize]
    }

    /// Read nucleotide sequence as reported in the BAM
    /// alignment. _Note_ that this is the complement of the read
    /// sequence itself for reverse-strand alignments.
    ///
    /// `None` is returned for read deletions.
    pub fn read_nt(&self) -> Option<u8> {
        self.pos.read_seq_pos().map(|pos| self.read_nt_at(pos))
    }

    /// Reference nucleotide sequence as reported in the BAM alignment.
    ///
    /// `None` is returned for read insertions and
    /// soft clipped positions.
    pub fn ref_nt(&self) -> Option<u8> {
        self.pos.ref_nt().or_else(|| {
            if self.ref_pos().is_some() {
                self.read_nt()
            } else {
                None
            }
        })
    }

    /// Read nucleotide quality score as reported in the BAM
    /// alignment.
    ///
    /// `None` is returned at read deletion positions.
    pub fn read_qual(&self) -> Option<u8> {
        self.pos.read_seq_pos().map(|pos| self.qual_at(pos))
    }

    /// Zero-based offset in the read sequence as reported in the BAM
    /// alignment. _Note_ that these positions will run from the last
    /// to the first base in the original input sequence for
    /// reverse-strand alignments.
    ///
    /// `None is returned for read deletions.
    pub fn read_seq_pos(&self) -> Option<u32> {
        self.pos.read_seq_pos()
    }

    /// Zero-based offset in the read sequence as reported in the BAM
    /// alignment. For read deletions, the offset of the next aligned
    /// read sequence nucleotide is reported. _Note_ that these
    /// positions will run from the last to the first base in the
    /// original input sequence for reverse-strand alignments.
    pub fn read_seq_pos_or_next(&self) -> u32 {
        self.pos.read_seq_pos_or_next()
    }

    /// Zero-based offset in the reference sequence.
    ///
    /// `None` is returned for read insertions and
    /// soft clipped positions.
    pub fn ref_pos(&self) -> Option<u32> {
        self.pos.ref_pos()
    }

    /// Zero-based offset in the reference sequence. For read
    /// insertions, the offset of the next aligned reference sequence
    /// nucleotide is reported.
    ///
    /// `None` is returned for soft clipped positions.
    pub fn ref_pos_or_next(&self) -> Option<u32> {
        self.pos.ref_pos_or_next()
    }

    /// Zero-based offset for this position in the original input read
    /// sequence. Position 0 is the first nucleotide from the input
    /// read sequence for either forward or reverse strand alignments.
    ///
    /// `None` is returned for read deletions.
    pub fn read_true_pos(&self) -> Option<u32> {
        self.read_seq_pos().map(|seq_pos| {
            if self.record.is_reverse() {
                self.record.qual().len() as u32 - (1 + seq_pos)
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
    pub fn read_true_nt(&self) -> Option<u8> {
        self.read_nt().map(|nt| {
            if self.record.is_reverse() {
                fast_compl(nt)
            } else {
                nt
            }
        })
    }

    /// Nucleotide for the reference sequence that is matched against
    /// `read_true_nt()`. _Note_ that this is the complement of the
    /// nucleotide in the reference sequence (as per `ref_nt`) for
    /// reverse-strand alignments.
    ///
    /// `None` is returned for read insertions and soft-clipped
    /// positions.
    pub fn ref_true_nt(&self) -> Option<u8> {
        self.ref_nt().map(|nt| {
            if self.record.is_reverse() {
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
    pub fn read_line_char(&self) -> char {
        if self.ref_pos_or_next().is_none() {
            self.read_nt()
                .map_or('-', |nt| nt.to_ascii_lowercase() as char)
        } else {
            self.read_nt().map_or('-', |nt| nt as char)
        }
    }

    /// Character for middle match line in pretty-printed alignment
    /// format. This is a vertical bar at match positions and a space
    /// otherwise.
    pub fn match_line_char(&self) -> char {
        let read = self.read_nt();
        if read.is_some() && read == self.ref_nt() {
            '|'
        } else {
            ' '
        }
    }

    /// Character for reference sequence line in pretty-printed
    /// alignment format. This character is the reference sequence
    /// base when present, a `-` for read insertions, and a space for
    /// soft-clipped positions.
    pub fn ref_line_char(&self) -> char {
        self.ref_nt().map_or_else(
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

/// Iterator that generates `MDAlignPos` alignment positions for a BAM
/// record.
///
/// Note that these positions are generated in reference sequence
/// order. For a reverse-strand alignment, they run from the last to
/// the first sequenced base.
pub struct MDAlignPosIter<'a> {
    align_iter: CigarMDPosIter,
    record: &'a bam::record::Record,
}

impl<'a> MDAlignPosIter<'a> {
    /// Create a new iterator for a BAM record.
    ///
    /// # Arguments
    ///
    /// * `record` is the BAM record whose alignment will be extracted
    ///
    /// # Errors
    ///
    /// An error variant is returned when `record` has no MD aux
    /// field, when the aux field cannot be parsed, or when there is a
    /// conflict between the MD, Cigar, and read sequence.
    pub fn new(record: &'a bam::record::Record) -> Result<Self, MDAlignError> {
        let align_iter = CigarMDPosIter::new_from_record(record)?;
        Ok(MDAlignPosIter {
            align_iter: align_iter,
            record: record,
        })
    }
}

impl<'a> Iterator for MDAlignPosIter<'a> {
    type Item = Result<MDAlignPos<'a>, MDAlignError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.align_iter
            .next()
            .map(|res| res.map(|mdpos| MDAlignPos::new(self.record, mdpos)))
    }
}

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

    /// Reference nucleotide, _when it differs from the read nucleotide_.
    ///
    /// `None` is returned for matches, read insertions and soft
    /// clipped positions.
    pub fn ref_nt(&self) -> Option<u8> {
        match self {
            CigarMDPos::Match {
                read_seq_pos,
                ref_pos,
            } => None,
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
}

/// Iterator over the `CigarMDPos` positions represented by a BAM
/// record.
///
/// Note that these positions are generated in reference sequence
/// order. For a reverse-strand alignment, they run from the last to
/// the first sequenced base.
#[derive(Debug)]
pub struct CigarMDPosIter {
    md_stack: Vec<MatchDesc>,
    cigar_stack: Vec<Cigar>,
    ref_pos_curr: u32,
    read_pos_curr: u32,
}

impl CigarMDPosIter {
    /// Create a new iterator for a BAM record.
    ///
    /// # Arguments
    ///
    /// * `record` is the BAM record whose alignment will be extracted
    pub fn new_from_record(record: &bam::record::Record) -> Result<Self, MDAlignError> {
        let mut md_stack = MatchDesc::parse_from_record(record)?;
        md_stack.reverse();

        let CigarString(mut cigar_stack) = (*record.cigar()).clone();
        cigar_stack.reverse();

        Ok(CigarMDPosIter {
            md_stack: md_stack,
            cigar_stack: cigar_stack,
            ref_pos_curr: record.pos() as u32,
            read_pos_curr: 0,
        })
    }

    // Utility function that yields the next CigarMDPos.
    // Requires the cigar stack is non-empty
    // Requires that 0-length matches and non-yielding cigar entries
    //   are cleared from the tops of those respective stacks.
    fn next_with_some(&mut self) -> Result<CigarMDPos, MDAlignError> {
        let pop_cigar;

        let res = match self.cigar_stack.last_mut().unwrap() {
            Cigar::Match(ref mut ciglen) => {
                let pop_md;
                let mmm = match self
                    .md_stack
                    .last_mut()
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
                        pop_md = *mdlen == 0;
                        *ciglen -= 1;
                        pop_cigar = *ciglen == 0;
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
                        pop_md = true;
                        *ciglen -= 1;
                        pop_cigar = *ciglen == 0;
                        pos
                    }
                    _ => {
                        return Err(MDAlignError::MDvsCIGAR);
                    }
                };
                if pop_md {
                    self.md_stack.pop();
                }
                mmm
            }
            Cigar::Equal(ref mut ciglen) => {
                let pop_md;
                let m = match self
                    .md_stack
                    .last_mut()
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
                        pop_md = *mdlen == 0;
                        *ciglen -= 1;
                        pop_cigar = *ciglen == 0;
                        pos
                    }
                    _ => return Err(MDAlignError::MDvsCIGAR),
                };
                if pop_md {
                    self.md_stack.pop();
                }
                m
            }
            Cigar::Diff(ref mut ciglen) => {
                let pop_md;
                let mm = match self
                    .md_stack
                    .last_mut()
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
                        pop_md = true;
                        *ciglen -= 1;
                        pop_cigar = *ciglen == 0;
                        pos
                    }
                    _ => return Err(MDAlignError::MDvsCIGAR),
                };
                if pop_md {
                    self.md_stack.pop();
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
                pop_cigar = *len == 0;
                pos
            }
            Cigar::Del(ref mut len) => {
                let pop_md;

                let del = match self.md_stack.last_mut() {
                    Some(MatchDesc::Deletion(ref mut ref_nts)) => {
                        if ref_nts.len() > 0 {
                            let ref_nt = ref_nts.remove(0);
                            pop_md = ref_nts.is_empty();
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
                pop_cigar = *len == 0;
                if pop_md {
                    self.md_stack.pop();
                }
                del
            }
            Cigar::SoftClip(ref mut len) => {
                let pos = CigarMDPos::SoftClip {
                    read_seq_pos: self.read_pos_curr,
                };
                self.read_pos_curr += 1;
                *len -= 1;
                pop_cigar = *len == 0;
                pos
            }
            Cigar::RefSkip(_) => panic!("RefSkip in next_with_some"),
            Cigar::HardClip(_) => panic!("HardClip in next_with_some"),
            Cigar::Pad(_) => panic!("Pad in next_with_some"),
        };

        if pop_cigar {
            self.cigar_stack.pop();
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
        while curr
            .ok_or_else(|| MDAlignError::EmptyAlign)?
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

impl Iterator for CigarMDPosIter {
    type Item = Result<CigarMDPos, MDAlignError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Process non-yielding cigar entries
        loop {
            let handled = match self.cigar_stack.last() {
                Some(Cigar::RefSkip(len)) => {
                    self.ref_pos_curr += *len;
                    true
                }
                Some(Cigar::HardClip(_len)) => true,
                Some(Cigar::Pad(_len)) => true,
                _ => false,
            };
            if handled {
                self.cigar_stack.pop();
            } else {
                break;
            }
        }

        // Consume 0-length match
        while self.md_stack.last() == Some(&MatchDesc::Matches(0)) {
            self.md_stack.pop();
        }

        if self.cigar_stack.is_empty() {
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

    /// Create a `Vec` of `MatchDesc` entries corresponding to the MD aux field on a record.
    ///
    /// # Arguments
    ///
    /// * `record` is the source of the MD aux field to parse
    ///
    /// # Errors
    ///
    /// If `record` has no MD aux field, or if the value of the aux field is malformed,
    /// an error variant will be returned.
    pub fn parse_from_record(record: &bam::record::Record) -> Result<Vec<MatchDesc>, MDAlignError> {
        Self::parse(record.aux_md().ok_or_else(|| MDAlignError::NoMD)?)
    }

    /// Parses a bytestring as an MD aux field description.
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
    /// use rust_htslib::bam::md_align::{MatchDesc,MDAlignError};
    /// # fn try_main() -> Result<(), MDAlignError> {
    /// assert_eq!(MatchDesc::parse(b"10A5^AC6")?,
    ///            [MatchDesc::Matches(10),
    ///             MatchDesc::Mismatch(b'A'),
    ///             MatchDesc::Matches(5),
    ///             MatchDesc::Deletion(b"AC".to_vec()),
    ///             MatchDesc::Matches(6)]);
    /// # Ok( () )
    /// # }
    /// # fn main() { try_main().unwrap(); }
    /// ```
    pub fn parse(mdstring: &[u8]) -> Result<Vec<MatchDesc>, MDAlignError> {
        MatchDescIter::new(mdstring).collect()
    }

    /// Convert a vector of `MatchDesc` entries into an MD aux field
    /// string. The string will be partly normalized by merging
    /// adjacent matches (but not deletions), which guarantees
    /// unambiguous parsing.
    ///
    /// # Arguments
    ///
    /// * `mds` are a vector of `MatchDesc` entries    
    /// * `strict_match0` controls whether 0-length matches are added
    /// whenever a mismatch or deletion is not followed immediately by a
    /// run of matches. Strictly adding a 0-length match between a mismatch
    /// and a subsequent mismatch, deletion, or end-of-string ensures
    /// that the resulting MD field string matches the regexp
    ///
    /// `[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*`
    ///
    /// _Note_ that 0-length matches are always added when a deletion is
    /// followed by a mismatch, which is required for unambiguous parsing.
    ///
    /// # Examples
    ///
    /// ```
    /// use rust_htslib::bam::md_align::{MatchDesc,MDAlignError};
    /// # fn try_main() -> Result<(), MDAlignError> {
    /// let mdvec = vec![MatchDesc::Matches(10),
    ///                  MatchDesc::Mismatch(b'A'),
    ///                  MatchDesc::Matches(5),
    ///                  MatchDesc::Deletion(b"AC".to_vec()),
    ///                  MatchDesc::Matches(6)];
    /// assert_eq!(MatchDesc::unparse(&mdvec, false), "10A5^AC6");
    /// assert_eq!(MatchDesc::unparse(&mdvec, true),  "10A5^AC6");
    /// let unnorm = vec![MatchDesc::Matches(10),
    ///                   MatchDesc::Mismatch(b'A'),
    ///                   MatchDesc::Mismatch(b'G'),
    ///                   MatchDesc::Matches(2),
    ///                   MatchDesc::Matches(2),
    ///                   MatchDesc::Deletion(b"AC".to_vec()),
    ///                   MatchDesc::Mismatch(b'T'),
    ///                   MatchDesc::Matches(5)];
    /// assert_eq!(MatchDesc::unparse(&unnorm, false), "10AG4^AC0T5");
    /// assert_eq!(MatchDesc::unparse(&unnorm, true), "10A0G4^AC0T5");
    /// # Ok( () )
    /// # }
    /// # fn main() { try_main().unwrap(); }
    /// ```
    pub fn unparse(mds: &[MatchDesc], strict_match0: bool) -> String {
        let mut mdstr = String::new();
        let mut match_accum = 0;

        let mut md_iter = mds.iter().peekable();
        while let Some(md) = md_iter.next() {
            match md {
                MatchDesc::Matches(mlen) => {
                    if md_iter.peek().map_or(false, |next| next.is_matches()) {
                        match_accum += mlen;
                    } else {
                        write!(&mut mdstr, "{}", match_accum + mlen).unwrap();
                        match_accum = 0;
                    }
                }
                MatchDesc::Mismatch(refnt) => {
                    if strict_match0 && md_iter.peek().map_or(true, |next| !next.is_matches()) {
                        write!(&mut mdstr, "{}0", *refnt as char).unwrap();
                    } else {
                        write!(&mut mdstr, "{}", *refnt as char).unwrap();
                    }
                }
                MatchDesc::Deletion(refnts) => {
                    if strict_match0 {
                        if md_iter.peek().map_or(true, |next| !next.is_matches()) {
                            write!(&mut mdstr, "^{}0", str::from_utf8(refnts).unwrap()).unwrap();
                        } else {
                            write!(&mut mdstr, "^{}", str::from_utf8(refnts).unwrap()).unwrap();
                        }
                    } else {
                        if md_iter.peek().map_or(false, |next| next.is_mismatch()) {
                            write!(&mut mdstr, "^{}0", str::from_utf8(refnts).unwrap()).unwrap();
                        } else {
                            write!(&mut mdstr, "^{}", str::from_utf8(refnts).unwrap()).unwrap();
                        }
                    }
                }
            };
        }

        mdstr
    }
}

/// Iterator over `MatchDesc` entries in an MD aux field
pub struct MatchDescIter<'a> {
    md: &'a [u8],
}

impl<'a> MatchDescIter<'a> {
    /// Create a new iterator over an MD aux field string
    ///
    /// # Arguments
    ///
    /// * `md` is a bytestring of an MD aux field entry
    pub fn new(md: &'a [u8]) -> Self {
        MatchDescIter { md: md }
    }

    /// Create a new iterator over the MD aux field from a BAM record
    ///
    /// # Arguments
    ///
    /// * `record` is the source of the MD aux field for the iterator
    ///
    /// # Errors
    ///
    /// If the record does not have a string-valued MD aux field, an
    /// error variant is returned.
    pub fn new_from_record(record: &'a bam::record::Record) -> Result<Self, MDAlignError> {
        Ok(Self::new(
            record.aux_md().ok_or_else(|| MDAlignError::NoMD)?,
        ))
    }

    // Requires self.md is non-empty, guaranteed to yield a MatchDesc
    fn next_nonempty(&mut self) -> Result<MatchDesc, MDAlignError> {
        let ch = self.md.first().unwrap();
        if ch.is_ascii_digit() {
            let endpos = self
                .md
                .iter()
                .position(|&c| !c.is_ascii_digit())
                .unwrap_or(self.md.len());
            let numstr = str::from_utf8(&self.md[0..endpos])?;
            self.md = &self.md[endpos..];
            Ok(MatchDesc::Matches(numstr.parse()?))
        } else if *ch == b'^' {
            let endpos = self
                .md
                .iter()
                .skip(1)
                .position(|&c| !c.is_ascii_uppercase())
                .unwrap_or(self.md.len() - 1) + 1;
            let del = MatchDesc::Deletion(self.md[1..endpos].to_vec());
            self.md = &self.md[endpos..];
            Ok(del)
        } else if ch.is_ascii_uppercase() {
            self.md = &self.md[1..];
            Ok(MatchDesc::Mismatch(*ch))
        } else {
            Err(MDAlignError::BadMD)
        }
    }
}

impl<'a> Iterator for MatchDescIter<'a> {
    type Item = Result<MatchDesc, MDAlignError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.md.is_empty() {
            None
        } else {
            Some(self.next_nonempty())
        }
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

            let rap_iter = MDAlignPosIter::new(&rec)
                .ok()
                .expect("Creating MDAlignPosIter");
            let rap_res: Result<Vec<MDAlignPos>, MDAlignError> = rap_iter.collect();
            let raps = rap_res.ok().expect("Unable to create Vec<MDAlignPos>");
            let rap_read_line: String = raps.iter().map(|rap| rap.read_line_char()).collect();
            assert_eq!(rap_read_line, TEST_READ_LINE[i]);
            let rap_ref_line: String = raps.iter().map(|rap| rap.ref_line_char()).collect();
            assert_eq!(rap_ref_line, TEST_REF_LINE[i]);

            let md_iter = CigarMDPosIter::new_from_record(&rec)
                .ok()
                .expect("Creating CigarMDPosIter");
            let alignment = md_iter
                .collect_into_alignment()
                .ok()
                .expect("collecting into alignment");
            let a_cigar = alignment.cigar(false);
            assert_eq!(a_cigar, TEST_ALIGN_CIGAR[i]);
        }
    }
}
