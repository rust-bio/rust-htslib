//// Copyright 2019 Johannes Köster and Florian Finkernagel.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Extensions for BAM records beyond htslib

use crate::bam;
use crate::bam::record::Cigar;
use crate::htslib;
use std::collections::HashMap;

pub struct IterAlignedBlockPairs {
    genome_pos: i64,
    read_pos: i64,
    cigar_index: usize,
    cigar: Vec<Cigar>,
}

impl Iterator for IterAlignedBlockPairs {
    type Item = ([i64; 2], [i64; 2]);
    fn next(&mut self) -> Option<Self::Item> {
        while self.cigar_index < self.cigar.len() {
            let entry = self.cigar[self.cigar_index];
            match entry {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    let qstart = self.read_pos;
                    let qend = qstart + len as i64;
                    let rstart = self.genome_pos;
                    let rend = self.genome_pos + len as i64;
                    self.read_pos += len as i64;
                    self.genome_pos += len as i64;
                    self.cigar_index += 1;
                    return Some(([qstart, qend], [rstart, rend]));
                }
                Cigar::Ins(len) | Cigar::SoftClip(len) => {
                    self.read_pos += len as i64;
                }
                Cigar::Del(len) | Cigar::RefSkip(len) => {
                    self.genome_pos += len as i64;
                }
                Cigar::HardClip(_) => {} // no advance
                Cigar::Pad(_) => panic!("Padding (Cigar::Pad) is not supported."), //padding is only used for multiple sequence alignment
            }
            self.cigar_index += 1;
        }
        None
    }
}

pub struct IterAlignedBlocks {
    pos: i64,
    cigar_index: usize,
    cigar: Vec<Cigar>,
}

impl Iterator for IterAlignedBlocks {
    type Item = [i64; 2];
    fn next(&mut self) -> Option<Self::Item> {
        while self.cigar_index < self.cigar.len() {
            let entry = self.cigar[self.cigar_index];
            match entry {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    let out_pos = self.pos;
                    //result.push([pos, pos + *len as i64]);
                    self.pos += len as i64;
                    self.cigar_index += 1;
                    return Some([out_pos, out_pos + len as i64]);
                }
                Cigar::Del(len) => self.pos += len as i64,
                Cigar::RefSkip(len) => self.pos += len as i64,
                _ => (),
            }
            self.cigar_index += 1;
        }
        None
    }
}

pub struct IterIntrons {
    pos: i64,
    cigar_index: usize,
    cigar: Vec<Cigar>,
}

impl Iterator for IterIntrons {
    type Item = [i64; 2];
    fn next(&mut self) -> Option<Self::Item> {
        while self.cigar_index < self.cigar.len() {
            let entry = self.cigar[self.cigar_index];
            match entry {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) | Cigar::Del(len) => {
                    self.pos += len as i64
                }
                Cigar::RefSkip(len) => {
                    let junc_start = self.pos;
                    self.pos += len as i64;
                    self.cigar_index += 1;
                    return Some([junc_start, self.pos]); //self.pos is  junc_start + len
                }
                _ => {}
            }
            self.cigar_index += 1;
        }
        None
    }
}

pub struct IterAlignedPairs {
    genome_pos: i64,
    read_pos: i64,
    cigar: Vec<Cigar>,
    remaining_match_bp: u32,
    cigar_index: usize,
}

impl Iterator for IterAlignedPairs {
    type Item = [i64; 2];
    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining_match_bp > 0 {
            self.remaining_match_bp -= 1;
            self.genome_pos += 1;
            self.read_pos += 1;
            return Some([self.read_pos - 1, self.genome_pos - 1]);
        }

        while self.cigar_index < self.cigar.len() {
            let entry = self.cigar[self.cigar_index];
            match entry {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    self.genome_pos += 1;
                    self.read_pos += 1;
                    self.remaining_match_bp = len - 1;
                    self.cigar_index += 1;
                    return Some([self.read_pos - 1, self.genome_pos - 1]);
                }
                Cigar::Ins(len) | Cigar::SoftClip(len) => {
                    self.read_pos += len as i64;
                }
                Cigar::Del(len) | Cigar::RefSkip(len) => {
                    self.genome_pos += len as i64;
                }
                Cigar::HardClip(_) => {} // no advance
                Cigar::Pad(_) => panic!("Padding (Cigar::Pad) is not supported."), //padding is only used for multiple sequence alignment
            }
            self.cigar_index += 1;
        }
        None
    }
}

pub struct IterAlignedPairsFull {
    genome_pos: i64,
    read_pos: i64,
    cigar: Vec<Cigar>,
    remaining_match_bp: u32,
    remaining_ins_bp: u32,
    remaining_del_bp: u32,
    cigar_index: usize,
}

impl Iterator for IterAlignedPairsFull {
    type Item = [Option<i64>; 2];
    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining_match_bp > 0 {
            self.remaining_match_bp -= 1;
            self.genome_pos += 1;
            self.read_pos += 1;
            return Some([Some(self.read_pos - 1), Some(self.genome_pos - 1)]);
        }
        if self.remaining_ins_bp > 0 {
            self.remaining_ins_bp -= 1;
            self.read_pos += 1;
            return Some([Some(self.read_pos - 1), None]);
        }
        if self.remaining_del_bp > 0 {
            self.remaining_del_bp -= 1;
            self.genome_pos += 1;
            return Some([None, Some(self.genome_pos - 1)]);
        }

        while self.cigar_index < self.cigar.len() {
            let entry = self.cigar[self.cigar_index];
            match entry {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    self.genome_pos += 1;
                    self.read_pos += 1;
                    self.remaining_match_bp = len - 1;
                    self.cigar_index += 1;
                    return Some([Some(self.read_pos - 1), Some(self.genome_pos - 1)]);
                }
                Cigar::Ins(len) | Cigar::SoftClip(len) => {
                    self.read_pos += 1;
                    self.remaining_ins_bp = len - 1;
                    self.cigar_index += 1;
                    return Some([Some(self.read_pos - 1), None]);
                }
                Cigar::Del(len) | Cigar::RefSkip(len) => {
                    self.genome_pos += 1;
                    self.remaining_del_bp = len - 1;
                    self.cigar_index += 1;
                    return Some([None, Some(self.genome_pos - 1)]);
                }
                Cigar::HardClip(_) => {
                    // no advance
                }
                Cigar::Pad(_) => panic!("Padding (Cigar::Pad) is not supported."), //padding is only used for multiple sequence alignment
            }
            self.cigar_index += 1;
        }
        None
    }
}

/// Extra functionality for BAM records
///
/// Inspired by pysam
pub trait BamRecordExtensions {
    /// iterator over start and end positions of aligned gapless blocks
    ///
    /// The start and end positions are in genomic coordinates.
    /// There is not necessarily a gap between blocks on the genome,
    /// this happens on insertions.
    ///
    /// pysam: blocks
    /// See also: [aligned_block_pairs](#tymethod.aligned_block_pairs) if you need
    /// the read coordinates as well.
    fn aligned_blocks(&self) -> IterAlignedBlocks;

    ///Iter over <([read_start, read_stop], [genome_start, genome_stop]) blocks
    ///of continously aligned reads.
    ///
    ///In contrast to [aligned_blocks](#tymethod.aligned_blocks), this returns
    ///read and genome coordinates.
    ///In contrast to aligned_pairs, this returns just the start-stop
    ///coordinates of each block.
    ///
    ///There is not necessarily a gap between blocks in either coordinate space
    ///(this happens in in-dels).
    fn aligned_block_pairs(&self) -> IterAlignedBlockPairs;

    /// This scans the CIGAR for reference skips
    /// and reports their positions.
    /// It does not inspect the reported regions
    /// for actual splice sites.
    /// pysam: get_introns
    fn introns(&self) -> IterIntrons;

    /// iter aligned read and reference positions on a basepair level
    ///
    /// No entry for insertions, deletions or skipped pairs
    ///
    /// pysam: get_aligned_pairs(matches_only = True)
    ///
    /// See also [aligned_block_pairs](#tymethod.aligned_block_pairs)
    /// if you just need start&end coordinates of each block.
    /// That way you can allocate less memory for the same
    /// informational content.
    fn aligned_pairs(&self) -> IterAlignedPairs;

    /// iter list of read and reference positions on a basepair level.
    ///
    /// Unlike `aligned_pairs` this returns None in
    /// either the read positions or the reference position
    /// for insertions, deletions or skipped pairs
    ///
    /// pysam: aligned_pairs(matches_only = False)
    fn aligned_pairs_full(&self) -> IterAlignedPairsFull;

    /// the number of nucleotides covered by each Cigar::* variant.
    ///
    /// Result is a Hashmap Cigar::*(0) => covered nucleotides
    ///
    /// pysam: first result from get_cigar_stats
    fn cigar_stats_nucleotides(&self) -> HashMap<Cigar, i32>;

    /// the number of occurrences of each each Cigar::* variant
    ///
    /// Result is a Hashmap Cigar::*(0) => number of times this Cigar::
    /// appeared
    ///
    /// pysam: second result from get_cigar_stats
    fn cigar_stats_blocks(&self) -> HashMap<Cigar, i32>;

    /// iter over  reference positions that this read aligns to
    ///
    /// only returns positions that are aligned, excluding any soft-clipped
    /// or unaligned positions within the read
    ///
    /// pysam: get_reference_positions(full_length=False)
    fn reference_positions(&self) -> Box<dyn Iterator<Item = i64>>;

    ///
    /// iter over reference positions that this read aligns to
    ///
    /// include soft-clipped or skipped positions as None
    ///
    /// pysam: get_reference_positions(full_length=True)
    fn reference_positions_full(&self) -> Box<dyn Iterator<Item = Option<i64>>>;

    /// left most aligned reference position of the read on the reference genome.
    fn reference_start(&self) -> i64;

    /// right most aligned reference position of the read on the reference genome.
    fn reference_end(&self) -> i64;

    /// infer the query length from the cigar string, optionally include hard clipped bases
    ///
    /// Contrast with record::seq_len which returns the length of the sequence stored
    /// in the BAM file, and as such is 0 if the BAM file omits sequences
    ///
    /// pysam: infer_query_length / infer_read_length
    fn seq_len_from_cigar(&self, include_hard_clip: bool) -> usize;
}

impl BamRecordExtensions for bam::Record {
    fn aligned_blocks(&self) -> IterAlignedBlocks {
        IterAlignedBlocks {
            pos: self.pos(),
            cigar: self.cigar().take().0,
            cigar_index: 0,
        }
    }

    fn introns(&self) -> IterIntrons {
        IterIntrons {
            pos: self.pos(),
            cigar: self.cigar().take().0,
            cigar_index: 0,
        }
    }

    fn aligned_block_pairs(&self) -> IterAlignedBlockPairs {
        IterAlignedBlockPairs {
            genome_pos: self.pos(),
            read_pos: 0,
            cigar: self.cigar().take().0,
            cigar_index: 0,
        }
    }

    fn aligned_pairs(&self) -> IterAlignedPairs {
        IterAlignedPairs {
            genome_pos: self.pos(),
            read_pos: 0,
            cigar: self.cigar().take().0,
            remaining_match_bp: 0,
            cigar_index: 0,
        }
    }

    fn aligned_pairs_full(&self) -> IterAlignedPairsFull {
        IterAlignedPairsFull {
            genome_pos: self.pos(),
            read_pos: 0,
            cigar: self.cigar().take().0,
            remaining_match_bp: 0,
            remaining_ins_bp: 0,
            remaining_del_bp: 0,
            cigar_index: 0,
        }
    }

    fn cigar_stats_nucleotides(&self) -> HashMap<Cigar, i32> {
        let mut result = HashMap::new();
        result.insert(Cigar::Match(0), 0); // M
        result.insert(Cigar::Ins(0), 0); // I
        result.insert(Cigar::Del(0), 0); // D
        result.insert(Cigar::RefSkip(0), 0); // N
        result.insert(Cigar::SoftClip(0), 0); // S
        result.insert(Cigar::HardClip(0), 0); // H
        result.insert(Cigar::Pad(0), 0); // P
        result.insert(Cigar::Equal(0), 0); // =
        result.insert(Cigar::Diff(0), 0); // X
        for entry in self.cigar().iter() {
            match entry {
                Cigar::Match(len) => *result.get_mut(&Cigar::Match(0)).unwrap() += *len as i32, // M
                Cigar::Ins(len) => *result.get_mut(&Cigar::Ins(0)).unwrap() += *len as i32,     // I
                Cigar::Del(len) => *result.get_mut(&Cigar::Del(0)).unwrap() += *len as i32,     // D
                Cigar::RefSkip(len) => *result.get_mut(&Cigar::RefSkip(0)).unwrap() += *len as i32, // N
                Cigar::SoftClip(len) => {
                    *result.get_mut(&Cigar::SoftClip(0)).unwrap() += *len as i32
                } // S
                Cigar::HardClip(len) => {
                    *result.get_mut(&Cigar::HardClip(0)).unwrap() += *len as i32
                } // H
                Cigar::Pad(len) => *result.get_mut(&Cigar::Pad(0)).unwrap() += *len as i32, // P
                Cigar::Equal(len) => *result.get_mut(&Cigar::Equal(0)).unwrap() += *len as i32, // =
                Cigar::Diff(len) => *result.get_mut(&Cigar::Diff(0)).unwrap() += *len as i32, // X
            }
        }
        result
    }

    fn cigar_stats_blocks(&self) -> HashMap<Cigar, i32> {
        let mut result = HashMap::new();
        result.insert(Cigar::Match(0), 0); // M
        result.insert(Cigar::Ins(0), 0); // I
        result.insert(Cigar::Del(0), 0); // D
        result.insert(Cigar::RefSkip(0), 0); // N
        result.insert(Cigar::SoftClip(0), 0); // S
        result.insert(Cigar::HardClip(0), 0); // H
        result.insert(Cigar::Pad(0), 0); // P
        result.insert(Cigar::Equal(0), 0); // =
        result.insert(Cigar::Diff(0), 0); // X
        for entry in self.cigar().iter() {
            match entry {
                Cigar::Match(_) => *result.get_mut(&Cigar::Match(0)).unwrap() += 1, // M
                Cigar::Ins(_) => *result.get_mut(&Cigar::Ins(0)).unwrap() += 1,     // I
                Cigar::Del(_) => *result.get_mut(&Cigar::Del(0)).unwrap() += 1,     // D
                Cigar::RefSkip(_) => *result.get_mut(&Cigar::RefSkip(0)).unwrap() += 1, // N
                Cigar::SoftClip(_) => *result.get_mut(&Cigar::SoftClip(0)).unwrap() += 1, // S
                Cigar::HardClip(_) => *result.get_mut(&Cigar::HardClip(0)).unwrap() += 1, // H
                Cigar::Pad(_) => *result.get_mut(&Cigar::Pad(0)).unwrap() += 1,     // P
                Cigar::Equal(_) => *result.get_mut(&Cigar::Equal(0)).unwrap() += 1, // =
                Cigar::Diff(_) => *result.get_mut(&Cigar::Diff(0)).unwrap() += 1,   // X
            }
        }
        result
    }

    fn reference_positions(&self) -> Box<dyn Iterator<Item = i64>> {
        Box::new(self.aligned_pairs().map(|x| x[1]))
    }

    fn reference_positions_full(&self) -> Box<dyn Iterator<Item = Option<i64>>> {
        Box::new(
            self.aligned_pairs_full()
                .filter(|x| x[0].is_some())
                .map(|x| x[1]),
        )
    }

    fn reference_start(&self) -> i64 {
        self.pos()
    }
    fn reference_end(&self) -> i64 {
        unsafe { htslib::bam_endpos(self.inner_ptr()) }
    }

    fn seq_len_from_cigar(&self, include_hard_clip: bool) -> usize {
        let mut result = 0;
        for entry in self.cigar().iter() {
            match entry {
                Cigar::Match(len)
                | Cigar::Ins(len)
                | Cigar::SoftClip(len)
                | Cigar::Equal(len)
                | Cigar::Diff(len) => {
                    result += len;
                }
                Cigar::HardClip(len) => {
                    if include_hard_clip {
                        result += len;
                    }
                }
                _ => {}
            }
        }
        result as usize
    }
}

#[cfg(test)]
mod tests {
    use crate::bam;
    use crate::bam::ext::BamRecordExtensions;
    use crate::bam::record::{Cigar, CigarString};
    use crate::bam::Read;
    use std::collections::HashMap;

    #[test]
    fn spliced_reads() {
        let mut bam = bam::Reader::from_path("./test/test_spliced_reads.bam").unwrap();
        let mut it = bam.records();
        let blocks: Vec<_> = it.next().expect("iter").unwrap().aligned_blocks().collect();
        //6S45M - 0
        assert!(blocks[0] == [16050676, 16050721]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //7M2D44M - 1
        assert!(blocks[0] == [16096878, 16096885]);
        //7M2D44M - 1
        assert!(blocks[1] == [16096887, 16096931]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //29M2D22M - 2
        assert!(blocks[0] == [16097145, 16097174]);
        //29M2D22M - 2
        assert!(blocks[1] == [16097176, 16097198]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //51M - 3
        assert!(blocks[0] == [16117350, 16117401]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //51M - 4
        assert!(blocks[0] == [16118483, 16118534]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //51M - 5
        assert!(blocks[0] == [16118499, 16118550]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //51M - 6
        assert!(blocks[0] == [16118499, 16118550]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //51M - 7
        assert!(blocks[0] == [16118499, 16118550]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //51M - 8
        assert!(blocks[0] == [16123411, 16123462]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //6S45M - 9
        assert!(blocks[0] == [16123417, 16123462]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //41M10S - 10
        assert!(blocks[0] == [16165860, 16165901]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //51M - 11
        assert!(blocks[0] == [16180871, 16180922]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //51M - 12
        assert!(blocks[0] == [16189705, 16189756]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //51M - 13
        assert!(blocks[0] == [16231271, 16231322]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //51M - 14
        assert!(blocks[0] == [16237657, 16237708]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //9S42M - 15
        assert!(blocks[0] == [16255012, 16255054]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //51M - 16
        assert!(blocks[0] == [16255391, 16255442]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //50M1S - 17
        assert!(blocks[0] == [16255392, 16255442]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //45M6S - 18
        assert!(blocks[0] == [16256084, 16256129]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //3S48M - 19
        assert!(blocks[0] == [16256224, 16256272]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //42M9S - 20
        assert!(blocks[0] == [16325199, 16325241]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //13S38M - 21
        assert!(blocks[0] == [16352865, 16352903]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //44M7S - 22
        assert!(blocks[0] == [16352968, 16353012]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //5S46M - 23
        assert!(blocks[0] == [16414998, 16415044]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //23M4I24M - 24
        assert!(blocks[0] == [17031591, 17031614]);
        //23M4I24M - 24
        assert!(blocks[1] == [17031614, 17031638]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //18M1I32M - 25
        assert!(blocks[0] == [17057382, 17057400]);
        //18M1I32M - 25
        assert!(blocks[1] == [17057400, 17057432]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //17M2183N34M - 26
        assert!(blocks[0] == [17092766, 17092783]);
        //17M2183N34M - 26
        assert!(blocks[1] == [17094966, 17095000]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //1M2183N50M - 27
        assert!(blocks[0] == [17092782, 17092783]);
        //1M2183N50M - 27
        assert!(blocks[1] == [17094966, 17095016]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //1M2183N50M - 28
        assert!(blocks[0] == [17092782, 17092783]);
        //1M2183N50M - 28
        assert!(blocks[1] == [17094966, 17095016]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //9S33M9S - 29
        assert!(blocks[0] == [17137287, 17137320]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //2S48M1S - 30
        assert!(blocks[0] == [17306238, 17306286]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //4S45M2S - 31
        assert!(blocks[0] == [17561868, 17561913]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //41M11832N10M - 32
        assert!(blocks[0] == [17566078, 17566119]);
        //41M11832N10M - 32
        assert!(blocks[1] == [17577951, 17577961]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //11M11832N25M710N15M - 33
        assert!(blocks[0] == [17566108, 17566119]);
        //11M11832N25M710N15M - 33
        assert!(blocks[1] == [17577951, 17577976]);
        //11M11832N25M710N15M - 33
        assert!(blocks[2] == [17578686, 17578701]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //8M11832N25M710N18M - 34
        assert!(blocks[0] == [17566111, 17566119]);
        //8M11832N25M710N18M - 34
        assert!(blocks[1] == [17577951, 17577976]);
        //8M11832N25M710N18M - 34
        assert!(blocks[2] == [17578686, 17578704]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //8M11832N25M710N18M - 35
        assert!(blocks[0] == [17566111, 17566119]);
        //8M11832N25M710N18M - 35
        assert!(blocks[1] == [17577951, 17577976]);
        //8M11832N25M710N18M - 35
        assert!(blocks[2] == [17578686, 17578704]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //8M11832N25M710N18M - 36
        assert!(blocks[0] == [17566111, 17566119]);
        //8M11832N25M710N18M - 36
        assert!(blocks[1] == [17577951, 17577976]);
        //8M11832N25M710N18M - 36
        assert!(blocks[2] == [17578686, 17578704]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //8M11832N25M710N18M - 37
        assert!(blocks[0] == [17566111, 17566119]);
        //8M11832N25M710N18M - 37
        assert!(blocks[1] == [17577951, 17577976]);
        //8M11832N25M710N18M - 37
        assert!(blocks[2] == [17578686, 17578704]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //7M11832N25M710N19M - 38
        assert!(blocks[0] == [17566112, 17566119]);
        //7M11832N25M710N19M - 38
        assert!(blocks[1] == [17577951, 17577976]);
        //7M11832N25M710N19M - 38
        assert!(blocks[2] == [17578686, 17578705]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //6M11832N25M710N20M - 39
        assert!(blocks[0] == [17566113, 17566119]);
        //6M11832N25M710N20M - 39
        assert!(blocks[1] == [17577951, 17577976]);
        //6M11832N25M710N20M - 39
        assert!(blocks[2] == [17578686, 17578706]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //6M11832N25M710N20M - 40
        assert!(blocks[0] == [17566113, 17566119]);
        //6M11832N25M710N20M - 40
        assert!(blocks[1] == [17577951, 17577976]);
        //6M11832N25M710N20M - 40
        assert!(blocks[2] == [17578686, 17578706]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //1S44M1467N6M - 41
        assert!(blocks[0] == [17579733, 17579777]);
        //1S44M1467N6M - 41
        assert!(blocks[1] == [17581244, 17581250]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //2M1514N48M95N1M - 42
        assert!(blocks[0] == [17581369, 17581371]);
        //2M1514N48M95N1M - 42
        assert!(blocks[1] == [17582885, 17582933]);
        //2M1514N48M95N1M - 42
        assert!(blocks[2] == [17583028, 17583029]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //1M1514N48M95N2M - 43
        assert!(blocks[0] == [17581370, 17581371]);
        //1M1514N48M95N2M - 43
        assert!(blocks[1] == [17582885, 17582933]);
        //1M1514N48M95N2M - 43
        assert!(blocks[2] == [17583028, 17583030]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //1M1514N48M95N2M - 44
        assert!(blocks[0] == [17581370, 17581371]);
        //1M1514N48M95N2M - 44
        assert!(blocks[1] == [17582885, 17582933]);
        //1M1514N48M95N2M - 44
        assert!(blocks[2] == [17583028, 17583030]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //1S22M95N28M - 45
        assert!(blocks[0] == [17582911, 17582933]);
        //1S22M95N28M - 45
        assert!(blocks[1] == [17583028, 17583056]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //37M538N13M1S - 46
        assert!(blocks[0] == [17588621, 17588658]);
        //37M538N13M1S - 46
        assert!(blocks[1] == [17589196, 17589209]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //37M538N13M1S - 47
        assert!(blocks[0] == [17588621, 17588658]);
        //37M538N13M1S - 47
        assert!(blocks[1] == [17589196, 17589209]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //37M538N13M1S - 48
        assert!(blocks[0] == [17588621, 17588658]);
        //37M538N13M1S - 48
        assert!(blocks[1] == [17589196, 17589209]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //1S25M1D25M - 49
        assert!(blocks[0] == [17591770, 17591795]);
        //1S25M1D25M - 49
        assert!(blocks[1] == [17591796, 17591821]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //24M1D24M3S - 50
        assert!(blocks[0] == [17593855, 17593879]);
        //24M1D24M3S - 50
        assert!(blocks[1] == [17593880, 17593904]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //16M1D28M7S - 51
        assert!(blocks[0] == [17593863, 17593879]);
        //16M1D28M7S - 51
        assert!(blocks[1] == [17593880, 17593908]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //11S7M1I32M - 52
        assert!(blocks[0] == [17596476, 17596483]);
        //11S7M1I32M - 52
        assert!(blocks[1] == [17596483, 17596515]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //5S9M1892N37M - 53
        assert!(blocks[0] == [17624012, 17624021]);
        //5S9M1892N37M - 53
        assert!(blocks[1] == [17625913, 17625950]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //2S9M1892N40M - 54
        assert!(blocks[0] == [17624012, 17624021]);
        //2S9M1892N40M - 54
        assert!(blocks[1] == [17625913, 17625953]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //1S7M3D19M2285N24M - 55
        assert!(blocks[0] == [31796700, 31796707]);
        //1S7M3D19M2285N24M - 55
        assert!(blocks[1] == [31796710, 31796729]);
        //1S7M3D19M2285N24M - 55
        assert!(blocks[2] == [31799014, 31799038]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //14M799N28M13881N7M2S - 56
        assert!(blocks[0] == [36722692, 36722706]);
        //14M799N28M13881N7M2S - 56
        assert!(blocks[1] == [36723505, 36723533]);
        //14M799N28M13881N7M2S - 56
        assert!(blocks[2] == [36737414, 36737421]);

        let blocks: Vec<_> = it.next().unwrap().unwrap().aligned_blocks().collect();
        //4S21M1696N23M2331N3M - 57
        assert!(blocks[0] == [44587963, 44587984]);
        //4S21M1696N23M2331N3M - 57
        assert!(blocks[1] == [44589680, 44589703]);
        //4S21M1696N23M2331N3M - 57
        assert!(blocks[2] == [44592034, 44592037]);
    }

    #[test]
    fn test_introns() {
        let mut bam = bam::Reader::from_path("./test/test_spliced_reads.bam").unwrap();
        let mut it = bam.records();

        //6S45M - 0
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //7M2D44M - 1
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //29M2D22M - 2
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //51M - 3
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //51M - 4
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //51M - 5
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //51M - 6
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //51M - 7
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //51M - 8
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //6S45M - 9
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //41M10S - 10
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //51M - 11
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //51M - 12
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //51M - 13
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //51M - 14
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //9S42M - 15
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //51M - 16
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //50M1S - 17
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //45M6S - 18
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //3S48M - 19
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //42M9S - 20
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //13S38M - 21
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //44M7S - 22
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //5S46M - 23
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //23M4I24M - 24
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //18M1I32M - 25
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //17M2183N34M - 26
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], [17092783, 17094966]);
        //1M2183N50M - 27
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], [17092783, 17094966]);
        //1M2183N50M - 28
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], [17092783, 17094966]);
        //9S33M9S - 29
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //2S48M1S - 30
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //4S45M2S - 31
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //41M11832N10M - 32
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], [17566119, 17577951]);
        //11M11832N25M710N15M - 33
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], [17566119, 17577951]);
        assert_eq!(introns[1], [17577976, 17578686]);
        //8M11832N25M710N18M - 34
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], [17566119, 17577951]);
        assert_eq!(introns[1], [17577976, 17578686]);
        //8M11832N25M710N18M - 35
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], [17566119, 17577951]);
        assert_eq!(introns[1], [17577976, 17578686]);
        //8M11832N25M710N18M - 36
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], [17566119, 17577951]);
        assert_eq!(introns[1], [17577976, 17578686]);
        //8M11832N25M710N18M - 37
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], [17566119, 17577951]);
        assert_eq!(introns[1], [17577976, 17578686]);
        //7M11832N25M710N19M - 38
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], [17566119, 17577951]);
        assert_eq!(introns[1], [17577976, 17578686]);
        //6M11832N25M710N20M - 39
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], [17566119, 17577951]);
        assert_eq!(introns[1], [17577976, 17578686]);
        //6M11832N25M710N20M - 40
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], [17566119, 17577951]);
        assert_eq!(introns[1], [17577976, 17578686]);
        //1S44M1467N6M - 41
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], [17579777, 17581244]);
        //2M1514N48M95N1M - 42
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], [17581371, 17582885]);
        assert_eq!(introns[1], [17582933, 17583028]);
        //1M1514N48M95N2M - 43
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], [17581371, 17582885]);
        assert_eq!(introns[1], [17582933, 17583028]);
        //1M1514N48M95N2M - 44
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], [17581371, 17582885]);
        assert_eq!(introns[1], [17582933, 17583028]);
        //1S22M95N28M - 45
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], [17582933, 17583028]);
        //37M538N13M1S - 46
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], [17588658, 17589196]);
        //37M538N13M1S - 47
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], [17588658, 17589196]);
        //37M538N13M1S - 48
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], [17588658, 17589196]);
        //1S25M1D25M - 49
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //24M1D24M3S - 50
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //16M1D28M7S - 51
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //11S7M1I32M - 52
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 0);
        //5S9M1892N37M - 53
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], [17624021, 17625913]);
        //2S9M1892N40M - 54
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], [17624021, 17625913]);
        //1S7M3D19M2285N24M - 55
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], [31796729, 31799014]);
        //14M799N28M13881N7M2S - 56
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], [36722706, 36723505]);
        assert_eq!(introns[1], [36723533, 36737414]);
        //4S21M1696N23M2331N3M - 57
        let introns: Vec<_> = it.next().unwrap().unwrap().introns().collect();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], [44587984, 44589680]);
        assert_eq!(introns[1], [44589703, 44592034]);
    }

    #[test]
    fn test_aligned_pairs() {
        let mut bam = bam::Reader::from_path("./test/test_spliced_reads.bam").unwrap();
        let mut it = bam.records();

        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [6, 16050676],
                [7, 16050677],
                [8, 16050678],
                [9, 16050679],
                [10, 16050680],
                [11, 16050681],
                [12, 16050682],
                [13, 16050683],
                [14, 16050684],
                [15, 16050685],
                [16, 16050686],
                [17, 16050687],
                [18, 16050688],
                [19, 16050689],
                [20, 16050690],
                [21, 16050691],
                [22, 16050692],
                [23, 16050693],
                [24, 16050694],
                [25, 16050695],
                [26, 16050696],
                [27, 16050697],
                [28, 16050698],
                [29, 16050699],
                [30, 16050700],
                [31, 16050701],
                [32, 16050702],
                [33, 16050703],
                [34, 16050704],
                [35, 16050705],
                [36, 16050706],
                [37, 16050707],
                [38, 16050708],
                [39, 16050709],
                [40, 16050710],
                [41, 16050711],
                [42, 16050712],
                [43, 16050713],
                [44, 16050714],
                [45, 16050715],
                [46, 16050716],
                [47, 16050717],
                [48, 16050718],
                [49, 16050719],
                [50, 16050720]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16096878],
                [1, 16096879],
                [2, 16096880],
                [3, 16096881],
                [4, 16096882],
                [5, 16096883],
                [6, 16096884],
                [7, 16096887],
                [8, 16096888],
                [9, 16096889],
                [10, 16096890],
                [11, 16096891],
                [12, 16096892],
                [13, 16096893],
                [14, 16096894],
                [15, 16096895],
                [16, 16096896],
                [17, 16096897],
                [18, 16096898],
                [19, 16096899],
                [20, 16096900],
                [21, 16096901],
                [22, 16096902],
                [23, 16096903],
                [24, 16096904],
                [25, 16096905],
                [26, 16096906],
                [27, 16096907],
                [28, 16096908],
                [29, 16096909],
                [30, 16096910],
                [31, 16096911],
                [32, 16096912],
                [33, 16096913],
                [34, 16096914],
                [35, 16096915],
                [36, 16096916],
                [37, 16096917],
                [38, 16096918],
                [39, 16096919],
                [40, 16096920],
                [41, 16096921],
                [42, 16096922],
                [43, 16096923],
                [44, 16096924],
                [45, 16096925],
                [46, 16096926],
                [47, 16096927],
                [48, 16096928],
                [49, 16096929],
                [50, 16096930]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16097145],
                [1, 16097146],
                [2, 16097147],
                [3, 16097148],
                [4, 16097149],
                [5, 16097150],
                [6, 16097151],
                [7, 16097152],
                [8, 16097153],
                [9, 16097154],
                [10, 16097155],
                [11, 16097156],
                [12, 16097157],
                [13, 16097158],
                [14, 16097159],
                [15, 16097160],
                [16, 16097161],
                [17, 16097162],
                [18, 16097163],
                [19, 16097164],
                [20, 16097165],
                [21, 16097166],
                [22, 16097167],
                [23, 16097168],
                [24, 16097169],
                [25, 16097170],
                [26, 16097171],
                [27, 16097172],
                [28, 16097173],
                [29, 16097176],
                [30, 16097177],
                [31, 16097178],
                [32, 16097179],
                [33, 16097180],
                [34, 16097181],
                [35, 16097182],
                [36, 16097183],
                [37, 16097184],
                [38, 16097185],
                [39, 16097186],
                [40, 16097187],
                [41, 16097188],
                [42, 16097189],
                [43, 16097190],
                [44, 16097191],
                [45, 16097192],
                [46, 16097193],
                [47, 16097194],
                [48, 16097195],
                [49, 16097196],
                [50, 16097197]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16117350],
                [1, 16117351],
                [2, 16117352],
                [3, 16117353],
                [4, 16117354],
                [5, 16117355],
                [6, 16117356],
                [7, 16117357],
                [8, 16117358],
                [9, 16117359],
                [10, 16117360],
                [11, 16117361],
                [12, 16117362],
                [13, 16117363],
                [14, 16117364],
                [15, 16117365],
                [16, 16117366],
                [17, 16117367],
                [18, 16117368],
                [19, 16117369],
                [20, 16117370],
                [21, 16117371],
                [22, 16117372],
                [23, 16117373],
                [24, 16117374],
                [25, 16117375],
                [26, 16117376],
                [27, 16117377],
                [28, 16117378],
                [29, 16117379],
                [30, 16117380],
                [31, 16117381],
                [32, 16117382],
                [33, 16117383],
                [34, 16117384],
                [35, 16117385],
                [36, 16117386],
                [37, 16117387],
                [38, 16117388],
                [39, 16117389],
                [40, 16117390],
                [41, 16117391],
                [42, 16117392],
                [43, 16117393],
                [44, 16117394],
                [45, 16117395],
                [46, 16117396],
                [47, 16117397],
                [48, 16117398],
                [49, 16117399],
                [50, 16117400]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16118483],
                [1, 16118484],
                [2, 16118485],
                [3, 16118486],
                [4, 16118487],
                [5, 16118488],
                [6, 16118489],
                [7, 16118490],
                [8, 16118491],
                [9, 16118492],
                [10, 16118493],
                [11, 16118494],
                [12, 16118495],
                [13, 16118496],
                [14, 16118497],
                [15, 16118498],
                [16, 16118499],
                [17, 16118500],
                [18, 16118501],
                [19, 16118502],
                [20, 16118503],
                [21, 16118504],
                [22, 16118505],
                [23, 16118506],
                [24, 16118507],
                [25, 16118508],
                [26, 16118509],
                [27, 16118510],
                [28, 16118511],
                [29, 16118512],
                [30, 16118513],
                [31, 16118514],
                [32, 16118515],
                [33, 16118516],
                [34, 16118517],
                [35, 16118518],
                [36, 16118519],
                [37, 16118520],
                [38, 16118521],
                [39, 16118522],
                [40, 16118523],
                [41, 16118524],
                [42, 16118525],
                [43, 16118526],
                [44, 16118527],
                [45, 16118528],
                [46, 16118529],
                [47, 16118530],
                [48, 16118531],
                [49, 16118532],
                [50, 16118533]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16118499],
                [1, 16118500],
                [2, 16118501],
                [3, 16118502],
                [4, 16118503],
                [5, 16118504],
                [6, 16118505],
                [7, 16118506],
                [8, 16118507],
                [9, 16118508],
                [10, 16118509],
                [11, 16118510],
                [12, 16118511],
                [13, 16118512],
                [14, 16118513],
                [15, 16118514],
                [16, 16118515],
                [17, 16118516],
                [18, 16118517],
                [19, 16118518],
                [20, 16118519],
                [21, 16118520],
                [22, 16118521],
                [23, 16118522],
                [24, 16118523],
                [25, 16118524],
                [26, 16118525],
                [27, 16118526],
                [28, 16118527],
                [29, 16118528],
                [30, 16118529],
                [31, 16118530],
                [32, 16118531],
                [33, 16118532],
                [34, 16118533],
                [35, 16118534],
                [36, 16118535],
                [37, 16118536],
                [38, 16118537],
                [39, 16118538],
                [40, 16118539],
                [41, 16118540],
                [42, 16118541],
                [43, 16118542],
                [44, 16118543],
                [45, 16118544],
                [46, 16118545],
                [47, 16118546],
                [48, 16118547],
                [49, 16118548],
                [50, 16118549]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16118499],
                [1, 16118500],
                [2, 16118501],
                [3, 16118502],
                [4, 16118503],
                [5, 16118504],
                [6, 16118505],
                [7, 16118506],
                [8, 16118507],
                [9, 16118508],
                [10, 16118509],
                [11, 16118510],
                [12, 16118511],
                [13, 16118512],
                [14, 16118513],
                [15, 16118514],
                [16, 16118515],
                [17, 16118516],
                [18, 16118517],
                [19, 16118518],
                [20, 16118519],
                [21, 16118520],
                [22, 16118521],
                [23, 16118522],
                [24, 16118523],
                [25, 16118524],
                [26, 16118525],
                [27, 16118526],
                [28, 16118527],
                [29, 16118528],
                [30, 16118529],
                [31, 16118530],
                [32, 16118531],
                [33, 16118532],
                [34, 16118533],
                [35, 16118534],
                [36, 16118535],
                [37, 16118536],
                [38, 16118537],
                [39, 16118538],
                [40, 16118539],
                [41, 16118540],
                [42, 16118541],
                [43, 16118542],
                [44, 16118543],
                [45, 16118544],
                [46, 16118545],
                [47, 16118546],
                [48, 16118547],
                [49, 16118548],
                [50, 16118549]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16118499],
                [1, 16118500],
                [2, 16118501],
                [3, 16118502],
                [4, 16118503],
                [5, 16118504],
                [6, 16118505],
                [7, 16118506],
                [8, 16118507],
                [9, 16118508],
                [10, 16118509],
                [11, 16118510],
                [12, 16118511],
                [13, 16118512],
                [14, 16118513],
                [15, 16118514],
                [16, 16118515],
                [17, 16118516],
                [18, 16118517],
                [19, 16118518],
                [20, 16118519],
                [21, 16118520],
                [22, 16118521],
                [23, 16118522],
                [24, 16118523],
                [25, 16118524],
                [26, 16118525],
                [27, 16118526],
                [28, 16118527],
                [29, 16118528],
                [30, 16118529],
                [31, 16118530],
                [32, 16118531],
                [33, 16118532],
                [34, 16118533],
                [35, 16118534],
                [36, 16118535],
                [37, 16118536],
                [38, 16118537],
                [39, 16118538],
                [40, 16118539],
                [41, 16118540],
                [42, 16118541],
                [43, 16118542],
                [44, 16118543],
                [45, 16118544],
                [46, 16118545],
                [47, 16118546],
                [48, 16118547],
                [49, 16118548],
                [50, 16118549]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16123411],
                [1, 16123412],
                [2, 16123413],
                [3, 16123414],
                [4, 16123415],
                [5, 16123416],
                [6, 16123417],
                [7, 16123418],
                [8, 16123419],
                [9, 16123420],
                [10, 16123421],
                [11, 16123422],
                [12, 16123423],
                [13, 16123424],
                [14, 16123425],
                [15, 16123426],
                [16, 16123427],
                [17, 16123428],
                [18, 16123429],
                [19, 16123430],
                [20, 16123431],
                [21, 16123432],
                [22, 16123433],
                [23, 16123434],
                [24, 16123435],
                [25, 16123436],
                [26, 16123437],
                [27, 16123438],
                [28, 16123439],
                [29, 16123440],
                [30, 16123441],
                [31, 16123442],
                [32, 16123443],
                [33, 16123444],
                [34, 16123445],
                [35, 16123446],
                [36, 16123447],
                [37, 16123448],
                [38, 16123449],
                [39, 16123450],
                [40, 16123451],
                [41, 16123452],
                [42, 16123453],
                [43, 16123454],
                [44, 16123455],
                [45, 16123456],
                [46, 16123457],
                [47, 16123458],
                [48, 16123459],
                [49, 16123460],
                [50, 16123461]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [6, 16123417],
                [7, 16123418],
                [8, 16123419],
                [9, 16123420],
                [10, 16123421],
                [11, 16123422],
                [12, 16123423],
                [13, 16123424],
                [14, 16123425],
                [15, 16123426],
                [16, 16123427],
                [17, 16123428],
                [18, 16123429],
                [19, 16123430],
                [20, 16123431],
                [21, 16123432],
                [22, 16123433],
                [23, 16123434],
                [24, 16123435],
                [25, 16123436],
                [26, 16123437],
                [27, 16123438],
                [28, 16123439],
                [29, 16123440],
                [30, 16123441],
                [31, 16123442],
                [32, 16123443],
                [33, 16123444],
                [34, 16123445],
                [35, 16123446],
                [36, 16123447],
                [37, 16123448],
                [38, 16123449],
                [39, 16123450],
                [40, 16123451],
                [41, 16123452],
                [42, 16123453],
                [43, 16123454],
                [44, 16123455],
                [45, 16123456],
                [46, 16123457],
                [47, 16123458],
                [48, 16123459],
                [49, 16123460],
                [50, 16123461]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16165860],
                [1, 16165861],
                [2, 16165862],
                [3, 16165863],
                [4, 16165864],
                [5, 16165865],
                [6, 16165866],
                [7, 16165867],
                [8, 16165868],
                [9, 16165869],
                [10, 16165870],
                [11, 16165871],
                [12, 16165872],
                [13, 16165873],
                [14, 16165874],
                [15, 16165875],
                [16, 16165876],
                [17, 16165877],
                [18, 16165878],
                [19, 16165879],
                [20, 16165880],
                [21, 16165881],
                [22, 16165882],
                [23, 16165883],
                [24, 16165884],
                [25, 16165885],
                [26, 16165886],
                [27, 16165887],
                [28, 16165888],
                [29, 16165889],
                [30, 16165890],
                [31, 16165891],
                [32, 16165892],
                [33, 16165893],
                [34, 16165894],
                [35, 16165895],
                [36, 16165896],
                [37, 16165897],
                [38, 16165898],
                [39, 16165899],
                [40, 16165900]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16180871],
                [1, 16180872],
                [2, 16180873],
                [3, 16180874],
                [4, 16180875],
                [5, 16180876],
                [6, 16180877],
                [7, 16180878],
                [8, 16180879],
                [9, 16180880],
                [10, 16180881],
                [11, 16180882],
                [12, 16180883],
                [13, 16180884],
                [14, 16180885],
                [15, 16180886],
                [16, 16180887],
                [17, 16180888],
                [18, 16180889],
                [19, 16180890],
                [20, 16180891],
                [21, 16180892],
                [22, 16180893],
                [23, 16180894],
                [24, 16180895],
                [25, 16180896],
                [26, 16180897],
                [27, 16180898],
                [28, 16180899],
                [29, 16180900],
                [30, 16180901],
                [31, 16180902],
                [32, 16180903],
                [33, 16180904],
                [34, 16180905],
                [35, 16180906],
                [36, 16180907],
                [37, 16180908],
                [38, 16180909],
                [39, 16180910],
                [40, 16180911],
                [41, 16180912],
                [42, 16180913],
                [43, 16180914],
                [44, 16180915],
                [45, 16180916],
                [46, 16180917],
                [47, 16180918],
                [48, 16180919],
                [49, 16180920],
                [50, 16180921]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16189705],
                [1, 16189706],
                [2, 16189707],
                [3, 16189708],
                [4, 16189709],
                [5, 16189710],
                [6, 16189711],
                [7, 16189712],
                [8, 16189713],
                [9, 16189714],
                [10, 16189715],
                [11, 16189716],
                [12, 16189717],
                [13, 16189718],
                [14, 16189719],
                [15, 16189720],
                [16, 16189721],
                [17, 16189722],
                [18, 16189723],
                [19, 16189724],
                [20, 16189725],
                [21, 16189726],
                [22, 16189727],
                [23, 16189728],
                [24, 16189729],
                [25, 16189730],
                [26, 16189731],
                [27, 16189732],
                [28, 16189733],
                [29, 16189734],
                [30, 16189735],
                [31, 16189736],
                [32, 16189737],
                [33, 16189738],
                [34, 16189739],
                [35, 16189740],
                [36, 16189741],
                [37, 16189742],
                [38, 16189743],
                [39, 16189744],
                [40, 16189745],
                [41, 16189746],
                [42, 16189747],
                [43, 16189748],
                [44, 16189749],
                [45, 16189750],
                [46, 16189751],
                [47, 16189752],
                [48, 16189753],
                [49, 16189754],
                [50, 16189755]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16231271],
                [1, 16231272],
                [2, 16231273],
                [3, 16231274],
                [4, 16231275],
                [5, 16231276],
                [6, 16231277],
                [7, 16231278],
                [8, 16231279],
                [9, 16231280],
                [10, 16231281],
                [11, 16231282],
                [12, 16231283],
                [13, 16231284],
                [14, 16231285],
                [15, 16231286],
                [16, 16231287],
                [17, 16231288],
                [18, 16231289],
                [19, 16231290],
                [20, 16231291],
                [21, 16231292],
                [22, 16231293],
                [23, 16231294],
                [24, 16231295],
                [25, 16231296],
                [26, 16231297],
                [27, 16231298],
                [28, 16231299],
                [29, 16231300],
                [30, 16231301],
                [31, 16231302],
                [32, 16231303],
                [33, 16231304],
                [34, 16231305],
                [35, 16231306],
                [36, 16231307],
                [37, 16231308],
                [38, 16231309],
                [39, 16231310],
                [40, 16231311],
                [41, 16231312],
                [42, 16231313],
                [43, 16231314],
                [44, 16231315],
                [45, 16231316],
                [46, 16231317],
                [47, 16231318],
                [48, 16231319],
                [49, 16231320],
                [50, 16231321]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16237657],
                [1, 16237658],
                [2, 16237659],
                [3, 16237660],
                [4, 16237661],
                [5, 16237662],
                [6, 16237663],
                [7, 16237664],
                [8, 16237665],
                [9, 16237666],
                [10, 16237667],
                [11, 16237668],
                [12, 16237669],
                [13, 16237670],
                [14, 16237671],
                [15, 16237672],
                [16, 16237673],
                [17, 16237674],
                [18, 16237675],
                [19, 16237676],
                [20, 16237677],
                [21, 16237678],
                [22, 16237679],
                [23, 16237680],
                [24, 16237681],
                [25, 16237682],
                [26, 16237683],
                [27, 16237684],
                [28, 16237685],
                [29, 16237686],
                [30, 16237687],
                [31, 16237688],
                [32, 16237689],
                [33, 16237690],
                [34, 16237691],
                [35, 16237692],
                [36, 16237693],
                [37, 16237694],
                [38, 16237695],
                [39, 16237696],
                [40, 16237697],
                [41, 16237698],
                [42, 16237699],
                [43, 16237700],
                [44, 16237701],
                [45, 16237702],
                [46, 16237703],
                [47, 16237704],
                [48, 16237705],
                [49, 16237706],
                [50, 16237707]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [9, 16255012],
                [10, 16255013],
                [11, 16255014],
                [12, 16255015],
                [13, 16255016],
                [14, 16255017],
                [15, 16255018],
                [16, 16255019],
                [17, 16255020],
                [18, 16255021],
                [19, 16255022],
                [20, 16255023],
                [21, 16255024],
                [22, 16255025],
                [23, 16255026],
                [24, 16255027],
                [25, 16255028],
                [26, 16255029],
                [27, 16255030],
                [28, 16255031],
                [29, 16255032],
                [30, 16255033],
                [31, 16255034],
                [32, 16255035],
                [33, 16255036],
                [34, 16255037],
                [35, 16255038],
                [36, 16255039],
                [37, 16255040],
                [38, 16255041],
                [39, 16255042],
                [40, 16255043],
                [41, 16255044],
                [42, 16255045],
                [43, 16255046],
                [44, 16255047],
                [45, 16255048],
                [46, 16255049],
                [47, 16255050],
                [48, 16255051],
                [49, 16255052],
                [50, 16255053]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16255391],
                [1, 16255392],
                [2, 16255393],
                [3, 16255394],
                [4, 16255395],
                [5, 16255396],
                [6, 16255397],
                [7, 16255398],
                [8, 16255399],
                [9, 16255400],
                [10, 16255401],
                [11, 16255402],
                [12, 16255403],
                [13, 16255404],
                [14, 16255405],
                [15, 16255406],
                [16, 16255407],
                [17, 16255408],
                [18, 16255409],
                [19, 16255410],
                [20, 16255411],
                [21, 16255412],
                [22, 16255413],
                [23, 16255414],
                [24, 16255415],
                [25, 16255416],
                [26, 16255417],
                [27, 16255418],
                [28, 16255419],
                [29, 16255420],
                [30, 16255421],
                [31, 16255422],
                [32, 16255423],
                [33, 16255424],
                [34, 16255425],
                [35, 16255426],
                [36, 16255427],
                [37, 16255428],
                [38, 16255429],
                [39, 16255430],
                [40, 16255431],
                [41, 16255432],
                [42, 16255433],
                [43, 16255434],
                [44, 16255435],
                [45, 16255436],
                [46, 16255437],
                [47, 16255438],
                [48, 16255439],
                [49, 16255440],
                [50, 16255441]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16255392],
                [1, 16255393],
                [2, 16255394],
                [3, 16255395],
                [4, 16255396],
                [5, 16255397],
                [6, 16255398],
                [7, 16255399],
                [8, 16255400],
                [9, 16255401],
                [10, 16255402],
                [11, 16255403],
                [12, 16255404],
                [13, 16255405],
                [14, 16255406],
                [15, 16255407],
                [16, 16255408],
                [17, 16255409],
                [18, 16255410],
                [19, 16255411],
                [20, 16255412],
                [21, 16255413],
                [22, 16255414],
                [23, 16255415],
                [24, 16255416],
                [25, 16255417],
                [26, 16255418],
                [27, 16255419],
                [28, 16255420],
                [29, 16255421],
                [30, 16255422],
                [31, 16255423],
                [32, 16255424],
                [33, 16255425],
                [34, 16255426],
                [35, 16255427],
                [36, 16255428],
                [37, 16255429],
                [38, 16255430],
                [39, 16255431],
                [40, 16255432],
                [41, 16255433],
                [42, 16255434],
                [43, 16255435],
                [44, 16255436],
                [45, 16255437],
                [46, 16255438],
                [47, 16255439],
                [48, 16255440],
                [49, 16255441]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16256084],
                [1, 16256085],
                [2, 16256086],
                [3, 16256087],
                [4, 16256088],
                [5, 16256089],
                [6, 16256090],
                [7, 16256091],
                [8, 16256092],
                [9, 16256093],
                [10, 16256094],
                [11, 16256095],
                [12, 16256096],
                [13, 16256097],
                [14, 16256098],
                [15, 16256099],
                [16, 16256100],
                [17, 16256101],
                [18, 16256102],
                [19, 16256103],
                [20, 16256104],
                [21, 16256105],
                [22, 16256106],
                [23, 16256107],
                [24, 16256108],
                [25, 16256109],
                [26, 16256110],
                [27, 16256111],
                [28, 16256112],
                [29, 16256113],
                [30, 16256114],
                [31, 16256115],
                [32, 16256116],
                [33, 16256117],
                [34, 16256118],
                [35, 16256119],
                [36, 16256120],
                [37, 16256121],
                [38, 16256122],
                [39, 16256123],
                [40, 16256124],
                [41, 16256125],
                [42, 16256126],
                [43, 16256127],
                [44, 16256128]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [3, 16256224],
                [4, 16256225],
                [5, 16256226],
                [6, 16256227],
                [7, 16256228],
                [8, 16256229],
                [9, 16256230],
                [10, 16256231],
                [11, 16256232],
                [12, 16256233],
                [13, 16256234],
                [14, 16256235],
                [15, 16256236],
                [16, 16256237],
                [17, 16256238],
                [18, 16256239],
                [19, 16256240],
                [20, 16256241],
                [21, 16256242],
                [22, 16256243],
                [23, 16256244],
                [24, 16256245],
                [25, 16256246],
                [26, 16256247],
                [27, 16256248],
                [28, 16256249],
                [29, 16256250],
                [30, 16256251],
                [31, 16256252],
                [32, 16256253],
                [33, 16256254],
                [34, 16256255],
                [35, 16256256],
                [36, 16256257],
                [37, 16256258],
                [38, 16256259],
                [39, 16256260],
                [40, 16256261],
                [41, 16256262],
                [42, 16256263],
                [43, 16256264],
                [44, 16256265],
                [45, 16256266],
                [46, 16256267],
                [47, 16256268],
                [48, 16256269],
                [49, 16256270],
                [50, 16256271]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16325199],
                [1, 16325200],
                [2, 16325201],
                [3, 16325202],
                [4, 16325203],
                [5, 16325204],
                [6, 16325205],
                [7, 16325206],
                [8, 16325207],
                [9, 16325208],
                [10, 16325209],
                [11, 16325210],
                [12, 16325211],
                [13, 16325212],
                [14, 16325213],
                [15, 16325214],
                [16, 16325215],
                [17, 16325216],
                [18, 16325217],
                [19, 16325218],
                [20, 16325219],
                [21, 16325220],
                [22, 16325221],
                [23, 16325222],
                [24, 16325223],
                [25, 16325224],
                [26, 16325225],
                [27, 16325226],
                [28, 16325227],
                [29, 16325228],
                [30, 16325229],
                [31, 16325230],
                [32, 16325231],
                [33, 16325232],
                [34, 16325233],
                [35, 16325234],
                [36, 16325235],
                [37, 16325236],
                [38, 16325237],
                [39, 16325238],
                [40, 16325239],
                [41, 16325240]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [13, 16352865],
                [14, 16352866],
                [15, 16352867],
                [16, 16352868],
                [17, 16352869],
                [18, 16352870],
                [19, 16352871],
                [20, 16352872],
                [21, 16352873],
                [22, 16352874],
                [23, 16352875],
                [24, 16352876],
                [25, 16352877],
                [26, 16352878],
                [27, 16352879],
                [28, 16352880],
                [29, 16352881],
                [30, 16352882],
                [31, 16352883],
                [32, 16352884],
                [33, 16352885],
                [34, 16352886],
                [35, 16352887],
                [36, 16352888],
                [37, 16352889],
                [38, 16352890],
                [39, 16352891],
                [40, 16352892],
                [41, 16352893],
                [42, 16352894],
                [43, 16352895],
                [44, 16352896],
                [45, 16352897],
                [46, 16352898],
                [47, 16352899],
                [48, 16352900],
                [49, 16352901],
                [50, 16352902]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 16352968],
                [1, 16352969],
                [2, 16352970],
                [3, 16352971],
                [4, 16352972],
                [5, 16352973],
                [6, 16352974],
                [7, 16352975],
                [8, 16352976],
                [9, 16352977],
                [10, 16352978],
                [11, 16352979],
                [12, 16352980],
                [13, 16352981],
                [14, 16352982],
                [15, 16352983],
                [16, 16352984],
                [17, 16352985],
                [18, 16352986],
                [19, 16352987],
                [20, 16352988],
                [21, 16352989],
                [22, 16352990],
                [23, 16352991],
                [24, 16352992],
                [25, 16352993],
                [26, 16352994],
                [27, 16352995],
                [28, 16352996],
                [29, 16352997],
                [30, 16352998],
                [31, 16352999],
                [32, 16353000],
                [33, 16353001],
                [34, 16353002],
                [35, 16353003],
                [36, 16353004],
                [37, 16353005],
                [38, 16353006],
                [39, 16353007],
                [40, 16353008],
                [41, 16353009],
                [42, 16353010],
                [43, 16353011]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [5, 16414998],
                [6, 16414999],
                [7, 16415000],
                [8, 16415001],
                [9, 16415002],
                [10, 16415003],
                [11, 16415004],
                [12, 16415005],
                [13, 16415006],
                [14, 16415007],
                [15, 16415008],
                [16, 16415009],
                [17, 16415010],
                [18, 16415011],
                [19, 16415012],
                [20, 16415013],
                [21, 16415014],
                [22, 16415015],
                [23, 16415016],
                [24, 16415017],
                [25, 16415018],
                [26, 16415019],
                [27, 16415020],
                [28, 16415021],
                [29, 16415022],
                [30, 16415023],
                [31, 16415024],
                [32, 16415025],
                [33, 16415026],
                [34, 16415027],
                [35, 16415028],
                [36, 16415029],
                [37, 16415030],
                [38, 16415031],
                [39, 16415032],
                [40, 16415033],
                [41, 16415034],
                [42, 16415035],
                [43, 16415036],
                [44, 16415037],
                [45, 16415038],
                [46, 16415039],
                [47, 16415040],
                [48, 16415041],
                [49, 16415042],
                [50, 16415043]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17031591],
                [1, 17031592],
                [2, 17031593],
                [3, 17031594],
                [4, 17031595],
                [5, 17031596],
                [6, 17031597],
                [7, 17031598],
                [8, 17031599],
                [9, 17031600],
                [10, 17031601],
                [11, 17031602],
                [12, 17031603],
                [13, 17031604],
                [14, 17031605],
                [15, 17031606],
                [16, 17031607],
                [17, 17031608],
                [18, 17031609],
                [19, 17031610],
                [20, 17031611],
                [21, 17031612],
                [22, 17031613],
                [27, 17031614],
                [28, 17031615],
                [29, 17031616],
                [30, 17031617],
                [31, 17031618],
                [32, 17031619],
                [33, 17031620],
                [34, 17031621],
                [35, 17031622],
                [36, 17031623],
                [37, 17031624],
                [38, 17031625],
                [39, 17031626],
                [40, 17031627],
                [41, 17031628],
                [42, 17031629],
                [43, 17031630],
                [44, 17031631],
                [45, 17031632],
                [46, 17031633],
                [47, 17031634],
                [48, 17031635],
                [49, 17031636],
                [50, 17031637]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17057382],
                [1, 17057383],
                [2, 17057384],
                [3, 17057385],
                [4, 17057386],
                [5, 17057387],
                [6, 17057388],
                [7, 17057389],
                [8, 17057390],
                [9, 17057391],
                [10, 17057392],
                [11, 17057393],
                [12, 17057394],
                [13, 17057395],
                [14, 17057396],
                [15, 17057397],
                [16, 17057398],
                [17, 17057399],
                [19, 17057400],
                [20, 17057401],
                [21, 17057402],
                [22, 17057403],
                [23, 17057404],
                [24, 17057405],
                [25, 17057406],
                [26, 17057407],
                [27, 17057408],
                [28, 17057409],
                [29, 17057410],
                [30, 17057411],
                [31, 17057412],
                [32, 17057413],
                [33, 17057414],
                [34, 17057415],
                [35, 17057416],
                [36, 17057417],
                [37, 17057418],
                [38, 17057419],
                [39, 17057420],
                [40, 17057421],
                [41, 17057422],
                [42, 17057423],
                [43, 17057424],
                [44, 17057425],
                [45, 17057426],
                [46, 17057427],
                [47, 17057428],
                [48, 17057429],
                [49, 17057430],
                [50, 17057431]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17092766],
                [1, 17092767],
                [2, 17092768],
                [3, 17092769],
                [4, 17092770],
                [5, 17092771],
                [6, 17092772],
                [7, 17092773],
                [8, 17092774],
                [9, 17092775],
                [10, 17092776],
                [11, 17092777],
                [12, 17092778],
                [13, 17092779],
                [14, 17092780],
                [15, 17092781],
                [16, 17092782],
                [17, 17094966],
                [18, 17094967],
                [19, 17094968],
                [20, 17094969],
                [21, 17094970],
                [22, 17094971],
                [23, 17094972],
                [24, 17094973],
                [25, 17094974],
                [26, 17094975],
                [27, 17094976],
                [28, 17094977],
                [29, 17094978],
                [30, 17094979],
                [31, 17094980],
                [32, 17094981],
                [33, 17094982],
                [34, 17094983],
                [35, 17094984],
                [36, 17094985],
                [37, 17094986],
                [38, 17094987],
                [39, 17094988],
                [40, 17094989],
                [41, 17094990],
                [42, 17094991],
                [43, 17094992],
                [44, 17094993],
                [45, 17094994],
                [46, 17094995],
                [47, 17094996],
                [48, 17094997],
                [49, 17094998],
                [50, 17094999]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17092782],
                [1, 17094966],
                [2, 17094967],
                [3, 17094968],
                [4, 17094969],
                [5, 17094970],
                [6, 17094971],
                [7, 17094972],
                [8, 17094973],
                [9, 17094974],
                [10, 17094975],
                [11, 17094976],
                [12, 17094977],
                [13, 17094978],
                [14, 17094979],
                [15, 17094980],
                [16, 17094981],
                [17, 17094982],
                [18, 17094983],
                [19, 17094984],
                [20, 17094985],
                [21, 17094986],
                [22, 17094987],
                [23, 17094988],
                [24, 17094989],
                [25, 17094990],
                [26, 17094991],
                [27, 17094992],
                [28, 17094993],
                [29, 17094994],
                [30, 17094995],
                [31, 17094996],
                [32, 17094997],
                [33, 17094998],
                [34, 17094999],
                [35, 17095000],
                [36, 17095001],
                [37, 17095002],
                [38, 17095003],
                [39, 17095004],
                [40, 17095005],
                [41, 17095006],
                [42, 17095007],
                [43, 17095008],
                [44, 17095009],
                [45, 17095010],
                [46, 17095011],
                [47, 17095012],
                [48, 17095013],
                [49, 17095014],
                [50, 17095015]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17092782],
                [1, 17094966],
                [2, 17094967],
                [3, 17094968],
                [4, 17094969],
                [5, 17094970],
                [6, 17094971],
                [7, 17094972],
                [8, 17094973],
                [9, 17094974],
                [10, 17094975],
                [11, 17094976],
                [12, 17094977],
                [13, 17094978],
                [14, 17094979],
                [15, 17094980],
                [16, 17094981],
                [17, 17094982],
                [18, 17094983],
                [19, 17094984],
                [20, 17094985],
                [21, 17094986],
                [22, 17094987],
                [23, 17094988],
                [24, 17094989],
                [25, 17094990],
                [26, 17094991],
                [27, 17094992],
                [28, 17094993],
                [29, 17094994],
                [30, 17094995],
                [31, 17094996],
                [32, 17094997],
                [33, 17094998],
                [34, 17094999],
                [35, 17095000],
                [36, 17095001],
                [37, 17095002],
                [38, 17095003],
                [39, 17095004],
                [40, 17095005],
                [41, 17095006],
                [42, 17095007],
                [43, 17095008],
                [44, 17095009],
                [45, 17095010],
                [46, 17095011],
                [47, 17095012],
                [48, 17095013],
                [49, 17095014],
                [50, 17095015]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [9, 17137287],
                [10, 17137288],
                [11, 17137289],
                [12, 17137290],
                [13, 17137291],
                [14, 17137292],
                [15, 17137293],
                [16, 17137294],
                [17, 17137295],
                [18, 17137296],
                [19, 17137297],
                [20, 17137298],
                [21, 17137299],
                [22, 17137300],
                [23, 17137301],
                [24, 17137302],
                [25, 17137303],
                [26, 17137304],
                [27, 17137305],
                [28, 17137306],
                [29, 17137307],
                [30, 17137308],
                [31, 17137309],
                [32, 17137310],
                [33, 17137311],
                [34, 17137312],
                [35, 17137313],
                [36, 17137314],
                [37, 17137315],
                [38, 17137316],
                [39, 17137317],
                [40, 17137318],
                [41, 17137319]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [2, 17306238],
                [3, 17306239],
                [4, 17306240],
                [5, 17306241],
                [6, 17306242],
                [7, 17306243],
                [8, 17306244],
                [9, 17306245],
                [10, 17306246],
                [11, 17306247],
                [12, 17306248],
                [13, 17306249],
                [14, 17306250],
                [15, 17306251],
                [16, 17306252],
                [17, 17306253],
                [18, 17306254],
                [19, 17306255],
                [20, 17306256],
                [21, 17306257],
                [22, 17306258],
                [23, 17306259],
                [24, 17306260],
                [25, 17306261],
                [26, 17306262],
                [27, 17306263],
                [28, 17306264],
                [29, 17306265],
                [30, 17306266],
                [31, 17306267],
                [32, 17306268],
                [33, 17306269],
                [34, 17306270],
                [35, 17306271],
                [36, 17306272],
                [37, 17306273],
                [38, 17306274],
                [39, 17306275],
                [40, 17306276],
                [41, 17306277],
                [42, 17306278],
                [43, 17306279],
                [44, 17306280],
                [45, 17306281],
                [46, 17306282],
                [47, 17306283],
                [48, 17306284],
                [49, 17306285]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [4, 17561868],
                [5, 17561869],
                [6, 17561870],
                [7, 17561871],
                [8, 17561872],
                [9, 17561873],
                [10, 17561874],
                [11, 17561875],
                [12, 17561876],
                [13, 17561877],
                [14, 17561878],
                [15, 17561879],
                [16, 17561880],
                [17, 17561881],
                [18, 17561882],
                [19, 17561883],
                [20, 17561884],
                [21, 17561885],
                [22, 17561886],
                [23, 17561887],
                [24, 17561888],
                [25, 17561889],
                [26, 17561890],
                [27, 17561891],
                [28, 17561892],
                [29, 17561893],
                [30, 17561894],
                [31, 17561895],
                [32, 17561896],
                [33, 17561897],
                [34, 17561898],
                [35, 17561899],
                [36, 17561900],
                [37, 17561901],
                [38, 17561902],
                [39, 17561903],
                [40, 17561904],
                [41, 17561905],
                [42, 17561906],
                [43, 17561907],
                [44, 17561908],
                [45, 17561909],
                [46, 17561910],
                [47, 17561911],
                [48, 17561912]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17566078],
                [1, 17566079],
                [2, 17566080],
                [3, 17566081],
                [4, 17566082],
                [5, 17566083],
                [6, 17566084],
                [7, 17566085],
                [8, 17566086],
                [9, 17566087],
                [10, 17566088],
                [11, 17566089],
                [12, 17566090],
                [13, 17566091],
                [14, 17566092],
                [15, 17566093],
                [16, 17566094],
                [17, 17566095],
                [18, 17566096],
                [19, 17566097],
                [20, 17566098],
                [21, 17566099],
                [22, 17566100],
                [23, 17566101],
                [24, 17566102],
                [25, 17566103],
                [26, 17566104],
                [27, 17566105],
                [28, 17566106],
                [29, 17566107],
                [30, 17566108],
                [31, 17566109],
                [32, 17566110],
                [33, 17566111],
                [34, 17566112],
                [35, 17566113],
                [36, 17566114],
                [37, 17566115],
                [38, 17566116],
                [39, 17566117],
                [40, 17566118],
                [41, 17577951],
                [42, 17577952],
                [43, 17577953],
                [44, 17577954],
                [45, 17577955],
                [46, 17577956],
                [47, 17577957],
                [48, 17577958],
                [49, 17577959],
                [50, 17577960]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17566108],
                [1, 17566109],
                [2, 17566110],
                [3, 17566111],
                [4, 17566112],
                [5, 17566113],
                [6, 17566114],
                [7, 17566115],
                [8, 17566116],
                [9, 17566117],
                [10, 17566118],
                [11, 17577951],
                [12, 17577952],
                [13, 17577953],
                [14, 17577954],
                [15, 17577955],
                [16, 17577956],
                [17, 17577957],
                [18, 17577958],
                [19, 17577959],
                [20, 17577960],
                [21, 17577961],
                [22, 17577962],
                [23, 17577963],
                [24, 17577964],
                [25, 17577965],
                [26, 17577966],
                [27, 17577967],
                [28, 17577968],
                [29, 17577969],
                [30, 17577970],
                [31, 17577971],
                [32, 17577972],
                [33, 17577973],
                [34, 17577974],
                [35, 17577975],
                [36, 17578686],
                [37, 17578687],
                [38, 17578688],
                [39, 17578689],
                [40, 17578690],
                [41, 17578691],
                [42, 17578692],
                [43, 17578693],
                [44, 17578694],
                [45, 17578695],
                [46, 17578696],
                [47, 17578697],
                [48, 17578698],
                [49, 17578699],
                [50, 17578700]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17566111],
                [1, 17566112],
                [2, 17566113],
                [3, 17566114],
                [4, 17566115],
                [5, 17566116],
                [6, 17566117],
                [7, 17566118],
                [8, 17577951],
                [9, 17577952],
                [10, 17577953],
                [11, 17577954],
                [12, 17577955],
                [13, 17577956],
                [14, 17577957],
                [15, 17577958],
                [16, 17577959],
                [17, 17577960],
                [18, 17577961],
                [19, 17577962],
                [20, 17577963],
                [21, 17577964],
                [22, 17577965],
                [23, 17577966],
                [24, 17577967],
                [25, 17577968],
                [26, 17577969],
                [27, 17577970],
                [28, 17577971],
                [29, 17577972],
                [30, 17577973],
                [31, 17577974],
                [32, 17577975],
                [33, 17578686],
                [34, 17578687],
                [35, 17578688],
                [36, 17578689],
                [37, 17578690],
                [38, 17578691],
                [39, 17578692],
                [40, 17578693],
                [41, 17578694],
                [42, 17578695],
                [43, 17578696],
                [44, 17578697],
                [45, 17578698],
                [46, 17578699],
                [47, 17578700],
                [48, 17578701],
                [49, 17578702],
                [50, 17578703]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17566111],
                [1, 17566112],
                [2, 17566113],
                [3, 17566114],
                [4, 17566115],
                [5, 17566116],
                [6, 17566117],
                [7, 17566118],
                [8, 17577951],
                [9, 17577952],
                [10, 17577953],
                [11, 17577954],
                [12, 17577955],
                [13, 17577956],
                [14, 17577957],
                [15, 17577958],
                [16, 17577959],
                [17, 17577960],
                [18, 17577961],
                [19, 17577962],
                [20, 17577963],
                [21, 17577964],
                [22, 17577965],
                [23, 17577966],
                [24, 17577967],
                [25, 17577968],
                [26, 17577969],
                [27, 17577970],
                [28, 17577971],
                [29, 17577972],
                [30, 17577973],
                [31, 17577974],
                [32, 17577975],
                [33, 17578686],
                [34, 17578687],
                [35, 17578688],
                [36, 17578689],
                [37, 17578690],
                [38, 17578691],
                [39, 17578692],
                [40, 17578693],
                [41, 17578694],
                [42, 17578695],
                [43, 17578696],
                [44, 17578697],
                [45, 17578698],
                [46, 17578699],
                [47, 17578700],
                [48, 17578701],
                [49, 17578702],
                [50, 17578703]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17566111],
                [1, 17566112],
                [2, 17566113],
                [3, 17566114],
                [4, 17566115],
                [5, 17566116],
                [6, 17566117],
                [7, 17566118],
                [8, 17577951],
                [9, 17577952],
                [10, 17577953],
                [11, 17577954],
                [12, 17577955],
                [13, 17577956],
                [14, 17577957],
                [15, 17577958],
                [16, 17577959],
                [17, 17577960],
                [18, 17577961],
                [19, 17577962],
                [20, 17577963],
                [21, 17577964],
                [22, 17577965],
                [23, 17577966],
                [24, 17577967],
                [25, 17577968],
                [26, 17577969],
                [27, 17577970],
                [28, 17577971],
                [29, 17577972],
                [30, 17577973],
                [31, 17577974],
                [32, 17577975],
                [33, 17578686],
                [34, 17578687],
                [35, 17578688],
                [36, 17578689],
                [37, 17578690],
                [38, 17578691],
                [39, 17578692],
                [40, 17578693],
                [41, 17578694],
                [42, 17578695],
                [43, 17578696],
                [44, 17578697],
                [45, 17578698],
                [46, 17578699],
                [47, 17578700],
                [48, 17578701],
                [49, 17578702],
                [50, 17578703]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17566111],
                [1, 17566112],
                [2, 17566113],
                [3, 17566114],
                [4, 17566115],
                [5, 17566116],
                [6, 17566117],
                [7, 17566118],
                [8, 17577951],
                [9, 17577952],
                [10, 17577953],
                [11, 17577954],
                [12, 17577955],
                [13, 17577956],
                [14, 17577957],
                [15, 17577958],
                [16, 17577959],
                [17, 17577960],
                [18, 17577961],
                [19, 17577962],
                [20, 17577963],
                [21, 17577964],
                [22, 17577965],
                [23, 17577966],
                [24, 17577967],
                [25, 17577968],
                [26, 17577969],
                [27, 17577970],
                [28, 17577971],
                [29, 17577972],
                [30, 17577973],
                [31, 17577974],
                [32, 17577975],
                [33, 17578686],
                [34, 17578687],
                [35, 17578688],
                [36, 17578689],
                [37, 17578690],
                [38, 17578691],
                [39, 17578692],
                [40, 17578693],
                [41, 17578694],
                [42, 17578695],
                [43, 17578696],
                [44, 17578697],
                [45, 17578698],
                [46, 17578699],
                [47, 17578700],
                [48, 17578701],
                [49, 17578702],
                [50, 17578703]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17566112],
                [1, 17566113],
                [2, 17566114],
                [3, 17566115],
                [4, 17566116],
                [5, 17566117],
                [6, 17566118],
                [7, 17577951],
                [8, 17577952],
                [9, 17577953],
                [10, 17577954],
                [11, 17577955],
                [12, 17577956],
                [13, 17577957],
                [14, 17577958],
                [15, 17577959],
                [16, 17577960],
                [17, 17577961],
                [18, 17577962],
                [19, 17577963],
                [20, 17577964],
                [21, 17577965],
                [22, 17577966],
                [23, 17577967],
                [24, 17577968],
                [25, 17577969],
                [26, 17577970],
                [27, 17577971],
                [28, 17577972],
                [29, 17577973],
                [30, 17577974],
                [31, 17577975],
                [32, 17578686],
                [33, 17578687],
                [34, 17578688],
                [35, 17578689],
                [36, 17578690],
                [37, 17578691],
                [38, 17578692],
                [39, 17578693],
                [40, 17578694],
                [41, 17578695],
                [42, 17578696],
                [43, 17578697],
                [44, 17578698],
                [45, 17578699],
                [46, 17578700],
                [47, 17578701],
                [48, 17578702],
                [49, 17578703],
                [50, 17578704]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17566113],
                [1, 17566114],
                [2, 17566115],
                [3, 17566116],
                [4, 17566117],
                [5, 17566118],
                [6, 17577951],
                [7, 17577952],
                [8, 17577953],
                [9, 17577954],
                [10, 17577955],
                [11, 17577956],
                [12, 17577957],
                [13, 17577958],
                [14, 17577959],
                [15, 17577960],
                [16, 17577961],
                [17, 17577962],
                [18, 17577963],
                [19, 17577964],
                [20, 17577965],
                [21, 17577966],
                [22, 17577967],
                [23, 17577968],
                [24, 17577969],
                [25, 17577970],
                [26, 17577971],
                [27, 17577972],
                [28, 17577973],
                [29, 17577974],
                [30, 17577975],
                [31, 17578686],
                [32, 17578687],
                [33, 17578688],
                [34, 17578689],
                [35, 17578690],
                [36, 17578691],
                [37, 17578692],
                [38, 17578693],
                [39, 17578694],
                [40, 17578695],
                [41, 17578696],
                [42, 17578697],
                [43, 17578698],
                [44, 17578699],
                [45, 17578700],
                [46, 17578701],
                [47, 17578702],
                [48, 17578703],
                [49, 17578704],
                [50, 17578705]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17566113],
                [1, 17566114],
                [2, 17566115],
                [3, 17566116],
                [4, 17566117],
                [5, 17566118],
                [6, 17577951],
                [7, 17577952],
                [8, 17577953],
                [9, 17577954],
                [10, 17577955],
                [11, 17577956],
                [12, 17577957],
                [13, 17577958],
                [14, 17577959],
                [15, 17577960],
                [16, 17577961],
                [17, 17577962],
                [18, 17577963],
                [19, 17577964],
                [20, 17577965],
                [21, 17577966],
                [22, 17577967],
                [23, 17577968],
                [24, 17577969],
                [25, 17577970],
                [26, 17577971],
                [27, 17577972],
                [28, 17577973],
                [29, 17577974],
                [30, 17577975],
                [31, 17578686],
                [32, 17578687],
                [33, 17578688],
                [34, 17578689],
                [35, 17578690],
                [36, 17578691],
                [37, 17578692],
                [38, 17578693],
                [39, 17578694],
                [40, 17578695],
                [41, 17578696],
                [42, 17578697],
                [43, 17578698],
                [44, 17578699],
                [45, 17578700],
                [46, 17578701],
                [47, 17578702],
                [48, 17578703],
                [49, 17578704],
                [50, 17578705]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [1, 17579733],
                [2, 17579734],
                [3, 17579735],
                [4, 17579736],
                [5, 17579737],
                [6, 17579738],
                [7, 17579739],
                [8, 17579740],
                [9, 17579741],
                [10, 17579742],
                [11, 17579743],
                [12, 17579744],
                [13, 17579745],
                [14, 17579746],
                [15, 17579747],
                [16, 17579748],
                [17, 17579749],
                [18, 17579750],
                [19, 17579751],
                [20, 17579752],
                [21, 17579753],
                [22, 17579754],
                [23, 17579755],
                [24, 17579756],
                [25, 17579757],
                [26, 17579758],
                [27, 17579759],
                [28, 17579760],
                [29, 17579761],
                [30, 17579762],
                [31, 17579763],
                [32, 17579764],
                [33, 17579765],
                [34, 17579766],
                [35, 17579767],
                [36, 17579768],
                [37, 17579769],
                [38, 17579770],
                [39, 17579771],
                [40, 17579772],
                [41, 17579773],
                [42, 17579774],
                [43, 17579775],
                [44, 17579776],
                [45, 17581244],
                [46, 17581245],
                [47, 17581246],
                [48, 17581247],
                [49, 17581248],
                [50, 17581249]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17581369],
                [1, 17581370],
                [2, 17582885],
                [3, 17582886],
                [4, 17582887],
                [5, 17582888],
                [6, 17582889],
                [7, 17582890],
                [8, 17582891],
                [9, 17582892],
                [10, 17582893],
                [11, 17582894],
                [12, 17582895],
                [13, 17582896],
                [14, 17582897],
                [15, 17582898],
                [16, 17582899],
                [17, 17582900],
                [18, 17582901],
                [19, 17582902],
                [20, 17582903],
                [21, 17582904],
                [22, 17582905],
                [23, 17582906],
                [24, 17582907],
                [25, 17582908],
                [26, 17582909],
                [27, 17582910],
                [28, 17582911],
                [29, 17582912],
                [30, 17582913],
                [31, 17582914],
                [32, 17582915],
                [33, 17582916],
                [34, 17582917],
                [35, 17582918],
                [36, 17582919],
                [37, 17582920],
                [38, 17582921],
                [39, 17582922],
                [40, 17582923],
                [41, 17582924],
                [42, 17582925],
                [43, 17582926],
                [44, 17582927],
                [45, 17582928],
                [46, 17582929],
                [47, 17582930],
                [48, 17582931],
                [49, 17582932],
                [50, 17583028]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17581370],
                [1, 17582885],
                [2, 17582886],
                [3, 17582887],
                [4, 17582888],
                [5, 17582889],
                [6, 17582890],
                [7, 17582891],
                [8, 17582892],
                [9, 17582893],
                [10, 17582894],
                [11, 17582895],
                [12, 17582896],
                [13, 17582897],
                [14, 17582898],
                [15, 17582899],
                [16, 17582900],
                [17, 17582901],
                [18, 17582902],
                [19, 17582903],
                [20, 17582904],
                [21, 17582905],
                [22, 17582906],
                [23, 17582907],
                [24, 17582908],
                [25, 17582909],
                [26, 17582910],
                [27, 17582911],
                [28, 17582912],
                [29, 17582913],
                [30, 17582914],
                [31, 17582915],
                [32, 17582916],
                [33, 17582917],
                [34, 17582918],
                [35, 17582919],
                [36, 17582920],
                [37, 17582921],
                [38, 17582922],
                [39, 17582923],
                [40, 17582924],
                [41, 17582925],
                [42, 17582926],
                [43, 17582927],
                [44, 17582928],
                [45, 17582929],
                [46, 17582930],
                [47, 17582931],
                [48, 17582932],
                [49, 17583028],
                [50, 17583029]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17581370],
                [1, 17582885],
                [2, 17582886],
                [3, 17582887],
                [4, 17582888],
                [5, 17582889],
                [6, 17582890],
                [7, 17582891],
                [8, 17582892],
                [9, 17582893],
                [10, 17582894],
                [11, 17582895],
                [12, 17582896],
                [13, 17582897],
                [14, 17582898],
                [15, 17582899],
                [16, 17582900],
                [17, 17582901],
                [18, 17582902],
                [19, 17582903],
                [20, 17582904],
                [21, 17582905],
                [22, 17582906],
                [23, 17582907],
                [24, 17582908],
                [25, 17582909],
                [26, 17582910],
                [27, 17582911],
                [28, 17582912],
                [29, 17582913],
                [30, 17582914],
                [31, 17582915],
                [32, 17582916],
                [33, 17582917],
                [34, 17582918],
                [35, 17582919],
                [36, 17582920],
                [37, 17582921],
                [38, 17582922],
                [39, 17582923],
                [40, 17582924],
                [41, 17582925],
                [42, 17582926],
                [43, 17582927],
                [44, 17582928],
                [45, 17582929],
                [46, 17582930],
                [47, 17582931],
                [48, 17582932],
                [49, 17583028],
                [50, 17583029]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [1, 17582911],
                [2, 17582912],
                [3, 17582913],
                [4, 17582914],
                [5, 17582915],
                [6, 17582916],
                [7, 17582917],
                [8, 17582918],
                [9, 17582919],
                [10, 17582920],
                [11, 17582921],
                [12, 17582922],
                [13, 17582923],
                [14, 17582924],
                [15, 17582925],
                [16, 17582926],
                [17, 17582927],
                [18, 17582928],
                [19, 17582929],
                [20, 17582930],
                [21, 17582931],
                [22, 17582932],
                [23, 17583028],
                [24, 17583029],
                [25, 17583030],
                [26, 17583031],
                [27, 17583032],
                [28, 17583033],
                [29, 17583034],
                [30, 17583035],
                [31, 17583036],
                [32, 17583037],
                [33, 17583038],
                [34, 17583039],
                [35, 17583040],
                [36, 17583041],
                [37, 17583042],
                [38, 17583043],
                [39, 17583044],
                [40, 17583045],
                [41, 17583046],
                [42, 17583047],
                [43, 17583048],
                [44, 17583049],
                [45, 17583050],
                [46, 17583051],
                [47, 17583052],
                [48, 17583053],
                [49, 17583054],
                [50, 17583055]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17588621],
                [1, 17588622],
                [2, 17588623],
                [3, 17588624],
                [4, 17588625],
                [5, 17588626],
                [6, 17588627],
                [7, 17588628],
                [8, 17588629],
                [9, 17588630],
                [10, 17588631],
                [11, 17588632],
                [12, 17588633],
                [13, 17588634],
                [14, 17588635],
                [15, 17588636],
                [16, 17588637],
                [17, 17588638],
                [18, 17588639],
                [19, 17588640],
                [20, 17588641],
                [21, 17588642],
                [22, 17588643],
                [23, 17588644],
                [24, 17588645],
                [25, 17588646],
                [26, 17588647],
                [27, 17588648],
                [28, 17588649],
                [29, 17588650],
                [30, 17588651],
                [31, 17588652],
                [32, 17588653],
                [33, 17588654],
                [34, 17588655],
                [35, 17588656],
                [36, 17588657],
                [37, 17589196],
                [38, 17589197],
                [39, 17589198],
                [40, 17589199],
                [41, 17589200],
                [42, 17589201],
                [43, 17589202],
                [44, 17589203],
                [45, 17589204],
                [46, 17589205],
                [47, 17589206],
                [48, 17589207],
                [49, 17589208]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17588621],
                [1, 17588622],
                [2, 17588623],
                [3, 17588624],
                [4, 17588625],
                [5, 17588626],
                [6, 17588627],
                [7, 17588628],
                [8, 17588629],
                [9, 17588630],
                [10, 17588631],
                [11, 17588632],
                [12, 17588633],
                [13, 17588634],
                [14, 17588635],
                [15, 17588636],
                [16, 17588637],
                [17, 17588638],
                [18, 17588639],
                [19, 17588640],
                [20, 17588641],
                [21, 17588642],
                [22, 17588643],
                [23, 17588644],
                [24, 17588645],
                [25, 17588646],
                [26, 17588647],
                [27, 17588648],
                [28, 17588649],
                [29, 17588650],
                [30, 17588651],
                [31, 17588652],
                [32, 17588653],
                [33, 17588654],
                [34, 17588655],
                [35, 17588656],
                [36, 17588657],
                [37, 17589196],
                [38, 17589197],
                [39, 17589198],
                [40, 17589199],
                [41, 17589200],
                [42, 17589201],
                [43, 17589202],
                [44, 17589203],
                [45, 17589204],
                [46, 17589205],
                [47, 17589206],
                [48, 17589207],
                [49, 17589208]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17588621],
                [1, 17588622],
                [2, 17588623],
                [3, 17588624],
                [4, 17588625],
                [5, 17588626],
                [6, 17588627],
                [7, 17588628],
                [8, 17588629],
                [9, 17588630],
                [10, 17588631],
                [11, 17588632],
                [12, 17588633],
                [13, 17588634],
                [14, 17588635],
                [15, 17588636],
                [16, 17588637],
                [17, 17588638],
                [18, 17588639],
                [19, 17588640],
                [20, 17588641],
                [21, 17588642],
                [22, 17588643],
                [23, 17588644],
                [24, 17588645],
                [25, 17588646],
                [26, 17588647],
                [27, 17588648],
                [28, 17588649],
                [29, 17588650],
                [30, 17588651],
                [31, 17588652],
                [32, 17588653],
                [33, 17588654],
                [34, 17588655],
                [35, 17588656],
                [36, 17588657],
                [37, 17589196],
                [38, 17589197],
                [39, 17589198],
                [40, 17589199],
                [41, 17589200],
                [42, 17589201],
                [43, 17589202],
                [44, 17589203],
                [45, 17589204],
                [46, 17589205],
                [47, 17589206],
                [48, 17589207],
                [49, 17589208]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [1, 17591770],
                [2, 17591771],
                [3, 17591772],
                [4, 17591773],
                [5, 17591774],
                [6, 17591775],
                [7, 17591776],
                [8, 17591777],
                [9, 17591778],
                [10, 17591779],
                [11, 17591780],
                [12, 17591781],
                [13, 17591782],
                [14, 17591783],
                [15, 17591784],
                [16, 17591785],
                [17, 17591786],
                [18, 17591787],
                [19, 17591788],
                [20, 17591789],
                [21, 17591790],
                [22, 17591791],
                [23, 17591792],
                [24, 17591793],
                [25, 17591794],
                [26, 17591796],
                [27, 17591797],
                [28, 17591798],
                [29, 17591799],
                [30, 17591800],
                [31, 17591801],
                [32, 17591802],
                [33, 17591803],
                [34, 17591804],
                [35, 17591805],
                [36, 17591806],
                [37, 17591807],
                [38, 17591808],
                [39, 17591809],
                [40, 17591810],
                [41, 17591811],
                [42, 17591812],
                [43, 17591813],
                [44, 17591814],
                [45, 17591815],
                [46, 17591816],
                [47, 17591817],
                [48, 17591818],
                [49, 17591819],
                [50, 17591820]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17593855],
                [1, 17593856],
                [2, 17593857],
                [3, 17593858],
                [4, 17593859],
                [5, 17593860],
                [6, 17593861],
                [7, 17593862],
                [8, 17593863],
                [9, 17593864],
                [10, 17593865],
                [11, 17593866],
                [12, 17593867],
                [13, 17593868],
                [14, 17593869],
                [15, 17593870],
                [16, 17593871],
                [17, 17593872],
                [18, 17593873],
                [19, 17593874],
                [20, 17593875],
                [21, 17593876],
                [22, 17593877],
                [23, 17593878],
                [24, 17593880],
                [25, 17593881],
                [26, 17593882],
                [27, 17593883],
                [28, 17593884],
                [29, 17593885],
                [30, 17593886],
                [31, 17593887],
                [32, 17593888],
                [33, 17593889],
                [34, 17593890],
                [35, 17593891],
                [36, 17593892],
                [37, 17593893],
                [38, 17593894],
                [39, 17593895],
                [40, 17593896],
                [41, 17593897],
                [42, 17593898],
                [43, 17593899],
                [44, 17593900],
                [45, 17593901],
                [46, 17593902],
                [47, 17593903]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 17593863],
                [1, 17593864],
                [2, 17593865],
                [3, 17593866],
                [4, 17593867],
                [5, 17593868],
                [6, 17593869],
                [7, 17593870],
                [8, 17593871],
                [9, 17593872],
                [10, 17593873],
                [11, 17593874],
                [12, 17593875],
                [13, 17593876],
                [14, 17593877],
                [15, 17593878],
                [16, 17593880],
                [17, 17593881],
                [18, 17593882],
                [19, 17593883],
                [20, 17593884],
                [21, 17593885],
                [22, 17593886],
                [23, 17593887],
                [24, 17593888],
                [25, 17593889],
                [26, 17593890],
                [27, 17593891],
                [28, 17593892],
                [29, 17593893],
                [30, 17593894],
                [31, 17593895],
                [32, 17593896],
                [33, 17593897],
                [34, 17593898],
                [35, 17593899],
                [36, 17593900],
                [37, 17593901],
                [38, 17593902],
                [39, 17593903],
                [40, 17593904],
                [41, 17593905],
                [42, 17593906],
                [43, 17593907]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [11, 17596476],
                [12, 17596477],
                [13, 17596478],
                [14, 17596479],
                [15, 17596480],
                [16, 17596481],
                [17, 17596482],
                [19, 17596483],
                [20, 17596484],
                [21, 17596485],
                [22, 17596486],
                [23, 17596487],
                [24, 17596488],
                [25, 17596489],
                [26, 17596490],
                [27, 17596491],
                [28, 17596492],
                [29, 17596493],
                [30, 17596494],
                [31, 17596495],
                [32, 17596496],
                [33, 17596497],
                [34, 17596498],
                [35, 17596499],
                [36, 17596500],
                [37, 17596501],
                [38, 17596502],
                [39, 17596503],
                [40, 17596504],
                [41, 17596505],
                [42, 17596506],
                [43, 17596507],
                [44, 17596508],
                [45, 17596509],
                [46, 17596510],
                [47, 17596511],
                [48, 17596512],
                [49, 17596513],
                [50, 17596514]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [5, 17624012],
                [6, 17624013],
                [7, 17624014],
                [8, 17624015],
                [9, 17624016],
                [10, 17624017],
                [11, 17624018],
                [12, 17624019],
                [13, 17624020],
                [14, 17625913],
                [15, 17625914],
                [16, 17625915],
                [17, 17625916],
                [18, 17625917],
                [19, 17625918],
                [20, 17625919],
                [21, 17625920],
                [22, 17625921],
                [23, 17625922],
                [24, 17625923],
                [25, 17625924],
                [26, 17625925],
                [27, 17625926],
                [28, 17625927],
                [29, 17625928],
                [30, 17625929],
                [31, 17625930],
                [32, 17625931],
                [33, 17625932],
                [34, 17625933],
                [35, 17625934],
                [36, 17625935],
                [37, 17625936],
                [38, 17625937],
                [39, 17625938],
                [40, 17625939],
                [41, 17625940],
                [42, 17625941],
                [43, 17625942],
                [44, 17625943],
                [45, 17625944],
                [46, 17625945],
                [47, 17625946],
                [48, 17625947],
                [49, 17625948],
                [50, 17625949]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [2, 17624012],
                [3, 17624013],
                [4, 17624014],
                [5, 17624015],
                [6, 17624016],
                [7, 17624017],
                [8, 17624018],
                [9, 17624019],
                [10, 17624020],
                [11, 17625913],
                [12, 17625914],
                [13, 17625915],
                [14, 17625916],
                [15, 17625917],
                [16, 17625918],
                [17, 17625919],
                [18, 17625920],
                [19, 17625921],
                [20, 17625922],
                [21, 17625923],
                [22, 17625924],
                [23, 17625925],
                [24, 17625926],
                [25, 17625927],
                [26, 17625928],
                [27, 17625929],
                [28, 17625930],
                [29, 17625931],
                [30, 17625932],
                [31, 17625933],
                [32, 17625934],
                [33, 17625935],
                [34, 17625936],
                [35, 17625937],
                [36, 17625938],
                [37, 17625939],
                [38, 17625940],
                [39, 17625941],
                [40, 17625942],
                [41, 17625943],
                [42, 17625944],
                [43, 17625945],
                [44, 17625946],
                [45, 17625947],
                [46, 17625948],
                [47, 17625949],
                [48, 17625950],
                [49, 17625951],
                [50, 17625952]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [1, 31796700],
                [2, 31796701],
                [3, 31796702],
                [4, 31796703],
                [5, 31796704],
                [6, 31796705],
                [7, 31796706],
                [8, 31796710],
                [9, 31796711],
                [10, 31796712],
                [11, 31796713],
                [12, 31796714],
                [13, 31796715],
                [14, 31796716],
                [15, 31796717],
                [16, 31796718],
                [17, 31796719],
                [18, 31796720],
                [19, 31796721],
                [20, 31796722],
                [21, 31796723],
                [22, 31796724],
                [23, 31796725],
                [24, 31796726],
                [25, 31796727],
                [26, 31796728],
                [27, 31799014],
                [28, 31799015],
                [29, 31799016],
                [30, 31799017],
                [31, 31799018],
                [32, 31799019],
                [33, 31799020],
                [34, 31799021],
                [35, 31799022],
                [36, 31799023],
                [37, 31799024],
                [38, 31799025],
                [39, 31799026],
                [40, 31799027],
                [41, 31799028],
                [42, 31799029],
                [43, 31799030],
                [44, 31799031],
                [45, 31799032],
                [46, 31799033],
                [47, 31799034],
                [48, 31799035],
                [49, 31799036],
                [50, 31799037]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [0, 36722692],
                [1, 36722693],
                [2, 36722694],
                [3, 36722695],
                [4, 36722696],
                [5, 36722697],
                [6, 36722698],
                [7, 36722699],
                [8, 36722700],
                [9, 36722701],
                [10, 36722702],
                [11, 36722703],
                [12, 36722704],
                [13, 36722705],
                [14, 36723505],
                [15, 36723506],
                [16, 36723507],
                [17, 36723508],
                [18, 36723509],
                [19, 36723510],
                [20, 36723511],
                [21, 36723512],
                [22, 36723513],
                [23, 36723514],
                [24, 36723515],
                [25, 36723516],
                [26, 36723517],
                [27, 36723518],
                [28, 36723519],
                [29, 36723520],
                [30, 36723521],
                [31, 36723522],
                [32, 36723523],
                [33, 36723524],
                [34, 36723525],
                [35, 36723526],
                [36, 36723527],
                [37, 36723528],
                [38, 36723529],
                [39, 36723530],
                [40, 36723531],
                [41, 36723532],
                [42, 36737414],
                [43, 36737415],
                [44, 36737416],
                [45, 36737417],
                [46, 36737418],
                [47, 36737419],
                [48, 36737420]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs().collect();
        assert_eq!(
            pairs,
            vec![
                [4, 44587963],
                [5, 44587964],
                [6, 44587965],
                [7, 44587966],
                [8, 44587967],
                [9, 44587968],
                [10, 44587969],
                [11, 44587970],
                [12, 44587971],
                [13, 44587972],
                [14, 44587973],
                [15, 44587974],
                [16, 44587975],
                [17, 44587976],
                [18, 44587977],
                [19, 44587978],
                [20, 44587979],
                [21, 44587980],
                [22, 44587981],
                [23, 44587982],
                [24, 44587983],
                [25, 44589680],
                [26, 44589681],
                [27, 44589682],
                [28, 44589683],
                [29, 44589684],
                [30, 44589685],
                [31, 44589686],
                [32, 44589687],
                [33, 44589688],
                [34, 44589689],
                [35, 44589690],
                [36, 44589691],
                [37, 44589692],
                [38, 44589693],
                [39, 44589694],
                [40, 44589695],
                [41, 44589696],
                [42, 44589697],
                [43, 44589698],
                [44, 44589699],
                [45, 44589700],
                [46, 44589701],
                [47, 44589702],
                [48, 44592034],
                [49, 44592035],
                [50, 44592036]
            ]
        );
    }

    #[test]
    fn test_aligned_pairs_full() {
        let mut bam = bam::Reader::from_path("./test/test_spliced_reads.bam").unwrap();
        let mut it = bam.records();

        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(
            pairs,
            vec![
                [Some(0), None],
                [Some(1), None],
                [Some(2), None],
                [Some(3), None],
                [Some(4), None],
                [Some(5), None],
                [Some(6), Some(16050676)],
                [Some(7), Some(16050677)],
                [Some(8), Some(16050678)],
                [Some(9), Some(16050679)],
                [Some(10), Some(16050680)],
                [Some(11), Some(16050681)],
                [Some(12), Some(16050682)],
                [Some(13), Some(16050683)],
                [Some(14), Some(16050684)],
                [Some(15), Some(16050685)],
                [Some(16), Some(16050686)],
                [Some(17), Some(16050687)],
                [Some(18), Some(16050688)],
                [Some(19), Some(16050689)],
                [Some(20), Some(16050690)],
                [Some(21), Some(16050691)],
                [Some(22), Some(16050692)],
                [Some(23), Some(16050693)],
                [Some(24), Some(16050694)],
                [Some(25), Some(16050695)],
                [Some(26), Some(16050696)],
                [Some(27), Some(16050697)],
                [Some(28), Some(16050698)],
                [Some(29), Some(16050699)],
                [Some(30), Some(16050700)],
                [Some(31), Some(16050701)],
                [Some(32), Some(16050702)],
                [Some(33), Some(16050703)],
                [Some(34), Some(16050704)],
                [Some(35), Some(16050705)],
                [Some(36), Some(16050706)],
                [Some(37), Some(16050707)],
                [Some(38), Some(16050708)],
                [Some(39), Some(16050709)],
                [Some(40), Some(16050710)],
                [Some(41), Some(16050711)],
                [Some(42), Some(16050712)],
                [Some(43), Some(16050713)],
                [Some(44), Some(16050714)],
                [Some(45), Some(16050715)],
                [Some(46), Some(16050716)],
                [Some(47), Some(16050717)],
                [Some(48), Some(16050718)],
                [Some(49), Some(16050719)],
                [Some(50), Some(16050720)]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(
            pairs,
            vec![
                [Some(0), Some(16096878)],
                [Some(1), Some(16096879)],
                [Some(2), Some(16096880)],
                [Some(3), Some(16096881)],
                [Some(4), Some(16096882)],
                [Some(5), Some(16096883)],
                [Some(6), Some(16096884)],
                [None, Some(16096885)],
                [None, Some(16096886)],
                [Some(7), Some(16096887)],
                [Some(8), Some(16096888)],
                [Some(9), Some(16096889)],
                [Some(10), Some(16096890)],
                [Some(11), Some(16096891)],
                [Some(12), Some(16096892)],
                [Some(13), Some(16096893)],
                [Some(14), Some(16096894)],
                [Some(15), Some(16096895)],
                [Some(16), Some(16096896)],
                [Some(17), Some(16096897)],
                [Some(18), Some(16096898)],
                [Some(19), Some(16096899)],
                [Some(20), Some(16096900)],
                [Some(21), Some(16096901)],
                [Some(22), Some(16096902)],
                [Some(23), Some(16096903)],
                [Some(24), Some(16096904)],
                [Some(25), Some(16096905)],
                [Some(26), Some(16096906)],
                [Some(27), Some(16096907)],
                [Some(28), Some(16096908)],
                [Some(29), Some(16096909)],
                [Some(30), Some(16096910)],
                [Some(31), Some(16096911)],
                [Some(32), Some(16096912)],
                [Some(33), Some(16096913)],
                [Some(34), Some(16096914)],
                [Some(35), Some(16096915)],
                [Some(36), Some(16096916)],
                [Some(37), Some(16096917)],
                [Some(38), Some(16096918)],
                [Some(39), Some(16096919)],
                [Some(40), Some(16096920)],
                [Some(41), Some(16096921)],
                [Some(42), Some(16096922)],
                [Some(43), Some(16096923)],
                [Some(44), Some(16096924)],
                [Some(45), Some(16096925)],
                [Some(46), Some(16096926)],
                [Some(47), Some(16096927)],
                [Some(48), Some(16096928)],
                [Some(49), Some(16096929)],
                [Some(50), Some(16096930)]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(
            pairs,
            vec![
                [Some(0), Some(16097145)],
                [Some(1), Some(16097146)],
                [Some(2), Some(16097147)],
                [Some(3), Some(16097148)],
                [Some(4), Some(16097149)],
                [Some(5), Some(16097150)],
                [Some(6), Some(16097151)],
                [Some(7), Some(16097152)],
                [Some(8), Some(16097153)],
                [Some(9), Some(16097154)],
                [Some(10), Some(16097155)],
                [Some(11), Some(16097156)],
                [Some(12), Some(16097157)],
                [Some(13), Some(16097158)],
                [Some(14), Some(16097159)],
                [Some(15), Some(16097160)],
                [Some(16), Some(16097161)],
                [Some(17), Some(16097162)],
                [Some(18), Some(16097163)],
                [Some(19), Some(16097164)],
                [Some(20), Some(16097165)],
                [Some(21), Some(16097166)],
                [Some(22), Some(16097167)],
                [Some(23), Some(16097168)],
                [Some(24), Some(16097169)],
                [Some(25), Some(16097170)],
                [Some(26), Some(16097171)],
                [Some(27), Some(16097172)],
                [Some(28), Some(16097173)],
                [None, Some(16097174)],
                [None, Some(16097175)],
                [Some(29), Some(16097176)],
                [Some(30), Some(16097177)],
                [Some(31), Some(16097178)],
                [Some(32), Some(16097179)],
                [Some(33), Some(16097180)],
                [Some(34), Some(16097181)],
                [Some(35), Some(16097182)],
                [Some(36), Some(16097183)],
                [Some(37), Some(16097184)],
                [Some(38), Some(16097185)],
                [Some(39), Some(16097186)],
                [Some(40), Some(16097187)],
                [Some(41), Some(16097188)],
                [Some(42), Some(16097189)],
                [Some(43), Some(16097190)],
                [Some(44), Some(16097191)],
                [Some(45), Some(16097192)],
                [Some(46), Some(16097193)],
                [Some(47), Some(16097194)],
                [Some(48), Some(16097195)],
                [Some(49), Some(16097196)],
                [Some(50), Some(16097197)]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(
            pairs,
            vec![
                [Some(0), Some(16117350)],
                [Some(1), Some(16117351)],
                [Some(2), Some(16117352)],
                [Some(3), Some(16117353)],
                [Some(4), Some(16117354)],
                [Some(5), Some(16117355)],
                [Some(6), Some(16117356)],
                [Some(7), Some(16117357)],
                [Some(8), Some(16117358)],
                [Some(9), Some(16117359)],
                [Some(10), Some(16117360)],
                [Some(11), Some(16117361)],
                [Some(12), Some(16117362)],
                [Some(13), Some(16117363)],
                [Some(14), Some(16117364)],
                [Some(15), Some(16117365)],
                [Some(16), Some(16117366)],
                [Some(17), Some(16117367)],
                [Some(18), Some(16117368)],
                [Some(19), Some(16117369)],
                [Some(20), Some(16117370)],
                [Some(21), Some(16117371)],
                [Some(22), Some(16117372)],
                [Some(23), Some(16117373)],
                [Some(24), Some(16117374)],
                [Some(25), Some(16117375)],
                [Some(26), Some(16117376)],
                [Some(27), Some(16117377)],
                [Some(28), Some(16117378)],
                [Some(29), Some(16117379)],
                [Some(30), Some(16117380)],
                [Some(31), Some(16117381)],
                [Some(32), Some(16117382)],
                [Some(33), Some(16117383)],
                [Some(34), Some(16117384)],
                [Some(35), Some(16117385)],
                [Some(36), Some(16117386)],
                [Some(37), Some(16117387)],
                [Some(38), Some(16117388)],
                [Some(39), Some(16117389)],
                [Some(40), Some(16117390)],
                [Some(41), Some(16117391)],
                [Some(42), Some(16117392)],
                [Some(43), Some(16117393)],
                [Some(44), Some(16117394)],
                [Some(45), Some(16117395)],
                [Some(46), Some(16117396)],
                [Some(47), Some(16117397)],
                [Some(48), Some(16117398)],
                [Some(49), Some(16117399)],
                [Some(50), Some(16117400)]
            ]
        );
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(
            pairs,
            vec![
                [Some(0), Some(16118483)],
                [Some(1), Some(16118484)],
                [Some(2), Some(16118485)],
                [Some(3), Some(16118486)],
                [Some(4), Some(16118487)],
                [Some(5), Some(16118488)],
                [Some(6), Some(16118489)],
                [Some(7), Some(16118490)],
                [Some(8), Some(16118491)],
                [Some(9), Some(16118492)],
                [Some(10), Some(16118493)],
                [Some(11), Some(16118494)],
                [Some(12), Some(16118495)],
                [Some(13), Some(16118496)],
                [Some(14), Some(16118497)],
                [Some(15), Some(16118498)],
                [Some(16), Some(16118499)],
                [Some(17), Some(16118500)],
                [Some(18), Some(16118501)],
                [Some(19), Some(16118502)],
                [Some(20), Some(16118503)],
                [Some(21), Some(16118504)],
                [Some(22), Some(16118505)],
                [Some(23), Some(16118506)],
                [Some(24), Some(16118507)],
                [Some(25), Some(16118508)],
                [Some(26), Some(16118509)],
                [Some(27), Some(16118510)],
                [Some(28), Some(16118511)],
                [Some(29), Some(16118512)],
                [Some(30), Some(16118513)],
                [Some(31), Some(16118514)],
                [Some(32), Some(16118515)],
                [Some(33), Some(16118516)],
                [Some(34), Some(16118517)],
                [Some(35), Some(16118518)],
                [Some(36), Some(16118519)],
                [Some(37), Some(16118520)],
                [Some(38), Some(16118521)],
                [Some(39), Some(16118522)],
                [Some(40), Some(16118523)],
                [Some(41), Some(16118524)],
                [Some(42), Some(16118525)],
                [Some(43), Some(16118526)],
                [Some(44), Some(16118527)],
                [Some(45), Some(16118528)],
                [Some(46), Some(16118529)],
                [Some(47), Some(16118530)],
                [Some(48), Some(16118531)],
                [Some(49), Some(16118532)],
                [Some(50), Some(16118533)]
            ]
        );

        //
        //for the rest, we just verify that they have the expected amount of None in each position
        fn some_count(pairs: &[[Option<i64>; 2]], pos: usize) -> i64 {
            pairs.iter().filter(|x| x[pos].is_some()).count() as i64
        }
        fn none_count(pairs: &[[Option<i64>; 2]], pos: usize) -> i64 {
            pairs.iter().filter(|x| x[pos].is_none()).count() as i64
        }

        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 51);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 51);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 51);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 51);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 45);
        assert_eq!(none_count(&pairs, 1), 6);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 41);
        assert_eq!(none_count(&pairs, 1), 10);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 51);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 51);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 51);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 51);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 42);
        assert_eq!(none_count(&pairs, 1), 9);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 51);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 50);
        assert_eq!(none_count(&pairs, 1), 1);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 45);
        assert_eq!(none_count(&pairs, 1), 6);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 48);
        assert_eq!(none_count(&pairs, 1), 3);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 42);
        assert_eq!(none_count(&pairs, 1), 9);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 38);
        assert_eq!(none_count(&pairs, 1), 13);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 44);
        assert_eq!(none_count(&pairs, 1), 7);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 46);
        assert_eq!(none_count(&pairs, 1), 5);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 47);
        assert_eq!(none_count(&pairs, 1), 4);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 50);
        assert_eq!(none_count(&pairs, 1), 1);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 2183);
        assert_eq!(some_count(&pairs, 1), 2234);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 2183);
        assert_eq!(some_count(&pairs, 1), 2234);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 2183);
        assert_eq!(some_count(&pairs, 1), 2234);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 33);
        assert_eq!(none_count(&pairs, 1), 18);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 48);
        assert_eq!(none_count(&pairs, 1), 3);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 45);
        assert_eq!(none_count(&pairs, 1), 6);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 11832);
        assert_eq!(some_count(&pairs, 1), 11883);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 12542);
        assert_eq!(some_count(&pairs, 1), 12593);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 12542);
        assert_eq!(some_count(&pairs, 1), 12593);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 12542);
        assert_eq!(some_count(&pairs, 1), 12593);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 12542);
        assert_eq!(some_count(&pairs, 1), 12593);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 12542);
        assert_eq!(some_count(&pairs, 1), 12593);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 12542);
        assert_eq!(some_count(&pairs, 1), 12593);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 12542);
        assert_eq!(some_count(&pairs, 1), 12593);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 12542);
        assert_eq!(some_count(&pairs, 1), 12593);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 1467);
        assert_eq!(some_count(&pairs, 1), 1517);
        assert_eq!(none_count(&pairs, 1), 1);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 1609);
        assert_eq!(some_count(&pairs, 1), 1660);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 1609);
        assert_eq!(some_count(&pairs, 1), 1660);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 1609);
        assert_eq!(some_count(&pairs, 1), 1660);
        assert_eq!(none_count(&pairs, 1), 0);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 95);
        assert_eq!(some_count(&pairs, 1), 145);
        assert_eq!(none_count(&pairs, 1), 1);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 538);
        assert_eq!(some_count(&pairs, 1), 588);
        assert_eq!(none_count(&pairs, 1), 1);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 538);
        assert_eq!(some_count(&pairs, 1), 588);
        assert_eq!(none_count(&pairs, 1), 1);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 538);
        assert_eq!(some_count(&pairs, 1), 588);
        assert_eq!(none_count(&pairs, 1), 1);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 1);
        assert_eq!(some_count(&pairs, 1), 51);
        assert_eq!(none_count(&pairs, 1), 1);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 1);
        assert_eq!(some_count(&pairs, 1), 49);
        assert_eq!(none_count(&pairs, 1), 3);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 1);
        assert_eq!(some_count(&pairs, 1), 45);
        assert_eq!(none_count(&pairs, 1), 7);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 0);
        assert_eq!(some_count(&pairs, 1), 39);
        assert_eq!(none_count(&pairs, 1), 12);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 1892);
        assert_eq!(some_count(&pairs, 1), 1938);
        assert_eq!(none_count(&pairs, 1), 5);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 1892);
        assert_eq!(some_count(&pairs, 1), 1941);
        assert_eq!(none_count(&pairs, 1), 2);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 2288);
        assert_eq!(some_count(&pairs, 1), 2338);
        assert_eq!(none_count(&pairs, 1), 1);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 14680);
        assert_eq!(some_count(&pairs, 1), 14729);
        assert_eq!(none_count(&pairs, 1), 2);
        let pairs: Vec<_> = it.next().unwrap().unwrap().aligned_pairs_full().collect();
        assert_eq!(some_count(&pairs, 0), 51);
        assert_eq!(none_count(&pairs, 0), 4027);
        assert_eq!(some_count(&pairs, 1), 4074);
        assert_eq!(none_count(&pairs, 1), 4);
    }

    #[test]
    fn test_aligned_block_pairs() {
        let mut bam = bam::Reader::from_path("./test/test_spliced_reads.bam").unwrap();
        let mut it = bam.records();

        let read = it.next().unwrap().unwrap();
        let pairs: Vec<_> = read.aligned_pairs().collect();
        let block_pairs: Vec<_> = read.aligned_block_pairs().collect();

        //first coordinates identical
        assert_eq!(pairs[0][0], block_pairs[0].0[0]); //read
        assert_eq!(pairs[0][1], block_pairs[0].1[0]); // genomic

        //end coordinates are + 1, so the ranges are the same...
        assert_eq!(
            pairs[pairs.len() - 1][0],
            block_pairs[block_pairs.len() - 1].0[1] - 1
        );
        assert_eq!(
            pairs[pairs.len() - 1][1],
            block_pairs[block_pairs.len() - 1].1[1] - 1
        );

        //let's see if they're really identical
        for read in it {
            let read = read.unwrap();
            let pairs: Vec<_> = read.aligned_pairs().collect();
            let block_pairs: Vec<_> = read.aligned_block_pairs().collect();
            let mut ii = 0;
            for ([read_start, read_stop], [genome_start, genome_stop]) in block_pairs {
                assert_eq!(read_stop - read_start, genome_stop - genome_start);
                for (read_pos, genome_pos) in (read_start..read_stop).zip(genome_start..genome_stop)
                {
                    assert_eq!(pairs[ii][0], read_pos);
                    assert_eq!(pairs[ii][1], genome_pos);
                    ii += 1;
                }
            }
        }
    }

    #[test]
    fn test_get_cigar_stats() {
        let mut bam = bam::Reader::from_path("./test/test_spliced_reads.bam").unwrap();
        let mut it = bam.records();

        fn to_arr(hm: HashMap<Cigar, i32>) -> [i32; 9] {
            [
                *hm.get(&Cigar::Match(0)).unwrap(),
                *hm.get(&Cigar::Ins(0)).unwrap(),
                *hm.get(&Cigar::Del(0)).unwrap(),
                *hm.get(&Cigar::RefSkip(0)).unwrap(),
                *hm.get(&Cigar::SoftClip(0)).unwrap(),
                *hm.get(&Cigar::HardClip(0)).unwrap(),
                *hm.get(&Cigar::Pad(0)).unwrap(),
                *hm.get(&Cigar::Equal(0)).unwrap(),
                *hm.get(&Cigar::Diff(0)).unwrap(),
            ]
        }

        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [45, 0, 0, 0, 6, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 2, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 1, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 2, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 1, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [45, 0, 0, 0, 6, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [41, 0, 0, 0, 10, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [42, 0, 0, 0, 9, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [50, 0, 0, 0, 1, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [45, 0, 0, 0, 6, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [48, 0, 0, 0, 3, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [42, 0, 0, 0, 9, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [38, 0, 0, 0, 13, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [44, 0, 0, 0, 7, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [46, 0, 0, 0, 5, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [47, 4, 0, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 1, 0, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [50, 1, 0, 0, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 1, 0, 0, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 2183, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 0, 1, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 2183, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 0, 1, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 2183, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 0, 1, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [33, 0, 0, 0, 18, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 2, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [48, 0, 0, 0, 3, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 2, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [45, 0, 0, 0, 6, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [1, 0, 0, 0, 2, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 11832, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 0, 1, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 12542, 0, 0, 0, 0, 0],);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [3, 0, 0, 2, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 12542, 0, 0, 0, 0, 0],);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [3, 0, 0, 2, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 12542, 0, 0, 0, 0, 0],);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [3, 0, 0, 2, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 12542, 0, 0, 0, 0, 0],);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [3, 0, 0, 2, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 12542, 0, 0, 0, 0, 0],);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [3, 0, 0, 2, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 12542, 0, 0, 0, 0, 0],);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [3, 0, 0, 2, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 12542, 0, 0, 0, 0, 0],);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [3, 0, 0, 2, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 12542, 0, 0, 0, 0, 0],);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [3, 0, 0, 2, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [50, 0, 0, 1467, 1, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 0, 1, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 1609, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [3, 0, 0, 2, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 1609, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [3, 0, 0, 2, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [51, 0, 0, 1609, 0, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [3, 0, 0, 2, 0, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [50, 0, 0, 95, 1, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 0, 1, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [50, 0, 0, 538, 1, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 0, 1, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [50, 0, 0, 538, 1, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 0, 1, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [50, 0, 0, 538, 1, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 0, 1, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [50, 0, 1, 0, 1, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 1, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [48, 0, 1, 0, 3, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 1, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [44, 0, 1, 0, 7, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 1, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [39, 1, 0, 0, 11, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 1, 0, 0, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [46, 0, 0, 1892, 5, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 0, 1, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [49, 0, 0, 1892, 2, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [2, 0, 0, 1, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [50, 0, 3, 2285, 1, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [3, 0, 1, 1, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [49, 0, 0, 14680, 2, 0, 0, 0, 0],);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [3, 0, 0, 2, 1, 0, 0, 0, 0]);
        let read = it.next().unwrap().unwrap();
        let cigar_nucleotides = read.cigar_stats_nucleotides();
        assert_eq!(to_arr(cigar_nucleotides), [47, 0, 0, 4027, 4, 0, 0, 0, 0]);
        let cigar_blocks = read.cigar_stats_blocks();
        assert_eq!(to_arr(cigar_blocks), [3, 0, 0, 2, 1, 0, 0, 0, 0]);
    }

    #[test]
    fn test_reference_positions() {
        let mut bam = bam::Reader::from_path("./test/test_spliced_reads.bam").unwrap();
        let mut it = bam.records();

        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16050676, 16050677, 16050678, 16050679, 16050680, 16050681, 16050682, 16050683,
                16050684, 16050685, 16050686, 16050687, 16050688, 16050689, 16050690, 16050691,
                16050692, 16050693, 16050694, 16050695, 16050696, 16050697, 16050698, 16050699,
                16050700, 16050701, 16050702, 16050703, 16050704, 16050705, 16050706, 16050707,
                16050708, 16050709, 16050710, 16050711, 16050712, 16050713, 16050714, 16050715,
                16050716, 16050717, 16050718, 16050719, 16050720
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16096878, 16096879, 16096880, 16096881, 16096882, 16096883, 16096884, 16096887,
                16096888, 16096889, 16096890, 16096891, 16096892, 16096893, 16096894, 16096895,
                16096896, 16096897, 16096898, 16096899, 16096900, 16096901, 16096902, 16096903,
                16096904, 16096905, 16096906, 16096907, 16096908, 16096909, 16096910, 16096911,
                16096912, 16096913, 16096914, 16096915, 16096916, 16096917, 16096918, 16096919,
                16096920, 16096921, 16096922, 16096923, 16096924, 16096925, 16096926, 16096927,
                16096928, 16096929, 16096930
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16097145, 16097146, 16097147, 16097148, 16097149, 16097150, 16097151, 16097152,
                16097153, 16097154, 16097155, 16097156, 16097157, 16097158, 16097159, 16097160,
                16097161, 16097162, 16097163, 16097164, 16097165, 16097166, 16097167, 16097168,
                16097169, 16097170, 16097171, 16097172, 16097173, 16097176, 16097177, 16097178,
                16097179, 16097180, 16097181, 16097182, 16097183, 16097184, 16097185, 16097186,
                16097187, 16097188, 16097189, 16097190, 16097191, 16097192, 16097193, 16097194,
                16097195, 16097196, 16097197
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16117350, 16117351, 16117352, 16117353, 16117354, 16117355, 16117356, 16117357,
                16117358, 16117359, 16117360, 16117361, 16117362, 16117363, 16117364, 16117365,
                16117366, 16117367, 16117368, 16117369, 16117370, 16117371, 16117372, 16117373,
                16117374, 16117375, 16117376, 16117377, 16117378, 16117379, 16117380, 16117381,
                16117382, 16117383, 16117384, 16117385, 16117386, 16117387, 16117388, 16117389,
                16117390, 16117391, 16117392, 16117393, 16117394, 16117395, 16117396, 16117397,
                16117398, 16117399, 16117400
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16118483, 16118484, 16118485, 16118486, 16118487, 16118488, 16118489, 16118490,
                16118491, 16118492, 16118493, 16118494, 16118495, 16118496, 16118497, 16118498,
                16118499, 16118500, 16118501, 16118502, 16118503, 16118504, 16118505, 16118506,
                16118507, 16118508, 16118509, 16118510, 16118511, 16118512, 16118513, 16118514,
                16118515, 16118516, 16118517, 16118518, 16118519, 16118520, 16118521, 16118522,
                16118523, 16118524, 16118525, 16118526, 16118527, 16118528, 16118529, 16118530,
                16118531, 16118532, 16118533
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16118499, 16118500, 16118501, 16118502, 16118503, 16118504, 16118505, 16118506,
                16118507, 16118508, 16118509, 16118510, 16118511, 16118512, 16118513, 16118514,
                16118515, 16118516, 16118517, 16118518, 16118519, 16118520, 16118521, 16118522,
                16118523, 16118524, 16118525, 16118526, 16118527, 16118528, 16118529, 16118530,
                16118531, 16118532, 16118533, 16118534, 16118535, 16118536, 16118537, 16118538,
                16118539, 16118540, 16118541, 16118542, 16118543, 16118544, 16118545, 16118546,
                16118547, 16118548, 16118549
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16118499, 16118500, 16118501, 16118502, 16118503, 16118504, 16118505, 16118506,
                16118507, 16118508, 16118509, 16118510, 16118511, 16118512, 16118513, 16118514,
                16118515, 16118516, 16118517, 16118518, 16118519, 16118520, 16118521, 16118522,
                16118523, 16118524, 16118525, 16118526, 16118527, 16118528, 16118529, 16118530,
                16118531, 16118532, 16118533, 16118534, 16118535, 16118536, 16118537, 16118538,
                16118539, 16118540, 16118541, 16118542, 16118543, 16118544, 16118545, 16118546,
                16118547, 16118548, 16118549
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16118499, 16118500, 16118501, 16118502, 16118503, 16118504, 16118505, 16118506,
                16118507, 16118508, 16118509, 16118510, 16118511, 16118512, 16118513, 16118514,
                16118515, 16118516, 16118517, 16118518, 16118519, 16118520, 16118521, 16118522,
                16118523, 16118524, 16118525, 16118526, 16118527, 16118528, 16118529, 16118530,
                16118531, 16118532, 16118533, 16118534, 16118535, 16118536, 16118537, 16118538,
                16118539, 16118540, 16118541, 16118542, 16118543, 16118544, 16118545, 16118546,
                16118547, 16118548, 16118549
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16123411, 16123412, 16123413, 16123414, 16123415, 16123416, 16123417, 16123418,
                16123419, 16123420, 16123421, 16123422, 16123423, 16123424, 16123425, 16123426,
                16123427, 16123428, 16123429, 16123430, 16123431, 16123432, 16123433, 16123434,
                16123435, 16123436, 16123437, 16123438, 16123439, 16123440, 16123441, 16123442,
                16123443, 16123444, 16123445, 16123446, 16123447, 16123448, 16123449, 16123450,
                16123451, 16123452, 16123453, 16123454, 16123455, 16123456, 16123457, 16123458,
                16123459, 16123460, 16123461
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16123417, 16123418, 16123419, 16123420, 16123421, 16123422, 16123423, 16123424,
                16123425, 16123426, 16123427, 16123428, 16123429, 16123430, 16123431, 16123432,
                16123433, 16123434, 16123435, 16123436, 16123437, 16123438, 16123439, 16123440,
                16123441, 16123442, 16123443, 16123444, 16123445, 16123446, 16123447, 16123448,
                16123449, 16123450, 16123451, 16123452, 16123453, 16123454, 16123455, 16123456,
                16123457, 16123458, 16123459, 16123460, 16123461
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16165860, 16165861, 16165862, 16165863, 16165864, 16165865, 16165866, 16165867,
                16165868, 16165869, 16165870, 16165871, 16165872, 16165873, 16165874, 16165875,
                16165876, 16165877, 16165878, 16165879, 16165880, 16165881, 16165882, 16165883,
                16165884, 16165885, 16165886, 16165887, 16165888, 16165889, 16165890, 16165891,
                16165892, 16165893, 16165894, 16165895, 16165896, 16165897, 16165898, 16165899,
                16165900
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16180871, 16180872, 16180873, 16180874, 16180875, 16180876, 16180877, 16180878,
                16180879, 16180880, 16180881, 16180882, 16180883, 16180884, 16180885, 16180886,
                16180887, 16180888, 16180889, 16180890, 16180891, 16180892, 16180893, 16180894,
                16180895, 16180896, 16180897, 16180898, 16180899, 16180900, 16180901, 16180902,
                16180903, 16180904, 16180905, 16180906, 16180907, 16180908, 16180909, 16180910,
                16180911, 16180912, 16180913, 16180914, 16180915, 16180916, 16180917, 16180918,
                16180919, 16180920, 16180921
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16189705, 16189706, 16189707, 16189708, 16189709, 16189710, 16189711, 16189712,
                16189713, 16189714, 16189715, 16189716, 16189717, 16189718, 16189719, 16189720,
                16189721, 16189722, 16189723, 16189724, 16189725, 16189726, 16189727, 16189728,
                16189729, 16189730, 16189731, 16189732, 16189733, 16189734, 16189735, 16189736,
                16189737, 16189738, 16189739, 16189740, 16189741, 16189742, 16189743, 16189744,
                16189745, 16189746, 16189747, 16189748, 16189749, 16189750, 16189751, 16189752,
                16189753, 16189754, 16189755
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16231271, 16231272, 16231273, 16231274, 16231275, 16231276, 16231277, 16231278,
                16231279, 16231280, 16231281, 16231282, 16231283, 16231284, 16231285, 16231286,
                16231287, 16231288, 16231289, 16231290, 16231291, 16231292, 16231293, 16231294,
                16231295, 16231296, 16231297, 16231298, 16231299, 16231300, 16231301, 16231302,
                16231303, 16231304, 16231305, 16231306, 16231307, 16231308, 16231309, 16231310,
                16231311, 16231312, 16231313, 16231314, 16231315, 16231316, 16231317, 16231318,
                16231319, 16231320, 16231321
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16237657, 16237658, 16237659, 16237660, 16237661, 16237662, 16237663, 16237664,
                16237665, 16237666, 16237667, 16237668, 16237669, 16237670, 16237671, 16237672,
                16237673, 16237674, 16237675, 16237676, 16237677, 16237678, 16237679, 16237680,
                16237681, 16237682, 16237683, 16237684, 16237685, 16237686, 16237687, 16237688,
                16237689, 16237690, 16237691, 16237692, 16237693, 16237694, 16237695, 16237696,
                16237697, 16237698, 16237699, 16237700, 16237701, 16237702, 16237703, 16237704,
                16237705, 16237706, 16237707
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16255012, 16255013, 16255014, 16255015, 16255016, 16255017, 16255018, 16255019,
                16255020, 16255021, 16255022, 16255023, 16255024, 16255025, 16255026, 16255027,
                16255028, 16255029, 16255030, 16255031, 16255032, 16255033, 16255034, 16255035,
                16255036, 16255037, 16255038, 16255039, 16255040, 16255041, 16255042, 16255043,
                16255044, 16255045, 16255046, 16255047, 16255048, 16255049, 16255050, 16255051,
                16255052, 16255053
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16255391, 16255392, 16255393, 16255394, 16255395, 16255396, 16255397, 16255398,
                16255399, 16255400, 16255401, 16255402, 16255403, 16255404, 16255405, 16255406,
                16255407, 16255408, 16255409, 16255410, 16255411, 16255412, 16255413, 16255414,
                16255415, 16255416, 16255417, 16255418, 16255419, 16255420, 16255421, 16255422,
                16255423, 16255424, 16255425, 16255426, 16255427, 16255428, 16255429, 16255430,
                16255431, 16255432, 16255433, 16255434, 16255435, 16255436, 16255437, 16255438,
                16255439, 16255440, 16255441
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16255392, 16255393, 16255394, 16255395, 16255396, 16255397, 16255398, 16255399,
                16255400, 16255401, 16255402, 16255403, 16255404, 16255405, 16255406, 16255407,
                16255408, 16255409, 16255410, 16255411, 16255412, 16255413, 16255414, 16255415,
                16255416, 16255417, 16255418, 16255419, 16255420, 16255421, 16255422, 16255423,
                16255424, 16255425, 16255426, 16255427, 16255428, 16255429, 16255430, 16255431,
                16255432, 16255433, 16255434, 16255435, 16255436, 16255437, 16255438, 16255439,
                16255440, 16255441
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16256084, 16256085, 16256086, 16256087, 16256088, 16256089, 16256090, 16256091,
                16256092, 16256093, 16256094, 16256095, 16256096, 16256097, 16256098, 16256099,
                16256100, 16256101, 16256102, 16256103, 16256104, 16256105, 16256106, 16256107,
                16256108, 16256109, 16256110, 16256111, 16256112, 16256113, 16256114, 16256115,
                16256116, 16256117, 16256118, 16256119, 16256120, 16256121, 16256122, 16256123,
                16256124, 16256125, 16256126, 16256127, 16256128
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16256224, 16256225, 16256226, 16256227, 16256228, 16256229, 16256230, 16256231,
                16256232, 16256233, 16256234, 16256235, 16256236, 16256237, 16256238, 16256239,
                16256240, 16256241, 16256242, 16256243, 16256244, 16256245, 16256246, 16256247,
                16256248, 16256249, 16256250, 16256251, 16256252, 16256253, 16256254, 16256255,
                16256256, 16256257, 16256258, 16256259, 16256260, 16256261, 16256262, 16256263,
                16256264, 16256265, 16256266, 16256267, 16256268, 16256269, 16256270, 16256271
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16325199, 16325200, 16325201, 16325202, 16325203, 16325204, 16325205, 16325206,
                16325207, 16325208, 16325209, 16325210, 16325211, 16325212, 16325213, 16325214,
                16325215, 16325216, 16325217, 16325218, 16325219, 16325220, 16325221, 16325222,
                16325223, 16325224, 16325225, 16325226, 16325227, 16325228, 16325229, 16325230,
                16325231, 16325232, 16325233, 16325234, 16325235, 16325236, 16325237, 16325238,
                16325239, 16325240
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16352865, 16352866, 16352867, 16352868, 16352869, 16352870, 16352871, 16352872,
                16352873, 16352874, 16352875, 16352876, 16352877, 16352878, 16352879, 16352880,
                16352881, 16352882, 16352883, 16352884, 16352885, 16352886, 16352887, 16352888,
                16352889, 16352890, 16352891, 16352892, 16352893, 16352894, 16352895, 16352896,
                16352897, 16352898, 16352899, 16352900, 16352901, 16352902
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16352968, 16352969, 16352970, 16352971, 16352972, 16352973, 16352974, 16352975,
                16352976, 16352977, 16352978, 16352979, 16352980, 16352981, 16352982, 16352983,
                16352984, 16352985, 16352986, 16352987, 16352988, 16352989, 16352990, 16352991,
                16352992, 16352993, 16352994, 16352995, 16352996, 16352997, 16352998, 16352999,
                16353000, 16353001, 16353002, 16353003, 16353004, 16353005, 16353006, 16353007,
                16353008, 16353009, 16353010, 16353011
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                16414998, 16414999, 16415000, 16415001, 16415002, 16415003, 16415004, 16415005,
                16415006, 16415007, 16415008, 16415009, 16415010, 16415011, 16415012, 16415013,
                16415014, 16415015, 16415016, 16415017, 16415018, 16415019, 16415020, 16415021,
                16415022, 16415023, 16415024, 16415025, 16415026, 16415027, 16415028, 16415029,
                16415030, 16415031, 16415032, 16415033, 16415034, 16415035, 16415036, 16415037,
                16415038, 16415039, 16415040, 16415041, 16415042, 16415043
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17031591, 17031592, 17031593, 17031594, 17031595, 17031596, 17031597, 17031598,
                17031599, 17031600, 17031601, 17031602, 17031603, 17031604, 17031605, 17031606,
                17031607, 17031608, 17031609, 17031610, 17031611, 17031612, 17031613, 17031614,
                17031615, 17031616, 17031617, 17031618, 17031619, 17031620, 17031621, 17031622,
                17031623, 17031624, 17031625, 17031626, 17031627, 17031628, 17031629, 17031630,
                17031631, 17031632, 17031633, 17031634, 17031635, 17031636, 17031637
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17057382, 17057383, 17057384, 17057385, 17057386, 17057387, 17057388, 17057389,
                17057390, 17057391, 17057392, 17057393, 17057394, 17057395, 17057396, 17057397,
                17057398, 17057399, 17057400, 17057401, 17057402, 17057403, 17057404, 17057405,
                17057406, 17057407, 17057408, 17057409, 17057410, 17057411, 17057412, 17057413,
                17057414, 17057415, 17057416, 17057417, 17057418, 17057419, 17057420, 17057421,
                17057422, 17057423, 17057424, 17057425, 17057426, 17057427, 17057428, 17057429,
                17057430, 17057431
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17092766, 17092767, 17092768, 17092769, 17092770, 17092771, 17092772, 17092773,
                17092774, 17092775, 17092776, 17092777, 17092778, 17092779, 17092780, 17092781,
                17092782, 17094966, 17094967, 17094968, 17094969, 17094970, 17094971, 17094972,
                17094973, 17094974, 17094975, 17094976, 17094977, 17094978, 17094979, 17094980,
                17094981, 17094982, 17094983, 17094984, 17094985, 17094986, 17094987, 17094988,
                17094989, 17094990, 17094991, 17094992, 17094993, 17094994, 17094995, 17094996,
                17094997, 17094998, 17094999
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17092782, 17094966, 17094967, 17094968, 17094969, 17094970, 17094971, 17094972,
                17094973, 17094974, 17094975, 17094976, 17094977, 17094978, 17094979, 17094980,
                17094981, 17094982, 17094983, 17094984, 17094985, 17094986, 17094987, 17094988,
                17094989, 17094990, 17094991, 17094992, 17094993, 17094994, 17094995, 17094996,
                17094997, 17094998, 17094999, 17095000, 17095001, 17095002, 17095003, 17095004,
                17095005, 17095006, 17095007, 17095008, 17095009, 17095010, 17095011, 17095012,
                17095013, 17095014, 17095015
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17092782, 17094966, 17094967, 17094968, 17094969, 17094970, 17094971, 17094972,
                17094973, 17094974, 17094975, 17094976, 17094977, 17094978, 17094979, 17094980,
                17094981, 17094982, 17094983, 17094984, 17094985, 17094986, 17094987, 17094988,
                17094989, 17094990, 17094991, 17094992, 17094993, 17094994, 17094995, 17094996,
                17094997, 17094998, 17094999, 17095000, 17095001, 17095002, 17095003, 17095004,
                17095005, 17095006, 17095007, 17095008, 17095009, 17095010, 17095011, 17095012,
                17095013, 17095014, 17095015
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17137287, 17137288, 17137289, 17137290, 17137291, 17137292, 17137293, 17137294,
                17137295, 17137296, 17137297, 17137298, 17137299, 17137300, 17137301, 17137302,
                17137303, 17137304, 17137305, 17137306, 17137307, 17137308, 17137309, 17137310,
                17137311, 17137312, 17137313, 17137314, 17137315, 17137316, 17137317, 17137318,
                17137319
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17306238, 17306239, 17306240, 17306241, 17306242, 17306243, 17306244, 17306245,
                17306246, 17306247, 17306248, 17306249, 17306250, 17306251, 17306252, 17306253,
                17306254, 17306255, 17306256, 17306257, 17306258, 17306259, 17306260, 17306261,
                17306262, 17306263, 17306264, 17306265, 17306266, 17306267, 17306268, 17306269,
                17306270, 17306271, 17306272, 17306273, 17306274, 17306275, 17306276, 17306277,
                17306278, 17306279, 17306280, 17306281, 17306282, 17306283, 17306284, 17306285
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17561868, 17561869, 17561870, 17561871, 17561872, 17561873, 17561874, 17561875,
                17561876, 17561877, 17561878, 17561879, 17561880, 17561881, 17561882, 17561883,
                17561884, 17561885, 17561886, 17561887, 17561888, 17561889, 17561890, 17561891,
                17561892, 17561893, 17561894, 17561895, 17561896, 17561897, 17561898, 17561899,
                17561900, 17561901, 17561902, 17561903, 17561904, 17561905, 17561906, 17561907,
                17561908, 17561909, 17561910, 17561911, 17561912
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17566078, 17566079, 17566080, 17566081, 17566082, 17566083, 17566084, 17566085,
                17566086, 17566087, 17566088, 17566089, 17566090, 17566091, 17566092, 17566093,
                17566094, 17566095, 17566096, 17566097, 17566098, 17566099, 17566100, 17566101,
                17566102, 17566103, 17566104, 17566105, 17566106, 17566107, 17566108, 17566109,
                17566110, 17566111, 17566112, 17566113, 17566114, 17566115, 17566116, 17566117,
                17566118, 17577951, 17577952, 17577953, 17577954, 17577955, 17577956, 17577957,
                17577958, 17577959, 17577960
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17566108, 17566109, 17566110, 17566111, 17566112, 17566113, 17566114, 17566115,
                17566116, 17566117, 17566118, 17577951, 17577952, 17577953, 17577954, 17577955,
                17577956, 17577957, 17577958, 17577959, 17577960, 17577961, 17577962, 17577963,
                17577964, 17577965, 17577966, 17577967, 17577968, 17577969, 17577970, 17577971,
                17577972, 17577973, 17577974, 17577975, 17578686, 17578687, 17578688, 17578689,
                17578690, 17578691, 17578692, 17578693, 17578694, 17578695, 17578696, 17578697,
                17578698, 17578699, 17578700
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17566111, 17566112, 17566113, 17566114, 17566115, 17566116, 17566117, 17566118,
                17577951, 17577952, 17577953, 17577954, 17577955, 17577956, 17577957, 17577958,
                17577959, 17577960, 17577961, 17577962, 17577963, 17577964, 17577965, 17577966,
                17577967, 17577968, 17577969, 17577970, 17577971, 17577972, 17577973, 17577974,
                17577975, 17578686, 17578687, 17578688, 17578689, 17578690, 17578691, 17578692,
                17578693, 17578694, 17578695, 17578696, 17578697, 17578698, 17578699, 17578700,
                17578701, 17578702, 17578703
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17566111, 17566112, 17566113, 17566114, 17566115, 17566116, 17566117, 17566118,
                17577951, 17577952, 17577953, 17577954, 17577955, 17577956, 17577957, 17577958,
                17577959, 17577960, 17577961, 17577962, 17577963, 17577964, 17577965, 17577966,
                17577967, 17577968, 17577969, 17577970, 17577971, 17577972, 17577973, 17577974,
                17577975, 17578686, 17578687, 17578688, 17578689, 17578690, 17578691, 17578692,
                17578693, 17578694, 17578695, 17578696, 17578697, 17578698, 17578699, 17578700,
                17578701, 17578702, 17578703
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17566111, 17566112, 17566113, 17566114, 17566115, 17566116, 17566117, 17566118,
                17577951, 17577952, 17577953, 17577954, 17577955, 17577956, 17577957, 17577958,
                17577959, 17577960, 17577961, 17577962, 17577963, 17577964, 17577965, 17577966,
                17577967, 17577968, 17577969, 17577970, 17577971, 17577972, 17577973, 17577974,
                17577975, 17578686, 17578687, 17578688, 17578689, 17578690, 17578691, 17578692,
                17578693, 17578694, 17578695, 17578696, 17578697, 17578698, 17578699, 17578700,
                17578701, 17578702, 17578703
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17566111, 17566112, 17566113, 17566114, 17566115, 17566116, 17566117, 17566118,
                17577951, 17577952, 17577953, 17577954, 17577955, 17577956, 17577957, 17577958,
                17577959, 17577960, 17577961, 17577962, 17577963, 17577964, 17577965, 17577966,
                17577967, 17577968, 17577969, 17577970, 17577971, 17577972, 17577973, 17577974,
                17577975, 17578686, 17578687, 17578688, 17578689, 17578690, 17578691, 17578692,
                17578693, 17578694, 17578695, 17578696, 17578697, 17578698, 17578699, 17578700,
                17578701, 17578702, 17578703
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17566112, 17566113, 17566114, 17566115, 17566116, 17566117, 17566118, 17577951,
                17577952, 17577953, 17577954, 17577955, 17577956, 17577957, 17577958, 17577959,
                17577960, 17577961, 17577962, 17577963, 17577964, 17577965, 17577966, 17577967,
                17577968, 17577969, 17577970, 17577971, 17577972, 17577973, 17577974, 17577975,
                17578686, 17578687, 17578688, 17578689, 17578690, 17578691, 17578692, 17578693,
                17578694, 17578695, 17578696, 17578697, 17578698, 17578699, 17578700, 17578701,
                17578702, 17578703, 17578704
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17566113, 17566114, 17566115, 17566116, 17566117, 17566118, 17577951, 17577952,
                17577953, 17577954, 17577955, 17577956, 17577957, 17577958, 17577959, 17577960,
                17577961, 17577962, 17577963, 17577964, 17577965, 17577966, 17577967, 17577968,
                17577969, 17577970, 17577971, 17577972, 17577973, 17577974, 17577975, 17578686,
                17578687, 17578688, 17578689, 17578690, 17578691, 17578692, 17578693, 17578694,
                17578695, 17578696, 17578697, 17578698, 17578699, 17578700, 17578701, 17578702,
                17578703, 17578704, 17578705
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17566113, 17566114, 17566115, 17566116, 17566117, 17566118, 17577951, 17577952,
                17577953, 17577954, 17577955, 17577956, 17577957, 17577958, 17577959, 17577960,
                17577961, 17577962, 17577963, 17577964, 17577965, 17577966, 17577967, 17577968,
                17577969, 17577970, 17577971, 17577972, 17577973, 17577974, 17577975, 17578686,
                17578687, 17578688, 17578689, 17578690, 17578691, 17578692, 17578693, 17578694,
                17578695, 17578696, 17578697, 17578698, 17578699, 17578700, 17578701, 17578702,
                17578703, 17578704, 17578705
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17579733, 17579734, 17579735, 17579736, 17579737, 17579738, 17579739, 17579740,
                17579741, 17579742, 17579743, 17579744, 17579745, 17579746, 17579747, 17579748,
                17579749, 17579750, 17579751, 17579752, 17579753, 17579754, 17579755, 17579756,
                17579757, 17579758, 17579759, 17579760, 17579761, 17579762, 17579763, 17579764,
                17579765, 17579766, 17579767, 17579768, 17579769, 17579770, 17579771, 17579772,
                17579773, 17579774, 17579775, 17579776, 17581244, 17581245, 17581246, 17581247,
                17581248, 17581249
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17581369, 17581370, 17582885, 17582886, 17582887, 17582888, 17582889, 17582890,
                17582891, 17582892, 17582893, 17582894, 17582895, 17582896, 17582897, 17582898,
                17582899, 17582900, 17582901, 17582902, 17582903, 17582904, 17582905, 17582906,
                17582907, 17582908, 17582909, 17582910, 17582911, 17582912, 17582913, 17582914,
                17582915, 17582916, 17582917, 17582918, 17582919, 17582920, 17582921, 17582922,
                17582923, 17582924, 17582925, 17582926, 17582927, 17582928, 17582929, 17582930,
                17582931, 17582932, 17583028
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17581370, 17582885, 17582886, 17582887, 17582888, 17582889, 17582890, 17582891,
                17582892, 17582893, 17582894, 17582895, 17582896, 17582897, 17582898, 17582899,
                17582900, 17582901, 17582902, 17582903, 17582904, 17582905, 17582906, 17582907,
                17582908, 17582909, 17582910, 17582911, 17582912, 17582913, 17582914, 17582915,
                17582916, 17582917, 17582918, 17582919, 17582920, 17582921, 17582922, 17582923,
                17582924, 17582925, 17582926, 17582927, 17582928, 17582929, 17582930, 17582931,
                17582932, 17583028, 17583029
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17581370, 17582885, 17582886, 17582887, 17582888, 17582889, 17582890, 17582891,
                17582892, 17582893, 17582894, 17582895, 17582896, 17582897, 17582898, 17582899,
                17582900, 17582901, 17582902, 17582903, 17582904, 17582905, 17582906, 17582907,
                17582908, 17582909, 17582910, 17582911, 17582912, 17582913, 17582914, 17582915,
                17582916, 17582917, 17582918, 17582919, 17582920, 17582921, 17582922, 17582923,
                17582924, 17582925, 17582926, 17582927, 17582928, 17582929, 17582930, 17582931,
                17582932, 17583028, 17583029
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17582911, 17582912, 17582913, 17582914, 17582915, 17582916, 17582917, 17582918,
                17582919, 17582920, 17582921, 17582922, 17582923, 17582924, 17582925, 17582926,
                17582927, 17582928, 17582929, 17582930, 17582931, 17582932, 17583028, 17583029,
                17583030, 17583031, 17583032, 17583033, 17583034, 17583035, 17583036, 17583037,
                17583038, 17583039, 17583040, 17583041, 17583042, 17583043, 17583044, 17583045,
                17583046, 17583047, 17583048, 17583049, 17583050, 17583051, 17583052, 17583053,
                17583054, 17583055
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17588621, 17588622, 17588623, 17588624, 17588625, 17588626, 17588627, 17588628,
                17588629, 17588630, 17588631, 17588632, 17588633, 17588634, 17588635, 17588636,
                17588637, 17588638, 17588639, 17588640, 17588641, 17588642, 17588643, 17588644,
                17588645, 17588646, 17588647, 17588648, 17588649, 17588650, 17588651, 17588652,
                17588653, 17588654, 17588655, 17588656, 17588657, 17589196, 17589197, 17589198,
                17589199, 17589200, 17589201, 17589202, 17589203, 17589204, 17589205, 17589206,
                17589207, 17589208
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17588621, 17588622, 17588623, 17588624, 17588625, 17588626, 17588627, 17588628,
                17588629, 17588630, 17588631, 17588632, 17588633, 17588634, 17588635, 17588636,
                17588637, 17588638, 17588639, 17588640, 17588641, 17588642, 17588643, 17588644,
                17588645, 17588646, 17588647, 17588648, 17588649, 17588650, 17588651, 17588652,
                17588653, 17588654, 17588655, 17588656, 17588657, 17589196, 17589197, 17589198,
                17589199, 17589200, 17589201, 17589202, 17589203, 17589204, 17589205, 17589206,
                17589207, 17589208
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17588621, 17588622, 17588623, 17588624, 17588625, 17588626, 17588627, 17588628,
                17588629, 17588630, 17588631, 17588632, 17588633, 17588634, 17588635, 17588636,
                17588637, 17588638, 17588639, 17588640, 17588641, 17588642, 17588643, 17588644,
                17588645, 17588646, 17588647, 17588648, 17588649, 17588650, 17588651, 17588652,
                17588653, 17588654, 17588655, 17588656, 17588657, 17589196, 17589197, 17589198,
                17589199, 17589200, 17589201, 17589202, 17589203, 17589204, 17589205, 17589206,
                17589207, 17589208
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17591770, 17591771, 17591772, 17591773, 17591774, 17591775, 17591776, 17591777,
                17591778, 17591779, 17591780, 17591781, 17591782, 17591783, 17591784, 17591785,
                17591786, 17591787, 17591788, 17591789, 17591790, 17591791, 17591792, 17591793,
                17591794, 17591796, 17591797, 17591798, 17591799, 17591800, 17591801, 17591802,
                17591803, 17591804, 17591805, 17591806, 17591807, 17591808, 17591809, 17591810,
                17591811, 17591812, 17591813, 17591814, 17591815, 17591816, 17591817, 17591818,
                17591819, 17591820
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17593855, 17593856, 17593857, 17593858, 17593859, 17593860, 17593861, 17593862,
                17593863, 17593864, 17593865, 17593866, 17593867, 17593868, 17593869, 17593870,
                17593871, 17593872, 17593873, 17593874, 17593875, 17593876, 17593877, 17593878,
                17593880, 17593881, 17593882, 17593883, 17593884, 17593885, 17593886, 17593887,
                17593888, 17593889, 17593890, 17593891, 17593892, 17593893, 17593894, 17593895,
                17593896, 17593897, 17593898, 17593899, 17593900, 17593901, 17593902, 17593903
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17593863, 17593864, 17593865, 17593866, 17593867, 17593868, 17593869, 17593870,
                17593871, 17593872, 17593873, 17593874, 17593875, 17593876, 17593877, 17593878,
                17593880, 17593881, 17593882, 17593883, 17593884, 17593885, 17593886, 17593887,
                17593888, 17593889, 17593890, 17593891, 17593892, 17593893, 17593894, 17593895,
                17593896, 17593897, 17593898, 17593899, 17593900, 17593901, 17593902, 17593903,
                17593904, 17593905, 17593906, 17593907
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17596476, 17596477, 17596478, 17596479, 17596480, 17596481, 17596482, 17596483,
                17596484, 17596485, 17596486, 17596487, 17596488, 17596489, 17596490, 17596491,
                17596492, 17596493, 17596494, 17596495, 17596496, 17596497, 17596498, 17596499,
                17596500, 17596501, 17596502, 17596503, 17596504, 17596505, 17596506, 17596507,
                17596508, 17596509, 17596510, 17596511, 17596512, 17596513, 17596514
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17624012, 17624013, 17624014, 17624015, 17624016, 17624017, 17624018, 17624019,
                17624020, 17625913, 17625914, 17625915, 17625916, 17625917, 17625918, 17625919,
                17625920, 17625921, 17625922, 17625923, 17625924, 17625925, 17625926, 17625927,
                17625928, 17625929, 17625930, 17625931, 17625932, 17625933, 17625934, 17625935,
                17625936, 17625937, 17625938, 17625939, 17625940, 17625941, 17625942, 17625943,
                17625944, 17625945, 17625946, 17625947, 17625948, 17625949
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                17624012, 17624013, 17624014, 17624015, 17624016, 17624017, 17624018, 17624019,
                17624020, 17625913, 17625914, 17625915, 17625916, 17625917, 17625918, 17625919,
                17625920, 17625921, 17625922, 17625923, 17625924, 17625925, 17625926, 17625927,
                17625928, 17625929, 17625930, 17625931, 17625932, 17625933, 17625934, 17625935,
                17625936, 17625937, 17625938, 17625939, 17625940, 17625941, 17625942, 17625943,
                17625944, 17625945, 17625946, 17625947, 17625948, 17625949, 17625950, 17625951,
                17625952
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                31796700, 31796701, 31796702, 31796703, 31796704, 31796705, 31796706, 31796710,
                31796711, 31796712, 31796713, 31796714, 31796715, 31796716, 31796717, 31796718,
                31796719, 31796720, 31796721, 31796722, 31796723, 31796724, 31796725, 31796726,
                31796727, 31796728, 31799014, 31799015, 31799016, 31799017, 31799018, 31799019,
                31799020, 31799021, 31799022, 31799023, 31799024, 31799025, 31799026, 31799027,
                31799028, 31799029, 31799030, 31799031, 31799032, 31799033, 31799034, 31799035,
                31799036, 31799037
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                36722692, 36722693, 36722694, 36722695, 36722696, 36722697, 36722698, 36722699,
                36722700, 36722701, 36722702, 36722703, 36722704, 36722705, 36723505, 36723506,
                36723507, 36723508, 36723509, 36723510, 36723511, 36723512, 36723513, 36723514,
                36723515, 36723516, 36723517, 36723518, 36723519, 36723520, 36723521, 36723522,
                36723523, 36723524, 36723525, 36723526, 36723527, 36723528, 36723529, 36723530,
                36723531, 36723532, 36737414, 36737415, 36737416, 36737417, 36737418, 36737419,
                36737420
            ]
        );
        let rp: Vec<i64> = it.next().unwrap().unwrap().reference_positions().collect();
        assert_eq!(
            rp,
            vec![
                44587963, 44587964, 44587965, 44587966, 44587967, 44587968, 44587969, 44587970,
                44587971, 44587972, 44587973, 44587974, 44587975, 44587976, 44587977, 44587978,
                44587979, 44587980, 44587981, 44587982, 44587983, 44589680, 44589681, 44589682,
                44589683, 44589684, 44589685, 44589686, 44589687, 44589688, 44589689, 44589690,
                44589691, 44589692, 44589693, 44589694, 44589695, 44589696, 44589697, 44589698,
                44589699, 44589700, 44589701, 44589702, 44592034, 44592035, 44592036
            ]
        );
    }

    #[test]
    fn test_reference_positions_full() {
        let mut bam = bam::Reader::from_path("./test/test_spliced_reads.bam").unwrap();
        let mut it = bam.records();

        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                None,
                None,
                None,
                None,
                None,
                Some(16050676),
                Some(16050677),
                Some(16050678),
                Some(16050679),
                Some(16050680),
                Some(16050681),
                Some(16050682),
                Some(16050683),
                Some(16050684),
                Some(16050685),
                Some(16050686),
                Some(16050687),
                Some(16050688),
                Some(16050689),
                Some(16050690),
                Some(16050691),
                Some(16050692),
                Some(16050693),
                Some(16050694),
                Some(16050695),
                Some(16050696),
                Some(16050697),
                Some(16050698),
                Some(16050699),
                Some(16050700),
                Some(16050701),
                Some(16050702),
                Some(16050703),
                Some(16050704),
                Some(16050705),
                Some(16050706),
                Some(16050707),
                Some(16050708),
                Some(16050709),
                Some(16050710),
                Some(16050711),
                Some(16050712),
                Some(16050713),
                Some(16050714),
                Some(16050715),
                Some(16050716),
                Some(16050717),
                Some(16050718),
                Some(16050719),
                Some(16050720)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16096878),
                Some(16096879),
                Some(16096880),
                Some(16096881),
                Some(16096882),
                Some(16096883),
                Some(16096884),
                Some(16096887),
                Some(16096888),
                Some(16096889),
                Some(16096890),
                Some(16096891),
                Some(16096892),
                Some(16096893),
                Some(16096894),
                Some(16096895),
                Some(16096896),
                Some(16096897),
                Some(16096898),
                Some(16096899),
                Some(16096900),
                Some(16096901),
                Some(16096902),
                Some(16096903),
                Some(16096904),
                Some(16096905),
                Some(16096906),
                Some(16096907),
                Some(16096908),
                Some(16096909),
                Some(16096910),
                Some(16096911),
                Some(16096912),
                Some(16096913),
                Some(16096914),
                Some(16096915),
                Some(16096916),
                Some(16096917),
                Some(16096918),
                Some(16096919),
                Some(16096920),
                Some(16096921),
                Some(16096922),
                Some(16096923),
                Some(16096924),
                Some(16096925),
                Some(16096926),
                Some(16096927),
                Some(16096928),
                Some(16096929),
                Some(16096930)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16097145),
                Some(16097146),
                Some(16097147),
                Some(16097148),
                Some(16097149),
                Some(16097150),
                Some(16097151),
                Some(16097152),
                Some(16097153),
                Some(16097154),
                Some(16097155),
                Some(16097156),
                Some(16097157),
                Some(16097158),
                Some(16097159),
                Some(16097160),
                Some(16097161),
                Some(16097162),
                Some(16097163),
                Some(16097164),
                Some(16097165),
                Some(16097166),
                Some(16097167),
                Some(16097168),
                Some(16097169),
                Some(16097170),
                Some(16097171),
                Some(16097172),
                Some(16097173),
                Some(16097176),
                Some(16097177),
                Some(16097178),
                Some(16097179),
                Some(16097180),
                Some(16097181),
                Some(16097182),
                Some(16097183),
                Some(16097184),
                Some(16097185),
                Some(16097186),
                Some(16097187),
                Some(16097188),
                Some(16097189),
                Some(16097190),
                Some(16097191),
                Some(16097192),
                Some(16097193),
                Some(16097194),
                Some(16097195),
                Some(16097196),
                Some(16097197)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16117350),
                Some(16117351),
                Some(16117352),
                Some(16117353),
                Some(16117354),
                Some(16117355),
                Some(16117356),
                Some(16117357),
                Some(16117358),
                Some(16117359),
                Some(16117360),
                Some(16117361),
                Some(16117362),
                Some(16117363),
                Some(16117364),
                Some(16117365),
                Some(16117366),
                Some(16117367),
                Some(16117368),
                Some(16117369),
                Some(16117370),
                Some(16117371),
                Some(16117372),
                Some(16117373),
                Some(16117374),
                Some(16117375),
                Some(16117376),
                Some(16117377),
                Some(16117378),
                Some(16117379),
                Some(16117380),
                Some(16117381),
                Some(16117382),
                Some(16117383),
                Some(16117384),
                Some(16117385),
                Some(16117386),
                Some(16117387),
                Some(16117388),
                Some(16117389),
                Some(16117390),
                Some(16117391),
                Some(16117392),
                Some(16117393),
                Some(16117394),
                Some(16117395),
                Some(16117396),
                Some(16117397),
                Some(16117398),
                Some(16117399),
                Some(16117400)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16118483),
                Some(16118484),
                Some(16118485),
                Some(16118486),
                Some(16118487),
                Some(16118488),
                Some(16118489),
                Some(16118490),
                Some(16118491),
                Some(16118492),
                Some(16118493),
                Some(16118494),
                Some(16118495),
                Some(16118496),
                Some(16118497),
                Some(16118498),
                Some(16118499),
                Some(16118500),
                Some(16118501),
                Some(16118502),
                Some(16118503),
                Some(16118504),
                Some(16118505),
                Some(16118506),
                Some(16118507),
                Some(16118508),
                Some(16118509),
                Some(16118510),
                Some(16118511),
                Some(16118512),
                Some(16118513),
                Some(16118514),
                Some(16118515),
                Some(16118516),
                Some(16118517),
                Some(16118518),
                Some(16118519),
                Some(16118520),
                Some(16118521),
                Some(16118522),
                Some(16118523),
                Some(16118524),
                Some(16118525),
                Some(16118526),
                Some(16118527),
                Some(16118528),
                Some(16118529),
                Some(16118530),
                Some(16118531),
                Some(16118532),
                Some(16118533)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16118499),
                Some(16118500),
                Some(16118501),
                Some(16118502),
                Some(16118503),
                Some(16118504),
                Some(16118505),
                Some(16118506),
                Some(16118507),
                Some(16118508),
                Some(16118509),
                Some(16118510),
                Some(16118511),
                Some(16118512),
                Some(16118513),
                Some(16118514),
                Some(16118515),
                Some(16118516),
                Some(16118517),
                Some(16118518),
                Some(16118519),
                Some(16118520),
                Some(16118521),
                Some(16118522),
                Some(16118523),
                Some(16118524),
                Some(16118525),
                Some(16118526),
                Some(16118527),
                Some(16118528),
                Some(16118529),
                Some(16118530),
                Some(16118531),
                Some(16118532),
                Some(16118533),
                Some(16118534),
                Some(16118535),
                Some(16118536),
                Some(16118537),
                Some(16118538),
                Some(16118539),
                Some(16118540),
                Some(16118541),
                Some(16118542),
                Some(16118543),
                Some(16118544),
                Some(16118545),
                Some(16118546),
                Some(16118547),
                Some(16118548),
                Some(16118549)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16118499),
                Some(16118500),
                Some(16118501),
                Some(16118502),
                Some(16118503),
                Some(16118504),
                Some(16118505),
                Some(16118506),
                Some(16118507),
                Some(16118508),
                Some(16118509),
                Some(16118510),
                Some(16118511),
                Some(16118512),
                Some(16118513),
                Some(16118514),
                Some(16118515),
                Some(16118516),
                Some(16118517),
                Some(16118518),
                Some(16118519),
                Some(16118520),
                Some(16118521),
                Some(16118522),
                Some(16118523),
                Some(16118524),
                Some(16118525),
                Some(16118526),
                Some(16118527),
                Some(16118528),
                Some(16118529),
                Some(16118530),
                Some(16118531),
                Some(16118532),
                Some(16118533),
                Some(16118534),
                Some(16118535),
                Some(16118536),
                Some(16118537),
                Some(16118538),
                Some(16118539),
                Some(16118540),
                Some(16118541),
                Some(16118542),
                Some(16118543),
                Some(16118544),
                Some(16118545),
                Some(16118546),
                Some(16118547),
                Some(16118548),
                Some(16118549)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16118499),
                Some(16118500),
                Some(16118501),
                Some(16118502),
                Some(16118503),
                Some(16118504),
                Some(16118505),
                Some(16118506),
                Some(16118507),
                Some(16118508),
                Some(16118509),
                Some(16118510),
                Some(16118511),
                Some(16118512),
                Some(16118513),
                Some(16118514),
                Some(16118515),
                Some(16118516),
                Some(16118517),
                Some(16118518),
                Some(16118519),
                Some(16118520),
                Some(16118521),
                Some(16118522),
                Some(16118523),
                Some(16118524),
                Some(16118525),
                Some(16118526),
                Some(16118527),
                Some(16118528),
                Some(16118529),
                Some(16118530),
                Some(16118531),
                Some(16118532),
                Some(16118533),
                Some(16118534),
                Some(16118535),
                Some(16118536),
                Some(16118537),
                Some(16118538),
                Some(16118539),
                Some(16118540),
                Some(16118541),
                Some(16118542),
                Some(16118543),
                Some(16118544),
                Some(16118545),
                Some(16118546),
                Some(16118547),
                Some(16118548),
                Some(16118549)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16123411),
                Some(16123412),
                Some(16123413),
                Some(16123414),
                Some(16123415),
                Some(16123416),
                Some(16123417),
                Some(16123418),
                Some(16123419),
                Some(16123420),
                Some(16123421),
                Some(16123422),
                Some(16123423),
                Some(16123424),
                Some(16123425),
                Some(16123426),
                Some(16123427),
                Some(16123428),
                Some(16123429),
                Some(16123430),
                Some(16123431),
                Some(16123432),
                Some(16123433),
                Some(16123434),
                Some(16123435),
                Some(16123436),
                Some(16123437),
                Some(16123438),
                Some(16123439),
                Some(16123440),
                Some(16123441),
                Some(16123442),
                Some(16123443),
                Some(16123444),
                Some(16123445),
                Some(16123446),
                Some(16123447),
                Some(16123448),
                Some(16123449),
                Some(16123450),
                Some(16123451),
                Some(16123452),
                Some(16123453),
                Some(16123454),
                Some(16123455),
                Some(16123456),
                Some(16123457),
                Some(16123458),
                Some(16123459),
                Some(16123460),
                Some(16123461)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                None,
                None,
                None,
                None,
                None,
                Some(16123417),
                Some(16123418),
                Some(16123419),
                Some(16123420),
                Some(16123421),
                Some(16123422),
                Some(16123423),
                Some(16123424),
                Some(16123425),
                Some(16123426),
                Some(16123427),
                Some(16123428),
                Some(16123429),
                Some(16123430),
                Some(16123431),
                Some(16123432),
                Some(16123433),
                Some(16123434),
                Some(16123435),
                Some(16123436),
                Some(16123437),
                Some(16123438),
                Some(16123439),
                Some(16123440),
                Some(16123441),
                Some(16123442),
                Some(16123443),
                Some(16123444),
                Some(16123445),
                Some(16123446),
                Some(16123447),
                Some(16123448),
                Some(16123449),
                Some(16123450),
                Some(16123451),
                Some(16123452),
                Some(16123453),
                Some(16123454),
                Some(16123455),
                Some(16123456),
                Some(16123457),
                Some(16123458),
                Some(16123459),
                Some(16123460),
                Some(16123461)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16165860),
                Some(16165861),
                Some(16165862),
                Some(16165863),
                Some(16165864),
                Some(16165865),
                Some(16165866),
                Some(16165867),
                Some(16165868),
                Some(16165869),
                Some(16165870),
                Some(16165871),
                Some(16165872),
                Some(16165873),
                Some(16165874),
                Some(16165875),
                Some(16165876),
                Some(16165877),
                Some(16165878),
                Some(16165879),
                Some(16165880),
                Some(16165881),
                Some(16165882),
                Some(16165883),
                Some(16165884),
                Some(16165885),
                Some(16165886),
                Some(16165887),
                Some(16165888),
                Some(16165889),
                Some(16165890),
                Some(16165891),
                Some(16165892),
                Some(16165893),
                Some(16165894),
                Some(16165895),
                Some(16165896),
                Some(16165897),
                Some(16165898),
                Some(16165899),
                Some(16165900),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16180871),
                Some(16180872),
                Some(16180873),
                Some(16180874),
                Some(16180875),
                Some(16180876),
                Some(16180877),
                Some(16180878),
                Some(16180879),
                Some(16180880),
                Some(16180881),
                Some(16180882),
                Some(16180883),
                Some(16180884),
                Some(16180885),
                Some(16180886),
                Some(16180887),
                Some(16180888),
                Some(16180889),
                Some(16180890),
                Some(16180891),
                Some(16180892),
                Some(16180893),
                Some(16180894),
                Some(16180895),
                Some(16180896),
                Some(16180897),
                Some(16180898),
                Some(16180899),
                Some(16180900),
                Some(16180901),
                Some(16180902),
                Some(16180903),
                Some(16180904),
                Some(16180905),
                Some(16180906),
                Some(16180907),
                Some(16180908),
                Some(16180909),
                Some(16180910),
                Some(16180911),
                Some(16180912),
                Some(16180913),
                Some(16180914),
                Some(16180915),
                Some(16180916),
                Some(16180917),
                Some(16180918),
                Some(16180919),
                Some(16180920),
                Some(16180921)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16189705),
                Some(16189706),
                Some(16189707),
                Some(16189708),
                Some(16189709),
                Some(16189710),
                Some(16189711),
                Some(16189712),
                Some(16189713),
                Some(16189714),
                Some(16189715),
                Some(16189716),
                Some(16189717),
                Some(16189718),
                Some(16189719),
                Some(16189720),
                Some(16189721),
                Some(16189722),
                Some(16189723),
                Some(16189724),
                Some(16189725),
                Some(16189726),
                Some(16189727),
                Some(16189728),
                Some(16189729),
                Some(16189730),
                Some(16189731),
                Some(16189732),
                Some(16189733),
                Some(16189734),
                Some(16189735),
                Some(16189736),
                Some(16189737),
                Some(16189738),
                Some(16189739),
                Some(16189740),
                Some(16189741),
                Some(16189742),
                Some(16189743),
                Some(16189744),
                Some(16189745),
                Some(16189746),
                Some(16189747),
                Some(16189748),
                Some(16189749),
                Some(16189750),
                Some(16189751),
                Some(16189752),
                Some(16189753),
                Some(16189754),
                Some(16189755)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16231271),
                Some(16231272),
                Some(16231273),
                Some(16231274),
                Some(16231275),
                Some(16231276),
                Some(16231277),
                Some(16231278),
                Some(16231279),
                Some(16231280),
                Some(16231281),
                Some(16231282),
                Some(16231283),
                Some(16231284),
                Some(16231285),
                Some(16231286),
                Some(16231287),
                Some(16231288),
                Some(16231289),
                Some(16231290),
                Some(16231291),
                Some(16231292),
                Some(16231293),
                Some(16231294),
                Some(16231295),
                Some(16231296),
                Some(16231297),
                Some(16231298),
                Some(16231299),
                Some(16231300),
                Some(16231301),
                Some(16231302),
                Some(16231303),
                Some(16231304),
                Some(16231305),
                Some(16231306),
                Some(16231307),
                Some(16231308),
                Some(16231309),
                Some(16231310),
                Some(16231311),
                Some(16231312),
                Some(16231313),
                Some(16231314),
                Some(16231315),
                Some(16231316),
                Some(16231317),
                Some(16231318),
                Some(16231319),
                Some(16231320),
                Some(16231321)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16237657),
                Some(16237658),
                Some(16237659),
                Some(16237660),
                Some(16237661),
                Some(16237662),
                Some(16237663),
                Some(16237664),
                Some(16237665),
                Some(16237666),
                Some(16237667),
                Some(16237668),
                Some(16237669),
                Some(16237670),
                Some(16237671),
                Some(16237672),
                Some(16237673),
                Some(16237674),
                Some(16237675),
                Some(16237676),
                Some(16237677),
                Some(16237678),
                Some(16237679),
                Some(16237680),
                Some(16237681),
                Some(16237682),
                Some(16237683),
                Some(16237684),
                Some(16237685),
                Some(16237686),
                Some(16237687),
                Some(16237688),
                Some(16237689),
                Some(16237690),
                Some(16237691),
                Some(16237692),
                Some(16237693),
                Some(16237694),
                Some(16237695),
                Some(16237696),
                Some(16237697),
                Some(16237698),
                Some(16237699),
                Some(16237700),
                Some(16237701),
                Some(16237702),
                Some(16237703),
                Some(16237704),
                Some(16237705),
                Some(16237706),
                Some(16237707)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                Some(16255012),
                Some(16255013),
                Some(16255014),
                Some(16255015),
                Some(16255016),
                Some(16255017),
                Some(16255018),
                Some(16255019),
                Some(16255020),
                Some(16255021),
                Some(16255022),
                Some(16255023),
                Some(16255024),
                Some(16255025),
                Some(16255026),
                Some(16255027),
                Some(16255028),
                Some(16255029),
                Some(16255030),
                Some(16255031),
                Some(16255032),
                Some(16255033),
                Some(16255034),
                Some(16255035),
                Some(16255036),
                Some(16255037),
                Some(16255038),
                Some(16255039),
                Some(16255040),
                Some(16255041),
                Some(16255042),
                Some(16255043),
                Some(16255044),
                Some(16255045),
                Some(16255046),
                Some(16255047),
                Some(16255048),
                Some(16255049),
                Some(16255050),
                Some(16255051),
                Some(16255052),
                Some(16255053)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16255391),
                Some(16255392),
                Some(16255393),
                Some(16255394),
                Some(16255395),
                Some(16255396),
                Some(16255397),
                Some(16255398),
                Some(16255399),
                Some(16255400),
                Some(16255401),
                Some(16255402),
                Some(16255403),
                Some(16255404),
                Some(16255405),
                Some(16255406),
                Some(16255407),
                Some(16255408),
                Some(16255409),
                Some(16255410),
                Some(16255411),
                Some(16255412),
                Some(16255413),
                Some(16255414),
                Some(16255415),
                Some(16255416),
                Some(16255417),
                Some(16255418),
                Some(16255419),
                Some(16255420),
                Some(16255421),
                Some(16255422),
                Some(16255423),
                Some(16255424),
                Some(16255425),
                Some(16255426),
                Some(16255427),
                Some(16255428),
                Some(16255429),
                Some(16255430),
                Some(16255431),
                Some(16255432),
                Some(16255433),
                Some(16255434),
                Some(16255435),
                Some(16255436),
                Some(16255437),
                Some(16255438),
                Some(16255439),
                Some(16255440),
                Some(16255441)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16255392),
                Some(16255393),
                Some(16255394),
                Some(16255395),
                Some(16255396),
                Some(16255397),
                Some(16255398),
                Some(16255399),
                Some(16255400),
                Some(16255401),
                Some(16255402),
                Some(16255403),
                Some(16255404),
                Some(16255405),
                Some(16255406),
                Some(16255407),
                Some(16255408),
                Some(16255409),
                Some(16255410),
                Some(16255411),
                Some(16255412),
                Some(16255413),
                Some(16255414),
                Some(16255415),
                Some(16255416),
                Some(16255417),
                Some(16255418),
                Some(16255419),
                Some(16255420),
                Some(16255421),
                Some(16255422),
                Some(16255423),
                Some(16255424),
                Some(16255425),
                Some(16255426),
                Some(16255427),
                Some(16255428),
                Some(16255429),
                Some(16255430),
                Some(16255431),
                Some(16255432),
                Some(16255433),
                Some(16255434),
                Some(16255435),
                Some(16255436),
                Some(16255437),
                Some(16255438),
                Some(16255439),
                Some(16255440),
                Some(16255441),
                None
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16256084),
                Some(16256085),
                Some(16256086),
                Some(16256087),
                Some(16256088),
                Some(16256089),
                Some(16256090),
                Some(16256091),
                Some(16256092),
                Some(16256093),
                Some(16256094),
                Some(16256095),
                Some(16256096),
                Some(16256097),
                Some(16256098),
                Some(16256099),
                Some(16256100),
                Some(16256101),
                Some(16256102),
                Some(16256103),
                Some(16256104),
                Some(16256105),
                Some(16256106),
                Some(16256107),
                Some(16256108),
                Some(16256109),
                Some(16256110),
                Some(16256111),
                Some(16256112),
                Some(16256113),
                Some(16256114),
                Some(16256115),
                Some(16256116),
                Some(16256117),
                Some(16256118),
                Some(16256119),
                Some(16256120),
                Some(16256121),
                Some(16256122),
                Some(16256123),
                Some(16256124),
                Some(16256125),
                Some(16256126),
                Some(16256127),
                Some(16256128),
                None,
                None,
                None,
                None,
                None,
                None
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                None,
                None,
                Some(16256224),
                Some(16256225),
                Some(16256226),
                Some(16256227),
                Some(16256228),
                Some(16256229),
                Some(16256230),
                Some(16256231),
                Some(16256232),
                Some(16256233),
                Some(16256234),
                Some(16256235),
                Some(16256236),
                Some(16256237),
                Some(16256238),
                Some(16256239),
                Some(16256240),
                Some(16256241),
                Some(16256242),
                Some(16256243),
                Some(16256244),
                Some(16256245),
                Some(16256246),
                Some(16256247),
                Some(16256248),
                Some(16256249),
                Some(16256250),
                Some(16256251),
                Some(16256252),
                Some(16256253),
                Some(16256254),
                Some(16256255),
                Some(16256256),
                Some(16256257),
                Some(16256258),
                Some(16256259),
                Some(16256260),
                Some(16256261),
                Some(16256262),
                Some(16256263),
                Some(16256264),
                Some(16256265),
                Some(16256266),
                Some(16256267),
                Some(16256268),
                Some(16256269),
                Some(16256270),
                Some(16256271)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16325199),
                Some(16325200),
                Some(16325201),
                Some(16325202),
                Some(16325203),
                Some(16325204),
                Some(16325205),
                Some(16325206),
                Some(16325207),
                Some(16325208),
                Some(16325209),
                Some(16325210),
                Some(16325211),
                Some(16325212),
                Some(16325213),
                Some(16325214),
                Some(16325215),
                Some(16325216),
                Some(16325217),
                Some(16325218),
                Some(16325219),
                Some(16325220),
                Some(16325221),
                Some(16325222),
                Some(16325223),
                Some(16325224),
                Some(16325225),
                Some(16325226),
                Some(16325227),
                Some(16325228),
                Some(16325229),
                Some(16325230),
                Some(16325231),
                Some(16325232),
                Some(16325233),
                Some(16325234),
                Some(16325235),
                Some(16325236),
                Some(16325237),
                Some(16325238),
                Some(16325239),
                Some(16325240),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                Some(16352865),
                Some(16352866),
                Some(16352867),
                Some(16352868),
                Some(16352869),
                Some(16352870),
                Some(16352871),
                Some(16352872),
                Some(16352873),
                Some(16352874),
                Some(16352875),
                Some(16352876),
                Some(16352877),
                Some(16352878),
                Some(16352879),
                Some(16352880),
                Some(16352881),
                Some(16352882),
                Some(16352883),
                Some(16352884),
                Some(16352885),
                Some(16352886),
                Some(16352887),
                Some(16352888),
                Some(16352889),
                Some(16352890),
                Some(16352891),
                Some(16352892),
                Some(16352893),
                Some(16352894),
                Some(16352895),
                Some(16352896),
                Some(16352897),
                Some(16352898),
                Some(16352899),
                Some(16352900),
                Some(16352901),
                Some(16352902)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(16352968),
                Some(16352969),
                Some(16352970),
                Some(16352971),
                Some(16352972),
                Some(16352973),
                Some(16352974),
                Some(16352975),
                Some(16352976),
                Some(16352977),
                Some(16352978),
                Some(16352979),
                Some(16352980),
                Some(16352981),
                Some(16352982),
                Some(16352983),
                Some(16352984),
                Some(16352985),
                Some(16352986),
                Some(16352987),
                Some(16352988),
                Some(16352989),
                Some(16352990),
                Some(16352991),
                Some(16352992),
                Some(16352993),
                Some(16352994),
                Some(16352995),
                Some(16352996),
                Some(16352997),
                Some(16352998),
                Some(16352999),
                Some(16353000),
                Some(16353001),
                Some(16353002),
                Some(16353003),
                Some(16353004),
                Some(16353005),
                Some(16353006),
                Some(16353007),
                Some(16353008),
                Some(16353009),
                Some(16353010),
                Some(16353011),
                None,
                None,
                None,
                None,
                None,
                None,
                None
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                None,
                None,
                None,
                None,
                Some(16414998),
                Some(16414999),
                Some(16415000),
                Some(16415001),
                Some(16415002),
                Some(16415003),
                Some(16415004),
                Some(16415005),
                Some(16415006),
                Some(16415007),
                Some(16415008),
                Some(16415009),
                Some(16415010),
                Some(16415011),
                Some(16415012),
                Some(16415013),
                Some(16415014),
                Some(16415015),
                Some(16415016),
                Some(16415017),
                Some(16415018),
                Some(16415019),
                Some(16415020),
                Some(16415021),
                Some(16415022),
                Some(16415023),
                Some(16415024),
                Some(16415025),
                Some(16415026),
                Some(16415027),
                Some(16415028),
                Some(16415029),
                Some(16415030),
                Some(16415031),
                Some(16415032),
                Some(16415033),
                Some(16415034),
                Some(16415035),
                Some(16415036),
                Some(16415037),
                Some(16415038),
                Some(16415039),
                Some(16415040),
                Some(16415041),
                Some(16415042),
                Some(16415043)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17031591),
                Some(17031592),
                Some(17031593),
                Some(17031594),
                Some(17031595),
                Some(17031596),
                Some(17031597),
                Some(17031598),
                Some(17031599),
                Some(17031600),
                Some(17031601),
                Some(17031602),
                Some(17031603),
                Some(17031604),
                Some(17031605),
                Some(17031606),
                Some(17031607),
                Some(17031608),
                Some(17031609),
                Some(17031610),
                Some(17031611),
                Some(17031612),
                Some(17031613),
                None,
                None,
                None,
                None,
                Some(17031614),
                Some(17031615),
                Some(17031616),
                Some(17031617),
                Some(17031618),
                Some(17031619),
                Some(17031620),
                Some(17031621),
                Some(17031622),
                Some(17031623),
                Some(17031624),
                Some(17031625),
                Some(17031626),
                Some(17031627),
                Some(17031628),
                Some(17031629),
                Some(17031630),
                Some(17031631),
                Some(17031632),
                Some(17031633),
                Some(17031634),
                Some(17031635),
                Some(17031636),
                Some(17031637)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17057382),
                Some(17057383),
                Some(17057384),
                Some(17057385),
                Some(17057386),
                Some(17057387),
                Some(17057388),
                Some(17057389),
                Some(17057390),
                Some(17057391),
                Some(17057392),
                Some(17057393),
                Some(17057394),
                Some(17057395),
                Some(17057396),
                Some(17057397),
                Some(17057398),
                Some(17057399),
                None,
                Some(17057400),
                Some(17057401),
                Some(17057402),
                Some(17057403),
                Some(17057404),
                Some(17057405),
                Some(17057406),
                Some(17057407),
                Some(17057408),
                Some(17057409),
                Some(17057410),
                Some(17057411),
                Some(17057412),
                Some(17057413),
                Some(17057414),
                Some(17057415),
                Some(17057416),
                Some(17057417),
                Some(17057418),
                Some(17057419),
                Some(17057420),
                Some(17057421),
                Some(17057422),
                Some(17057423),
                Some(17057424),
                Some(17057425),
                Some(17057426),
                Some(17057427),
                Some(17057428),
                Some(17057429),
                Some(17057430),
                Some(17057431)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17092766),
                Some(17092767),
                Some(17092768),
                Some(17092769),
                Some(17092770),
                Some(17092771),
                Some(17092772),
                Some(17092773),
                Some(17092774),
                Some(17092775),
                Some(17092776),
                Some(17092777),
                Some(17092778),
                Some(17092779),
                Some(17092780),
                Some(17092781),
                Some(17092782),
                Some(17094966),
                Some(17094967),
                Some(17094968),
                Some(17094969),
                Some(17094970),
                Some(17094971),
                Some(17094972),
                Some(17094973),
                Some(17094974),
                Some(17094975),
                Some(17094976),
                Some(17094977),
                Some(17094978),
                Some(17094979),
                Some(17094980),
                Some(17094981),
                Some(17094982),
                Some(17094983),
                Some(17094984),
                Some(17094985),
                Some(17094986),
                Some(17094987),
                Some(17094988),
                Some(17094989),
                Some(17094990),
                Some(17094991),
                Some(17094992),
                Some(17094993),
                Some(17094994),
                Some(17094995),
                Some(17094996),
                Some(17094997),
                Some(17094998),
                Some(17094999)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17092782),
                Some(17094966),
                Some(17094967),
                Some(17094968),
                Some(17094969),
                Some(17094970),
                Some(17094971),
                Some(17094972),
                Some(17094973),
                Some(17094974),
                Some(17094975),
                Some(17094976),
                Some(17094977),
                Some(17094978),
                Some(17094979),
                Some(17094980),
                Some(17094981),
                Some(17094982),
                Some(17094983),
                Some(17094984),
                Some(17094985),
                Some(17094986),
                Some(17094987),
                Some(17094988),
                Some(17094989),
                Some(17094990),
                Some(17094991),
                Some(17094992),
                Some(17094993),
                Some(17094994),
                Some(17094995),
                Some(17094996),
                Some(17094997),
                Some(17094998),
                Some(17094999),
                Some(17095000),
                Some(17095001),
                Some(17095002),
                Some(17095003),
                Some(17095004),
                Some(17095005),
                Some(17095006),
                Some(17095007),
                Some(17095008),
                Some(17095009),
                Some(17095010),
                Some(17095011),
                Some(17095012),
                Some(17095013),
                Some(17095014),
                Some(17095015)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17092782),
                Some(17094966),
                Some(17094967),
                Some(17094968),
                Some(17094969),
                Some(17094970),
                Some(17094971),
                Some(17094972),
                Some(17094973),
                Some(17094974),
                Some(17094975),
                Some(17094976),
                Some(17094977),
                Some(17094978),
                Some(17094979),
                Some(17094980),
                Some(17094981),
                Some(17094982),
                Some(17094983),
                Some(17094984),
                Some(17094985),
                Some(17094986),
                Some(17094987),
                Some(17094988),
                Some(17094989),
                Some(17094990),
                Some(17094991),
                Some(17094992),
                Some(17094993),
                Some(17094994),
                Some(17094995),
                Some(17094996),
                Some(17094997),
                Some(17094998),
                Some(17094999),
                Some(17095000),
                Some(17095001),
                Some(17095002),
                Some(17095003),
                Some(17095004),
                Some(17095005),
                Some(17095006),
                Some(17095007),
                Some(17095008),
                Some(17095009),
                Some(17095010),
                Some(17095011),
                Some(17095012),
                Some(17095013),
                Some(17095014),
                Some(17095015)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                Some(17137287),
                Some(17137288),
                Some(17137289),
                Some(17137290),
                Some(17137291),
                Some(17137292),
                Some(17137293),
                Some(17137294),
                Some(17137295),
                Some(17137296),
                Some(17137297),
                Some(17137298),
                Some(17137299),
                Some(17137300),
                Some(17137301),
                Some(17137302),
                Some(17137303),
                Some(17137304),
                Some(17137305),
                Some(17137306),
                Some(17137307),
                Some(17137308),
                Some(17137309),
                Some(17137310),
                Some(17137311),
                Some(17137312),
                Some(17137313),
                Some(17137314),
                Some(17137315),
                Some(17137316),
                Some(17137317),
                Some(17137318),
                Some(17137319),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                None,
                Some(17306238),
                Some(17306239),
                Some(17306240),
                Some(17306241),
                Some(17306242),
                Some(17306243),
                Some(17306244),
                Some(17306245),
                Some(17306246),
                Some(17306247),
                Some(17306248),
                Some(17306249),
                Some(17306250),
                Some(17306251),
                Some(17306252),
                Some(17306253),
                Some(17306254),
                Some(17306255),
                Some(17306256),
                Some(17306257),
                Some(17306258),
                Some(17306259),
                Some(17306260),
                Some(17306261),
                Some(17306262),
                Some(17306263),
                Some(17306264),
                Some(17306265),
                Some(17306266),
                Some(17306267),
                Some(17306268),
                Some(17306269),
                Some(17306270),
                Some(17306271),
                Some(17306272),
                Some(17306273),
                Some(17306274),
                Some(17306275),
                Some(17306276),
                Some(17306277),
                Some(17306278),
                Some(17306279),
                Some(17306280),
                Some(17306281),
                Some(17306282),
                Some(17306283),
                Some(17306284),
                Some(17306285),
                None
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                None,
                None,
                None,
                Some(17561868),
                Some(17561869),
                Some(17561870),
                Some(17561871),
                Some(17561872),
                Some(17561873),
                Some(17561874),
                Some(17561875),
                Some(17561876),
                Some(17561877),
                Some(17561878),
                Some(17561879),
                Some(17561880),
                Some(17561881),
                Some(17561882),
                Some(17561883),
                Some(17561884),
                Some(17561885),
                Some(17561886),
                Some(17561887),
                Some(17561888),
                Some(17561889),
                Some(17561890),
                Some(17561891),
                Some(17561892),
                Some(17561893),
                Some(17561894),
                Some(17561895),
                Some(17561896),
                Some(17561897),
                Some(17561898),
                Some(17561899),
                Some(17561900),
                Some(17561901),
                Some(17561902),
                Some(17561903),
                Some(17561904),
                Some(17561905),
                Some(17561906),
                Some(17561907),
                Some(17561908),
                Some(17561909),
                Some(17561910),
                Some(17561911),
                Some(17561912),
                None,
                None
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17566078),
                Some(17566079),
                Some(17566080),
                Some(17566081),
                Some(17566082),
                Some(17566083),
                Some(17566084),
                Some(17566085),
                Some(17566086),
                Some(17566087),
                Some(17566088),
                Some(17566089),
                Some(17566090),
                Some(17566091),
                Some(17566092),
                Some(17566093),
                Some(17566094),
                Some(17566095),
                Some(17566096),
                Some(17566097),
                Some(17566098),
                Some(17566099),
                Some(17566100),
                Some(17566101),
                Some(17566102),
                Some(17566103),
                Some(17566104),
                Some(17566105),
                Some(17566106),
                Some(17566107),
                Some(17566108),
                Some(17566109),
                Some(17566110),
                Some(17566111),
                Some(17566112),
                Some(17566113),
                Some(17566114),
                Some(17566115),
                Some(17566116),
                Some(17566117),
                Some(17566118),
                Some(17577951),
                Some(17577952),
                Some(17577953),
                Some(17577954),
                Some(17577955),
                Some(17577956),
                Some(17577957),
                Some(17577958),
                Some(17577959),
                Some(17577960)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17566108),
                Some(17566109),
                Some(17566110),
                Some(17566111),
                Some(17566112),
                Some(17566113),
                Some(17566114),
                Some(17566115),
                Some(17566116),
                Some(17566117),
                Some(17566118),
                Some(17577951),
                Some(17577952),
                Some(17577953),
                Some(17577954),
                Some(17577955),
                Some(17577956),
                Some(17577957),
                Some(17577958),
                Some(17577959),
                Some(17577960),
                Some(17577961),
                Some(17577962),
                Some(17577963),
                Some(17577964),
                Some(17577965),
                Some(17577966),
                Some(17577967),
                Some(17577968),
                Some(17577969),
                Some(17577970),
                Some(17577971),
                Some(17577972),
                Some(17577973),
                Some(17577974),
                Some(17577975),
                Some(17578686),
                Some(17578687),
                Some(17578688),
                Some(17578689),
                Some(17578690),
                Some(17578691),
                Some(17578692),
                Some(17578693),
                Some(17578694),
                Some(17578695),
                Some(17578696),
                Some(17578697),
                Some(17578698),
                Some(17578699),
                Some(17578700)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17566111),
                Some(17566112),
                Some(17566113),
                Some(17566114),
                Some(17566115),
                Some(17566116),
                Some(17566117),
                Some(17566118),
                Some(17577951),
                Some(17577952),
                Some(17577953),
                Some(17577954),
                Some(17577955),
                Some(17577956),
                Some(17577957),
                Some(17577958),
                Some(17577959),
                Some(17577960),
                Some(17577961),
                Some(17577962),
                Some(17577963),
                Some(17577964),
                Some(17577965),
                Some(17577966),
                Some(17577967),
                Some(17577968),
                Some(17577969),
                Some(17577970),
                Some(17577971),
                Some(17577972),
                Some(17577973),
                Some(17577974),
                Some(17577975),
                Some(17578686),
                Some(17578687),
                Some(17578688),
                Some(17578689),
                Some(17578690),
                Some(17578691),
                Some(17578692),
                Some(17578693),
                Some(17578694),
                Some(17578695),
                Some(17578696),
                Some(17578697),
                Some(17578698),
                Some(17578699),
                Some(17578700),
                Some(17578701),
                Some(17578702),
                Some(17578703)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17566111),
                Some(17566112),
                Some(17566113),
                Some(17566114),
                Some(17566115),
                Some(17566116),
                Some(17566117),
                Some(17566118),
                Some(17577951),
                Some(17577952),
                Some(17577953),
                Some(17577954),
                Some(17577955),
                Some(17577956),
                Some(17577957),
                Some(17577958),
                Some(17577959),
                Some(17577960),
                Some(17577961),
                Some(17577962),
                Some(17577963),
                Some(17577964),
                Some(17577965),
                Some(17577966),
                Some(17577967),
                Some(17577968),
                Some(17577969),
                Some(17577970),
                Some(17577971),
                Some(17577972),
                Some(17577973),
                Some(17577974),
                Some(17577975),
                Some(17578686),
                Some(17578687),
                Some(17578688),
                Some(17578689),
                Some(17578690),
                Some(17578691),
                Some(17578692),
                Some(17578693),
                Some(17578694),
                Some(17578695),
                Some(17578696),
                Some(17578697),
                Some(17578698),
                Some(17578699),
                Some(17578700),
                Some(17578701),
                Some(17578702),
                Some(17578703)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17566111),
                Some(17566112),
                Some(17566113),
                Some(17566114),
                Some(17566115),
                Some(17566116),
                Some(17566117),
                Some(17566118),
                Some(17577951),
                Some(17577952),
                Some(17577953),
                Some(17577954),
                Some(17577955),
                Some(17577956),
                Some(17577957),
                Some(17577958),
                Some(17577959),
                Some(17577960),
                Some(17577961),
                Some(17577962),
                Some(17577963),
                Some(17577964),
                Some(17577965),
                Some(17577966),
                Some(17577967),
                Some(17577968),
                Some(17577969),
                Some(17577970),
                Some(17577971),
                Some(17577972),
                Some(17577973),
                Some(17577974),
                Some(17577975),
                Some(17578686),
                Some(17578687),
                Some(17578688),
                Some(17578689),
                Some(17578690),
                Some(17578691),
                Some(17578692),
                Some(17578693),
                Some(17578694),
                Some(17578695),
                Some(17578696),
                Some(17578697),
                Some(17578698),
                Some(17578699),
                Some(17578700),
                Some(17578701),
                Some(17578702),
                Some(17578703)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17566111),
                Some(17566112),
                Some(17566113),
                Some(17566114),
                Some(17566115),
                Some(17566116),
                Some(17566117),
                Some(17566118),
                Some(17577951),
                Some(17577952),
                Some(17577953),
                Some(17577954),
                Some(17577955),
                Some(17577956),
                Some(17577957),
                Some(17577958),
                Some(17577959),
                Some(17577960),
                Some(17577961),
                Some(17577962),
                Some(17577963),
                Some(17577964),
                Some(17577965),
                Some(17577966),
                Some(17577967),
                Some(17577968),
                Some(17577969),
                Some(17577970),
                Some(17577971),
                Some(17577972),
                Some(17577973),
                Some(17577974),
                Some(17577975),
                Some(17578686),
                Some(17578687),
                Some(17578688),
                Some(17578689),
                Some(17578690),
                Some(17578691),
                Some(17578692),
                Some(17578693),
                Some(17578694),
                Some(17578695),
                Some(17578696),
                Some(17578697),
                Some(17578698),
                Some(17578699),
                Some(17578700),
                Some(17578701),
                Some(17578702),
                Some(17578703)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17566112),
                Some(17566113),
                Some(17566114),
                Some(17566115),
                Some(17566116),
                Some(17566117),
                Some(17566118),
                Some(17577951),
                Some(17577952),
                Some(17577953),
                Some(17577954),
                Some(17577955),
                Some(17577956),
                Some(17577957),
                Some(17577958),
                Some(17577959),
                Some(17577960),
                Some(17577961),
                Some(17577962),
                Some(17577963),
                Some(17577964),
                Some(17577965),
                Some(17577966),
                Some(17577967),
                Some(17577968),
                Some(17577969),
                Some(17577970),
                Some(17577971),
                Some(17577972),
                Some(17577973),
                Some(17577974),
                Some(17577975),
                Some(17578686),
                Some(17578687),
                Some(17578688),
                Some(17578689),
                Some(17578690),
                Some(17578691),
                Some(17578692),
                Some(17578693),
                Some(17578694),
                Some(17578695),
                Some(17578696),
                Some(17578697),
                Some(17578698),
                Some(17578699),
                Some(17578700),
                Some(17578701),
                Some(17578702),
                Some(17578703),
                Some(17578704)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17566113),
                Some(17566114),
                Some(17566115),
                Some(17566116),
                Some(17566117),
                Some(17566118),
                Some(17577951),
                Some(17577952),
                Some(17577953),
                Some(17577954),
                Some(17577955),
                Some(17577956),
                Some(17577957),
                Some(17577958),
                Some(17577959),
                Some(17577960),
                Some(17577961),
                Some(17577962),
                Some(17577963),
                Some(17577964),
                Some(17577965),
                Some(17577966),
                Some(17577967),
                Some(17577968),
                Some(17577969),
                Some(17577970),
                Some(17577971),
                Some(17577972),
                Some(17577973),
                Some(17577974),
                Some(17577975),
                Some(17578686),
                Some(17578687),
                Some(17578688),
                Some(17578689),
                Some(17578690),
                Some(17578691),
                Some(17578692),
                Some(17578693),
                Some(17578694),
                Some(17578695),
                Some(17578696),
                Some(17578697),
                Some(17578698),
                Some(17578699),
                Some(17578700),
                Some(17578701),
                Some(17578702),
                Some(17578703),
                Some(17578704),
                Some(17578705)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17566113),
                Some(17566114),
                Some(17566115),
                Some(17566116),
                Some(17566117),
                Some(17566118),
                Some(17577951),
                Some(17577952),
                Some(17577953),
                Some(17577954),
                Some(17577955),
                Some(17577956),
                Some(17577957),
                Some(17577958),
                Some(17577959),
                Some(17577960),
                Some(17577961),
                Some(17577962),
                Some(17577963),
                Some(17577964),
                Some(17577965),
                Some(17577966),
                Some(17577967),
                Some(17577968),
                Some(17577969),
                Some(17577970),
                Some(17577971),
                Some(17577972),
                Some(17577973),
                Some(17577974),
                Some(17577975),
                Some(17578686),
                Some(17578687),
                Some(17578688),
                Some(17578689),
                Some(17578690),
                Some(17578691),
                Some(17578692),
                Some(17578693),
                Some(17578694),
                Some(17578695),
                Some(17578696),
                Some(17578697),
                Some(17578698),
                Some(17578699),
                Some(17578700),
                Some(17578701),
                Some(17578702),
                Some(17578703),
                Some(17578704),
                Some(17578705)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                Some(17579733),
                Some(17579734),
                Some(17579735),
                Some(17579736),
                Some(17579737),
                Some(17579738),
                Some(17579739),
                Some(17579740),
                Some(17579741),
                Some(17579742),
                Some(17579743),
                Some(17579744),
                Some(17579745),
                Some(17579746),
                Some(17579747),
                Some(17579748),
                Some(17579749),
                Some(17579750),
                Some(17579751),
                Some(17579752),
                Some(17579753),
                Some(17579754),
                Some(17579755),
                Some(17579756),
                Some(17579757),
                Some(17579758),
                Some(17579759),
                Some(17579760),
                Some(17579761),
                Some(17579762),
                Some(17579763),
                Some(17579764),
                Some(17579765),
                Some(17579766),
                Some(17579767),
                Some(17579768),
                Some(17579769),
                Some(17579770),
                Some(17579771),
                Some(17579772),
                Some(17579773),
                Some(17579774),
                Some(17579775),
                Some(17579776),
                Some(17581244),
                Some(17581245),
                Some(17581246),
                Some(17581247),
                Some(17581248),
                Some(17581249)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17581369),
                Some(17581370),
                Some(17582885),
                Some(17582886),
                Some(17582887),
                Some(17582888),
                Some(17582889),
                Some(17582890),
                Some(17582891),
                Some(17582892),
                Some(17582893),
                Some(17582894),
                Some(17582895),
                Some(17582896),
                Some(17582897),
                Some(17582898),
                Some(17582899),
                Some(17582900),
                Some(17582901),
                Some(17582902),
                Some(17582903),
                Some(17582904),
                Some(17582905),
                Some(17582906),
                Some(17582907),
                Some(17582908),
                Some(17582909),
                Some(17582910),
                Some(17582911),
                Some(17582912),
                Some(17582913),
                Some(17582914),
                Some(17582915),
                Some(17582916),
                Some(17582917),
                Some(17582918),
                Some(17582919),
                Some(17582920),
                Some(17582921),
                Some(17582922),
                Some(17582923),
                Some(17582924),
                Some(17582925),
                Some(17582926),
                Some(17582927),
                Some(17582928),
                Some(17582929),
                Some(17582930),
                Some(17582931),
                Some(17582932),
                Some(17583028)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17581370),
                Some(17582885),
                Some(17582886),
                Some(17582887),
                Some(17582888),
                Some(17582889),
                Some(17582890),
                Some(17582891),
                Some(17582892),
                Some(17582893),
                Some(17582894),
                Some(17582895),
                Some(17582896),
                Some(17582897),
                Some(17582898),
                Some(17582899),
                Some(17582900),
                Some(17582901),
                Some(17582902),
                Some(17582903),
                Some(17582904),
                Some(17582905),
                Some(17582906),
                Some(17582907),
                Some(17582908),
                Some(17582909),
                Some(17582910),
                Some(17582911),
                Some(17582912),
                Some(17582913),
                Some(17582914),
                Some(17582915),
                Some(17582916),
                Some(17582917),
                Some(17582918),
                Some(17582919),
                Some(17582920),
                Some(17582921),
                Some(17582922),
                Some(17582923),
                Some(17582924),
                Some(17582925),
                Some(17582926),
                Some(17582927),
                Some(17582928),
                Some(17582929),
                Some(17582930),
                Some(17582931),
                Some(17582932),
                Some(17583028),
                Some(17583029)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17581370),
                Some(17582885),
                Some(17582886),
                Some(17582887),
                Some(17582888),
                Some(17582889),
                Some(17582890),
                Some(17582891),
                Some(17582892),
                Some(17582893),
                Some(17582894),
                Some(17582895),
                Some(17582896),
                Some(17582897),
                Some(17582898),
                Some(17582899),
                Some(17582900),
                Some(17582901),
                Some(17582902),
                Some(17582903),
                Some(17582904),
                Some(17582905),
                Some(17582906),
                Some(17582907),
                Some(17582908),
                Some(17582909),
                Some(17582910),
                Some(17582911),
                Some(17582912),
                Some(17582913),
                Some(17582914),
                Some(17582915),
                Some(17582916),
                Some(17582917),
                Some(17582918),
                Some(17582919),
                Some(17582920),
                Some(17582921),
                Some(17582922),
                Some(17582923),
                Some(17582924),
                Some(17582925),
                Some(17582926),
                Some(17582927),
                Some(17582928),
                Some(17582929),
                Some(17582930),
                Some(17582931),
                Some(17582932),
                Some(17583028),
                Some(17583029)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                Some(17582911),
                Some(17582912),
                Some(17582913),
                Some(17582914),
                Some(17582915),
                Some(17582916),
                Some(17582917),
                Some(17582918),
                Some(17582919),
                Some(17582920),
                Some(17582921),
                Some(17582922),
                Some(17582923),
                Some(17582924),
                Some(17582925),
                Some(17582926),
                Some(17582927),
                Some(17582928),
                Some(17582929),
                Some(17582930),
                Some(17582931),
                Some(17582932),
                Some(17583028),
                Some(17583029),
                Some(17583030),
                Some(17583031),
                Some(17583032),
                Some(17583033),
                Some(17583034),
                Some(17583035),
                Some(17583036),
                Some(17583037),
                Some(17583038),
                Some(17583039),
                Some(17583040),
                Some(17583041),
                Some(17583042),
                Some(17583043),
                Some(17583044),
                Some(17583045),
                Some(17583046),
                Some(17583047),
                Some(17583048),
                Some(17583049),
                Some(17583050),
                Some(17583051),
                Some(17583052),
                Some(17583053),
                Some(17583054),
                Some(17583055)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17588621),
                Some(17588622),
                Some(17588623),
                Some(17588624),
                Some(17588625),
                Some(17588626),
                Some(17588627),
                Some(17588628),
                Some(17588629),
                Some(17588630),
                Some(17588631),
                Some(17588632),
                Some(17588633),
                Some(17588634),
                Some(17588635),
                Some(17588636),
                Some(17588637),
                Some(17588638),
                Some(17588639),
                Some(17588640),
                Some(17588641),
                Some(17588642),
                Some(17588643),
                Some(17588644),
                Some(17588645),
                Some(17588646),
                Some(17588647),
                Some(17588648),
                Some(17588649),
                Some(17588650),
                Some(17588651),
                Some(17588652),
                Some(17588653),
                Some(17588654),
                Some(17588655),
                Some(17588656),
                Some(17588657),
                Some(17589196),
                Some(17589197),
                Some(17589198),
                Some(17589199),
                Some(17589200),
                Some(17589201),
                Some(17589202),
                Some(17589203),
                Some(17589204),
                Some(17589205),
                Some(17589206),
                Some(17589207),
                Some(17589208),
                None
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17588621),
                Some(17588622),
                Some(17588623),
                Some(17588624),
                Some(17588625),
                Some(17588626),
                Some(17588627),
                Some(17588628),
                Some(17588629),
                Some(17588630),
                Some(17588631),
                Some(17588632),
                Some(17588633),
                Some(17588634),
                Some(17588635),
                Some(17588636),
                Some(17588637),
                Some(17588638),
                Some(17588639),
                Some(17588640),
                Some(17588641),
                Some(17588642),
                Some(17588643),
                Some(17588644),
                Some(17588645),
                Some(17588646),
                Some(17588647),
                Some(17588648),
                Some(17588649),
                Some(17588650),
                Some(17588651),
                Some(17588652),
                Some(17588653),
                Some(17588654),
                Some(17588655),
                Some(17588656),
                Some(17588657),
                Some(17589196),
                Some(17589197),
                Some(17589198),
                Some(17589199),
                Some(17589200),
                Some(17589201),
                Some(17589202),
                Some(17589203),
                Some(17589204),
                Some(17589205),
                Some(17589206),
                Some(17589207),
                Some(17589208),
                None
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17588621),
                Some(17588622),
                Some(17588623),
                Some(17588624),
                Some(17588625),
                Some(17588626),
                Some(17588627),
                Some(17588628),
                Some(17588629),
                Some(17588630),
                Some(17588631),
                Some(17588632),
                Some(17588633),
                Some(17588634),
                Some(17588635),
                Some(17588636),
                Some(17588637),
                Some(17588638),
                Some(17588639),
                Some(17588640),
                Some(17588641),
                Some(17588642),
                Some(17588643),
                Some(17588644),
                Some(17588645),
                Some(17588646),
                Some(17588647),
                Some(17588648),
                Some(17588649),
                Some(17588650),
                Some(17588651),
                Some(17588652),
                Some(17588653),
                Some(17588654),
                Some(17588655),
                Some(17588656),
                Some(17588657),
                Some(17589196),
                Some(17589197),
                Some(17589198),
                Some(17589199),
                Some(17589200),
                Some(17589201),
                Some(17589202),
                Some(17589203),
                Some(17589204),
                Some(17589205),
                Some(17589206),
                Some(17589207),
                Some(17589208),
                None
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                Some(17591770),
                Some(17591771),
                Some(17591772),
                Some(17591773),
                Some(17591774),
                Some(17591775),
                Some(17591776),
                Some(17591777),
                Some(17591778),
                Some(17591779),
                Some(17591780),
                Some(17591781),
                Some(17591782),
                Some(17591783),
                Some(17591784),
                Some(17591785),
                Some(17591786),
                Some(17591787),
                Some(17591788),
                Some(17591789),
                Some(17591790),
                Some(17591791),
                Some(17591792),
                Some(17591793),
                Some(17591794),
                Some(17591796),
                Some(17591797),
                Some(17591798),
                Some(17591799),
                Some(17591800),
                Some(17591801),
                Some(17591802),
                Some(17591803),
                Some(17591804),
                Some(17591805),
                Some(17591806),
                Some(17591807),
                Some(17591808),
                Some(17591809),
                Some(17591810),
                Some(17591811),
                Some(17591812),
                Some(17591813),
                Some(17591814),
                Some(17591815),
                Some(17591816),
                Some(17591817),
                Some(17591818),
                Some(17591819),
                Some(17591820)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17593855),
                Some(17593856),
                Some(17593857),
                Some(17593858),
                Some(17593859),
                Some(17593860),
                Some(17593861),
                Some(17593862),
                Some(17593863),
                Some(17593864),
                Some(17593865),
                Some(17593866),
                Some(17593867),
                Some(17593868),
                Some(17593869),
                Some(17593870),
                Some(17593871),
                Some(17593872),
                Some(17593873),
                Some(17593874),
                Some(17593875),
                Some(17593876),
                Some(17593877),
                Some(17593878),
                Some(17593880),
                Some(17593881),
                Some(17593882),
                Some(17593883),
                Some(17593884),
                Some(17593885),
                Some(17593886),
                Some(17593887),
                Some(17593888),
                Some(17593889),
                Some(17593890),
                Some(17593891),
                Some(17593892),
                Some(17593893),
                Some(17593894),
                Some(17593895),
                Some(17593896),
                Some(17593897),
                Some(17593898),
                Some(17593899),
                Some(17593900),
                Some(17593901),
                Some(17593902),
                Some(17593903),
                None,
                None,
                None
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(17593863),
                Some(17593864),
                Some(17593865),
                Some(17593866),
                Some(17593867),
                Some(17593868),
                Some(17593869),
                Some(17593870),
                Some(17593871),
                Some(17593872),
                Some(17593873),
                Some(17593874),
                Some(17593875),
                Some(17593876),
                Some(17593877),
                Some(17593878),
                Some(17593880),
                Some(17593881),
                Some(17593882),
                Some(17593883),
                Some(17593884),
                Some(17593885),
                Some(17593886),
                Some(17593887),
                Some(17593888),
                Some(17593889),
                Some(17593890),
                Some(17593891),
                Some(17593892),
                Some(17593893),
                Some(17593894),
                Some(17593895),
                Some(17593896),
                Some(17593897),
                Some(17593898),
                Some(17593899),
                Some(17593900),
                Some(17593901),
                Some(17593902),
                Some(17593903),
                Some(17593904),
                Some(17593905),
                Some(17593906),
                Some(17593907),
                None,
                None,
                None,
                None,
                None,
                None,
                None
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                Some(17596476),
                Some(17596477),
                Some(17596478),
                Some(17596479),
                Some(17596480),
                Some(17596481),
                Some(17596482),
                None,
                Some(17596483),
                Some(17596484),
                Some(17596485),
                Some(17596486),
                Some(17596487),
                Some(17596488),
                Some(17596489),
                Some(17596490),
                Some(17596491),
                Some(17596492),
                Some(17596493),
                Some(17596494),
                Some(17596495),
                Some(17596496),
                Some(17596497),
                Some(17596498),
                Some(17596499),
                Some(17596500),
                Some(17596501),
                Some(17596502),
                Some(17596503),
                Some(17596504),
                Some(17596505),
                Some(17596506),
                Some(17596507),
                Some(17596508),
                Some(17596509),
                Some(17596510),
                Some(17596511),
                Some(17596512),
                Some(17596513),
                Some(17596514)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                None,
                None,
                None,
                None,
                Some(17624012),
                Some(17624013),
                Some(17624014),
                Some(17624015),
                Some(17624016),
                Some(17624017),
                Some(17624018),
                Some(17624019),
                Some(17624020),
                Some(17625913),
                Some(17625914),
                Some(17625915),
                Some(17625916),
                Some(17625917),
                Some(17625918),
                Some(17625919),
                Some(17625920),
                Some(17625921),
                Some(17625922),
                Some(17625923),
                Some(17625924),
                Some(17625925),
                Some(17625926),
                Some(17625927),
                Some(17625928),
                Some(17625929),
                Some(17625930),
                Some(17625931),
                Some(17625932),
                Some(17625933),
                Some(17625934),
                Some(17625935),
                Some(17625936),
                Some(17625937),
                Some(17625938),
                Some(17625939),
                Some(17625940),
                Some(17625941),
                Some(17625942),
                Some(17625943),
                Some(17625944),
                Some(17625945),
                Some(17625946),
                Some(17625947),
                Some(17625948),
                Some(17625949)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                None,
                Some(17624012),
                Some(17624013),
                Some(17624014),
                Some(17624015),
                Some(17624016),
                Some(17624017),
                Some(17624018),
                Some(17624019),
                Some(17624020),
                Some(17625913),
                Some(17625914),
                Some(17625915),
                Some(17625916),
                Some(17625917),
                Some(17625918),
                Some(17625919),
                Some(17625920),
                Some(17625921),
                Some(17625922),
                Some(17625923),
                Some(17625924),
                Some(17625925),
                Some(17625926),
                Some(17625927),
                Some(17625928),
                Some(17625929),
                Some(17625930),
                Some(17625931),
                Some(17625932),
                Some(17625933),
                Some(17625934),
                Some(17625935),
                Some(17625936),
                Some(17625937),
                Some(17625938),
                Some(17625939),
                Some(17625940),
                Some(17625941),
                Some(17625942),
                Some(17625943),
                Some(17625944),
                Some(17625945),
                Some(17625946),
                Some(17625947),
                Some(17625948),
                Some(17625949),
                Some(17625950),
                Some(17625951),
                Some(17625952)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                Some(31796700),
                Some(31796701),
                Some(31796702),
                Some(31796703),
                Some(31796704),
                Some(31796705),
                Some(31796706),
                Some(31796710),
                Some(31796711),
                Some(31796712),
                Some(31796713),
                Some(31796714),
                Some(31796715),
                Some(31796716),
                Some(31796717),
                Some(31796718),
                Some(31796719),
                Some(31796720),
                Some(31796721),
                Some(31796722),
                Some(31796723),
                Some(31796724),
                Some(31796725),
                Some(31796726),
                Some(31796727),
                Some(31796728),
                Some(31799014),
                Some(31799015),
                Some(31799016),
                Some(31799017),
                Some(31799018),
                Some(31799019),
                Some(31799020),
                Some(31799021),
                Some(31799022),
                Some(31799023),
                Some(31799024),
                Some(31799025),
                Some(31799026),
                Some(31799027),
                Some(31799028),
                Some(31799029),
                Some(31799030),
                Some(31799031),
                Some(31799032),
                Some(31799033),
                Some(31799034),
                Some(31799035),
                Some(31799036),
                Some(31799037)
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                Some(36722692),
                Some(36722693),
                Some(36722694),
                Some(36722695),
                Some(36722696),
                Some(36722697),
                Some(36722698),
                Some(36722699),
                Some(36722700),
                Some(36722701),
                Some(36722702),
                Some(36722703),
                Some(36722704),
                Some(36722705),
                Some(36723505),
                Some(36723506),
                Some(36723507),
                Some(36723508),
                Some(36723509),
                Some(36723510),
                Some(36723511),
                Some(36723512),
                Some(36723513),
                Some(36723514),
                Some(36723515),
                Some(36723516),
                Some(36723517),
                Some(36723518),
                Some(36723519),
                Some(36723520),
                Some(36723521),
                Some(36723522),
                Some(36723523),
                Some(36723524),
                Some(36723525),
                Some(36723526),
                Some(36723527),
                Some(36723528),
                Some(36723529),
                Some(36723530),
                Some(36723531),
                Some(36723532),
                Some(36737414),
                Some(36737415),
                Some(36737416),
                Some(36737417),
                Some(36737418),
                Some(36737419),
                Some(36737420),
                None,
                None
            ]
        );
        let rp: Vec<_> = it
            .next()
            .unwrap()
            .unwrap()
            .reference_positions_full()
            .collect();
        assert_eq!(
            rp,
            vec![
                None,
                None,
                None,
                None,
                Some(44587963),
                Some(44587964),
                Some(44587965),
                Some(44587966),
                Some(44587967),
                Some(44587968),
                Some(44587969),
                Some(44587970),
                Some(44587971),
                Some(44587972),
                Some(44587973),
                Some(44587974),
                Some(44587975),
                Some(44587976),
                Some(44587977),
                Some(44587978),
                Some(44587979),
                Some(44587980),
                Some(44587981),
                Some(44587982),
                Some(44587983),
                Some(44589680),
                Some(44589681),
                Some(44589682),
                Some(44589683),
                Some(44589684),
                Some(44589685),
                Some(44589686),
                Some(44589687),
                Some(44589688),
                Some(44589689),
                Some(44589690),
                Some(44589691),
                Some(44589692),
                Some(44589693),
                Some(44589694),
                Some(44589695),
                Some(44589696),
                Some(44589697),
                Some(44589698),
                Some(44589699),
                Some(44589700),
                Some(44589701),
                Some(44589702),
                Some(44592034),
                Some(44592035),
                Some(44592036)
            ]
        );
    }

    #[test]
    fn test_reference_start_end() {
        let mut bam = bam::Reader::from_path("./test/test_spliced_reads.bam").unwrap();
        let mut it = bam.records();
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16050676);
        assert_eq!(read.reference_end(), 16050721);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16096878);
        assert_eq!(read.reference_end(), 16096931);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16097145);
        assert_eq!(read.reference_end(), 16097198);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16117350);
        assert_eq!(read.reference_end(), 16117401);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16118483);
        assert_eq!(read.reference_end(), 16118534);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16118499);
        assert_eq!(read.reference_end(), 16118550);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16118499);
        assert_eq!(read.reference_end(), 16118550);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16118499);
        assert_eq!(read.reference_end(), 16118550);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16123411);
        assert_eq!(read.reference_end(), 16123462);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16123417);
        assert_eq!(read.reference_end(), 16123462);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16165860);
        assert_eq!(read.reference_end(), 16165901);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16180871);
        assert_eq!(read.reference_end(), 16180922);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16189705);
        assert_eq!(read.reference_end(), 16189756);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16231271);
        assert_eq!(read.reference_end(), 16231322);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16237657);
        assert_eq!(read.reference_end(), 16237708);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16255012);
        assert_eq!(read.reference_end(), 16255054);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16255391);
        assert_eq!(read.reference_end(), 16255442);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16255392);
        assert_eq!(read.reference_end(), 16255442);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16256084);
        assert_eq!(read.reference_end(), 16256129);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16256224);
        assert_eq!(read.reference_end(), 16256272);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16325199);
        assert_eq!(read.reference_end(), 16325241);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16352865);
        assert_eq!(read.reference_end(), 16352903);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16352968);
        assert_eq!(read.reference_end(), 16353012);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 16414998);
        assert_eq!(read.reference_end(), 16415044);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17031591);
        assert_eq!(read.reference_end(), 17031638);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17057382);
        assert_eq!(read.reference_end(), 17057432);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17092766);
        assert_eq!(read.reference_end(), 17095000);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17092782);
        assert_eq!(read.reference_end(), 17095016);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17092782);
        assert_eq!(read.reference_end(), 17095016);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17137287);
        assert_eq!(read.reference_end(), 17137320);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17306238);
        assert_eq!(read.reference_end(), 17306286);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17561868);
        assert_eq!(read.reference_end(), 17561913);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17566078);
        assert_eq!(read.reference_end(), 17577961);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17566108);
        assert_eq!(read.reference_end(), 17578701);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17566111);
        assert_eq!(read.reference_end(), 17578704);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17566111);
        assert_eq!(read.reference_end(), 17578704);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17566111);
        assert_eq!(read.reference_end(), 17578704);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17566111);
        assert_eq!(read.reference_end(), 17578704);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17566112);
        assert_eq!(read.reference_end(), 17578705);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17566113);
        assert_eq!(read.reference_end(), 17578706);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17566113);
        assert_eq!(read.reference_end(), 17578706);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17579733);
        assert_eq!(read.reference_end(), 17581250);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17581369);
        assert_eq!(read.reference_end(), 17583029);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17581370);
        assert_eq!(read.reference_end(), 17583030);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17581370);
        assert_eq!(read.reference_end(), 17583030);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17582911);
        assert_eq!(read.reference_end(), 17583056);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17588621);
        assert_eq!(read.reference_end(), 17589209);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17588621);
        assert_eq!(read.reference_end(), 17589209);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17588621);
        assert_eq!(read.reference_end(), 17589209);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17591770);
        assert_eq!(read.reference_end(), 17591821);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17593855);
        assert_eq!(read.reference_end(), 17593904);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17593863);
        assert_eq!(read.reference_end(), 17593908);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17596476);
        assert_eq!(read.reference_end(), 17596515);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17624012);
        assert_eq!(read.reference_end(), 17625950);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 17624012);
        assert_eq!(read.reference_end(), 17625953);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 31796700);
        assert_eq!(read.reference_end(), 31799038);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 36722692);
        assert_eq!(read.reference_end(), 36737421);
        let read = it.next().unwrap().unwrap();
        assert_eq!(read.reference_start(), 44587963);
        assert_eq!(read.reference_end(), 44592037);
    }

    #[test]
    fn test_infer_seq_len() {
        use std::convert::TryFrom;
        let mut bam = bam::Reader::from_path("./test/test_spliced_reads.bam").unwrap();
        for read in bam.records() {
            let read = read.unwrap();
            assert_eq!(read.seq_len_from_cigar(false), 51);
            assert_eq!(read.seq_len_from_cigar(true), 51);
        }

        let mut read = bam::Record::new();
        for (input_cigar, supposed_length_wo_hard_clip, supposed_length_with_hard_clip) in &[
            ("40M", 40, 40),
            ("40=", 40, 40),
            ("40X", 40, 40),
            ("20M5I20M", 45, 45),
            ("20M5D20M", 40, 40),
            ("5H35M", 35, 40),
            ("5S35M", 40, 40),
            ("35M5H", 35, 40),
            ("35M5S", 40, 40),
            ("", 0, 0),
        ] {
            read.set(
                b"test",
                Some(&CigarString::try_from(*input_cigar).unwrap()),
                b"agtc",
                b"BBBB",
            );
            assert_eq!(
                &read.seq_len_from_cigar(false),
                supposed_length_wo_hard_clip
            );
            assert_eq!(
                &read.seq_len_from_cigar(true),
                supposed_length_with_hard_clip
            );
        }
    }
}
