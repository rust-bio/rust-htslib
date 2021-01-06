// Copyright 2021 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use strum_macros::AsRefStr;

use crate::bam;

/// Representation of read orientation
/// (e.g. F1R2 means that the forward read comes first on the reference contig,
/// followed by the reverse read, on the same contig).
#[derive(Debug, Clone, Copy, PartialEq, Eq, AsRefStr)]
#[cfg_attr(feature = "serde_feature", derive(Serialize, Deserialize))]
pub(crate) enum ReadOrientation {
    F1R2,
    F2R1,
    R1F2,
    R2F1,
    F1F2,
    R1R2,
    F2F1,
    R2R1,
    None,
}

impl<'a> From<&'a bam::Record> for ReadOrientation {
    /// Infer read orientation from given BAM record. Returns `ReadOrientation::None` if record
    /// is not properly paired, mates are not mapping to the same contig, or mates start at the
    /// same position.
    fn from(record: &bam::Record) -> Self {
        if record.is_paired() && record.is_proper_pair() && record.tid() == record.mtid() {
            if record.pos() == record.mpos() {
                // both reads start at the same position, we cannot decide on the orientation.
                return ReadOrientation::None;
            }

            let (is_reverse, is_first_in_template, is_mate_reverse) =
                if record.pos() < record.mpos() {
                    // given record is the left one
                    (
                        record.is_reverse(),
                        record.is_first_in_template(),
                        record.is_mate_reverse(),
                    )
                } else {
                    // given record is the right one
                    (
                        record.is_mate_reverse(),
                        record.is_last_in_template(),
                        record.is_reverse(),
                    )
                };
            match (is_reverse, is_first_in_template, is_mate_reverse) {
                (false, false, false) => ReadOrientation::F2F1,
                (false, false, true) => ReadOrientation::F2R1,
                (false, true, false) => ReadOrientation::F1F2,
                (true, false, false) => ReadOrientation::R2F1,
                (false, true, true) => ReadOrientation::F1R2,
                (true, false, true) => ReadOrientation::R2R1,
                (true, true, false) => ReadOrientation::R1F2,
                (true, true, true) => ReadOrientation::R1R2,
            }
        } else {
            ReadOrientation::None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::bam::{self, Read};

    #[test]
    fn test_read_orientation_f1r2() {
        let mut bam = bam::Reader::from_path(&"test/test_paired.sam").unwrap();
        let mut record = bam::Record::new();
        bam.read(&mut record).unwrap().unwrap();

        assert_eq!(ReadOrientation::from(&record), ReadOrientation::F1R2);
    }
}
