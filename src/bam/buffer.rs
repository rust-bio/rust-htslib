// Copyright 2017 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::{vec_deque, VecDeque};
use std::error::Error;
use std::str;

use crate::bam;
use crate::bam::Read;

/// A buffer for BAM records. This allows access regions in a sorted BAM file while iterating
/// over it in a single pass.
/// The buffer is implemented as a ringbuffer, such that extension or movement to the right has
/// linear complexity. The buffer makes use of indexed random access. Hence, when fetching a
/// region at the very end of the BAM, everything before is omitted without cost.
#[derive(Debug)]
pub struct RecordBuffer {
    reader: bam::IndexedReader,
    inner: VecDeque<bam::Record>,
    overflow: Option<bam::Record>,
}

unsafe impl Sync for RecordBuffer {}
unsafe impl Send for RecordBuffer {}

impl RecordBuffer {
    /// Create a new `RecordBuffer`.
    pub fn new(bam: bam::IndexedReader) -> Self {
        RecordBuffer {
            reader: bam,
            inner: VecDeque::new(),
            overflow: None,
        }
    }

    /// Return start position of buffer
    fn start(&self) -> Option<u32> {
        self.inner.front().map(|rec| rec.pos() as u32)
    }

    /// Return end position of buffer.
    fn end(&self) -> Option<u32> {
        self.inner.back().map(|rec| rec.pos() as u32)
    }

    fn tid(&self) -> Option<i32> {
        self.inner.back().map(|rec| rec.tid())
    }

    /// Fill buffer at the given interval. If the start coordinate is left of
    /// the previous start coordinate, this will use an additional BAM fetch IO operation.
    /// Coordinates are 0-based, and end is exclusive.
    /// Returns tuple with numbers of added and deleted records since the previous fetch.
    #[allow(unused_assignments)] // TODO this is needed because rustc thinks that deleted is unused
    pub fn fetch(
        &mut self,
        chrom: &[u8],
        start: u32,
        end: u32,
    ) -> Result<(usize, usize), Box<Error>> {
        let mut added = 0;
        // move overflow from last fetch into ringbuffer
        if self.overflow.is_some() {
            added += 1;
            self.inner.push_back(self.overflow.take().unwrap());
        }

        if let Some(tid) = self.reader.header.tid(chrom) {
            let mut deleted = 0;
            let window_start = start;
            if self.inner.is_empty()
                || self.end().unwrap() < window_start
                || self.tid().unwrap() != tid as i32
                || self.start().unwrap() > window_start
            {
                let end = self.reader.header.target_len(tid).unwrap();
                self.reader.fetch(tid, window_start, end)?;
                deleted = self.inner.len();
                self.inner.clear();
            } else {
                // remove records too far left
                let to_remove = self
                    .inner
                    .iter()
                    .take_while(|rec| rec.pos() < window_start as i32)
                    .count();
                for _ in 0..to_remove {
                    self.inner.pop_front();
                }
                deleted = to_remove;
            }

            // extend to the right
            loop {
                let mut record = bam::Record::new();
                if let Err(e) = self.reader.read(&mut record) {
                    if e.is_eof() {
                        break;
                    }
                    return Err(Box::new(e));
                }

                if record.is_unmapped() {
                    continue;
                }

                let pos = record.pos();
                if pos >= end as i32 {
                    self.overflow = Some(record);
                    break;
                } else {
                    self.inner.push_back(record);
                    added += 1;
                }
            }

            Ok((added, deleted))
        } else {
            Err(Box::new(RecordBufferError::UnknownSequence(
                str::from_utf8(chrom).unwrap().to_owned(),
            )))
        }
    }

    /// Iterate over records that have been fetched with `fetch`.
    pub fn iter(&self) -> vec_deque::Iter<bam::Record> {
        self.inner.iter()
    }

    /// Iterate over mutable references to records that have been fetched with `fetch`.
    pub fn iter_mut(&mut self) -> vec_deque::IterMut<bam::Record> {
        self.inner.iter_mut()
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum RecordBufferError {
        UnknownSequence(chrom: String) {
            description("unknown sequence")
            display("sequence {} cannot be found in BAM", chrom)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam;
    use itertools::Itertools;

    #[test]
    fn test_buffer() {
        let reader = bam::IndexedReader::from_path(&"test/test.bam").unwrap();

        let mut buffer = RecordBuffer::new(reader);

        buffer.fetch(b"CHROMOSOME_I", 1, 5).unwrap();
        {
            let records = buffer.iter().collect_vec();
            assert_eq!(records.len(), 6);
            assert_eq!(records[0].pos(), 1);
            assert_eq!(records[1].pos(), 1);
            assert_eq!(records[2].pos(), 1);
            assert_eq!(records[3].pos(), 1);
            assert_eq!(records[4].pos(), 1);
            assert_eq!(records[5].pos(), 1);
        }
    }
}
