// Copyright 2017 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::{vec_deque, VecDeque};
use std::mem;
use std::rc::Rc;
use std::str;

use crate::bam;
use crate::bam::Read;
use crate::errors::{Error, Result};

/// A buffer for BAM records. This allows access regions in a sorted BAM file while iterating
/// over it in a single pass.
/// The buffer is implemented as a ringbuffer, such that extension or movement to the right has
/// linear complexity. The buffer makes use of indexed random access. Hence, when fetching a
/// region at the very end of the BAM, everything before is omitted without cost.
#[derive(Debug)]
pub struct RecordBuffer {
    reader: bam::IndexedReader,
    inner: VecDeque<Rc<bam::Record>>,
    overflow: Option<Rc<bam::Record>>,
    cache_cigar: bool,
    min_refetch_distance: u64,
    buffer_record: Rc<bam::Record>,
}

unsafe impl Sync for RecordBuffer {}
unsafe impl Send for RecordBuffer {}

impl RecordBuffer {
    /// Create a new `RecordBuffer`.
    ///
    /// # Arguments
    ///
    /// * `bam` - BAM reader
    /// * `cache_cigar` - whether to call `bam::Record::cache_cigar()` for each record.
    pub fn new(bam: bam::IndexedReader, cache_cigar: bool) -> Self {
        RecordBuffer {
            reader: bam,
            inner: VecDeque::new(),
            overflow: None,
            cache_cigar,
            min_refetch_distance: 1,
            buffer_record: Rc::new(bam::Record::new()),
        }
    }

    /// maximum distance to previous fetch window such that a
    /// new fetch operation is performed. If the distance is smaller, buffer will simply
    /// read through until the start of the new fetch window (probably saving some time
    /// by avoiding the random access).
    pub fn set_min_refetch_distance(&mut self, min_refetch_distance: u64) {
        self.min_refetch_distance = min_refetch_distance;
    }

    /// Return start position of buffer
    fn start(&self) -> Option<u64> {
        self.inner.front().map(|rec| rec.pos() as u64)
    }

    /// Return end position of buffer.
    fn end(&self) -> Option<u64> {
        self.inner.back().map(|rec| rec.pos() as u64)
    }

    fn tid(&self) -> Option<i32> {
        self.inner.back().map(|rec| rec.tid())
    }

    /// Fill buffer at the given interval. If the start coordinate is left of
    /// the previous start coordinate, this will use an additional BAM fetch IO operation.
    /// Coordinates are 0-based, and end is exclusive.
    /// Returns tuple with numbers of added and deleted records since the previous fetch.
    #[allow(unused_assignments)] // TODO this is needed because rustc thinks that deleted is unused
    pub fn fetch(&mut self, chrom: &[u8], start: u64, end: u64) -> Result<(usize, usize)> {
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
                || window_start.saturating_sub(self.end().unwrap()) >= self.min_refetch_distance
                || self.tid().unwrap() != tid as i32
                || self.start().unwrap() > window_start
            {
                let end = self.reader.header.target_len(tid).unwrap();
                self.reader.fetch((tid, window_start, end))?;
                deleted = self.inner.len();
                self.inner.clear();
            } else {
                // remove records too far left
                let to_remove = self
                    .inner
                    .iter()
                    .take_while(|rec| rec.pos() < window_start as i64)
                    .count();
                for _ in 0..to_remove {
                    self.inner.pop_front();
                }
                deleted = to_remove;
            }

            // extend to the right
            loop {
                match self
                    .reader
                    .read(Rc::get_mut(&mut self.buffer_record).unwrap())
                {
                    None => break,
                    Some(res) => res?,
                }

                if self.buffer_record.is_unmapped() {
                    continue;
                }

                let pos = self.buffer_record.pos();

                // skip records before the start
                if pos < start as i64 {
                    continue;
                }

                // Record is kept, do not reuse it for next iteration
                // and thus create a new one.
                let mut record = mem::replace(&mut self.buffer_record, Rc::new(bam::Record::new()));

                if self.cache_cigar {
                    Rc::get_mut(&mut record).unwrap().cache_cigar();
                }

                if pos >= end as i64 {
                    self.overflow = Some(record);
                    break;
                } else {
                    self.inner.push_back(record);
                    added += 1;
                }
            }

            Ok((added, deleted))
        } else {
            Err(Error::UnknownSequence {
                sequence: str::from_utf8(chrom).unwrap().to_owned(),
            })
        }
    }

    /// Iterate over records that have been fetched with `fetch`.
    pub fn iter(&self) -> vec_deque::Iter<Rc<bam::Record>> {
        self.inner.iter()
    }

    /// Iterate over mutable references to records that have been fetched with `fetch`.
    pub fn iter_mut(&mut self) -> vec_deque::IterMut<Rc<bam::Record>> {
        self.inner.iter_mut()
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam;

    #[test]
    fn test_buffer() {
        let reader = bam::IndexedReader::from_path(&"test/test.bam").unwrap();
        let mut buffer = RecordBuffer::new(reader, false);

        buffer.fetch(b"CHROMOSOME_I", 1, 5).unwrap();
        {
            let records: Vec<_> = buffer.iter().collect();
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
