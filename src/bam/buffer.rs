// Copyright 2017 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


use std::collections::{vec_deque, VecDeque};
use std::error::Error;
use std::str;

use bam;
use bam::Read;


/// A buffer for BAM records. This allows access regions in a sorted BAM file while iterating
/// over it in a single pass.
/// The buffer is implemented as a ringbuffer, such that extension or movement to the right has
/// linear complexity. The buffer makes use of indexed random access. Hence, when fetching a
/// region at the very end of the BAM, everything before is omitted without cost.
pub struct RecordBuffer {
    reader: bam::IndexedReader,
    inner: VecDeque<bam::Record>,
    overflow: Option<bam::Record>
}


unsafe impl Sync for RecordBuffer {}
unsafe impl Send for RecordBuffer {}


impl RecordBuffer {
    /// Create a new `RecordBuffer`.
    pub fn new(bam: bam::IndexedReader) -> Self {
        RecordBuffer {
            reader: bam,
            inner: VecDeque::new(),
            overflow: None
        }
    }

    /// Return end position of buffer.
    fn end(&self) -> Option<u32> {
        self.inner.back().map(|rec| rec.pos() as u32)
    }

    fn tid(&self) -> Option<i32> {
        self.inner.back().map(|rec| rec.tid())
    }

    /// Fill buffer at the given interval. The start coordinate has to be left of
    /// the start coordinate of any previous `fill` operation.
    /// Coordinates are 0-based, and end is exclusive.
    pub fn fetch(&mut self, chrom: &[u8], start: u32, end: u32) -> Result<(), Box<Error>> {
        // move overflow from last fetch into ringbuffer
        if self.overflow.is_some() {
            self.inner.push_back(self.overflow.take().unwrap());
        }

        if let Some(tid) = self.reader.header.tid(chrom) {
            let window_start = start;
            if self.inner.is_empty() || self.end().unwrap() < window_start || self.tid().unwrap() != tid as i32 {
                let end = self.reader.header.target_len(tid).unwrap();
                try!(self.reader.fetch(tid, window_start, end));
                self.inner.clear();
            } else {
                // remove records too far left
                let to_remove = self.inner.iter().take_while(|rec| rec.pos() < window_start as i32).count();
                for _ in 0..to_remove {
                    self.inner.pop_front();
                }
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
                }
            }

            Ok(())
        } else {
            Err(Box::new(RecordBufferError::UnknownSequence(str::from_utf8(chrom).unwrap().to_owned())))
        }
    }

    /// Iterate over records that have been fetched with `fetch`.
    pub fn iter(&self) -> vec_deque::Iter<bam::Record> {
        self.inner.iter()
    }
}


quick_error! {
    #[derive(Debug)]
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
    use itertools::Itertools;
    use bam;

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
