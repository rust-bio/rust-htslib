// Copyright 2017 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::{vec_deque, VecDeque};
use std::mem;
use std::error::Error;

use bcf;


/// A buffer for BCF records. This allows access regions in a sorted BCF file while iterating
/// over it in a single pass.
/// The buffer is implemented as a ringbuffer, such that extension or movement to the right has
/// linear complexity. The buffer does not use any indexed random access. Hence, for getting a
/// region at the very end of the BCF, you will have to wait until all records before have
/// been read.
#[derive(Debug)]
pub struct RecordBuffer {
    reader: bcf::Reader,
    ringbuffer: VecDeque<bcf::Record>,
    ringbuffer2: VecDeque<bcf::Record>,
    overflow: Option<bcf::Record>
}


unsafe impl Sync for RecordBuffer {}
unsafe impl Send for RecordBuffer {}


impl RecordBuffer {
    /// Create new buffer.
    pub fn new(reader: bcf::Reader) -> Self {
        RecordBuffer {
            reader: reader,
            ringbuffer: VecDeque::new(),
            ringbuffer2: VecDeque::new(),
            overflow: None
        }
    }

    fn last_rid(&self) -> Option<u32> {
        self.ringbuffer.back().map(|rec| rec.rid().unwrap())
    }

    fn next_rid(&self) -> Option<u32> {
        self.ringbuffer2.back().map(|rec| rec.rid().unwrap())
    }

    fn swap_buffers(&mut self) {
        // swap with buffer for next rid
        mem::swap(&mut self.ringbuffer2, &mut self.ringbuffer);
        // clear second buffer
        self.ringbuffer2.clear();
    }

    fn drain_left(&mut self, rid: u32, window_start: u32) -> usize {
        // remove records too far left or from wrong rid
        // rec.rid() will always yield Some(), because otherwise we won't put the rec into the
        // buffer.
        let to_remove = self.ringbuffer.iter().take_while(
            |rec| rec.pos() < window_start || rec.rid().unwrap() != rid
        ).count();
        self.ringbuffer.drain(..to_remove);
        to_remove
    }

    /// Fill the buffer with variants in the given window. The start coordinate has to be right of
    /// the start coordinate of any previous `fill` operation.
    /// Coordinates are 0-based, and end is exclusive.
    /// Returns tuple with numbers of added and deleted records compared to previous fetch.
    pub fn fetch(&mut self, chrom: &[u8], start: u32, end: u32) -> Result<(usize, usize), Box<Error>> {
        // TODO panic if start is left of previous start or we have moved past the given chrom
        // before.
        let rid = try!(self.reader.header.name2rid(chrom));
        let mut added = 0;
        let mut deleted = 0;

        // shrink and swap
        match (self.last_rid(), self.next_rid()) {
            (Some(last_rid), _) => {
                if last_rid != rid {
                    deleted = self.ringbuffer.len();
                    self.swap_buffers();
                    added = self.ringbuffer.len();
                    // TODO drain left?
                } else {
                    deleted = self.drain_left(rid, start);
                }
            },
            (_, Some(_)) => {
                // TODO is this really necessary? If there was no fetch before, there is nothing
                // to delete.
                deleted = self.ringbuffer.len();
                self.swap_buffers();
                deleted += self.drain_left(rid, start);
                added = self.ringbuffer.len();
            },
            _ => ()
        }

        if !self.ringbuffer2.is_empty() {
            // We have already read beyond the current rid. Hence we can't extend to the right for
            // this rid.
            return Ok((added, deleted))
        }

        // move overflow from last fill into ringbuffer
        if self.overflow.is_some() {
            let pos = self.overflow.as_ref().unwrap().pos();
            if pos >= start {
                if pos <= end {
                    self.ringbuffer.push_back(self.overflow.take().unwrap());
                    added += 1;
                } else {
                    return Ok((added, deleted));
                }
            } else {
                // discard overflow
                self.overflow.take();
            }
        }

        // extend to the right
        loop {
            let mut rec = self.reader.empty_record();
            if let Err(e) = self.reader.read(&mut rec) {
                if e.is_eof() {
                    break;
                }
                return Err(Box::new(e));
            }
            let pos = rec.pos();
            if let Some(rec_rid) = rec.rid() {
                if rec_rid == rid {
                    if pos >= end {
                        // Record is beyond our window. Store it anyways but stop.
                        self.overflow = Some(rec);
                        break;
                    } else if pos >= start {
                        // Record is within our window.
                        self.ringbuffer.push_back(rec);
                        added += 1;
                    } else {
                        // Record is upstream of our window, ignore it
                        continue
                    }
                } else if rec_rid > rid {
                    // record comes from next rid. Store it in second buffer but stop filling.
                    self.ringbuffer2.push_back(rec);
                    break;
                } else {
                    // Record comes from previous rid. Ignore it.
                    continue;
                }
            } else {
                // skip records without proper rid
                continue;
            }
        }

        Ok((added, deleted))
    }

    /// Iterate over records that have been fetched with `fetch`.
    pub fn iter(&self) -> vec_deque::Iter<bcf::Record> {
        self.ringbuffer.iter()
    }

    /// Iterate over mutable references to records that have been fetched with `fetch`.
    pub fn iter_mut(&mut self) -> vec_deque::IterMut<bcf::Record> {
        self.ringbuffer.iter_mut()
    }

    pub fn len(&self) -> usize {
        self.ringbuffer.len()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use bcf;
    use itertools::Itertools;

    #[test]
    fn test_buffer() {
        let reader = bcf::Reader::from_path(&"test/test.bcf").unwrap();
        let mut buffer = RecordBuffer::new(reader);

        buffer.fetch(b"1", 100, 10023).unwrap();
        {
            let records = buffer.iter().collect_vec();
            assert_eq!(records.len(), 2);
            assert_eq!(records[0].pos(), 10021);
            assert_eq!(records[1].pos(), 10022);
        }

        buffer.fetch(b"1", 10023, 10024).unwrap();
        {
            let records = buffer.iter().collect_vec();
            assert_eq!(records.len(), 1);
            assert_eq!(records[0].pos(), 10023);
        }
    }
}
