// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use bam::HeaderView;


/// A BAM header.
pub struct Header {
    records: Vec<Vec<u8>>
}


impl Header {
    /// Create a new header.
    pub fn new() -> Self {
        Header { records: Vec::new() }
    }

    pub fn from_template(header: &HeaderView) -> Self {
        let mut record = header.as_bytes().to_owned();
        // Strip off any trailing newline character.
        // Otherwise there could be a blank line in the
        // header which samtools (<=1.6) will complain
        // about
        while let Some(&last_char) = record.last() {
            if last_char == b'\n' {
                record.pop();
            } else {
                break;
            }
        }
        Header { records: vec![record] }
    }

    /// Add a record to the header.
    pub fn push_record(&mut self, record: &HeaderRecord) -> &mut Self {
        self.records.push(record.to_bytes());
        self
    }

    /// Add a comment to the header.
    pub fn push_comment(&mut self, comment: &[u8]) -> &mut Self {
        self.records.push([&b"@CO"[..], comment].join(&b'\t'));
        self
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        self.records.join(&b'\n')
    }
}


/// Header record.
pub struct HeaderRecord<'a> {
    rec_type: Vec<u8>,
    tags: Vec<(&'a [u8], Vec<u8>)>,
}


impl<'a> HeaderRecord<'a> {
    /// Create a new header record.
    /// See SAM format specification for possible record types.
    pub fn new(rec_type: &'a [u8]) -> Self {
        HeaderRecord { rec_type: [&b"@"[..], rec_type].concat(), tags: Vec::new() }
    }

    /// Add a new tag to the record.
    ///
    /// # Arguments
    ///
    /// * `tag` - the tag identifier
    /// * `value` - the value. Can be any type convertible into a string. Preferably numbers or strings.
    pub fn push_tag<V: ToString>(&mut self, tag: &'a [u8], value: &V) -> &mut Self {
        self.tags.push((tag, value.to_string().into_bytes()));
        self
    }

    fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend(self.rec_type.iter());
        for &(tag, ref value) in self.tags.iter() {
            out.push(b'\t');
            out.extend(tag.iter());
            out.push(b':');
            out.extend(value.iter());
        }
        out
    }
}
