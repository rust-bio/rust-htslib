// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


pub struct HeaderRecord<'a> {
    rec_type: &'a [u8],
    tags: Vec<(&'a [u8], Vec<u8>)>,
}


impl<'a> HeaderRecord<'a> {
    pub fn new(rec_type: &'a [u8]) -> Self {
        HeaderRecord { rec_type: rec_type, tags: Vec::new() }
    }

    pub fn push_tag<V: ToString>(&mut self, tag: &'a [u8], value: &V) -> &mut Self {
        self.tags.push((tag, value.to_string().into_bytes()));
        self
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.push_all(self.rec_type);
        for &(tag, ref value) in self.tags.iter() {
            out.push(b'\t');
            out.push_all(tag);
            out.push_all(&value);
        }
        out
    }
}


pub struct Header {
    records: Vec<Vec<u8>>
}


impl Header {
    pub fn new() -> Self {
        Header { records: Vec::new() }
    }

    pub fn push_record(&mut self, record: HeaderRecord) {
        self.records.push(record.to_bytes());
    }

    pub fn push_comment(&mut self, comment: &[u8]) {
        self.records.push([b"@CO", comment].connect(&b'\t'));
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        self.records.connect(&b'\n')
    }
}
