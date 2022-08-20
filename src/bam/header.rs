// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::bam::HeaderView;
use lazy_static::lazy_static;
use linear_map::LinearMap;
use regex::Regex;
use std::collections::HashMap;

/// A BAM header.
#[derive(Debug, Clone)]
pub struct Header {
    records: Vec<Vec<u8>>,
}

impl Default for Header {
    fn default() -> Self {
        Self::new()
    }
}

impl Header {
    /// Create a new header.
    pub fn new() -> Self {
        Header {
            records: Vec::new(),
        }
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
        Header {
            records: vec![record],
        }
    }

    /// Add a record to the header.
    pub fn push_record(&mut self, record: &HeaderRecord<'_>) -> &mut Self {
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

    /// This returns a header as a HashMap.
    /// Comment lines starting with "@CO" will NOT be included in the HashMap.
    /// Comment lines can be obtained by the `get_comments` function.
    pub fn to_hashmap(&self) -> HashMap<String, Vec<LinearMap<String, String>>> {
        let mut header_map = HashMap::default();

        lazy_static! {
            static ref REC_TYPE_RE: Regex = Regex::new(r"@([A-Z][A-Z])").unwrap();
            static ref TAG_RE: Regex = Regex::new(r"([A-Za-z][A-Za-z0-9]):([ -~]+)").unwrap();
        }

        let header_string = String::from_utf8(self.to_bytes()).unwrap();

        for line in header_string.split('\n').filter(|x| !x.is_empty()) {
            let parts: Vec<_> = line.split('\t').filter(|x| !x.is_empty()).collect();
            // assert!(rec_type_re.is_match(parts[0]));
            let record_type = REC_TYPE_RE
                .captures(parts[0])
                .unwrap()
                .get(1)
                .unwrap()
                .as_str()
                .to_owned();
            if record_type.eq("CO") {
                continue;
            }
            let mut field = LinearMap::default();
            for part in parts.iter().skip(1) {
                let cap = TAG_RE.captures(part).unwrap();
                let tag = cap.get(1).unwrap().as_str().to_owned();
                let value = cap.get(2).unwrap().as_str().to_owned();
                field.insert(tag, value);
            }
            header_map
                .entry(record_type)
                .or_insert_with(Vec::new)
                .push(field);
        }
        header_map
    }
    
    /// This returns comments in a header as a vector.
    /// If a header does not have any comment line, this returns `None`.
    pub fn get_comments(&self) -> Option<Vec<String>> {
        let mut vec: Vec<String> = Vec::new();
        let header_string = String::from_utf8(self.to_bytes()).unwrap();
        for line in header_string.split('\n').filter(|x| !x.is_empty()) {
            if line.starts_with("@CO\t") {
                match line.split_once("\t") {
                    Some((tag, value)) => vec.push(value.to_string()),
                    None => (),  // This should never happen.
                }
            }
        }
        if vec.len() >= 1 {
            Some(vec)
        } else {
            None
        }
    }
}

/// Header record.
#[derive(Debug, Clone)]
pub struct HeaderRecord<'a> {
    rec_type: Vec<u8>,
    tags: Vec<(&'a [u8], Vec<u8>)>,
}

impl<'a> HeaderRecord<'a> {
    /// Create a new header record.
    /// See SAM format specification for possible record types.
    pub fn new(rec_type: &'a [u8]) -> Self {
        HeaderRecord {
            rec_type: [&b"@"[..], rec_type].concat(),
            tags: Vec::new(),
        }
    }

    /// Add a new tag to the record.
    ///
    /// # Arguments
    ///
    /// * `tag` - the tag identifier
    /// * `value` - the value. Can be any type convertible into a string. Preferably numbers or
    ///   strings.
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
