// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::bam::HeaderView;
use crate::errors::{Error, Result};
use lazy_static::lazy_static;
use linear_map::LinearMap;
use regex::Regex;
use std::borrow::Cow;
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
    /// Comment lines can be obtained by the `comments` function.
    pub fn to_hashmap(&self) -> Result<HashMap<String, Vec<LinearMap<String, String>>>> {
        let mut header_map = HashMap::default();

        lazy_static! {
            static ref REC_TYPE_RE: Regex = Regex::new(r"@([A-Z][A-Z])").unwrap();
            static ref TAG_RE: Regex = Regex::new(r"([A-Za-z][A-Za-z0-9]):([ -~]*)").unwrap();
        }
        if let Ok(header_string) = String::from_utf8(self.to_bytes()) {
            for line in header_string.split('\n').filter(|x| !x.is_empty()) {
                let parts: Vec<_> = line.split('\t').filter(|x| !x.is_empty()).collect();
                if parts.is_empty() {
                    continue;
                }
                let record_type = REC_TYPE_RE
                    .captures(parts[0])
                    .and_then(|captures| captures.get(1))
                    .map(|m| m.as_str().to_owned());

                if let Some(record_type) = record_type {
                    if record_type == "CO" {
                        continue;
                    }
                    let mut field = LinearMap::default();
                    for part in parts.iter().skip(1) {
                        if let Some(cap) = TAG_RE.captures(part) {
                            let tag = cap.get(1).unwrap().as_str().to_owned();
                            let value = cap.get(2).unwrap().as_str().to_owned();
                            field.insert(tag, value);
                        } else {
                            return Err(Error::HeaderParse);
                        }
                    }
                    header_map
                        .entry(record_type)
                        .or_insert_with(Vec::new)
                        .push(field);
                } else {
                    return Err(Error::HeaderParse);
                }
            }
            Ok(header_map)
        } else {
            Err(Error::HeaderParse)
        }
    }

    /// Returns an iterator of comment lines.
    pub fn comments(&'_ self) -> impl Iterator<Item = Cow<'_, str>> {
        self.records.iter().flat_map(|r| {
            r.split(|x| x == &b'\n')
                .filter(|x| x.starts_with(b"@CO\t"))
                .map(|x| String::from_utf8_lossy(&x[4..]))
        })
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
    pub fn push_tag<V: ToString>(&mut self, tag: &'a [u8], value: V) -> &mut Self {
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

#[cfg(test)]
mod tests {
    use super::HeaderRecord;
    use crate::bam::Header;

    #[test]
    fn test_push_tag() {
        let mut record = HeaderRecord::new(b"HD");
        record.push_tag(b"X1", 0);
        record.push_tag(b"X2", 0);

        let x = "x".to_string();
        record.push_tag(b"X3", x.as_str());
        record.push_tag(b"X4", &x);
        record.push_tag(b"X5", x);

        assert_eq!(record.to_bytes(), b"@HD\tX1:0\tX2:0\tX3:x\tX4:x\tX5:x");
    }

    #[test]
    fn test_header_hash_map() {
        let mut records = Vec::new();
        let mut record = HeaderRecord::new(b"HD");
        record.push_tag(b"X1", 0);
        records.push(record);
        let mut record = HeaderRecord::new(b"PG");
        record.push_tag(b"ID", "mytool");
        records.push(record);
        let mut record = HeaderRecord::new(b"PG");
        record.push_tag(b"ID", "other_tool");
        records.push(record);
        let header = Header {
            records: records.into_iter().map(|rec| rec.to_bytes()).collect(),
        };
        let hm = header.to_hashmap().unwrap();
        assert!(hm.contains_key("HD"));
        assert!(hm.contains_key("PG"));
        assert_eq!(hm.get("HD").unwrap().len(), 1);
        assert_eq!(hm.get("PG").unwrap().len(), 2);
    }
}
