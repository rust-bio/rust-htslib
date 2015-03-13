// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::slice;
use std::iter;

use htslib;

use bam::record;


pub type Alignments<'a> = iter::Map<
    slice::Iter<'a, htslib::bam_pileup1_t>,
    fn(&'a htslib::bam_pileup1_t) -> Alignment<'a>
>;


pub struct Pileup {
    inner: Vec<htslib::bam_pileup1_t>,
    pub tid: u32,
    pub pos: u32,
}


impl Pileup {
    pub fn alignments(&self) -> Alignments {
        self.inner.iter().map(Alignment::new)
    }
}


pub struct Alignment<'a> {
    inner: &'a htslib::bam_pileup1_t,
}


impl<'a> Alignment<'a> {
    pub fn new(inner: &'a htslib::bam_pileup1_t) -> Self {
        Alignment { inner: inner }
    }

    /// Position within the read.
    pub fn qpos(&self) -> u32 {
        self.inner.qpos as u32
    }

    /// Insertion, deletion (with length) or None.
    pub fn indel(&self) -> Indel {
        match self.inner.indel {
            len if len < 0 => Indel::Del(-len as u32),
            len if len > 0 => Indel::Ins(len as u32),
            _              => Indel::None
        }
    }

    pub fn record(&self) -> record::Record {
        record::Record::from_inner(self.inner.b)
    }
}


pub enum Indel {
    Ins(u32),
    Del(u32),
    None
}


pub struct Pileups {
    itr: htslib::bam_plp_t,
}


impl Pileups {
    pub fn new(itr: htslib::bam_plp_t) -> Self {
        Pileups { itr: itr }
    }
}


impl Iterator for Pileups {
    type Item = Pileup;

    fn next(&mut self) -> Option<Pileup> {
        let (mut tid, mut pos, mut len) = (0i32, 0i32, 0i32);
        let inner = unsafe {
            htslib::bam_plp_auto(self.itr, &mut tid, &mut pos, &mut len)
        };

        if inner.is_null() {
            None
        }
        else {
            Some(Pileup {
                inner: unsafe {
                    Vec::from_raw_parts(inner as *mut htslib::bam_pileup1_t, len as usize, len as usize)
                },
                tid: tid as u32,
                pos: pos as u32,
            })
        }
    }
}


impl Drop for Pileups {
    fn drop(&mut self) {
        unsafe {
            htslib::bam_plp_destroy(self.itr);
        }
    }
}
