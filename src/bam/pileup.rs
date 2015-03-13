// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::slice;
use std::iter;

use htslib;

use bam::record;


/// Iterator over alignments of a pileup.
pub type Alignments<'a> = iter::Map<
    slice::Iter<'a, htslib::bam_pileup1_t>,
    fn(&'a htslib::bam_pileup1_t) -> Alignment<'a>
>;


/// A pileup over one genomic position.
pub struct Pileup {
    inner: *const htslib::bam_pileup1_t,
    len: usize,
    pub tid: u32,
    pub pos: u32,
}


impl Pileup {
    pub fn alignments(&self) -> Alignments {
        self.inner().iter().map(Alignment::new)
    }

    fn inner(&self) -> &[htslib::bam_pileup1_t] {
        unsafe { slice::from_raw_parts(self.inner as *mut htslib::bam_pileup1_t, self.len) }
    }
}


impl Drop for Pileup {
    fn drop(&mut self) {
        // TODO think about what to drop here, the following causes a double free
        //for &a in self.inner().iter() {
        //    unsafe { htslib::bam_destroy1(a.b); }
        //}
    }
}


/// An aligned read in a pileup.
pub struct Alignment<'a> {
    inner: &'a htslib::bam_pileup1_t,
}


impl<'a> Alignment<'a> {
    pub fn new(inner: &'a htslib::bam_pileup1_t) -> Self {
        Alignment { inner: inner }
    }

    /// Position within the read.
    pub fn qpos(&self) -> usize {
        self.inner.qpos as usize
    }

    /// Insertion, deletion (with length) or None if no indel.
    pub fn indel(&self) -> Indel {
        match self.inner.indel {
            len if len < 0 => Indel::Del(-len as u32),
            len if len > 0 => Indel::Ins(len as u32),
            _              => Indel::None
        }
    }

    /// The corresponding record.
    pub fn record(&self) -> record::Record {
        record::Record::from_inner(self.inner.b)
    }
}


#[derive(PartialEq)]
#[derive(Debug)]
pub enum Indel {
    Ins(u32),
    Del(u32),
    None
}


/// Iterator over pileups.
pub struct Pileups {
    itr: htslib::bam_plp_t,
}


impl Pileups {
    pub fn new(itr: htslib::bam_plp_t) -> Self {
        Pileups { itr: itr }
    }
}


impl Iterator for Pileups {
    type Item = Result<Pileup, ()>;

    fn next(&mut self) -> Option<Result<Pileup, ()>> {
        let (mut tid, mut pos, mut len) = (0i32, 0i32, 0i32);
        let inner = unsafe {
            htslib::bam_plp_auto(self.itr, &mut tid, &mut pos, &mut len)
        };
        //let x = unsafe {
        //    slice::from_raw_parts(inner as *mut htslib::bam_pileup1_t, len as usize)
        //};
        //println!("{:?} {} {:?}", inner, len, x.len());

        //return None;
        match inner.is_null() {
            true if len == -1 => Some(Err(())),
            true              => None,
            false             => Some(Ok(
                    Pileup {
                        inner: inner,
                        len: len as usize,
                        tid: tid as u32,
                        pos: pos as u32,
                    }
            ))
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
