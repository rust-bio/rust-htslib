// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


use htslib;


pub struct Pileup {
    inner: Vec<htslib::bam_pileup1_t>,
    pub tid: u32,
    pub pos: u32,
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
