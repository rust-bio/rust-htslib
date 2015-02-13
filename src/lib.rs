#![feature(libc)]
extern crate libc;
pub mod htslib;
pub mod bam;

#[cfg(test)]
mod tests {
    use htslib;
    use bam;
    use std::ffi::c_str_to_bytes;

    #[test]
    fn test1() {
        let version = unsafe { *htslib::hts_version() };    
        let f = unsafe { htslib::bgzf_open(b"/vol/huge/exomate/pipeline2/bams/M46539TCS.bam\0".as_ptr() as *const i8, b"r\0".as_ptr() as *const i8) };
        let header = unsafe { htslib::bam_hdr_read(f) };
        let b = unsafe { htslib::bam_init1() };
        let mut status = unsafe { htslib::bam_read1(f, b) };

        let x = unsafe { (*b).data as *const i8 };
        let qname = unsafe { c_str_to_bytes(&x) };
        unsafe { println!("{:?}", qname) };


        //let aux = unsafe { htslib::bam_aux_get(b, b"RG".as_ptr() as *mut i8 ) };
        //unsafe { println!("{:?}", from_raw_parts(aux, 5)) };

        unsafe { htslib::bam_destroy1(b) };
        //unsafe { hts_close(f) };

        //let idx = unsafe { hts_idx_load(b"/vol/home/schroeder/smaller.bam.bai\0".as_ptr() as *const i8, 1) };
        //let h = unsafe { sam_hdr_read(fq) };
        //let iter = unsafe { sam_itr_queryi(idx, 0, 100000, 200000) };
        //let mut null = 0us;
        //let r =  unsafe { hts_itr_next(*(*f).fp.bgzf(), iter, b as *mut libc::types::common::c95::c_void, null as *mut libc::types::common::c95::c_void) };


        //println!("{}", "hallo");

        //unsafe { hts_itr_destroy(iter) };
        //unsafe { bam_destroy1(b) };
        //unsafe { hts_close(f) };
    }

    #[test]
    fn test_record() {
        let f = bam::Samfile::new(b"/vol/huge/exomate/pipeline2/bams/M46539TCS.bam");
        for record in f.take(10) {
            println!("{:?}", String::from_utf8_lossy(record.qname()));
        }


//        println!("{}", String::from_utf8_lossy(record.aux(b"MD").ok().unwrap().string()).as_slice());
//        println!("{}", record.aux(b"SM").ok().unwrap().integer());

        assert!(false);
    }
}
