use htslib;
use std::ffi::c_str_to_bytes;
use std::slice::from_raw_parts;
use std::ffi::CString;
use std::mem::copy_lifetime;


#[derive(Debug)]
pub enum Aux<'a> {
    Integer(i32),
    String(&'a [u8]),
    Float(f64),
    Char(u8),
}

impl<'a> Aux<'a> {
    pub fn string(&self) -> &'a [u8] {
        match *self {
            Aux::String(x) => x,
            _ => panic!("not a string"),
        }
    }

    pub fn float(&self) -> f64 {
        match *self {
            Aux::Float(x) => x,
            _ => panic!("not a float"),
        }
    }

    pub fn integer(&self) -> i32 {
        match *self {
            Aux::Integer(x) => x,
            _ => panic!("not an integer"),
        }
    }

    pub fn char(&self) -> u8 {
        match *self {
            Aux::Char(x) => x,
            _ => panic!("not a character"),
        }
    }
}

pub struct Record<'a> {
    b: &'a htslib::bam1_t,
    data: &'a [u8],
}

impl<'a> Record<'a>{
    pub fn new(b: &'a htslib::bam1_t) -> Record {
        Record { b: b, data: unsafe { from_raw_parts((*b).data, b.l_data as usize) } }
    }

    fn tid(&self) -> i32 {
        self.b.core.tid
    }

    fn pos(&self) -> i32 {
        self.b.core.pos
    }

    fn bin(&self) -> u16 {
        self.b.core.bin
    }

    fn map_qual(&self) -> u8 {
        self.b.core.qual
    }

    fn flag(&self) -> u16 {
        self.b.core.flag
    }

    fn mtid(&self) -> i32 {
        self.b.core.mtid
    }

    fn mpos(&self) -> i32 {
        self.b.core.mpos
    }

    fn isize(&self) -> i32 {
        self.b.core.isize
    }

    fn qname_len(&self) -> usize {
        self.b.core.l_qname as usize
    }

    pub fn qname(&self) -> &[u8] {
        self.data[0..self.qname_len()].as_slice()
    }

    fn cigar_len(&self) -> usize {
        self.b.core.n_cigar as usize
    }

    pub fn bam_get_cigar(&self) -> &[u32] {
        let x = self.b.core.l_qname as usize;
        unsafe { from_raw_parts(self.data[self.qname_len()..].as_ptr() as *const u32, self.cigar_len()) }
    }

    fn seq_len(&self) -> usize {
        self.b.core.l_qseq as usize
    }

    pub fn seq(&self) -> &[u8] {
        self.data[self.qname_len() + self.cigar_len()*4..][..self.seq_len()].as_slice()
    }

    pub fn qual(&self) -> &[u8] {
        self.data[self.qname_len() + self.cigar_len()*4 + (self.seq_len()+1)/2..][..self.seq_len()].as_slice()
    }

    pub fn aux(&self, name: &[u8]) -> Result<Aux, &str> {
        let aux = unsafe { htslib::bam_aux_get(self.b, name.as_ptr() as *mut i8 ) };
        unsafe { println!("{}", *aux); }

        unsafe {
            match *aux {
                b'c'|b'C'|b's'|b'S'|b'i'|b'I' => Ok(Aux::Integer(htslib::bam_aux2i(aux))),
                b'f'|b'd' => Ok(Aux::Float(htslib::bam_aux2f(aux))),
                b'A' => Ok(Aux::Char(htslib::bam_aux2A(aux) as u8)),
                b'Z'|b'H' => {
                    let f = (aux.offset(1) as *const i8);
                    let x = c_str_to_bytes(&f);
                    Ok(Aux::String(copy_lifetime(self, x)))
                },
                _ => Err("unexpected aux type"),
            }
        }
    }
}

pub struct Samfile<'a> {
    f: *mut htslib::Struct_BGZF,
    header: *mut htslib::bam_hdr_t,
}

impl<'a> Samfile<'a>{
     pub fn new(filename: &[u8]) -> Samfile {
        let f = unsafe { htslib::bgzf_open(filename.as_ptr() as *const i8, b"r\0".as_ptr() as *const i8) };
        let header = unsafe { htslib::bam_hdr_read(f) };
        Samfile { f : f, header : header }
    }

    pub fn read(&self) -> Record {
        let b = unsafe { htslib::bam_init1() };
        let mut status = unsafe { htslib::bam_read1(self.f, b) };
        Record::new(unsafe { &(*b) })
    }
}
