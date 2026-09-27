use super::header::HeaderRecord;
use super::record::{Aux, Cigar, CigarString};
use super::*;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::str;

type GoldType = (
    [&'static [u8]; 6],
    [u16; 6],
    [&'static [u8]; 6],
    [&'static [u8]; 6],
    [CigarString; 6],
);
fn gold() -> GoldType {
    let names = [
        &b"I"[..],
        &b"II.14978392"[..],
        &b"III"[..],
        &b"IV"[..],
        &b"V"[..],
        &b"VI"[..],
    ];
    let flags = [16u16, 16u16, 16u16, 16u16, 16u16, 2048u16];
    let seqs = [
        &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCC\
TAAGCCTAAGCCTAAGCCTAA"[..],
        &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCC\
TAAGCCTAAGCCTAAGCCTAA"[..],
        &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCC\
TAAGCCTAAGCCTAAGCCTAA"[..],
        &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCC\
TAAGCCTAAGCCTAAGCCTAA"[..],
        &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCC\
TAAGCCTAAGCCTAAGCCTAA"[..],
        &b"ACTAAGCCTAAGCCTAAGCCTAAGCCAATTATCGATTTCTGAAAAAATTATCGAATTTTCTAGAAATTTTGCAAATTTT\
TTCATAAAATTATCGATTTTA"[..],
    ];
    let quals = [
        &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
CCCCCCCCCCCCCCCCCCC"[..],
        &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
CCCCCCCCCCCCCCCCCCC"[..],
        &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
CCCCCCCCCCCCCCCCCCC"[..],
        &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
CCCCCCCCCCCCCCCCCCC"[..],
        &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
CCCCCCCCCCCCCCCCCCC"[..],
        &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
CCCCCCCCCCCCCCCCCCC"[..],
    ];
    let cigars = [
        CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
        CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
        CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
        CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
        CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
        CigarString(vec![Cigar::Match(27), Cigar::Del(100000), Cigar::Match(73)]),
    ];
    (names, flags, seqs, quals, cigars)
}

fn compare_inner_bam_cram_records(cram_records: &[Record], bam_records: &[Record]) {
    // Selectively compares bam1_t struct fields from BAM and CRAM
    for (c1, b1) in cram_records.iter().zip(bam_records.iter()) {
        // CRAM vs BAM l_data is off by 3, see: https://github.com/rust-bio/rust-htslib/pull/184#issuecomment-590133544
        // The rest of the fields should be identical:
        assert_eq!(c1.cigar(), b1.cigar());
        assert_eq!(c1.inner().core.pos, b1.inner().core.pos);
        assert_eq!(c1.inner().core.mpos, b1.inner().core.mpos);
        assert_eq!(c1.inner().core.mtid, b1.inner().core.mtid);
        assert_eq!(c1.inner().core.tid, b1.inner().core.tid);
        assert_eq!(c1.inner().core.bin, b1.inner().core.bin);
        assert_eq!(c1.inner().core.qual, b1.inner().core.qual);
        assert_eq!(c1.inner().core.l_extranul, b1.inner().core.l_extranul);
        assert_eq!(c1.inner().core.flag, b1.inner().core.flag);
        assert_eq!(c1.inner().core.l_qname, b1.inner().core.l_qname);
        assert_eq!(c1.inner().core.n_cigar, b1.inner().core.n_cigar);
        assert_eq!(c1.inner().core.l_qseq, b1.inner().core.l_qseq);
        assert_eq!(c1.inner().core.isize_, b1.inner().core.isize_);
        //... except m_data
    }
}

#[test]
fn test_read() {
    let (names, flags, seqs, quals, cigars) = gold();
    let mut bam = Reader::from_path(Path::new("test/test.bam")).expect("Error opening file.");
    let del_len = [1, 1, 1, 1, 1, 100000];

    for (i, record) in bam.records().enumerate() {
        let rec = record.expect("Expected valid record");
        assert_eq!(rec.qname(), names[i]);
        assert_eq!(rec.flags(), flags[i]);
        assert_eq!(rec.seq().as_bytes(), seqs[i]);

        let cigar = rec.cigar();
        assert_eq!(*cigar, cigars[i]);

        let end_pos = cigar.end_pos();
        assert_eq!(end_pos, rec.pos() + 100 + del_len[i]);
        assert_eq!(
            cigar
                .read_pos(end_pos as u32 - 10, false, false)
                .unwrap()
                .unwrap(),
            90
        );
        assert_eq!(
            cigar
                .read_pos(rec.pos() as u32 + 20, false, false)
                .unwrap()
                .unwrap(),
            20
        );
        assert_eq!(cigar.read_pos(4000000, false, false).unwrap(), None);
        // fix qual offset
        let qual: Vec<u8> = quals[i].iter().map(|&q| q - 33).collect();
        assert_eq!(rec.qual(), &qual[..]);
    }
}

#[test]
fn test_seek() {
    let mut bam = Reader::from_path(Path::new("test/test.bam")).expect("Error opening file.");

    let mut names_by_voffset = HashMap::new();

    let mut offset = bam.tell();
    let mut rec = Record::new();
    while let Some(r) = bam.read(&mut rec) {
        r.expect("error reading bam");
        let qname = str::from_utf8(rec.qname()).unwrap().to_string();
        println!("{} {}", offset, qname);
        names_by_voffset.insert(offset, qname);
        offset = bam.tell();
    }

    for (offset, qname) in names_by_voffset.iter() {
        println!("{} {}", offset, qname);
        bam.seek(*offset).unwrap();
        if let Some(r) = bam.read(&mut rec) {
            r.unwrap();
        };
        let rec_qname = str::from_utf8(rec.qname()).unwrap().to_string();
        assert_eq!(qname, &rec_qname);
    }
}

#[test]
fn test_read_sam_header() {
    let bam = Reader::from_path("test/test.bam").expect("Error opening file.");

    let true_header = "@SQ\tSN:CHROMOSOME_I\tLN:15072423\n@SQ\tSN:CHROMOSOME_II\tLN:15279345\
             \n@SQ\tSN:CHROMOSOME_III\tLN:13783700\n@SQ\tSN:CHROMOSOME_IV\tLN:17493793\n@SQ\t\
             SN:CHROMOSOME_V\tLN:20924149\n"
        .to_string();
    let header_text = String::from_utf8(bam.header.as_bytes().to_owned()).unwrap();
    assert_eq!(header_text, true_header);
}

#[test]
fn test_read_against_sam() {
    let mut bam = Reader::from_path("./test/bam2sam_out.sam").unwrap();
    for read in bam.records() {
        let _read = read.unwrap();
    }
}

fn _test_read_indexed_common(mut bam: IndexedReader) {
    let (names, flags, seqs, quals, cigars) = gold();
    let sq_1 = b"CHROMOSOME_I";
    let sq_2 = b"CHROMOSOME_II";
    let tid_1 = bam.header.tid(sq_1).expect("Expected tid.");
    let tid_2 = bam.header.tid(sq_2).expect("Expected tid.");
    assert!(bam.header.target_len(tid_1).expect("Expected target len.") == 15072423);

    // fetch to position containing reads
    bam.fetch((tid_1, 0, 2))
        .expect("Expected successful fetch.");
    assert!(bam.records().count() == 6);

    // compare reads
    bam.fetch((tid_1, 0, 2))
        .expect("Expected successful fetch.");
    for (i, record) in bam.records().enumerate() {
        let rec = record.expect("Expected valid record");

        println!("{}", str::from_utf8(rec.qname()).unwrap());
        assert_eq!(rec.qname(), names[i]);
        assert_eq!(rec.flags(), flags[i]);
        assert_eq!(rec.seq().as_bytes(), seqs[i]);
        assert_eq!(*rec.cigar(), cigars[i]);
        // fix qual offset
        let qual: Vec<u8> = quals[i].iter().map(|&q| q - 33).collect();
        assert_eq!(rec.qual(), &qual[..]);
        assert_eq!(rec.aux(b"X"), Err(Error::BamAuxStringError));
        assert_eq!(rec.aux(b"NotAvailableAux"), Err(Error::BamAuxTagNotFound));
    }

    // fetch to empty position
    bam.fetch((tid_2, 1, 1))
        .expect("Expected successful fetch.");
    assert!(bam.records().count() == 0);

    // repeat with byte-string based fetch

    // fetch to position containing reads
    // using coordinate-string chr:start-stop
    bam.fetch(format!("{}:{}-{}", str::from_utf8(sq_1).unwrap(), 0, 2).as_bytes())
        .expect("Expected successful fetch.");
    assert!(bam.records().count() == 6);
    // using &str and exercising some of the coordinate conversion funcs
    bam.fetch((str::from_utf8(sq_1).unwrap(), 0_u32, 2_u64))
        .expect("Expected successful fetch.");
    assert!(bam.records().count() == 6);
    // using a slice
    bam.fetch((&sq_1[..], 0, 2))
        .expect("Expected successful fetch.");
    assert!(bam.records().count() == 6);
    // using a literal
    bam.fetch((sq_1, 0, 2)).expect("Expected successful fetch.");
    assert!(bam.records().count() == 6);

    // using a tid
    bam.fetch((0i32, 0u32, 2i64))
        .expect("Expected successful fetch.");
    assert!(bam.records().count() == 6);
    // using a tid:u32
    bam.fetch((0u32, 0u32, 2i64))
        .expect("Expected successful fetch.");
    assert!(bam.records().count() == 6);

    // compare reads
    bam.fetch(format!("{}:{}-{}", str::from_utf8(sq_1).unwrap(), 0, 2).as_bytes())
        .expect("Expected successful fetch.");
    for (i, record) in bam.records().enumerate() {
        let rec = record.expect("Expected valid record");

        println!("{}", str::from_utf8(rec.qname()).unwrap());
        assert_eq!(rec.qname(), names[i]);
        assert_eq!(rec.flags(), flags[i]);
        assert_eq!(rec.seq().as_bytes(), seqs[i]);
        assert_eq!(*rec.cigar(), cigars[i]);
        // fix qual offset
        let qual: Vec<u8> = quals[i].iter().map(|&q| q - 33).collect();
        assert_eq!(rec.qual(), &qual[..]);
        assert_eq!(rec.aux(b"NotAvailableAux"), Err(Error::BamAuxTagNotFound));
    }

    // fetch to empty position
    bam.fetch(format!("{}:{}-{}", str::from_utf8(sq_2).unwrap(), 1, 1).as_bytes())
        .expect("Expected successful fetch.");
    assert!(bam.records().count() == 0);

    //all on a tid
    bam.fetch(0).expect("Expected successful fetch.");
    assert!(bam.records().count() == 6);
    //all on a tid:u32
    bam.fetch(0u32).expect("Expected successful fetch.");
    assert!(bam.records().count() == 6);

    //all on a tid - by &[u8]
    bam.fetch(sq_1).expect("Expected successful fetch.");
    assert!(bam.records().count() == 6);
    //all on a tid - by str
    bam.fetch(str::from_utf8(sq_1).unwrap())
        .expect("Expected successful fetch.");
    assert!(bam.records().count() == 6);

    //all reads
    bam.fetch(FetchDefinition::All)
        .expect("Expected successful fetch.");
    assert!(bam.records().count() == 6);

    //all reads
    bam.fetch(".").expect("Expected successful fetch.");
    assert!(bam.records().count() == 6);

    //all unmapped
    bam.fetch(FetchDefinition::Unmapped)
        .expect("Expected successful fetch.");
    assert_eq!(bam.records().count(), 1); // expect one 'truncade record' Record.

    bam.fetch("*").expect("Expected successful fetch.");
    assert_eq!(bam.records().count(), 1); // expect one 'truncade record' Record.
}

#[test]
fn test_read_indexed() {
    let bam = IndexedReader::from_path("test/test.bam").expect("Expected valid index.");
    _test_read_indexed_common(bam);
}

#[test]
fn test_read_indexed_cram() {
    let mut reader = IndexedReader::from_path("test/test_cram.cram").unwrap();
    reader.set_reference("test/test_cram.fa").unwrap();
    reader.fetch(("chr1", 0, 120)).unwrap();

    let mut record = Record::new();
    reader.read(&mut record).unwrap().unwrap();
    assert_eq!(record.qname(), b"chr1.1");

    drop(reader);
}

#[test]
fn test_read_indexed_different_index_name() {
    let bam = IndexedReader::from_path_and_index(
        &"test/test_different_index_name.bam",
        &"test/test.bam.bai",
    )
    .expect("Expected valid index.");
    _test_read_indexed_common(bam);
}

#[test]
fn test_set_record() {
    let (names, _, seqs, quals, cigars) = gold();

    let mut rec = record::Record::new();
    rec.set_reverse();
    rec.set(names[0], Some(&cigars[0]), seqs[0], quals[0]);
    // note: this segfaults if you push_aux() before set()
    //       because set() obliterates aux
    rec.push_aux(b"NM", Aux::I32(15)).unwrap();

    assert_eq!(rec.qname(), names[0]);
    assert_eq!(*rec.cigar(), cigars[0]);
    assert_eq!(rec.seq().as_bytes(), seqs[0]);
    assert_eq!(rec.qual(), quals[0]);
    assert!(rec.is_reverse());
    assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));
}

#[test]
fn test_set_repeated() {
    let mut rec = Record::new();
    rec.set(
        b"123",
        Some(&CigarString(vec![Cigar::Match(3)])),
        b"AAA",
        b"III",
    );
    rec.push_aux(b"AS", Aux::I32(12345)).unwrap();
    assert_eq!(rec.qname(), b"123");
    assert_eq!(rec.seq().as_bytes(), b"AAA");
    assert_eq!(rec.qual(), b"III");
    assert_eq!(rec.aux(b"AS").unwrap(), Aux::I32(12345));

    rec.set(
        b"1234",
        Some(&CigarString(vec![Cigar::SoftClip(1), Cigar::Match(3)])),
        b"AAAA",
        b"IIII",
    );
    assert_eq!(rec.qname(), b"1234");
    assert_eq!(rec.seq().as_bytes(), b"AAAA");
    assert_eq!(rec.qual(), b"IIII");
    assert_eq!(rec.aux(b"AS").unwrap(), Aux::I32(12345));

    rec.set(
        b"12",
        Some(&CigarString(vec![Cigar::Match(2)])),
        b"AA",
        b"II",
    );
    assert_eq!(rec.qname(), b"12");
    assert_eq!(rec.seq().as_bytes(), b"AA");
    assert_eq!(rec.qual(), b"II");
    assert_eq!(rec.aux(b"AS").unwrap(), Aux::I32(12345));
}

#[test]
fn test_set_qname() {
    let (names, _, seqs, quals, cigars) = gold();

    assert!(names[0] != names[1]);

    for i in 0..names.len() {
        let mut rec = record::Record::new();
        rec.set(names[i], Some(&cigars[i]), seqs[i], quals[i]);
        rec.push_aux(b"NM", Aux::I32(15)).unwrap();

        assert_eq!(rec.qname(), names[i]);
        assert_eq!(*rec.cigar(), cigars[i]);
        assert_eq!(rec.seq().as_bytes(), seqs[i]);
        assert_eq!(rec.qual(), quals[i]);
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));

        // Equal length qname
        assert!(rec.qname()[0] != b'X');
        rec.set_qname(b"X");
        assert_eq!(rec.qname(), b"X");

        // Longer qname
        let mut longer_name = names[i].to_owned().clone();
        let extension = b"BuffaloBUffaloBUFFaloBUFFAloBUFFALoBUFFALO";
        longer_name.extend(extension.iter());
        rec.set_qname(&longer_name);

        assert_eq!(rec.qname(), longer_name.as_slice());
        assert_eq!(*rec.cigar(), cigars[i]);
        assert_eq!(rec.seq().as_bytes(), seqs[i]);
        assert_eq!(rec.qual(), quals[i]);
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));

        // Shorter qname
        let shorter_name = b"42";
        rec.set_qname(shorter_name);

        assert_eq!(rec.qname(), shorter_name);
        assert_eq!(*rec.cigar(), cigars[i]);
        assert_eq!(rec.seq().as_bytes(), seqs[i]);
        assert_eq!(rec.qual(), quals[i]);
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));

        // Zero-length qname
        rec.set_qname(b"");

        assert_eq!(rec.qname(), b"");
        assert_eq!(*rec.cigar(), cigars[i]);
        assert_eq!(rec.seq().as_bytes(), seqs[i]);
        assert_eq!(rec.qual(), quals[i]);
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));
    }
}

#[test]
fn test_set_qname2() {
    let mut _header = Header::new();
    _header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "1")
            .push_tag(b"LN", 10000000),
    );
    let header = HeaderView::from_header(&_header);

    let line =
            b"blah1	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1\tli:i:0\ttf:Z:cC";

    let mut rec = Record::from_sam(&header, line).unwrap();
    assert_eq!(rec.qname(), b"blah1");
    rec.set_qname(b"r0");
    assert_eq!(rec.qname(), b"r0");
}

#[test]
fn test_set_cigar() {
    let (names, _, seqs, quals, cigars) = gold();

    assert!(names[0] != names[1]);

    for i in 0..names.len() {
        let mut rec = record::Record::new();
        rec.set(names[i], Some(&cigars[i]), seqs[i], quals[i]);
        rec.push_aux(b"NM", Aux::I32(15)).unwrap();

        assert_eq!(rec.qname(), names[i]);
        assert_eq!(*rec.cigar(), cigars[i]);
        assert_eq!(rec.seq().as_bytes(), seqs[i]);
        assert_eq!(rec.qual(), quals[i]);
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));

        // boring cigar
        let new_cigar = CigarString(vec![Cigar::Match(rec.seq_len() as u32)]);
        assert_ne!(*rec.cigar(), new_cigar);
        rec.set_cigar(Some(&new_cigar));
        assert_eq!(*rec.cigar(), new_cigar);

        assert_eq!(rec.qname(), names[i]);
        assert_eq!(rec.seq().as_bytes(), seqs[i]);
        assert_eq!(rec.qual(), quals[i]);
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));

        // bizarre cigar
        let new_cigar = (0..rec.seq_len())
            .map(|i| {
                if i % 2 == 0 {
                    Cigar::Match(1)
                } else {
                    Cigar::Ins(1)
                }
            })
            .collect::<Vec<_>>();
        let new_cigar = CigarString(new_cigar);
        assert_ne!(*rec.cigar(), new_cigar);
        rec.set_cigar(Some(&new_cigar));
        assert_eq!(*rec.cigar(), new_cigar);

        assert_eq!(rec.qname(), names[i]);
        assert_eq!(rec.seq().as_bytes(), seqs[i]);
        assert_eq!(rec.qual(), quals[i]);
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));

        // empty cigar
        let new_cigar = CigarString(Vec::new());
        assert_ne!(*rec.cigar(), new_cigar);
        rec.set_cigar(None);
        assert_eq!(*rec.cigar(), new_cigar);

        assert_eq!(rec.qname(), names[i]);
        assert_eq!(rec.seq().as_bytes(), seqs[i]);
        assert_eq!(rec.qual(), quals[i]);
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));
    }
}

#[test]
fn test_remove_aux() {
    let mut bam = Reader::from_path(Path::new("test/test.bam")).expect("Error opening file.");

    for record in bam.records() {
        let mut rec = record.expect("Expected valid record");

        if rec.aux(b"XS").is_ok() {
            rec.remove_aux(b"XS").unwrap();
        }

        if rec.aux(b"YT").is_ok() {
            rec.remove_aux(b"YT").unwrap();
        }

        assert_eq!(rec.remove_aux(b"X"), Err(Error::BamAuxStringError));
        assert_eq!(rec.remove_aux(b"ab"), Err(Error::BamAuxTagNotFound));

        assert_eq!(rec.aux(b"XS"), Err(Error::BamAuxTagNotFound));
        assert_eq!(rec.aux(b"YT"), Err(Error::BamAuxTagNotFound));
    }
}

#[test]
fn test_write() {
    let (names, _, seqs, quals, cigars) = gold();

    let tmp = tempfile::Builder::new()
        .prefix("rust-htslib")
        .tempdir()
        .expect("Cannot create temp dir");
    let bampath = tmp.path().join("test.bam");
    println!("{:?}", bampath);
    {
        let mut bam = Writer::from_path(
            &bampath,
            Header::new().push_record(
                HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", "chr1")
                    .push_tag(b"LN", 15072423),
            ),
            Format::Bam,
        )
        .expect("Error opening file.");

        for i in 0..names.len() {
            let mut rec = record::Record::new();
            rec.set(names[i], Some(&cigars[i]), seqs[i], quals[i]);
            rec.push_aux(b"NM", Aux::I32(15)).unwrap();

            bam.write(&rec).expect("Failed to write record.");
        }
    }

    {
        let mut bam = Reader::from_path(bampath).expect("Error opening file.");

        for i in 0..names.len() {
            let mut rec = record::Record::new();
            if let Some(r) = bam.read(&mut rec) {
                r.expect("Failed to read record.");
            };

            assert_eq!(rec.qname(), names[i]);
            assert_eq!(*rec.cigar(), cigars[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);
            assert_eq!(rec.qual(), quals[i]);
            assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));
        }
    }

    tmp.close().expect("Failed to delete temp dir");
}

#[test]
fn test_write_threaded() {
    let (names, _, seqs, quals, cigars) = gold();

    let tmp = tempfile::Builder::new()
        .prefix("rust-htslib")
        .tempdir()
        .expect("Cannot create temp dir");
    let bampath = tmp.path().join("test.bam");
    println!("{:?}", bampath);
    {
        let mut bam = Writer::from_path(
            &bampath,
            Header::new().push_record(
                HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", "chr1")
                    .push_tag(b"LN", 15072423),
            ),
            Format::Bam,
        )
        .expect("Error opening file.");
        bam.set_threads(4).unwrap();

        for i in 0..10000 {
            let mut rec = record::Record::new();
            let idx = i % names.len();
            rec.set(names[idx], Some(&cigars[idx]), seqs[idx], quals[idx]);
            rec.push_aux(b"NM", Aux::I32(15)).unwrap();
            rec.set_pos(i as i64);

            bam.write(&rec).expect("Failed to write record.");
        }
    }

    {
        let mut bam = Reader::from_path(bampath).expect("Error opening file.");

        for (i, _rec) in bam.records().enumerate() {
            let idx = i % names.len();

            let rec = _rec.expect("Failed to read record.");

            assert_eq!(rec.pos(), i as i64);
            assert_eq!(rec.qname(), names[idx]);
            assert_eq!(*rec.cigar(), cigars[idx]);
            assert_eq!(rec.seq().as_bytes(), seqs[idx]);
            assert_eq!(rec.qual(), quals[idx]);
            assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));
        }
    }

    tmp.close().expect("Failed to delete temp dir");
}

#[test]
fn test_write_shared_tpool() {
    let (names, _, seqs, quals, cigars) = gold();

    let tmp = tempfile::Builder::new()
        .prefix("rust-htslib")
        .tempdir()
        .expect("Cannot create temp dir");
    let bampath1 = tmp.path().join("test1.bam");
    let bampath2 = tmp.path().join("test2.bam");

    {
        let (mut bam1, mut bam2) = {
            let pool = crate::tpool::ThreadPool::new(4).unwrap();

            let mut bam1 = Writer::from_path(
                &bampath1,
                Header::new().push_record(
                    HeaderRecord::new(b"SQ")
                        .push_tag(b"SN", "chr1")
                        .push_tag(b"LN", 15072423),
                ),
                Format::Bam,
            )
            .expect("Error opening file.");

            let mut bam2 = Writer::from_path(
                &bampath2,
                Header::new().push_record(
                    HeaderRecord::new(b"SQ")
                        .push_tag(b"SN", "chr1")
                        .push_tag(b"LN", 15072423),
                ),
                Format::Bam,
            )
            .expect("Error opening file.");

            bam1.set_thread_pool(&pool).unwrap();
            bam2.set_thread_pool(&pool).unwrap();
            (bam1, bam2)
        };

        for i in 0..10000 {
            let mut rec = record::Record::new();
            let idx = i % names.len();
            rec.set(names[idx], Some(&cigars[idx]), seqs[idx], quals[idx]);
            rec.push_aux(b"NM", Aux::I32(15)).unwrap();
            rec.set_pos(i as i64);

            bam1.write(&rec).expect("Failed to write record.");
            bam2.write(&rec).expect("Failed to write record.");
        }
    }

    {
        let pool = crate::tpool::ThreadPool::new(2).unwrap();

        for p in [bampath1, bampath2] {
            let mut bam = Reader::from_path(p).expect("Error opening file.");
            bam.set_thread_pool(&pool).unwrap();

            for (i, _rec) in bam.iter_chunk(None, None).enumerate() {
                let idx = i % names.len();

                let rec = _rec.expect("Failed to read record.");

                assert_eq!(rec.pos(), i as i64);
                assert_eq!(rec.qname(), names[idx]);
                assert_eq!(*rec.cigar(), cigars[idx]);
                assert_eq!(rec.seq().as_bytes(), seqs[idx]);
                assert_eq!(rec.qual(), quals[idx]);
                assert_eq!(rec.aux(b"NM").unwrap(), Aux::I32(15));
            }
        }
    }

    tmp.close().expect("Failed to delete temp dir");
}

#[test]
fn test_copy_template() {
    // Verify that BAM headers are transmitted correctly when using an existing BAM as a
    // template for headers.

    let tmp = tempfile::Builder::new()
        .prefix("rust-htslib")
        .tempdir()
        .expect("Cannot create temp dir");
    let bampath = tmp.path().join("test.bam");
    println!("{:?}", bampath);

    let mut input_bam = Reader::from_path("test/test.bam").expect("Error opening file.");

    {
        let mut bam = Writer::from_path(
            &bampath,
            &Header::from_template(input_bam.header()),
            Format::Bam,
        )
        .expect("Error opening file.");

        for rec in input_bam.records() {
            bam.write(&rec.unwrap()).expect("Failed to write record.");
        }
    }

    {
        let copy_bam = Reader::from_path(bampath).expect("Error opening file.");

        // Verify that the header came across correctly
        assert_eq!(input_bam.header().as_bytes(), copy_bam.header().as_bytes());
    }

    tmp.close().expect("Failed to delete temp dir");
}

#[test]
fn test_pileup() {
    let (_, _, seqs, quals, _) = gold();

    let mut bam = Reader::from_path("test/test.bam").expect("Error opening file.");
    let pileups = bam.pileup();
    for pileup in pileups.take(26) {
        let _pileup = pileup.expect("Expected successful pileup.");
        let pos = _pileup.pos() as usize;
        assert_eq!(_pileup.depth(), 6);
        assert!(_pileup.tid() == 0);
        for (i, a) in _pileup.alignments().enumerate() {
            assert_eq!(a.indel(), pileup::Indel::None);
            let qpos = a.qpos().unwrap();
            assert_eq!(qpos, pos - 1);
            assert_eq!(a.record().seq()[qpos], seqs[i][qpos]);
            assert_eq!(a.record().qual()[qpos], quals[i][qpos] - 33);
        }
    }
}

#[test]
fn test_idx_pileup() {
    let mut bam = IndexedReader::from_path("test/test.bam").expect("Error opening file.");
    // read without fetch
    for pileup in bam.pileup() {
        pileup.unwrap();
    }
    // go back again
    let tid = bam.header().tid(b"CHROMOSOME_I").unwrap();
    bam.fetch((tid, 0, 5)).unwrap();
    for p in bam.pileup() {
        println!("{}", p.unwrap().pos())
    }
}

#[test]
fn parse_from_sam() {
    use std::fs::File;
    use std::io::Read;

    let bamfile = "./test/bam2sam_test.bam";
    let samfile = "./test/bam2sam_expected.sam";

    // Load BAM file:
    let mut rdr = Reader::from_path(bamfile).unwrap();
    let bam_recs: Vec<Record> = rdr.records().map(|v| v.unwrap()).collect();

    let mut sam = Vec::new();
    assert!(File::open(samfile).unwrap().read_to_end(&mut sam).is_ok());

    let sam_recs: Vec<Record> = sam
        .split(|x| *x == b'\n')
        .filter(|x| !x.is_empty() && x[0] != b'@')
        .map(|line| Record::from_sam(rdr.header(), line).unwrap())
        .collect();

    for (b1, s1) in bam_recs.iter().zip(sam_recs.iter()) {
        assert!(b1 == s1);
    }
}

#[test]
fn test_cigar_modes() {
    // test the cached and uncached ways of getting the cigar string.

    let (_, _, _, _, cigars) = gold();
    let mut bam = Reader::from_path(Path::new("test/test.bam")).expect("Error opening file.");

    for (i, record) in bam.records().enumerate() {
        let rec = record.expect("Expected valid record");

        let cigar = rec.cigar();
        assert_eq!(*cigar, cigars[i]);
    }

    for (i, record) in bam.records().enumerate() {
        let mut rec = record.expect("Expected valid record");
        rec.cache_cigar();

        let cigar = rec.cigar_cached().unwrap();
        assert_eq!(**cigar, cigars[i]);

        let cigar = rec.cigar();
        assert_eq!(*cigar, cigars[i]);
    }
}

#[test]
fn test_read_cram() {
    let cram_path = "./test/test_cram.cram";
    let bam_path = "./test/test_cram.bam";
    let ref_path = "./test/test_cram.fa";

    // Load CRAM file, records
    let mut cram_reader = Reader::from_path(cram_path).unwrap();
    cram_reader.set_reference(ref_path).unwrap();
    let cram_records: Vec<Record> = cram_reader.records().map(|v| v.unwrap()).collect();

    // Load BAM file, records
    let mut bam_reader = Reader::from_path(bam_path).unwrap();
    let bam_records: Vec<Record> = bam_reader.records().map(|v| v.unwrap()).collect();

    compare_inner_bam_cram_records(&cram_records, &bam_records);
}

#[test]
fn test_write_cram() {
    // BAM file, records
    let bam_path = "./test/test_cram.bam";
    let ref_path = "./test/test_cram.fa";
    let mut bam_reader = Reader::from_path(bam_path).unwrap();
    let bam_records: Vec<Record> = bam_reader.records().map(|v| v.unwrap()).collect();

    // New CRAM file
    let tmp = tempfile::Builder::new()
        .prefix("rust-htslib")
        .tempdir()
        .expect("Cannot create temp dir");
    let cram_path = tmp.path().join("test.cram");

    // Write BAM records to new CRAM file
    {
        let mut header = Header::new();
        header.push_record(
            HeaderRecord::new(b"HD")
                .push_tag(b"VN", "1.5")
                .push_tag(b"SO", "coordinate"),
        );
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 120)
                .push_tag(b"M5", "20a9a0fb770814e6c5e49946750f9724")
                .push_tag(b"UR", "test/test_cram.fa"),
        );
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr2")
                .push_tag(b"LN", 120)
                .push_tag(b"M5", "7a2006ccca94ea92b6dae5997e1b0d70")
                .push_tag(b"UR", "test/test_cram.fa"),
        );
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr3")
                .push_tag(b"LN", 120)
                .push_tag(b"M5", "a66b336bfe3ee8801c744c9545c87e24")
                .push_tag(b"UR", "test/test_cram.fa"),
        );

        let mut cram_writer =
            Writer::from_path(&cram_path, &header, Format::Cram).expect("Error opening CRAM file.");
        cram_writer.set_reference(ref_path).unwrap();

        // Write BAM records to CRAM file
        for rec in bam_records.iter() {
            cram_writer
                .write(rec)
                .expect("Faied to write record to CRAM.");
        }
    }

    // Compare written CRAM records with BAM records
    {
        // Load written CRAM file
        let mut cram_reader = Reader::from_path(cram_path).unwrap();
        cram_reader.set_reference(ref_path).unwrap();
        let cram_records: Vec<Record> = cram_reader.records().map(|v| v.unwrap()).collect();

        // Compare CRAM records to BAM records
        compare_inner_bam_cram_records(&cram_records, &bam_records);
    }

    tmp.close().expect("Failed to delete temp dir");
}

#[test]
fn test_compression_level_conversion() {
    // predefined compression levels
    assert_eq!(CompressionLevel::Uncompressed.convert().unwrap(), 0);
    assert_eq!(CompressionLevel::Fastest.convert().unwrap(), 1);
    assert_eq!(CompressionLevel::Maximum.convert().unwrap(), 9);

    // numeric compression levels
    for level in 0..=9 {
        assert_eq!(CompressionLevel::Level(level).convert().unwrap(), level);
    }
    // invalid levels
    assert!(CompressionLevel::Level(10).convert().is_err());
}

#[test]
fn test_write_compression() {
    let tmp = tempfile::Builder::new()
        .prefix("rust-htslib")
        .tempdir()
        .expect("Cannot create temp dir");
    let input_bam_path = "test/test.bam";

    // test levels with decreasing compression factor
    let levels_to_test = vec![
        CompressionLevel::Maximum,
        CompressionLevel::Level(6),
        CompressionLevel::Fastest,
        CompressionLevel::Uncompressed,
    ];
    let file_sizes: Vec<_> = levels_to_test
        .iter()
        .map(|level| {
            let output_bam_path = tmp.path().join("test.bam");
            {
                let mut reader = Reader::from_path(input_bam_path).unwrap();
                let header = Header::from_template(reader.header());
                let mut writer = Writer::from_path(&output_bam_path, &header, Format::Bam).unwrap();
                writer.set_compression_level(*level).unwrap();
                for record in reader.records() {
                    let r = record.unwrap();
                    writer.write(&r).unwrap();
                }
            }
            fs::metadata(output_bam_path).unwrap().len()
        })
        .collect();

    // check that out BAM file sizes are in decreasing order, in line with the expected compression factor
    println!("testing compression leves: {:?}", levels_to_test);
    println!("got compressed sizes: {:?}", file_sizes);

    // libdeflate comes out with a slightly bigger file on Max compression
    // than on Level(6), so skip that check
    #[cfg(feature = "libdeflate")]
    assert!(file_sizes[1..].windows(2).all(|size| size[0] <= size[1]));

    #[cfg(not(feature = "libdeflate"))]
    assert!(file_sizes.windows(2).all(|size| size[0] <= size[1]));

    tmp.close().expect("Failed to delete temp dir");
}

#[test]
fn test_bam_fails_on_vcf() {
    let bam_path = "./test/test_left.vcf";
    let bam_reader = Reader::from_path(bam_path);
    assert!(bam_reader.is_err());
}

#[test]
fn test_indexde_bam_fails_on_vcf() {
    let bam_path = "./test/test_left.vcf";
    let bam_reader = IndexedReader::from_path(bam_path);
    assert!(bam_reader.is_err());
}

#[test]
fn test_bam_fails_on_toml() {
    let bam_path = "./Cargo.toml";
    let bam_reader = Reader::from_path(bam_path);
    assert!(bam_reader.is_err());
}

#[test]
fn test_sam_writer_example() {
    fn from_bam_with_filter<F>(bamfile: &str, samfile: &str, f: F) -> bool
    where
        F: Fn(&record::Record) -> Option<bool>,
    {
        let mut bam_reader = Reader::from_path(bamfile).unwrap(); // internal functions, just unwrap
        let header = header::Header::from_template(bam_reader.header());
        let mut sam_writer = Writer::from_path(samfile, &header, Format::Sam).unwrap();
        for record in bam_reader.records() {
            if record.is_err() {
                return false;
            }
            let parsed = record.unwrap();
            match f(&parsed) {
                None => return true,
                Some(false) => {}
                Some(true) => {
                    if sam_writer.write(&parsed).is_err() {
                        return false;
                    }
                }
            }
        }
        true
    }
    use std::fs::File;
    use std::io::Read;
    let bamfile = "./test/bam2sam_test.bam";
    let samfile = "./test/bam2sam_out.sam";
    let expectedfile = "./test/bam2sam_expected.sam";
    let result = from_bam_with_filter(bamfile, samfile, |_| Some(true));
    assert!(result);
    let mut expected = Vec::new();
    let mut written = Vec::new();
    assert!(File::open(expectedfile)
        .unwrap()
        .read_to_end(&mut expected)
        .is_ok());
    assert!(File::open(samfile)
        .unwrap()
        .read_to_end(&mut written)
        .is_ok());
    assert_eq!(expected, written);
}

// #[cfg(feature = "curl")]
// #[test]
// fn test_http_connect() {
//     let url: Url = Url::parse(
//         "https://raw.githubusercontent.com/brainstorm/tiny-test-data/master/wgs/mt.bam",
//     )
//     .unwrap();
//     let r = Reader::from_url(&url);
//     println!("{:#?}", r);
//     let r = r.unwrap();

//     assert_eq!(r.header().target_names()[0], b"chr1");
// }

#[test]
fn test_rc_records() {
    let (names, flags, seqs, quals, cigars) = gold();
    let mut bam = Reader::from_path(Path::new("test/test.bam")).expect("Error opening file.");
    let del_len = [1, 1, 1, 1, 1, 100000];

    for (i, record) in bam.rc_records().enumerate() {
        //let rec = record.expect("Expected valid record");
        let rec = record.unwrap();
        println!("{}", str::from_utf8(rec.qname()).ok().unwrap());
        assert_eq!(rec.qname(), names[i]);
        assert_eq!(rec.flags(), flags[i]);
        assert_eq!(rec.seq().as_bytes(), seqs[i]);

        let cigar = rec.cigar();
        assert_eq!(*cigar, cigars[i]);

        let end_pos = cigar.end_pos();
        assert_eq!(end_pos, rec.pos() + 100 + del_len[i]);
        assert_eq!(
            cigar
                .read_pos(end_pos as u32 - 10, false, false)
                .unwrap()
                .unwrap(),
            90
        );
        assert_eq!(
            cigar
                .read_pos(rec.pos() as u32 + 20, false, false)
                .unwrap()
                .unwrap(),
            20
        );
        assert_eq!(cigar.read_pos(4000000, false, false).unwrap(), None);
        // fix qual offset
        let qual: Vec<u8> = quals[i].iter().map(|&q| q - 33).collect();
        assert_eq!(rec.qual(), &qual[..]);
    }
}

#[test]
fn test_aux_arrays() {
    let bam_header = Header::new();
    let mut test_record = Record::from_sam(
        &HeaderView::from_header(&bam_header),
        "ali1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tFFFF".as_bytes(),
    )
    .unwrap();

    let array_i8: Vec<i8> = vec![i8::MIN, -1, 0, 1, i8::MAX];
    let array_u8: Vec<u8> = vec![u8::MIN, 0, 1, u8::MAX];
    let array_i16: Vec<i16> = vec![i16::MIN, -1, 0, 1, i16::MAX];
    let array_u16: Vec<u16> = vec![u16::MIN, 0, 1, u16::MAX];
    let array_i32: Vec<i32> = vec![i32::MIN, -1, 0, 1, i32::MAX];
    let array_u32: Vec<u32> = vec![u32::MIN, 0, 1, u32::MAX];
    let array_f32: Vec<f32> = vec![f32::MIN, 0.0, -0.0, 0.1, 0.99, f32::MAX];

    test_record
        .push_aux(b"XA", Aux::ArrayI8((&array_i8).into()))
        .unwrap();
    test_record
        .push_aux(b"XB", Aux::ArrayU8((&array_u8).into()))
        .unwrap();
    test_record
        .push_aux(b"XC", Aux::ArrayI16((&array_i16).into()))
        .unwrap();
    test_record
        .push_aux(b"XD", Aux::ArrayU16((&array_u16).into()))
        .unwrap();
    test_record
        .push_aux(b"XE", Aux::ArrayI32((&array_i32).into()))
        .unwrap();
    test_record
        .push_aux(b"XF", Aux::ArrayU32((&array_u32).into()))
        .unwrap();
    test_record
        .push_aux(b"XG", Aux::ArrayFloat((&array_f32).into()))
        .unwrap();

    {
        let tag = b"XA";
        if let Ok(Aux::ArrayI8(array)) = test_record.aux(tag) {
            // Retrieve aux array
            let aux_array_content = array.iter().collect::<Vec<_>>();
            assert_eq!(aux_array_content, array_i8);

            // Copy the stored aux array to another record
            {
                let mut copy_test_record = test_record.clone();

                // Pushing a field with an existing tag should fail
                assert!(copy_test_record.push_aux(tag, Aux::I8(3)).is_err());

                // Remove aux array from target record
                copy_test_record.remove_aux(tag).unwrap();
                assert!(copy_test_record.aux(tag).is_err());

                // Copy array to target record
                let src_aux = test_record.aux(tag).unwrap();
                assert!(copy_test_record.push_aux(tag, src_aux).is_ok());
                if let Ok(Aux::ArrayI8(array)) = copy_test_record.aux(tag) {
                    let aux_array_content_copied = array.iter().collect::<Vec<_>>();
                    assert_eq!(aux_array_content_copied, array_i8);
                } else {
                    panic!("Aux tag not found");
                }
            }
        } else {
            panic!("Aux tag not found");
        }
    }

    {
        let tag = b"XB";
        if let Ok(Aux::ArrayU8(array)) = test_record.aux(tag) {
            // Retrieve aux array
            let aux_array_content = array.iter().collect::<Vec<_>>();
            assert_eq!(aux_array_content, array_u8);

            // Copy the stored aux array to another record
            {
                let mut copy_test_record = test_record.clone();

                // Pushing a field with an existing tag should fail
                assert!(copy_test_record.push_aux(tag, Aux::U8(3)).is_err());

                // Remove aux array from target record
                copy_test_record.remove_aux(tag).unwrap();
                assert!(copy_test_record.aux(tag).is_err());

                // Copy array to target record
                let src_aux = test_record.aux(tag).unwrap();
                assert!(copy_test_record.push_aux(tag, src_aux).is_ok());
                if let Ok(Aux::ArrayU8(array)) = copy_test_record.aux(tag) {
                    let aux_array_content_copied = array.iter().collect::<Vec<_>>();
                    assert_eq!(aux_array_content_copied, array_u8);
                } else {
                    panic!("Aux tag not found");
                }
            }
        } else {
            panic!("Aux tag not found");
        }
    }

    {
        let tag = b"XC";
        if let Ok(Aux::ArrayI16(array)) = test_record.aux(tag) {
            // Retrieve aux array
            let aux_array_content = array.iter().collect::<Vec<_>>();
            assert_eq!(aux_array_content, array_i16);

            // Copy the stored aux array to another record
            {
                let mut copy_test_record = test_record.clone();

                // Pushing a field with an existing tag should fail
                assert!(copy_test_record.push_aux(tag, Aux::I16(3)).is_err());

                // Remove aux array from target record
                copy_test_record.remove_aux(tag).unwrap();
                assert!(copy_test_record.aux(tag).is_err());

                // Copy array to target record
                let src_aux = test_record.aux(tag).unwrap();
                assert!(copy_test_record.push_aux(tag, src_aux).is_ok());
                if let Ok(Aux::ArrayI16(array)) = copy_test_record.aux(tag) {
                    let aux_array_content_copied = array.iter().collect::<Vec<_>>();
                    assert_eq!(aux_array_content_copied, array_i16);
                } else {
                    panic!("Aux tag not found");
                }
            }
        } else {
            panic!("Aux tag not found");
        }
    }

    {
        let tag = b"XD";
        if let Ok(Aux::ArrayU16(array)) = test_record.aux(tag) {
            // Retrieve aux array
            let aux_array_content = array.iter().collect::<Vec<_>>();
            assert_eq!(aux_array_content, array_u16);

            // Copy the stored aux array to another record
            {
                let mut copy_test_record = test_record.clone();

                // Pushing a field with an existing tag should fail
                assert!(copy_test_record.push_aux(tag, Aux::U16(3)).is_err());

                // Remove aux array from target record
                copy_test_record.remove_aux(tag).unwrap();
                assert!(copy_test_record.aux(tag).is_err());

                // Copy array to target record
                let src_aux = test_record.aux(tag).unwrap();
                assert!(copy_test_record.push_aux(tag, src_aux).is_ok());
                if let Ok(Aux::ArrayU16(array)) = copy_test_record.aux(tag) {
                    let aux_array_content_copied = array.iter().collect::<Vec<_>>();
                    assert_eq!(aux_array_content_copied, array_u16);
                } else {
                    panic!("Aux tag not found");
                }
            }
        } else {
            panic!("Aux tag not found");
        }
    }

    {
        let tag = b"XE";
        if let Ok(Aux::ArrayI32(array)) = test_record.aux(tag) {
            // Retrieve aux array
            let aux_array_content = array.iter().collect::<Vec<_>>();
            assert_eq!(aux_array_content, array_i32);

            // Copy the stored aux array to another record
            {
                let mut copy_test_record = test_record.clone();

                // Pushing a field with an existing tag should fail
                assert!(copy_test_record.push_aux(tag, Aux::I32(3)).is_err());

                // Remove aux array from target record
                copy_test_record.remove_aux(tag).unwrap();
                assert!(copy_test_record.aux(tag).is_err());

                // Copy array to target record
                let src_aux = test_record.aux(tag).unwrap();
                assert!(copy_test_record.push_aux(tag, src_aux).is_ok());
                if let Ok(Aux::ArrayI32(array)) = copy_test_record.aux(tag) {
                    let aux_array_content_copied = array.iter().collect::<Vec<_>>();
                    assert_eq!(aux_array_content_copied, array_i32);
                } else {
                    panic!("Aux tag not found");
                }
            }
        } else {
            panic!("Aux tag not found");
        }
    }

    {
        let tag = b"XF";
        if let Ok(Aux::ArrayU32(array)) = test_record.aux(tag) {
            // Retrieve aux array
            let aux_array_content = array.iter().collect::<Vec<_>>();
            assert_eq!(aux_array_content, array_u32);

            // Copy the stored aux array to another record
            {
                let mut copy_test_record = test_record.clone();

                // Pushing a field with an existing tag should fail
                assert!(copy_test_record.push_aux(tag, Aux::U32(3)).is_err());

                // Remove aux array from target record
                copy_test_record.remove_aux(tag).unwrap();
                assert!(copy_test_record.aux(tag).is_err());

                // Copy array to target record
                let src_aux = test_record.aux(tag).unwrap();
                assert!(copy_test_record.push_aux(tag, src_aux).is_ok());
                if let Ok(Aux::ArrayU32(array)) = copy_test_record.aux(tag) {
                    let aux_array_content_copied = array.iter().collect::<Vec<_>>();
                    assert_eq!(aux_array_content_copied, array_u32);
                } else {
                    panic!("Aux tag not found");
                }
            }
        } else {
            panic!("Aux tag not found");
        }
    }

    {
        let tag = b"XG";
        if let Ok(Aux::ArrayFloat(array)) = test_record.aux(tag) {
            // Retrieve aux array
            let aux_array_content = array.iter().collect::<Vec<_>>();
            assert_eq!(aux_array_content, array_f32);

            // Copy the stored aux array to another record
            {
                let mut copy_test_record = test_record.clone();

                // Pushing a field with an existing tag should fail
                assert!(copy_test_record.push_aux(tag, Aux::Float(3.0)).is_err());

                // Remove aux array from target record
                copy_test_record.remove_aux(tag).unwrap();
                assert!(copy_test_record.aux(tag).is_err());

                // Copy array to target record
                let src_aux = test_record.aux(tag).unwrap();
                assert!(copy_test_record.push_aux(tag, src_aux).is_ok());
                if let Ok(Aux::ArrayFloat(array)) = copy_test_record.aux(tag) {
                    let aux_array_content_copied = array.iter().collect::<Vec<_>>();
                    assert_eq!(aux_array_content_copied, array_f32);
                } else {
                    panic!("Aux tag not found");
                }
            }
        } else {
            panic!("Aux tag not found");
        }
    }

    // Test via `Iterator` impl
    for item in test_record.aux_iter() {
        match item.unwrap() {
            (b"XA", Aux::ArrayI8(array)) => {
                assert_eq!(&array.iter().collect::<Vec<_>>(), &array_i8);
            }
            (b"XB", Aux::ArrayU8(array)) => {
                assert_eq!(&array.iter().collect::<Vec<_>>(), &array_u8);
            }
            (b"XC", Aux::ArrayI16(array)) => {
                assert_eq!(&array.iter().collect::<Vec<_>>(), &array_i16);
            }
            (b"XD", Aux::ArrayU16(array)) => {
                assert_eq!(&array.iter().collect::<Vec<_>>(), &array_u16);
            }
            (b"XE", Aux::ArrayI32(array)) => {
                assert_eq!(&array.iter().collect::<Vec<_>>(), &array_i32);
            }
            (b"XF", Aux::ArrayU32(array)) => {
                assert_eq!(&array.iter().collect::<Vec<_>>(), &array_u32);
            }
            (b"XG", Aux::ArrayFloat(array)) => {
                assert_eq!(&array.iter().collect::<Vec<_>>(), &array_f32);
            }
            _ => {
                panic!();
            }
        }
    }

    // Test via `PartialEq` impl
    assert_eq!(
        test_record.aux(b"XA").unwrap(),
        Aux::ArrayI8((&array_i8).into())
    );
    assert_eq!(
        test_record.aux(b"XB").unwrap(),
        Aux::ArrayU8((&array_u8).into())
    );
    assert_eq!(
        test_record.aux(b"XC").unwrap(),
        Aux::ArrayI16((&array_i16).into())
    );
    assert_eq!(
        test_record.aux(b"XD").unwrap(),
        Aux::ArrayU16((&array_u16).into())
    );
    assert_eq!(
        test_record.aux(b"XE").unwrap(),
        Aux::ArrayI32((&array_i32).into())
    );
    assert_eq!(
        test_record.aux(b"XF").unwrap(),
        Aux::ArrayU32((&array_u32).into())
    );
    assert_eq!(
        test_record.aux(b"XG").unwrap(),
        Aux::ArrayFloat((&array_f32).into())
    );
}

#[test]
fn test_aux_scalars() {
    let bam_header = Header::new();
    let mut test_record = Record::from_sam(
        &HeaderView::from_header(&bam_header),
        "ali1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tFFFF".as_bytes(),
    )
    .unwrap();

    test_record.push_aux(b"XA", Aux::I8(i8::MIN)).unwrap();
    test_record.push_aux(b"XB", Aux::I8(i8::MAX)).unwrap();
    test_record.push_aux(b"XC", Aux::U8(u8::MIN)).unwrap();
    test_record.push_aux(b"XD", Aux::U8(u8::MAX)).unwrap();
    test_record.push_aux(b"XE", Aux::I16(i16::MIN)).unwrap();
    test_record.push_aux(b"XF", Aux::I16(i16::MAX)).unwrap();
    test_record.push_aux(b"XG", Aux::U16(u16::MIN)).unwrap();
    test_record.push_aux(b"XH", Aux::U16(u16::MAX)).unwrap();
    test_record.push_aux(b"XI", Aux::I32(i32::MIN)).unwrap();
    test_record.push_aux(b"XJ", Aux::I32(i32::MAX)).unwrap();
    test_record.push_aux(b"XK", Aux::U32(u32::MIN)).unwrap();
    test_record.push_aux(b"XL", Aux::U32(u32::MAX)).unwrap();
    test_record
        .push_aux(b"XM", Aux::Float(std::f32::consts::PI))
        .unwrap();
    test_record
        .push_aux(b"XN", Aux::Double(std::f64::consts::PI))
        .unwrap();
    test_record
        .push_aux(b"XO", Aux::String("Test str"))
        .unwrap();
    test_record.push_aux(b"XP", Aux::I8(0)).unwrap();

    let collected_aux_fields = test_record.aux_iter().collect::<Result<Vec<_>>>().unwrap();
    assert_eq!(
        collected_aux_fields,
        vec![
            (&b"XA"[..], Aux::I8(i8::MIN)),
            (&b"XB"[..], Aux::I8(i8::MAX)),
            (&b"XC"[..], Aux::U8(u8::MIN)),
            (&b"XD"[..], Aux::U8(u8::MAX)),
            (&b"XE"[..], Aux::I16(i16::MIN)),
            (&b"XF"[..], Aux::I16(i16::MAX)),
            (&b"XG"[..], Aux::U16(u16::MIN)),
            (&b"XH"[..], Aux::U16(u16::MAX)),
            (&b"XI"[..], Aux::I32(i32::MIN)),
            (&b"XJ"[..], Aux::I32(i32::MAX)),
            (&b"XK"[..], Aux::U32(u32::MIN)),
            (&b"XL"[..], Aux::U32(u32::MAX)),
            (&b"XM"[..], Aux::Float(std::f32::consts::PI)),
            (&b"XN"[..], Aux::Double(std::f64::consts::PI)),
            (&b"XO"[..], Aux::String("Test str")),
            (&b"XP"[..], Aux::I8(0)),
        ]
    );
}

#[test]
fn test_aux_array_partial_eq() {
    use record::AuxArray;

    // Target types
    let one_data: Vec<i8> = vec![0, 1, 2, 3, 4, 5, 6];
    let one_aux_array = AuxArray::from(&one_data);

    let two_data: Vec<i8> = vec![0, 1, 2, 3, 4, 5];
    let two_aux_array = AuxArray::from(&two_data);

    assert_ne!(&one_data, &two_data);
    assert_ne!(&one_aux_array, &two_aux_array);

    let one_aux = Aux::ArrayI8(one_aux_array);
    let two_aux = Aux::ArrayI8(two_aux_array);
    assert_ne!(&one_aux, &two_aux);

    // Raw bytes
    let bam_header = Header::new();
    let mut test_record = Record::from_sam(
        &HeaderView::from_header(&bam_header),
        "ali1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tFFFF".as_bytes(),
    )
    .unwrap();

    test_record.push_aux(b"XA", one_aux).unwrap();
    test_record.push_aux(b"XB", two_aux).unwrap();

    // RawLeBytes == RawLeBytes
    assert_eq!(
        test_record.aux(b"XA").unwrap(),
        test_record.aux(b"XA").unwrap()
    );
    // RawLeBytes != RawLeBytes
    assert_ne!(
        test_record.aux(b"XA").unwrap(),
        test_record.aux(b"XB").unwrap()
    );

    // RawLeBytes == TargetType
    assert_eq!(
        test_record.aux(b"XA").unwrap(),
        Aux::ArrayI8((&one_data).into())
    );
    assert_eq!(
        test_record.aux(b"XB").unwrap(),
        Aux::ArrayI8((&two_data).into())
    );
    // RawLeBytes != TargetType
    assert_ne!(
        test_record.aux(b"XA").unwrap(),
        Aux::ArrayI8((&two_data).into())
    );
    assert_ne!(
        test_record.aux(b"XB").unwrap(),
        Aux::ArrayI8((&one_data).into())
    );
}

/// Test if both text and binary representations of a BAM header are in sync (#156)
#[test]
fn test_bam_header_sync() {
    let reader = Reader::from_path("test/test_issue_156_no_text.bam").unwrap();
    let header_hashmap = Header::from_template(reader.header()).to_hashmap().unwrap();
    let header_refseqs = header_hashmap.get("SQ").unwrap();

    assert_eq!(header_refseqs[0].get("SN").unwrap(), "ref_1",);
    assert_eq!(header_refseqs[0].get("LN").unwrap(), "10000000",);
}

#[test]
fn test_bam_new() {
    // Create the path to write the tmp test BAM
    let tmp = tempfile::Builder::new()
        .prefix("rust-htslib")
        .tempdir()
        .expect("Cannot create temp dir");
    let bampath = tmp.path().join("test.bam");

    // write an unmapped BAM record (uBAM)
    {
        // Build the header
        let mut header = Header::new();

        // Add the version
        header.push_record(
            HeaderRecord::new(b"HD")
                .push_tag(b"VN", "1.6")
                .push_tag(b"SO", "unsorted"),
        );

        // Build the writer
        let mut writer = Writer::from_path(&bampath, &header, Format::Bam).unwrap();

        // Build an empty record
        let record = Record::new();

        // Write the record (this previously seg-faulted)
        assert!(writer.write(&record).is_ok());
    }

    // Read the record
    {
        // Build th reader
        let mut reader = Reader::from_path(bampath).expect("Error opening file.");

        // Read the record
        let mut rec = Record::new();
        match reader.read(&mut rec) {
            Some(r) => r.expect("Failed to read record."),
            None => panic!("No record read."),
        };

        // Check a few things
        assert!(rec.is_unmapped());
        assert_eq!(rec.tid(), -1);
        assert_eq!(rec.pos(), -1);
        assert_eq!(rec.mtid(), -1);
        assert_eq!(rec.mpos(), -1);
    }
}

#[test]
fn test_idxstats_bam() {
    let mut reader = IndexedReader::from_path("test/test.bam").unwrap();
    let expected = vec![
        (0, 15072423, 6, 0),
        (1, 15279345, 0, 0),
        (2, 13783700, 0, 0),
        (3, 17493793, 0, 0),
        (4, 20924149, 0, 0),
        (-1, 0, 0, 0),
    ];
    let actual = reader.index_stats().unwrap();
    assert_eq!(expected, actual);
}

#[test]
fn test_number_mapped_and_unmapped_bam() {
    let reader = IndexedReader::from_path("test/test.bam").unwrap();
    let expected = (6, 0);
    let actual = reader.index().number_mapped_unmapped(0);
    assert_eq!(expected, actual);
}

#[test]
fn test_number_unmapped_global_bam() {
    let reader = IndexedReader::from_path("test/test_unmapped.bam").unwrap();
    let expected = 8;
    let actual = reader.index().number_unmapped();
    assert_eq!(expected, actual);
}

#[test]
fn test_idxstats_cram() {
    let mut reader = IndexedReader::from_path("test/test_cram.cram").unwrap();
    reader.set_reference("test/test_cram.fa").unwrap();
    let expected = vec![
        (0, 120, 2, 0),
        (1, 120, 2, 0),
        (2, 120, 2, 0),
        (-1, 0, 0, 0),
    ];
    let actual = reader.index_stats().unwrap();
    assert_eq!(expected, actual);
}

#[test]
fn test_slow_idxstats_cram() {
    let mut reader = IndexedReader::from_path("test/test_cram.cram").unwrap();
    reader.set_reference("test/test_cram.fa").unwrap();
    let expected = vec![
        (0, 120, 2, 0),
        (1, 120, 2, 0),
        (2, 120, 2, 0),
        (-1, 0, 0, 0),
    ];
    let actual = reader.index_stats().unwrap();
    assert_eq!(expected, actual);
}

#[test]
fn test_slow_idxstats_cram_unmapped() {
    let mut reader = IndexedReader::from_path("test/test_cram_unmapped.cram").unwrap();
    reader.set_reference("test/test_cram.fa").unwrap();
    let expected = vec![
        (0, 120, 2, 0),
        (1, 120, 2, 0),
        (2, 120, 2, 0),
        (-1, 0, 0, 2),
    ];
    let actual = reader.index_stats().unwrap();
    assert_eq!(expected, actual);
}

#[test]
fn test_nonexistent_tidname() {
    let header = Header::new();
    let header_view = HeaderView::from_header(&header);
    assert_eq!(b"", header_view.tid2name(0));
}

// #[test]
// fn test_number_mapped_and_unmapped_cram() {
//     let mut reader = IndexedReader::from_path("test/test_cram.cram").unwrap();
//     reader.set_reference("test/test_cram.fa").unwrap();
//     let expected = (2, 0);
//     let actual = reader.index().number_mapped_unmapped(0);
//     assert_eq!(expected, actual);
// }
//
// #[test]
// fn test_number_unmapped_global_cram() {
//     let mut reader = IndexedReader::from_path("test/test_unmapped.cram").unwrap();
//     let expected = 8;
//     let actual = reader.index().number_unmapped();
//     assert_eq!(expected, actual);
// }
