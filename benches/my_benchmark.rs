use criterion::Criterion;
use criterion::{criterion_group, criterion_main, BenchmarkId};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use streaming_iterator::StreamingIterator;

fn count_reads_records(reader: &mut bam::IndexedReader, should: usize) {
    let mut count = 0;
    let it = 0..reader.header().target_count();
    for tid in it {
        reader
            .fetch(tid, 0, reader.header().target_len(tid).unwrap())
            .unwrap();
        count += reader.records().count();
    }
    assert!(count == should);
}

fn count_reads_rc_records(reader: &mut bam::IndexedReader, should: usize) {
    let mut count = 0;
    let it = 0..reader.header().target_count();
    for tid in it {
        reader
            .fetch(tid, 0, reader.header().target_len(tid).unwrap())
            .unwrap();
        count += reader.rc_records().count();
    }
    assert!(count == should);
}

fn count_reads_fetch_iter(reader: &mut bam::IndexedReader, should: usize) {
    let mut count = 0;
    let it = 0..reader.header().target_count();
    for tid in it {
        count += reader.fetch_iter(tid).unwrap().count();
    }
    assert!(count == should);
}

fn count_reads_read(reader: &mut bam::IndexedReader, should: usize) {
    let mut count = 0;
    let it = 0..reader.header().target_count();
    let mut record = bam::record::Record::new();

    for tid in it {
        reader
            .fetch(tid, 0, reader.header().target_len(tid).unwrap())
            .unwrap();
        while let Ok(true) = reader.read(&mut record) {
            count += 1;
        }
    }
    assert!(count == should);
}

fn count_reads_stream(reader: &mut bam::IndexedReader, should: usize) {
    let mut count = 0;
    let it = 0..reader.header().target_count();
    for tid in it {
        reader
            .fetch(tid, 0, reader.header().target_len(tid).unwrap())
            .unwrap();
        count += reader.stream().count();
    }
    assert!(count > 0);
}

fn filter_reads_records(reader: &mut bam::IndexedReader, should: usize) {
    let mut count = 0;
    let it = 0..reader.header().target_count();
    for tid in it {
        reader
            .fetch(tid, 0, reader.header().target_len(tid).unwrap())
            .unwrap();
        count = reader
            .records()
            .filter(|x| x.as_ref().unwrap().is_reverse())
            .count();
        //count += reader.records().count();
    }
    assert!(count > 0);
}

fn filter_reads_rc_records(reader: &mut bam::IndexedReader, should: usize) {
    let mut count = 0;
    let it = 0..reader.header().target_count();
    for tid in it {
        reader
            .fetch(tid, 0, reader.header().target_len(tid).unwrap())
            .unwrap();
        count = reader
            .rc_records()
            .filter(|x| x.as_ref().unwrap().is_reverse())
            .count();
    }
    assert!(count > 0);
}

fn filter_reads_read(reader: &mut bam::IndexedReader, should: usize) {
    let mut count = 0;
    let it = 0..reader.header().target_count();
    let mut record = bam::record::Record::new();

    for tid in it {
        reader
            .fetch(tid, 0, reader.header().target_len(tid).unwrap())
            .unwrap();
        while let Ok(true) = reader.read(&mut record) {
            if record.is_reverse() {
                count += 1;
            }
        }
    }
    assert!(count > 0);
}

fn filter_reads_stream(reader: &mut bam::IndexedReader, should: usize) {
    let mut count = 0;
    let it = 0..reader.header().target_count();
    for tid in it {
        reader
            .fetch(tid, 0, reader.header().target_len(tid).unwrap())
            .unwrap();
        count = reader.stream().filter(|x| x.is_reverse()).count();
    }
    assert!(count > 0);
}

fn criterion_benchmarks_from_input(
    c: &mut Criterion,
    group_name: &str,
    inputs: Vec<(&str, usize)>,
) {
    let mut group = c.benchmark_group(group_name);
    group.sample_size(100);
    for i in inputs.iter() {
        group.bench_with_input(BenchmarkId::new("records_count", i.1), i, move |b, i| {
            let mut reader = bam::IndexedReader::from_path(i.0).unwrap();
            b.iter(move || count_reads_records(&mut reader, i.1))
        });
        group.bench_with_input(BenchmarkId::new("read_count", i.1), i, move |b, i| {
            let mut reader = bam::IndexedReader::from_path(i.0).unwrap();
            b.iter(move || count_reads_read(&mut reader, i.1))
        });
        group.bench_with_input(BenchmarkId::new("stream_count", i.1), i, move |b, i| {
            let mut reader = bam::IndexedReader::from_path(i.0).unwrap();
            b.iter(move || count_reads_stream(&mut reader, i.1))
        });
        group.bench_with_input(BenchmarkId::new("rcrecords_count", i.1), i, move |b, i| {
            let mut reader = bam::IndexedReader::from_path(i.0).unwrap();
            b.iter(move || count_reads_rc_records(&mut reader, i.1))
        });
        group.bench_with_input(BenchmarkId::new("fetch_iter_count", i.1), i, move |b, i| {
            let mut reader = bam::IndexedReader::from_path(i.0).unwrap();
            b.iter(move || count_reads_fetch_iter(&mut reader, i.1))
        });
        /*
            group.bench_with_input(BenchmarkId::new("records_filter", i.1), i, move |b, i| {
                let mut reader = bam::IndexedReader::from_path(i.0).unwrap();
                b.iter(move || filter_reads_records(&mut reader, i.1))
            });
            group.bench_with_input(BenchmarkId::new("read_filter", i.1), i, move |b, i| {
                let mut reader = bam::IndexedReader::from_path(i.0).unwrap();
                b.iter(move || filter_reads_read(&mut reader, i.1))
            });
            group.bench_with_input(BenchmarkId::new("stream_filter", i.1), i, move |b, i| {
                let mut reader = bam::IndexedReader::from_path(i.0).unwrap();
                b.iter(move || filter_reads_stream(&mut reader, i.1))
            });
            group.bench_with_input(BenchmarkId::new("rcrecords_filter", i.1), i, move |b, i| {
                let mut reader = bam::IndexedReader::from_path(i.0).unwrap();
                b.iter(move || filter_reads_rc_records(&mut reader, i.1))
            });
        */
    }
}

fn criterion_benchmarks_no_seek(c: &mut Criterion) {
    let inputs = vec![
        ("test/bench_1k.bam", 1000),
        ("test/bench_50k.bam", 50000),
        ("test/bench_100k.bam", 100000),
    ];
    criterion_benchmarks_from_input(c, "No_seek", inputs);
}

fn criterion_benchmarks_seek(c: &mut Criterion) {
    let inputs = vec![
        ("test/bench_seek_50k.bam", 50000),
        ("test/bench_seek_10k.bam", 10000),
    ];
    criterion_benchmarks_from_input(c, "seek", inputs);
}

fn criterion_benchmarks(c: &mut Criterion) {
    criterion_benchmarks_no_seek(c);
    //criterion_benchmarks_seek(c);
}

criterion_group!(benches, criterion_benchmarks);
criterion_main!(benches);
