[![Crates.io](https://img.shields.io/crates/d/rust-htslib.svg?style=flat-square)](https://crates.io/crates/rust-htslib)
[![Crates.io](https://img.shields.io/crates/v/rust-htslib.svg?style=flat-square)](https://crates.io/crates/rust-htslib)
[![Crates.io](https://img.shields.io/crates/l/rust-htslib.svg?style=flat-square)](https://crates.io/crates/rust-htslib)
[![Travis](https://img.shields.io/travis/rust-bio/rust-htslib.svg?style=flat-square)](https://travis-ci.org/rust-bio/rust-htslib)

# HTSlib bindings for Rust

This library provides HTSlib bindings and a high level Rust API for reading and writing BAM files.

To clone this repository, issue

```
git clone --recursive https://github.com/rust-bio/rust-htslib.git
```

ensuring that the HTSlib submodule is fetched, too.
If you only want to use the library, there is no need to clone the repository. Go on to the **Usage** section in this case.

# Resources

* Homepage: https://github.com/rust-bio/rust-htslib
* API documentation: http://rust-bio.github.io/rust-htslib
* Continuous integration tests: https://travis-ci.org/rust-bio/rust-htslib

# Usage

To use Rust-HTSlib in your Rust project, import the crate from your source code:

```rust
extern crate htslib;
use htslib::bam;
```

Rust-HTSlib provides a high level BAM API.
Reading and writing BAM files is as easy as
```rust
let bam = bam::Reader::new("path/to/some.bam").ok().expect("Error opening bam.");
let out = bam::Writer::with_template("path/to/some.bam", "reversereads.bam").ok().expect("Error opening bam.");

// copy reverse reads to new BAM file
for r in bam.records() {
    let record = r.ok().expect("Error reading BAM file.");
    if record.is_reverse() {
        out.write(record);
    }
}
```

Pileups can be performed with
```rust
let bam = bam::Reader::new("path/to/some.bam").ok().expect("Error opening bam.");

// pileup over all covered sites
for p in bam.pileup() {
    let pileup = p.ok().expect("Error reading BAM file.");
    println!("{}:{} depth {}", bam.pileup.tid(), pileup.pos(), pileup.depth());

    for alignment in pileup.alignments() {
        match alignment.indel() {
            bam::pileup::Indel::Ins(len) => println!("Insertion of length {}", len),
            bam::pileup::Indel::Del(len) => println!("Deletion of length {}", len),
            _ => println!("Base {}", alignment.record.seq()[alignment.qpos()])
        }
    }
}
```
In both cases, indexed BAM files can be seeked for specific regions, constraining either the record iterator or the pileups:

```rust
let bam = bam::IndexedReader::new("path/to/some.bam");

// seek to chr1:50000-100000
bam.seek(bam.header.tid(b"chr1"), 50000, 100000).ok().expect("Error seeking BAM file.");
// afterwards, read or pileup in this region
```

# Authors

* [Christopher Schröder](https://github.com/christopher-schroeder)
* [Johannes Köster](https://github.com/johanneskoester)

# Contributors

* [Fedor Gusev](https://github.com/gusevfe)


## License

Licensed under the MIT license http://opensource.org/licenses/MIT. This project may not be copied, modified, or distributed except according to those terms.
