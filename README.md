# HTSlib bindings for Rust

This library provides HTSlib bindings and a high level Rust API for reading and writing SAM/BAM files.

To clone this repository, issue

```
git clone --recursive https://github.com/christopher-schroeder/rust-htslib.git
```

ensuring that the HTSlib submodule is fetched, too.

# Usage

To use Rust-HTSlib in your Rust project, import the crate from your source code:

```rust
extern crate htslib;
use htslib::bam;
```

Rust-HTSlib provides a high level BAM API.
Reading and writing BAM files is as easy as
```rust
// copy reverse reads to new BAM file
let bam = bam::Reader::new("path/to/some.bam");
let out = bam::Writer::with_template("path/to/some.bam", "reversereads.bam");

for record in bam.records() {
    if record.is_reverse() {
        out.write(record);
    }
}
```

Pileups can be performed with
```rust
let bam = bam::Reader::new("path/to/some.bam");

for pileup in bam.pileup() {
    for alignment in pileup.alignments() {
        match alignment.indel() {
            bam::pileup::Indel::Ins(len) => println!("Insertion of length {}", len),
            bam::pileup::Indel::Del(len) => println!("Deletion of length {}", len),
            _ => println!("Base {}", alignment.record.seq()[alignment.qpos()])
        }
    }
}
```


# Authors

* Christopher Schröder (christopher.schroeder@tu-dortmund.de)
* Johannes Köster (johannes.koester@tu-dortmund.de)


## License

Licensed under the MIT license http://opensource.org/licenses/MIT. This project may not be copied, modified, or distributed except according to those terms.
