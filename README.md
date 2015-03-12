# HTSlib bindings for Rust

This library provides HTSlib bindings and a high level Rust API for reading and writing SAM/BAM files.

To clone this repository, issue

```
git clone --recursive https://github.com/christopher-schroeder/rust-htslib.git
```

ensuring that the HTSlib submodule is fetched, too.

# Usage

To use rust-bio in your Rust project, import the crate from your source code:
```rust
extern create htslib;
use htslib::bam;

// copy reverse reads to new BAM file
let bam = bam::Reader::new("path/to/some.bam");
let out = bam::Writer::with_template("path/to/some.bam", "reversereads.bam");

for record in bam.records() {
    if record.is_reverse() {
        out.write(record);
    }
}
```


# Authors

* Christopher Schröder (christopher.schroeder@tu-dortmund.de)
* Johannes Köster (johannes.koester@tu-dortmund.de)


## License

Licensed under the MIT license http://opensource.org/licenses/MIT. This project may not be copied, modified, or distributed except according to those terms.
