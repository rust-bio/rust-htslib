// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


use std::process::Command;
use std::env;

fn main() {
    let root = env::var("CARGO_MANIFEST_DIR").unwrap();

    Command::new("make").current_dir("htslib")
                        .arg("CFLAGS=-g -Wall -O2 -fPIC")
                        .arg("lib-static")
                        .status().ok().expect("Failed to build htslib");

    println!("cargo:rustc-flags=-L {}/htslib -l static=hts -l z", root);
}
