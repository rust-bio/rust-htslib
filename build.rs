// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

extern crate bindgen;
extern crate fs_utils;

use fs_utils::copy::copy_directory;

use std::process::Command;
use std::env;
use std::path::PathBuf;

fn main() {
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    if !out.join("htslib").exists() {
        copy_directory("htslib", &out).unwrap();
    }

    if Command::new("make")
        .current_dir(out.join("htslib"))
        .arg("CFLAGS=-g -Wall -O2 -fPIC")
        .arg("lib-static")
        .arg("-B")
        .status()
        .unwrap()
        .success() != true
    {
        panic!("failed to build htslib");
    }

    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
        .blacklist_type("max_align_t")
        .generate()
        .expect("Unable to generate bindings.");
    bindings
        .write_to_file(out.join("bindings.rs"))
        .expect("Could not write bindings.");

    println!(
        "cargo:rustc-flags=-L {out}/htslib -L {out} -l static=hts -l z -l lzma -l bz2",
        out = out.to_str().unwrap()
    );
}
