// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

#![feature(rust_2018_preview)]

extern crate bindgen;
extern crate fs_utils;

use fs_utils::copy::copy_directory;

use std::env;
use std::path::PathBuf;
use std::process::Command;

fn sed_htslib_makefile(out: &PathBuf, patterns: &Vec<&str>, feature: &str) {
    for pattern in patterns {
        if Command::new("sed")
            .current_dir(out.join("htslib"))
            .arg("-i")
            .arg("-e")
            .arg(pattern)
            .arg("Makefile")
            .status()
            .unwrap()
            .success() != true
        {
            panic!("failed to strip {} support", feature);
        }
    }
}

fn main() {
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    if !out.join("htslib").exists() {
        copy_directory("htslib", &out).unwrap();
    }

    let use_bzip2 = env::var("CARGO_FEATURE_BZIP2").is_ok();
    if !use_bzip2 {
        let bzip2_patterns = vec!["s/ -lbz2//", "/#define HAVE_LIBBZ2/d"];
        sed_htslib_makefile(&out, &bzip2_patterns, "bzip2");
    }

    let use_lzma = env::var("CARGO_FEATURE_LZMA").is_ok();
    if !use_lzma {
        let lzma_patterns = vec!["s/ -llzma//", "/#define HAVE_LIBLZMA/d"];
        sed_htslib_makefile(&out, &lzma_patterns, "lzma");
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
        .generate_comments(false)
        .blacklist_type("max_align_t")
        .generate()
        .expect("Unable to generate bindings.");
    bindings
        .write_to_file(out.join("bindings.rs"))
        .expect("Could not write bindings.");

    println!(
        "cargo:rustc-flags=-L {out}/htslib -L {out} -l static=hts -l z {lzma} {bzip2}",
        out = out.to_str().unwrap(),
        lzma = if use_lzma { "-l lzma" } else { "" },
        bzip2 = if use_bzip2 { "-l bz2" } else { "" },
    );
}
